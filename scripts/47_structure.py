#!/usr/bin/env python3
"""Structural hardening and coherence of the dark proteome, without
predicted coordinates.

    python scripts/47_structure.py                  # discover and plan
    python scripts/47_structure.py --go --stage db
    python scripts/47_structure.py --go --stage search
    python scripts/47_structure.py --go --stage cluster

WHY THIS IS NOW POSSIBLE

  The manuscript records, as its largest caveat, that the homology-free
  subset was hardened against sequence search only, because structural
  hardening needed ESMFold and ESMFold was not available. That was true
  of predicted coordinates and false of the thing we actually wanted.
  Foldseek takes a ProstT5 model directly (`--prostt5-model`) and
  predicts the 3Di structural alphabet from amino-acid sequence, which
  is enough for monomer search and clustering. No coordinates, no
  ESMFold, and ProstT5 is already cached for this project.

TWO SEPARATE QUESTIONS, TWO SEPARATE DESIGNS

  1. HARDENING — do the proteins with no sequence hit have a structural
     relative that sequence search missed?

     The search must run against THE SAME LOOKUP SET that DIAMOND and
     jackhmmer searched. Searching PDB instead would change the method
     and the database together, and any hit found could be explained by
     either. The lookup set is the control; the modality is the variable.

  2. COHERENCE — do the dark proteins share folds with EACH OTHER more
     than comparable proteins do?

     A clustered fraction on its own means nothing: short proteins,
     low-complexity proteins and compositionally biased proteins cluster
     more readily whatever their biology, and the dark set is enriched
     for all three. So the dark set is compared against a background of
     characterised proteins sampled to match its length distribution,
     clustered under identical settings. The comparison is the result;
     the raw number is not.

WHAT THIS CANNOT SETTLE

  A 3Di string from ProstT5 is a prediction, not a structure, and a
  structural hit resting on it is weaker evidence than one resting on
  coordinates. It is also not independent of the annotation work: the
  same ProstT5 model is one of the three used for transfer. If ProstT5
  represents a protein badly, its embedding transfer and its 3Di
  prediction fail together, and this check would not notice. The
  headline model is ESM-2, so the hardening is independent of THAT
  arm, and the report states the dependence rather than burying it.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod                              # noqa: E402

TRANSFER_EVALUE = 1e-3        # the same cut the sequence baselines used


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def human(n: float) -> str:
    for u in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or u == "TB":
            return f"{n:,.1f} {u}" if u != "B" else f"{int(n)} B"
        n /= 1024
    return f"{n:.1f} TB"


# ---------------------------------------------------------------------------
# discovery
# ---------------------------------------------------------------------------

def find_foldseek() -> str | None:
    p = shutil.which("foldseek")
    if p:
        return p
    cfg = cfgmod.load()
    work = cfg.p("work")
    # the static tarball, which is where this ends up when conda refuses
    for c in (work / "foldseek" / "bin" / "foldseek",
              Path.home() / "pmdark_work" / "foldseek" / "bin" / "foldseek",
              Path.home() / "foldseek" / "bin" / "foldseek"):
        if c.is_file() and os.access(c, os.X_OK):
            return str(c)
    for env in ("pmbio", "busco", "serprot", "esmprobe"):
        for base in (work / "envs" / env,
                     Path.home() / "miniconda3" / "envs" / env):
            c = base / "bin" / "foldseek"
            if c.is_file():
                return str(c)
    return None


def read_fasta_ids(path: Path) -> list[str]:
    ids = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def fasta_lengths(path: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    cur, n = None, 0
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    out[cur] = n
                cur, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if cur:
        out[cur] = n
    return out


def write_subset(src: Path, want: set[str], dst: Path) -> int:
    """Copy the named records. Streams, so a large FASTA costs nothing."""
    n, keep = 0, False
    with src.open() as fh, dst.open("w") as out:
        for line in fh:
            if line.startswith(">"):
                keep = line[1:].split()[0] in want
                n += keep
            if keep:
                out.write(line)
    return n


ITER_SUFFIX = re.compile(r"-i\d+$")


def parse_tbl(path: Path, max_e: float, drop_self: bool = True
              ) -> tuple[dict[str, float], dict]:
    """jackhmmer --tblout -> ({query: best E-value}, diagnostics).

    Two traps, both silent, both hit in this project:

    1. jackhmmer appends an iteration suffix to the query name (`-i2`,
       `-i3`). Left on, every lookup against the orphan list misses,
       the hardened set comes back equal to the orphan set, and the
       check passes while testing nothing.

    2. A query hits ITSELF at E=0 whenever the search database contains
       the queries. Counting that, every orphan looks rescued and the
       hardened set empties. parse_m8 dropped self-hits from the first
       version; this did not, which is how `hardened set : 0` happened.

    The diagnostics are returned rather than logged so the caller can
    print them: a derivation that collapses to 0 or to n must show its
    working, not just its answer.
    """
    best: dict[str, float] = {}
    d = {"rows": 0, "cols": None, "first": None, "self_rows": 0,
         "parse_fail": 0, "queries": set(), "below_cut": 0}
    if not path.exists():
        return best, d
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if d["cols"] is None:
                d["cols"], d["first"] = len(f), line.rstrip()[:160]
            if len(f) < 5:
                continue
            d["rows"] += 1
            target = f[0]
            q = ITER_SUFFIX.sub("", f[2])
            d["queries"].add(q)
            try:
                e = float(f[4])
            except ValueError:
                d["parse_fail"] += 1
                continue
            if ITER_SUFFIX.sub("", target) == q:
                d["self_rows"] += 1
                if drop_self:
                    continue
            if e <= max_e:
                d["below_cut"] += 1
                if q not in best or e < best[q]:
                    best[q] = e
    d["queries"] = len(d["queries"])
    return best, d


def find_tbl(root: Path) -> Path | None:
    cands = [p for p in root.rglob("*")
             if p.is_file() and any(p.name.endswith(x) for x in
                                    (".tbl", ".tblout", ".tbl.txt",
                                     "_tbl.txt", ".domtbl", ".domtblout"))]
    return max(cands, key=lambda p: p.stat().st_size) if cands else None


def discover(cfg, res: Path) -> dict:
    """Locate every input, and say plainly what is missing.

    Nothing here guesses. A set that cannot be built from a file on disk
    is reported as absent, because a structural claim resting on a set
    assembled from the wrong identifiers would be worse than no claim.
    """
    found: dict = {"fasta": {}, "ids": {}, "notes": []}

    prot = cfg.p("proteins")
    cands = sorted(
        [p for p in prot.rglob("*.f*a") if p.is_file()] +
        [p for p in prot.rglob("*.fasta") if p.is_file()],
        key=lambda p: -p.stat().st_size)
    seen = set()
    cands = [p for p in cands if not (p in seen or seen.add(p))]
    found["fasta"]["candidates"] = [(str(p), p.stat().st_size) for p in cands[:12]]

    # Exact paths first. Hint matching picked proteins_novel_exonic.faa
    # as the query set, which is a different population entirely — the
    # benchmark was built from benchmark/query.faa, and scoring the
    # wrong set would have been invisible in the output.
    exact = {"query": prot / "benchmark" / "query.faa",
             "lookup": prot / "benchmark" / "reference.faa",
             "orphan": prot / "benchmark" / "orphan_check" / "orphans.faa",
             "targets": prot / "annotate" / "targets.faa"}
    for key, p in exact.items():
        found["fasta"][key] = str(p) if p.is_file() else None
    for key, hints in (("query", ("query", "transdecoder", "pep")),
                       ("lookup", ("reference", "lookup", "swiss"))):
        if found["fasta"].get(key):
            continue
        for p in cands:
            if any(h in p.name.lower() for h in hints):
                found["fasta"][key] = str(p)
                break

    # the hardened set is derived, not stored: orphans with no profile hit
    found["tbl"] = None
    oc = prot / "benchmark" / "orphan_check"
    if oc.is_dir():
        t = find_tbl(oc)
        found["tbl"] = str(t) if t else None

    # identifier lists for the three sets we need
    for name, files, keys in (
            ("hardened", ("orphan_hardening.json",),
             ("strata.HARDENED.ids", "hardened_ids", "ids")),
            ("orphan", ("orphan_hardening.json", "bootstrap_ci.json"),
             ("n_orphan_diamond_ids", "orphan_ids", "ids")),
            ("targets", ("annotation_summary.json", "target_annotations.json"),
             ("target_ids", "ids", "targets"))):
        got = None
        for f in files:
            p = res / f
            if not p.exists():
                continue
            try:
                blob = json.loads(p.read_text())
            except ValueError:
                continue
            for k in keys:
                node = blob
                for part in k.split("."):
                    node = node.get(part) if isinstance(node, dict) else None
                    if node is None:
                        break
                if isinstance(node, list) and node and isinstance(node[0], str):
                    got = (f, k, len(node))
                    break
            if got:
                break
        found["ids"][name] = got

    # a plain text list is the fallback, and the most reliable route
    idsdir = res / "id_lists"
    found["ids"]["text_files"] = (
        sorted(str(p) for p in idsdir.glob("*.txt")) if idsdir.is_dir() else [])
    return found


# ---------------------------------------------------------------------------
# matched background
# ---------------------------------------------------------------------------

def match_by_length(dark: list[str], pool: list[str],
                    lengths: dict[str, int], seed: int,
                    per: int = 1) -> list[str]:
    """Sample a background matched on length, by decile.

    Short proteins cluster more readily than long ones, and the dark set
    skews short in a fragmented assembly, so an unmatched background
    would manufacture the coherence result it is meant to test.
    """
    rng = random.Random(seed)
    dl = sorted(lengths[i] for i in dark if i in lengths)
    if not dl:
        return []
    edges = [dl[int(q * (len(dl) - 1) / 10)] for q in range(11)]
    buckets: dict[int, list[str]] = defaultdict(list)
    for p in pool:
        L = lengths.get(p)
        if L is None:
            continue
        for b in range(10):
            if edges[b] <= L <= edges[b + 1]:
                buckets[b].append(p)
                break
    want: dict[int, int] = defaultdict(int)
    for i in dark:
        L = lengths.get(i)
        if L is None:
            continue
        for b in range(10):
            if edges[b] <= L <= edges[b + 1]:
                want[b] += 1
                break
    out: list[str] = []
    short = 0
    for b, k in want.items():
        pool_b = buckets.get(b, [])
        k = k * per
        if len(pool_b) <= k:
            out += pool_b
            short += k - len(pool_b)
        else:
            out += rng.sample(pool_b, k)
    if short:
        log(f"  background short by {short} in sparse length bins")
    return out


# ---------------------------------------------------------------------------
# foldseek
# ---------------------------------------------------------------------------

def run(cmd: list[str], cwd: Path | None = None, timeout: int = 36000):
    log("  $ " + " ".join(str(c) for c in cmd[:9]) +
        (" ..." if len(cmd) > 9 else ""))
    r = subprocess.run([str(c) for c in cmd], cwd=cwd, text=True,
                       capture_output=True, timeout=timeout)
    if r.returncode != 0:
        sys.stderr.write((r.stderr or r.stdout)[-4000:] + "\n")
    return r.returncode


def parse_m8(path: Path, max_e: float) -> dict[str, list[tuple[str, float]]]:
    """query -> [(target, evalue)], self-hits dropped."""
    out: dict[str, list[tuple[str, float]]] = defaultdict(list)
    if not path.exists():
        return out
    with path.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            q, t = f[0], f[1]
            if q == t:
                continue
            try:
                e = float(f[10])
            except ValueError:
                continue
            if e <= max_e:
                out[q].append((t, e))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--go", action="store_true")
    ap.add_argument("--stage", default="all",
                    choices=["db", "search", "cluster", "all"])
    ap.add_argument("--query-fasta", default=None,
                    help="protein FASTA holding the query set")
    ap.add_argument("--lookup-fasta", default=None,
                    help="the annotated lookup set DIAMOND and jackhmmer "
                         "searched; the hardening search uses this and "
                         "not PDB, so only the modality changes")
    ap.add_argument("--hardened-ids", default=None,
                    help="text file, one identifier per line")
    ap.add_argument("--target-ids", default=None)
    ap.add_argument("--orphan-fasta", default=None,
                    help="FASTA of the DIAMOND-orphan set; the hardened "
                         "subset is derived from it and the profile table")
    ap.add_argument("--jackhmmer-tbl", default=None,
                    help="jackhmmer --tblout from 37_harden_orphans.py")
    ap.add_argument("--characterised-ids", default=None,
                    help="pool the matched background is drawn from")
    ap.add_argument("--evalue", type=float, default=TRANSFER_EVALUE)
    ap.add_argument("--cluster-c", type=float, default=0.8,
                    help="foldseek coverage for clustering (default 0.8)")
    ap.add_argument("--gpu", type=int, default=1)
    ap.add_argument("--rebuild-db", action="store_true",
                    help="rebuild the lookup database even if it exists")
    ap.add_argument("--seed", type=int, default=None)
    args = ap.parse_args()

    cfg = cfgmod.load()
    res = cfg.p("results").resolve()
    out = cfg.p("structures").resolve()
    seed = args.seed if args.seed is not None else cfg.raw.get("seed", 1)

    fs = find_foldseek()
    print("=" * 78)
    print("STRUCTURAL HARDENING AND COHERENCE — plan")
    print("=" * 78)
    print(f"  foldseek     : {fs or 'NOT FOUND'}")
    if not fs:
        print("\n  install it into the pmbio environment:\n")
        print("    conda install -n pmbio -c conda-forge -c bioconda "
              "foldseek\n")
        print("  or take the static binary, which needs no conda:\n")
        print("    wget https://mmseqs.com/foldseek/foldseek-linux-gpu.tar.gz")
        print("    tar xzf foldseek-linux-gpu.tar.gz")
        print("    export PATH=$PWD/foldseek/bin:$PATH\n")

    d = discover(cfg, res)
    print(f"\n  protein FASTA candidates under {cfg.p('proteins')}:")
    if not d["fasta"]["candidates"]:
        print("    none found")
    for p, sz in d["fasta"]["candidates"]:
        print(f"    {human(sz):>10}  {p}")
    print("\n  resolved inputs:")
    for k in ("query", "lookup", "orphan", "targets"):
        v = d["fasta"].get(k)
        print(f"    {k:<9}{v or '-- none --'}")
    print(f"    tblout   {d.get('tbl') or '-- none --'}")

    print("\n  identifier lists:")
    for k in ("hardened", "orphan", "targets"):
        g = d["ids"].get(k)
        print(f"    {k:<10}" + (f"{g[2]} ids from {g[0]} [{g[1]}]" if g
                                else "NOT FOUND in the result JSONs"))
    if d["ids"]["text_files"]:
        print("    text lists found:")
        for t in d["ids"]["text_files"]:
            print(f"      {t}")

    qf = Path(args.query_fasta) if args.query_fasta else (
        Path(d["fasta"]["query"]) if d["fasta"]["query"] else None)
    lf = Path(args.lookup_fasta) if args.lookup_fasta else (
        Path(d["fasta"]["lookup"]) if d["fasta"]["lookup"] else None)

    print(f"\n  output dir   : {out}")
    print(f"  E-value      : {args.evalue:g}  (same cut as the sequence "
          f"baselines)")
    print(f"  cluster cov  : {args.cluster_c}")
    print(f"  seed         : {seed}")
    print("=" * 78)

    missing = []
    if not fs:
        missing.append("foldseek")
    if not qf or not qf.is_file():
        missing.append("--query-fasta")
    if not lf or not lf.is_file():
        missing.append("--lookup-fasta")
    # ---- the three sets -------------------------------------------------
    orphan_fa = Path(args.orphan_fasta) if args.orphan_fasta else (
        Path(d["fasta"]["orphan"]) if d["fasta"].get("orphan") else None)
    tbl = Path(args.jackhmmer_tbl) if args.jackhmmer_tbl else (
        Path(d["tbl"]) if d.get("tbl") else None)

    hard_ids, orphan_ids, derivation = None, None, ""
    if args.hardened_ids and Path(args.hardened_ids).is_file():
        hard_ids = [l.strip() for l in
                    Path(args.hardened_ids).read_text().split() if l.strip()]
        derivation = f"read from {args.hardened_ids}"
    elif orphan_fa and orphan_fa.is_file():
        orphan_ids = read_fasta_ids(orphan_fa)
        if tbl and tbl.is_file():
            hit, diag = parse_tbl(tbl, args.evalue)
            print(f"\n  profile table: {tbl.name}")
            print(f"    {diag['rows']} rows, {diag['cols']} columns, "
                  f"{diag['queries']} distinct queries")
            print(f"    self-hits dropped : {diag['self_rows']}")
            print(f"    hits kept at E<={args.evalue:g} : "
                  f"{diag['below_cut']}  ({len(hit)} queries)")
            if diag["parse_fail"]:
                print(f"    column 5 unparseable on {diag['parse_fail']} "
                      f"rows — this may be a domtblout, whose E-value is "
                      f"in column 7; check the first row below")
            if diag["first"]:
                print(f"    first row: {diag['first']}")
            hard_ids = [i for i in orphan_ids if i not in hit]
            derivation = (f"{len(orphan_ids)} orphans minus "
                          f"{len(orphan_ids) - len(hard_ids)} with a "
                          f"profile hit at E<={args.evalue:g}, from "
                          f"{tbl.name}")
            if len(hard_ids) == len(orphan_ids):
                print("\n  WARNING the profile table removed nothing. "
                      "Either no orphan has a hit, or the identifiers in "
                      "the table do not match the FASTA — the second is "
                      "the usual cause and it fails silently. Check a "
                      "name from each before trusting this.")
        else:
            missing.append("jackhmmer tblout (pass --jackhmmer-tbl)")
    if not hard_ids and "jackhmmer tblout (pass --jackhmmer-tbl)" not in missing:
        missing.append("--hardened-ids or --orphan-fasta")

    if hard_ids is not None:
        print(f"\n  hardened set : {len(hard_ids)}  ({derivation})")

    if missing:
        print(f"\n  CANNOT RUN — missing: {', '.join(missing)}")
        print("\n  Send me this output. The identifier lists are the part I")
        print("  cannot guess from here: if they are not inside the result")
        print("  JSONs, they can be regenerated from the DIAMOND and")
        print("  jackhmmer outputs, and I will add that step rather than")
        print("  have this script assemble a set from the wrong column.")
        return 1

    # Background pool: the benchmark's OWN non-orphan queries, not the
    # curated reference set. Same assembly, same ORF caller, same length
    # and fragmentation profile — so the only systematic difference left
    # between the two arms is whether homology finds them.
    pool: list[str] = []
    if args.characterised_ids and Path(args.characterised_ids).is_file():
        pool = [l.strip() for l in
                Path(args.characterised_ids).read_text().split() if l.strip()]
        pool_src = args.characterised_ids
    elif qf and qf.is_file():
        allq = read_fasta_ids(qf)
        drop = set(orphan_ids or hard_ids or [])
        pool = [i for i in allq if i not in drop]
        pool_src = f"{qf.name} minus the {len(drop)} orphans"

    print(f"  background   : {len(pool)} to draw from ({pool_src})")

    if not args.go:
        print("\nPLAN ONLY — nothing written. Re-run with --go.")
        return 0

    out.mkdir(parents=True, exist_ok=True)
    tmp = out / "tmp"
    tmp.mkdir(exist_ok=True)
    weights = out / "prostt5_weights"

    # ---- stage: db ------------------------------------------------------
    if args.stage in ("db", "all"):
        if not weights.exists():
            log("fetching the ProstT5 weights foldseek wants")
            if run([fs, "databases", "ProstT5", weights, tmp]) != 0:
                log("weight download failed")
                return 1
        lengths = fasta_lengths(qf)
        sets = {"hardened": set(hard_ids)}
        # The 328 application targets live in their own FASTA, not in the
        # benchmark query set, so subsetting query.faa for them finds
        # nothing. Copy that file straight across instead. Without this
        # the targets are silently never searched, which is how the first
        # run finished with "targetsdb_ss does not exist".
        tsrc = (Path(args.target_ids) if args.target_ids and
                Path(args.target_ids).is_file() and
                Path(args.target_ids).suffix in (".faa", ".fasta", ".fa")
                else (Path(d["fasta"]["targets"])
                      if d["fasta"].get("targets") else None))
        if tsrc and tsrc.is_file():
            shutil.copy2(tsrc, out / "targets.fasta")
            log(f"targets: copied {tsrc.name} "
                f"({len(read_fasta_ids(tsrc))} proteins)")
        elif args.target_ids and Path(args.target_ids).is_file():
            sets["targets"] = set(Path(args.target_ids).read_text().split())
        if pool:
            bg = match_by_length(list(sets["hardened"]), pool, lengths, seed)
            sets["background"] = set(bg)
            log(f"length-matched background: {len(bg)} proteins")
        else:
            log("NO BACKGROUND: --characterised-ids not given. The "
                "coherence test needs it; hardening does not.")
        for name, ids in sets.items():
            fa = out / f"{name}.fasta"
            n = write_subset(qf, ids, fa)
            log(f"{name}: {n}/{len(ids)} records written to {fa.name}")
            if n < len(ids):
                log(f"  WARNING {len(ids) - n} identifiers absent from "
                    f"{qf.name} — check the identifier namespace matches")
        built = (out / "lookupdb_ss.0").exists() or \
            (out / "lookupdb_ss").exists()
        if built and not args.rebuild_db:
            # Naming the stage should not destroy 15-35 minutes of work.
            # Nothing about this database changes between runs, and
            # rebuilding it because "db" appeared in a command line is a
            # trap, not a feature.
            log("lookup database already built — reusing it "
                "(--rebuild-db to force)")
        else:
            log("building the lookup database (3Di predicted by ProstT5)")
            cmd = [fs, "createdb", lf, out / "lookupdb",
                   "--prostt5-model", weights]
            if args.gpu:
                cmd += ["--gpu", "1"]
            if run(cmd) != 0:
                return 1
            log("lookup database built")

    # ---- stage: search --------------------------------------------------
    # createdb + search + convertalis rather than easy-search.
    #
    # easy-search died with "No datafile could be found for .../query_ca".
    # A ProstT5 database carries 3Di only and no C-alpha coordinates, and
    # easy-search builds its query database internally where we cannot see
    # or control it. Doing the three steps by hand means the query db is
    # built with the same ProstT5 flag as the target db, the alignment
    # type is stated rather than defaulted, and each step gets its own tmp
    # directory so a stale parameter-hash folder cannot be picked up.
    if args.stage in ("search", "all"):
        for name in ("hardened", "targets", "background"):
            fa = out / f"{name}.fasta"
            if not fa.is_file():
                continue
            qdb = out / f"{name}db"
            if not qdb.exists():
                cmd = [fs, "createdb", fa, qdb, "--prostt5-model", weights]
                if args.gpu:
                    cmd += ["--gpu", "1"]
                if run(cmd) != 0:
                    return 1
            if name == "background":
                continue            # built for clustering, not searched
            aln = out / f"aln_{name}"
            t = out / f"tmp_search_{name}"
            t.mkdir(exist_ok=True)
            if run([fs, "search", qdb, out / "lookupdb", aln, t,
                    "-e", str(args.evalue), "--alignment-type", "2"]) != 0:
                return 1
            m8 = out / f"{name}_vs_lookup.m8"
            if run([fs, "convertalis", qdb, out / "lookupdb", aln, m8,
                    "--format-output",
                    "query,target,fident,alnlen,mismatch,gapopen,qstart,"
                    "qend,tstart,tend,evalue,bits"]) != 0:
                return 1
            hits = parse_m8(m8, args.evalue)
            tot = len(read_fasta_ids(fa))
            log(f"{name}: {len(hits)}/{tot} have a structural hit at "
                f"E<={args.evalue:g}")

    # ---- stage: cluster -------------------------------------------------
    coherence = {}
    if args.stage in ("cluster", "all"):
        for name in ("hardened", "background", "targets"):
            fa = out / f"{name}.fasta"
            if not fa.is_file():
                continue
            qdb = out / f"{name}db"
            if not qdb.exists():
                cmd = [fs, "createdb", fa, qdb, "--prostt5-model", weights]
                if args.gpu:
                    cmd += ["--gpu", "1"]
                if run(cmd) != 0:
                    return 1
            cdb = out / f"clu_{name}"
            t = out / f"tmp_clu_{name}"
            t.mkdir(exist_ok=True)
            if run([fs, "cluster", qdb, cdb, t, "-c", str(args.cluster_c),
                    "-e", str(args.evalue), "--alignment-type", "2"]) != 0:
                return 1
            tsv = out / f"clu_{name}.tsv"
            if run([fs, "createtsv", qdb, qdb, cdb, tsv]) != 0:
                return 1
            members: dict[str, set[str]] = defaultdict(set)
            if tsv.exists():
                for line in tsv.read_text().splitlines():
                    f = line.split("\t")
                    if len(f) >= 2:
                        members[f[0]].add(f[1])
            # Count only members that are actually in this FASTA. A
            # cluster table built against the wrong database yields a
            # fraction above 1, which is impossible and must not be
            # reported as a result; restricting the count makes the
            # mismatch visible instead of plausible.
            want = set(read_fasta_ids(fa))
            n = len(want)
            stray = {m for v in members.values() for m in v} - want
            if stray:
                log(f"  WARNING {len(stray)} identifiers in "
                    f"{tsv.name} are not in {fa.name} — wrong database? "
                    f"excluded from the count")
            multi = {k: (v & want) for k, v in members.items()
                     if len(v & want) > 1}
            inmulti = sum(len(v) for v in multi.values())
            if inmulti > n:
                log(f"  REFUSING to report {name}: {inmulti} clustered "
                    f"of {n} present. Cluster table and FASTA disagree.")
                return 1
            coherence[name] = {
                "n": n, "clusters": len(members),
                "multi_member_clusters": len(multi),
                "proteins_in_multi": inmulti,
                "clustered_fraction": inmulti / n if n else None,
                "largest_cluster": max((len(v) for v in members.values()),
                                       default=0),
            }
            log(f"{name}: {inmulti}/{n} "
                f"({100 * inmulti / n if n else 0:.1f}%) fall in a "
                f"multi-member cluster; largest "
                f"{coherence[name]['largest_cluster']}")

    # ---- report ---------------------------------------------------------
    report = {
        "evalue": args.evalue, "cluster_c": args.cluster_c, "seed": seed,
        "query_fasta": str(qf), "lookup_fasta": str(lf),
        "n_hardened": len(hard_ids), "coherence": coherence,
        "caveat": ("3Di predicted by ProstT5, not derived from "
                   "coordinates. ProstT5 is also one of the three models "
                   "used for annotation transfer, so this check is not "
                   "independent of that arm; it is independent of the "
                   "ESM-2 result, which is the headline."),
    }
    for name in ("hardened", "targets"):
        m8 = out / f"{name}_vs_lookup.m8"
        if m8.exists():
            hits = parse_m8(m8, args.evalue)
            fa = out / f"{name}.fasta"
            tot = len(read_fasta_ids(fa)) if fa.is_file() else 0
            report[f"{name}_structural_hits"] = {
                "n": tot, "with_hit": len(hits),
                "fraction": len(hits) / tot if tot else None,
            }
    (res / "structural_hardening.json").write_text(
        json.dumps(report, indent=2) + "\n")
    log(f"wrote {res / 'structural_hardening.json'}")

    print("\n" + "=" * 78)
    h = report.get("hardened_structural_hits")
    if h and h["n"]:
        print(f"  HARDENING  {h['with_hit']} of {h['n']} sequence-orphans "
              f"have a structural relative in the same lookup set")
        print(f"             ({100 * h['fraction']:.1f}%) at "
              f"E<={args.evalue:g}")
        print(f"             still orphan after structure: "
              f"{h['n'] - h['with_hit']}")
    if "hardened" in coherence and "background" in coherence:
        a, b = coherence["hardened"], coherence["background"]
        print(f"\n  COHERENCE  dark {100 * a['clustered_fraction']:.1f}% "
              f"clustered vs length-matched background "
              f"{100 * b['clustered_fraction']:.1f}%")
        print("             (n is small; the interval comes next, in 48)")
    elif "hardened" in coherence:
        print("\n  COHERENCE  no matched background — rerun with "
              "--characterised-ids for a number that means anything")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
