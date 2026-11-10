#!/usr/bin/env python3
"""Phylostratigraphy: how far from Penaeus monodon can each protein be traced?

    python scripts/52_lineage.py                      # plan, writes nothing
    python scripts/52_lineage.py --go --stage panel   # download (network)
    python scripts/52_lineage.py --go --stage seq     # DIAMOND + strata

WHAT "LINEAGE-RESTRICTED" WOULD MEAN

  A protein whose most distant detectable relative lies inside Arthropoda
  has no traceable ancestry beyond the phylum. The claim in the working
  title was "lineage-restricted host factors", and nothing in the study
  so far tests it. This does.

  Each query is searched against a panel of species chosen to sit at
  known distances from P. monodon, and assigned the rank of the most
  distant clade it hits. The ranks are a ladder, not a tree: it is the
  FURTHEST hit that matters, because one hit in E. coli makes a protein
  ancient however many arthropod hits it also has.

THE OBJECTION THIS INVITES, AND WHY IT MATTERS HERE MORE THAN USUAL

  Phylostratigraphy infers age from the absence of detectable homologs,
  and absence of detection is not absence of homolog. Fast-evolving
  genes look young because the search fails, which is the standard
  criticism of the method.

  This paper's entire argument is that homology search fails on exactly
  these proteins. So inferring "lineage-restricted" from a failed
  sequence search would be inferring from the failure mode we spent the
  rest of the study documenting — a reviewer would notice, and would be
  right to.

  The answer is to repeat the search structurally, with the ProstT5 and
  Foldseek pipeline already built, and call a protein lineage-restricted
  only when NEITHER sequence nor structure reaches beyond the phylum.
  That stage is deliberately separate and not run by default: the
  structural database costs hours of GPU, and it is only worth building
  if the sequence arm finds restriction to confirm.

NO HEADER REWRITING, NO TAXONOMY FILES

  Each downloaded proteome has its clade rank prefixed onto every
  header, so a hit's rank is read from the hit itself. That removes the
  accession-to-taxid mapping step, which is the part of this analysis
  that usually goes quietly wrong.
"""

from __future__ import annotations

import argparse
import json
import math
import gzip
import random
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod                              # noqa: E402

# rank, clade label, [(taxon id, species)]
PANEL: list[tuple[int, str, list[tuple[int, str]]]] = [
    (1, "Decapoda",            [(6689, "Penaeus vannamei")]),
    (2, "Crustacea",           [(6669, "Daphnia pulex")]),
    (3, "Arthropoda",          [(7227, "Drosophila melanogaster"),
                                (7070, "Tribolium castaneum")]),
    (4, "Ecdysozoa",           [(6239, "Caenorhabditis elegans")]),
    (5, "Protostomia",         [(29159, "Crassostrea gigas")]),
    (6, "Bilateria",           [(9606, "Homo sapiens"),
                                (7955, "Danio rerio")]),
    (7, "Metazoa",             [(45351, "Nematostella vectensis"),
                                (400682, "Amphimedon queenslandica")]),
    (8, "Eukaryota",           [(559292, "Saccharomyces cerevisiae"),
                                (3702, "Arabidopsis thaliana")]),
    (9, "Cellular organisms",  [(83333, "Escherichia coli K-12")]),
]
RANK_NAME = {r: n for r, n, _ in PANEL}
RANK_NAME[0] = "P. monodon only"
RESTRICTED_MAX = 3          # Arthropoda or narrower
STREAM = "https://rest.uniprot.org/uniprotkb/stream"
MIN_SEQS = 200              # a proteome smaller than this is a failed fetch


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def human(n: float) -> str:
    for u in ("B", "KB", "MB", "GB"):
        if n < 1024 or u == "GB":
            return f"{n:,.1f} {u}" if u != "B" else f"{int(n)} B"
        n /= 1024
    return f"{n:.1f} GB"


LOCRE = re.compile(r"\bLOC(\d+)\b")


def refseq_by_gene(path: Path) -> dict[str, str]:
    """{LOC id: longest protein} from a RefSeq protein FASTA.

    The file is keyed by protein accession (`XP_037789123.1`) while the
    targets and the DE table are keyed by gene (`gene-LOC119568046`).
    The link is the LOC number in the header DESCRIPTION, not in the
    accession, so matching on the first token finds nothing — which is
    exactly what happened: 0 of 1567 control genes resolved.

    Genes with several isoforms keep the longest, so each gene is one
    sequence and the search is not weighted by isoform count.
    """
    best: dict[str, str] = {}
    gene, buf = None, []

    def flush():
        if gene and buf:
            seq = "".join(buf)
            if len(seq) > len(best.get(gene, "")):
                best[gene] = seq

    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                m = LOCRE.search(line)
                gene, buf = (f"LOC{m.group(1)}" if m else None), []
            else:
                buf.append(line.strip())
    flush()
    return best


def gene_map_from_gff(gff: Path) -> dict[str, str]:
    """{protein accession: gene id} from an NCBI GFF.

    Needed because the protein FASTA names the gene only in the
    description, and only for LOC-named genes: 32,816 proteins but just
    8,669 headers mentioning LOC. Keying off the FASTA therefore built a
    background from the LOC-named subset, which is precisely the
    uncharacterised fraction — so the "proteome background" was an
    uncharacterised background, and the comparison was rigged without
    anyone choosing to rig it.

    CDS features carry protein_id and gene together, which closes it.
    """
    out: dict[str, str] = {}
    opener = gzip.open if gff.suffix == ".gz" else open
    with opener(gff, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            attrs = {}
            for kv in f[8].split(";"):
                if "=" in kv:
                    k, _, v = kv.partition("=")
                    attrs[k] = v
            pid, gene = attrs.get("protein_id"), attrs.get("gene")
            if pid and gene:
                out[pid] = f"gene-{gene}"
    return out


def proteins_by_gene(faa: Path, gff: Path | None) -> tuple[dict, str]:
    """{gene id: longest protein}. Uses the GFF when available."""
    if gff and gff.exists():
        acc2gene = gene_map_from_gff(gff)
        if acc2gene:
            best: dict[str, str] = {}
            gene, buf = None, []

            def flush():
                if gene and buf:
                    sq = "".join(buf)
                    if len(sq) > len(best.get(gene, "")):
                        best[gene] = sq

            with faa.open() as fh:
                for line in fh:
                    if line.startswith(">"):
                        flush()
                        acc = line[1:].split()[0]
                        gene, buf = acc2gene.get(acc), []
                    else:
                        buf.append(line.strip())
            flush()
            return best, f"GFF ({len(acc2gene):,} CDS mapped)"
    loc = refseq_by_gene(faa)
    return ({f"gene-{k}": v for k, v in loc.items()},
            "FASTA headers only — LOC-named genes, a biased subset")


def as_loc(x: str) -> str | None:
    m = LOCRE.search(x)
    return f"LOC{m.group(1)}" if m else None


def query_stamp(fa: Path) -> str:
    import hashlib
    return hashlib.md5(fa.read_bytes()).hexdigest()


def cached_ok(m8: Path, fa: Path) -> bool:
    """A result file is only reusable if it is about the SAME queries."""
    stamp = m8.with_suffix(".queries")
    return (m8.exists() and stamp.exists()
            and stamp.read_text().strip() == query_stamp(fa))


def gene_keyed(path: Path, probe: int = 200) -> float:
    """Fraction of the first headers that look like RefSeq gene ids.

    The targets and the differential-expression table are keyed
    `gene-LOC119568046`. The TransDecoder output is keyed
    `MSTRG.1001.1.p1` — de novo ORFs from a fragmented assembly, a
    different population of proteins entirely. Comparing one against the
    other measures the assembly as much as the biology, so the
    background has to come from the same annotation as the targets.
    """
    n = hit = 0
    try:
        with path.open() as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                n += 1
                h = line[1:].split()[0]
                if h.startswith("gene-") or "LOC" in h or \
                        h.startswith(("NP_", "XP_")):
                    hit += 1
                if n >= probe:
                    break
    except OSError:
        return 0.0
    return hit / n if n else 0.0


def find_tool(cfg, name: str) -> str | None:
    p = shutil.which(name)
    if p:
        return p
    work = cfg.p("work")
    for env in ("busco", "pmbio", "serprot", "esmprobe"):
        for base in (work / "envs" / env,
                     Path.home() / "miniconda3" / "envs" / env):
            c = base / "bin" / name
            if c.is_file():
                return str(c)
    return None


def read_fasta(path: Path) -> dict[str, str]:
    out, cur, buf = {}, None, []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    out[cur] = "".join(buf)
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if cur:
        out[cur] = "".join(buf)
    return out


def count_seqs(path: Path) -> int:
    n = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def fetch(taxid: int, dest: Path, rank: int) -> tuple[bool, str]:
    """Download one proteome, prefixing the clade rank onto every header."""
    url = f"{STREAM}?query=taxonomy_id:{taxid}&format=fasta"
    tmp = dest.with_suffix(".part")
    r = subprocess.run(["curl", "-fsSL", "--retry", "3", "--max-time",
                        "1800", "-o", str(tmp), url],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return False, f"curl exit {r.returncode}: {r.stderr.strip()[:120]}"
    n = count_seqs(tmp)
    if n < MIN_SEQS:
        return False, (f"only {n} sequences — wrong taxon id, or the "
                       f"query syntax has changed")
    with tmp.open() as fin, dest.open("w") as fout:
        for line in fin:
            fout.write(f">r{rank}|{line[1:]}" if line.startswith(">")
                       else line)
    tmp.unlink()
    return True, f"{n:,} sequences, {human(dest.stat().st_size)}"


def deep_structural(path: Path, max_e: float) -> set[str]:
    """Query identifiers with a structural hit in the far panel.

    The far panel holds only ranks at or beyond Ecdysozoa, so any hit
    here is a relative outside Arthropoda — which is precisely what a
    sequence-based call of "lineage-restricted" denies. Self-hits are
    dropped: the query set is included in no panel file, but foldseek
    reports a query against itself whenever the two databases share a
    sequence, and a self-hit counted here would manufacture a
    structural relative for every protein.
    """
    out: set[str] = set()
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if len(f) < 11 or f[0] == f[1]:
            continue
        try:
            if float(f[10]) <= max_e:
                out.add(f[0])
        except ValueError:
            continue
    return out


def strata_from_m8(path: Path, max_e: float) -> dict[str, int]:
    """query -> highest clade rank hit. The FURTHEST hit decides."""
    best: dict[str, int] = {}
    if not path.exists():
        return best
    with path.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11 or not f[1].startswith("r"):
                continue
            try:
                if float(f[10]) > max_e:
                    continue
                rank = int(f[1].split("|")[0][1:])
            except ValueError:
                continue
            if rank > best.get(f[0], 0):
                best[f[0]] = rank
    return best


def fisher(a, b, c, d):
    def lc(n, k):
        return (math.lgamma(n + 1) - math.lgamma(k + 1)
                - math.lgamma(n - k + 1))
    r1, r2, c1, nn = a + b, c + d, a + c, a + b + c + d

    def pr(x):
        return math.exp(lc(r1, x) + lc(r2, c1 - x) - lc(nn, c1))
    p0 = pr(a) * (1 + 1e-9)
    lo, hi = max(0, c1 - r2), min(r1, c1)
    return min(1.0, sum(pr(x) for x in range(lo, hi + 1) if pr(x) <= p0))


def wilson(k, n):
    if not n:
        return None, None
    z, p = 1.959964, k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, c - h), min(1.0, c + h)


def match_length(dark: dict[str, int], pool: dict[str, int],
                 seed: int, per: int = 2) -> list[str]:
    rng = random.Random(seed)
    dl = sorted(dark.values())
    if not dl:
        return []
    edges = [dl[int(q * (len(dl) - 1) / 10)] for q in range(11)]
    buckets: dict[int, list[str]] = defaultdict(list)
    for p, L in pool.items():
        for b in range(10):
            if edges[b] <= L <= edges[b + 1]:
                buckets[b].append(p)
                break
    want: Counter = Counter()
    for L in dark.values():
        for b in range(10):
            if edges[b] <= L <= edges[b + 1]:
                want[b] += 1
                break
    out: list[str] = []
    for b, k in want.items():
        cand = buckets.get(b, [])
        k *= per
        out += cand if len(cand) <= k else rng.sample(cand, k)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--go", action="store_true")
    ap.add_argument("--stage", default="all",
                    choices=["panel", "seq", "control", "struct",
                             "report", "all"])
    ap.add_argument("--evalue", type=float, default=1e-3)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--background-fasta", default=None,
                    help="all P. monodon proteins, for the matched control")
    ap.add_argument("--control-ids", default=None,
                    help="gene identifiers of UNCHARACTERISED genes that "
                         "are NOT STAT-dependent. Holds the selection "
                         "criterion constant, which the whole-proteome "
                         "background does not.")
    ap.add_argument("--control-fasta", default=None,
                    help="protein FASTA for those identifiers, if they "
                         "are not in the background file")
    ap.add_argument("--max-struct-rank", type=int, default=4,
                    help="structural database covers ranks >= this; the "
                         "question is only whether restriction is broken "
                         "BEYOND Arthropoda (default 4)")
    ap.add_argument("--linclust-id", type=float, default=0.9,
                    help="redundancy reduction before ProstT5 (default "
                         "0.9); 327k human sequences are mostly isoforms")
    ap.add_argument("--seed", type=int, default=None)
    args = ap.parse_args()

    cfg = cfgmod.load()
    res = cfg.p("results").resolve()
    prot = cfg.p("proteins").resolve()
    out = (cfg.p("dbs") if cfg.p("dbs").parent.exists()
           else cfg.p("share") / "03_dbs").resolve() / "lineage"
    seed = args.seed if args.seed is not None else cfg.raw.get("seed", 1)
    diamond = find_tool(cfg, "diamond")

    targets_fa = prot / "annotate" / "targets.faa"
    # The full protein set almost certainly ends in .pep, which earlier
    # globs in this project did not match — which is why the 28,423
    # proteins never appeared in any discovery listing.
    roots = [prot]
    for key in ("genome", "interim", "dbs"):
        try:
            d = cfg.p(key)
            if d.is_dir():
                roots.append(d)
        except Exception:
            pass
    cands = []
    for root in roots:
        cands += [p for p in root.rglob("*")
                  if p.is_file() and p.suffix in
                  (".pep", ".faa", ".fasta", ".fa")]
    cands = sorted(set(cands), key=lambda p: -p.stat().st_size)
    bg_fa = Path(args.background_fasta) if args.background_fasta else None
    if bg_fa is None:
        # `*.transdecoder.pep` is the FILTERED protein set. The larger
        # `*_dir/longest_orfs.pep` is the candidate set before filtering
        # — every ORF in every frame, most of them spurious. Spurious
        # ORFs have homologs nowhere, so a background drawn from them
        # would be maximally "lineage-restricted" for a reason that has
        # nothing to do with biology. Size is the wrong way to choose.
        # Prefer a GENE-KEYED protein set: the targets are RefSeq genes,
        # so a TransDecoder background would compare two different
        # populations and measure the assembly as much as the biology.
        clean = [p for p in cands if "decoys" not in str(p).lower()
                 and "longest_orfs" not in p.name.lower()]
        gk = [p for p in clean if gene_keyed(p) > 0.8]
        if gk:
            bg_fa = max(gk, key=lambda p: p.stat().st_size)
        else:
            ranked = [p for p in clean
                      if p.name.endswith(".transdecoder.pep")]
            bg_fa = ranked[0] if ranked else (
                clean[0] if clean else None)

    print("=" * 78)
    print("PHYLOSTRATIGRAPHY — plan")
    print("=" * 78)
    print(f"  diamond    : {diamond or 'NOT FOUND'}")
    print(f"  targets    : {targets_fa if targets_fa.exists() else 'MISSING'}")
    print(f"  background : {bg_fa or 'NOT FOUND — pass --background-fasta'}")
    if bg_fa:
        gk = gene_keyed(bg_fa)
        print(f"               {count_seqs(bg_fa):,} proteins, "
              f"{100 * gk:.0f}% gene-keyed headers")
        if gk <= 0.8:
            print("               WARNING not gene-keyed. The targets "
                  "are RefSeq genes, so this\n"
                  "               compares two different protein "
                  "populations. Pass\n"
                  "               --background-fasta with the RefSeq "
                  "protein set.")
    print(f"\n  protein files, largest first (gene-keyed fraction):")
    for p in cands[:14]:
        print(f"    {human(p.stat().st_size):>10}  "
              f"{100 * gene_keyed(p):>3.0f}%  {p}")
    print(f"  panel dir  : {out}")
    print(f"\n  {'rank':<6}{'clade':<22}{'species':<40}")
    for rank, clade, spp in PANEL:
        for taxid, name in spp:
            f = out / f"r{rank}_{taxid}.faa"
            mark = f"  [{count_seqs(f):,}]" if f.exists() else ""
            print(f"  {rank:<6}{clade:<22}{name + ' (' + str(taxid) + ')':<40}"
                  f"{mark}")
    print(f"\n  lineage-restricted := furthest hit at rank <= "
          f"{RESTRICTED_MAX} ({RANK_NAME[RESTRICTED_MAX]})")
    print(f"  E-value    : {args.evalue:g}   seed {seed}")
    print("=" * 78)

    if not args.go:
        print("\nPLAN ONLY — nothing written. Re-run with --go --stage panel")
        return 0

    out.mkdir(parents=True, exist_ok=True)

    # ---- panel ----------------------------------------------------------
    if args.stage in ("panel", "all"):
        failed = []
        for rank, clade, spp in PANEL:
            for taxid, name in spp:
                dest = out / f"r{rank}_{taxid}.faa"
                if dest.exists() and count_seqs(dest) >= MIN_SEQS:
                    log(f"have  r{rank} {name}")
                    continue
                log(f"fetch r{rank} {name} (taxid {taxid})")
                ok, msg = fetch(taxid, dest, rank)
                log(f"      {'ok' if ok else 'FAILED'}: {msg}")
                if not ok:
                    failed.append((rank, taxid, name, msg))
        if failed:
            print("\n  these did not download:")
            for rank, taxid, name, msg in failed:
                print(f"    r{rank} {name} ({taxid}): {msg}")
            print("\n  A clade with no proteome makes every protein look "
                  "younger than\n  it is, because the rung it would have "
                  "hit is missing. Fix these\n  before running the "
                  "search.")
            return 1
        log("panel complete")
        if args.stage == "panel":
            return 0

    panel_files = sorted(out.glob("r*_*.faa"))
    if not panel_files:
        print("no panel — run --stage panel", file=sys.stderr)
        return 1

    # ---- search ---------------------------------------------------------
    if args.stage in ("seq", "all"):
        if not diamond:
            print("diamond not found; it lives in the busco environment",
                  file=sys.stderr)
            return 1
        combined = out / "panel.faa"
        with combined.open("w") as fh:
            for p in panel_files:
                fh.write(p.read_text())
        log(f"panel.faa: {count_seqs(combined):,} sequences, "
            f"{human(combined.stat().st_size)}")
        db = out / "panel"
        if not Path(f"{db}.dmnd").exists():
            log("diamond makedb")
            if subprocess.run([diamond, "makedb", "--in", str(combined),
                               "-d", str(db), "--quiet"]).returncode != 0:
                return 1

        queries: list[tuple[str, Path]] = []
        if targets_fa.exists():
            queries.append(("targets", targets_fa))
        if bg_fa and bg_fa.exists():
            gff = next((q for q in bg_fa.parent.glob("*.gff*")), None)
            by_gene, how = proteins_by_gene(bg_fa, gff)
            log(f"background source: {len(by_gene):,} genes via {how}")
            if "biased subset" in how:
                print("  WARNING no GFF — the background is the "
                      "LOC-named subset, which\n"
                      "  is the uncharacterised fraction, not the "
                      "proteome.")
            tgt = set(read_fasta(targets_fa)) if targets_fa.exists() \
                else set()
            hit = len(tgt & set(by_gene))
            log(f"  {hit}/{len(tgt)} target genes resolve in it")
            if hit < 0.5 * len(tgt):
                print("  REFUSING: the targets do not resolve against "
                      "this protein set, so a background drawn from it "
                      "is a different population.")
                return 1
            dark_len = {g: len(by_gene[g]) for g in tgt if g in by_gene}
            log(f"  background pool: {len(by_gene) - len(tgt & set(by_gene)):,} "
                f"genes after removing the targets")
            pool_len = {g: len(v) for g, v in by_gene.items()
                        if g not in tgt}
            picked = match_length(dark_len, pool_len, seed)
            bgsub = out / "background.faa"
            bgsub.write_text("".join(f">{g}\n{by_gene[g]}\n"
                                     for g in picked))
            log(f"length-matched background: {len(picked)} genes")
            queries.append(("background", bgsub))

        for name, q in queries:
            m8 = out / f"{name}_vs_panel.m8"
            if cached_ok(m8, q):
                log(f"{name}: reusing {m8.name} (same queries)")
                continue
            if m8.exists():
                log(f"{name}: {m8.name} is about different queries — "
                    f"recomputing")
            log(f"diamond blastp {name} ({count_seqs(q):,} queries)")
            cmd = [diamond, "blastp", "-q", str(q), "-d", str(db),
                   "-o", str(m8), "--ultra-sensitive", "-e",
                   str(args.evalue), "-p", str(args.threads),
                   "--max-target-seqs", "500", "--quiet"]
            if subprocess.run(cmd).returncode != 0:
                return 1
            m8.with_suffix(".queries").write_text(query_stamp(q) + "\n")

    # ---- control: uncharacterised but NOT STAT-dependent ----------------
    # The whole-proteome background conflates lineage restriction with
    # the selection criterion: these targets were chosen for being
    # uncharacterised, and a gene is uncharacterised largely because
    # homology search failed to name it. Comparing against other
    # uncharacterised genes holds that constant.
    if args.stage in ("control", "all"):
        ctrl_ids = None
        if args.control_ids and Path(args.control_ids).is_file():
            ctrl_ids = [l.strip() for l in
                        Path(args.control_ids).read_text().split()
                        if l.strip()]
            log(f"control: {len(ctrl_ids)} ids from {args.control_ids}")
        elif (res / "de_stat_dependence.csv").exists():
            import csv as _csv
            de = res / "de_stat_dependence.csv"
            tgt_ids = set(read_fasta(targets_fa)) if targets_fa.exists() \
                else set()
            TRUE = {"true", "1", "yes", "t", "y"}
            with de.open() as fh:
                rd = _csv.DictReader(fh)
                gcol = rd.fieldnames[0]
                rows = [r for r in rd]
            unch = [r[gcol] for r in rows
                    if (r.get("uncharacterised") or "").strip().lower()
                    in TRUE
                    and (r.get("biotype") or "").strip() in
                    ("", "protein_coding")
                    and r[gcol] not in tgt_ids]
            ctrl_ids = unch
            log(f"control: {len(ctrl_ids)} uncharacterised genes that "
                f"are not among the {len(tgt_ids)} targets "
                f"(from {de.name})")
            if not ctrl_ids:
                print("    the 'uncharacterised' column matched nothing "
                      "— check its values")
        else:
            de = res / "de_stat_dependence.csv"
            sd = res / "stat_dependent_uncharacterised.txt"
            print("\n  CONTROL SET — cannot build it yet, reporting what "
                  "is available:")
            if de.exists():
                import csv as _csv
                with de.open() as fh:
                    rd = _csv.reader(fh)
                    hdr = next(rd, [])
                    row = next(rd, [])
                print(f"    {de.name} columns : {hdr}")
                print(f"    first row         : {row}")
            else:
                print(f"    {de.name}: not found")
            if sd.exists():
                lines = sd.read_text().split()
                print(f"    {sd.name}: {len(lines)} ids, "
                      f"first {lines[0] if lines else '-'}")
            tg = read_fasta(targets_fa) if targets_fa.exists() else {}
            if tg:
                print(f"    targets.faa first header : "
                      f"{list(tg)[0]}")
            if bg_fa and bg_fa.exists():
                bgf = read_fasta(bg_fa)
                print(f"    background first header  : {list(bgf)[0]}")
            print("    -> pass --control-ids once the namespaces are "
                  "known; skipping")
        if ctrl_ids:
            src = (Path(args.control_fasta) if args.control_fasta
                   else bg_fa)
            gff = next((q for q in src.parent.glob("*.gff*")), None) \
                if src else None
            by_gene, how = proteins_by_gene(src, gff) if src and \
                src.exists() else ({}, "none")
            log(f"control source: {len(by_gene):,} genes via {how}")
            keep = [g for g in ctrl_ids if g in by_gene]
            if len(keep) < 0.5 * len(ctrl_ids):
                print(f"    only {len(keep)}/{len(ctrl_ids)} control "
                      f"genes resolve in {src.name if src else '-'}; "
                      f"refusing")
            else:
                log(f"control: {len(keep)}/{len(ctrl_ids)} genes resolved")
                cfa = out / "control.faa"
                cfa.write_text("".join(f">{g}\n{by_gene[g]}\n"
                                       for g in keep))
                m8 = out / "control_vs_panel.m8"
                if not cached_ok(m8, cfa):
                    log(f"diamond blastp control ({len(keep):,} queries)")
                    subprocess.run([diamond, "blastp", "-q", str(cfa),
                                    "-d", str(out / "panel"), "-o",
                                    str(m8), "--ultra-sensitive", "-e",
                                    str(args.evalue), "-p",
                                    str(args.threads), "--max-target-seqs",
                                    "500", "--quiet"])
                    m8.with_suffix(".queries").write_text(
                        query_stamp(cfa) + "\n")

    # ---- structural arm -------------------------------------------------
    # A protein restricted by SEQUENCE search might simply have diverged
    # past detection — the criticism that sinks most phylostratigraphy,
    # and the exact failure mode the rest of this paper documents.
    # Searching structurally for a relative beyond Arthropoda tests
    # restriction with evidence the sequence search cannot supply.
    if args.stage in ("struct", "all"):
        fsk = find_tool(cfg, "foldseek")
        mms = find_tool(cfg, "mmseqs")
        weights = cfg.p("structures") / "prostt5_weights"
        far = [p for p in panel_files
               if int(p.name[1:p.name.index("_")]) >= args.max_struct_rank]
        print(f"\n  STRUCTURAL ARM")
        print(f"    foldseek {fsk or 'NOT FOUND'}")
        print(f"    mmseqs   {mms or 'NOT FOUND'}")
        print(f"    weights  {weights if weights.exists() else 'MISSING'}")
        print(f"    panel    ranks >= {args.max_struct_rank}: "
              f"{len(far)} files, "
              f"{sum(count_seqs(p) for p in far):,} sequences")
        if not (fsk and mms and weights.exists()):
            print("    cannot run; mmseqs lives in the serprot "
                  "environment")
        else:
            big = out / "panel_far.faa"
            if not big.exists():
                with big.open("w") as fh:
                    for p in far:
                        fh.write(p.read_text())
            rep_fa = out / f"panel_far_rep_{args.linclust_id}.fasta"
            if not rep_fa.exists():
                log(f"mmseqs linclust at {args.linclust_id} "
                    f"(removes isoform redundancy before ProstT5)")
                pref = out / "lc"
                t = out / "tmp_lc"
                t.mkdir(exist_ok=True)
                if subprocess.run([mms, "easy-linclust", str(big),
                                   str(pref), str(t), "--min-seq-id",
                                   str(args.linclust_id), "-c", "0.8",
                                   "--threads", str(args.threads)]
                                  ).returncode != 0:
                    return 1
                shutil.move(f"{pref}_rep_seq.fasta", rep_fa)
            n_rep = count_seqs(rep_fa)
            log(f"representatives: {n_rep:,} "
                f"(from {sum(count_seqs(p) for p in far):,})")
            est = n_rep / 25000 * 14
            log(f"ProstT5 estimate: ~{est:.0f} min on an L40, "
                f"~{est * 2.5:.0f} min on an A4000")
            db = out / "fardb"
            if not (Path(f"{db}_ss.0").exists() or
                    Path(f"{db}_ss").exists()):
                cmd = [fsk, "createdb", str(rep_fa), str(db),
                       "--prostt5-model", str(weights), "--gpu", "1"]
                if subprocess.run(cmd).returncode != 0:
                    return 1
            qdb = out / "tgtdb"
            if not Path(f"{qdb}.dbtype").exists():
                if subprocess.run([fsk, "createdb", str(targets_fa),
                                   str(qdb), "--prostt5-model",
                                   str(weights), "--gpu", "1"]
                                  ).returncode != 0:
                    return 1
            aln, t2 = out / "aln_far", out / "tmp_far"
            t2.mkdir(exist_ok=True)
            if subprocess.run([fsk, "search", str(qdb), str(db), str(aln),
                               str(t2), "-e", str(args.evalue),
                               "--alignment-type", "2"]).returncode != 0:
                return 1
            sm8 = out / "targets_vs_far.m8"
            if subprocess.run([fsk, "convertalis", str(qdb), str(db),
                               str(aln), str(sm8), "--format-output",
                               "query,target,fident,alnlen,mismatch,"
                               "gapopen,qstart,qend,tstart,tend,evalue,"
                               "bits"]).returncode != 0:
                return 1
            deep = deep_structural(sm8, args.evalue)
            log(f"{len(deep)} targets have a STRUCTURAL relative beyond "
                f"Arthropoda")
            log("run --stage report to cross-tabulate this against the "
                "sequence strata and record it")

    # ---- report ---------------------------------------------------------
    # Guarded: without this, `--stage seq` fell through to a report with
    # no inputs, printed "no search output" and returned 1 — which broke
    # the `&&` chain and stopped the next stage from running at all.
    if args.stage not in ("report", "seq", "all"):
        return 0

    report: dict = {"evalue": args.evalue, "seed": seed,
                    "restricted_max_rank": RESTRICTED_MAX,
                    "panel": {f"r{r}": n for r, n, _ in PANEL}}
    sets = {}
    for name in ("targets", "control", "background"):
        m8 = out / f"{name}_vs_panel.m8"
        fa = (targets_fa if name == "targets" else out / f"{name}.faa")
        if not m8.exists() or not fa.exists():
            continue
        if not cached_ok(m8, fa):
            print(f"  SKIPPING {name}: {m8.name} is about different "
                  f"queries than {fa.name}. Re-run --stage seq.")
            continue
        ids = list(read_fasta(fa))
        st = strata_from_m8(m8, args.evalue)
        strata = {i: st.get(i, 0) for i in ids}
        sets[name] = strata

    if not sets:
        print("\n  no search output to report on")
        return 1

    print("\n" + "=" * 78)
    print("PHYLOSTRATA")
    print("=" * 78)
    print(f"  {'rank':<5}{'clade':<22}" +
          "".join(f"{n:>14}" for n in sets))
    for rank in range(0, 10):
        row = [sum(1 for v in s.values() if v == rank) for s in sets.values()]
        if not any(row):
            continue
        cells = "".join(f"{c:>8} {100 * c / len(list(sets.values())[i]):>4.1f}%"
                        for i, c in enumerate(row))
        print(f"  {rank:<5}{RANK_NAME.get(rank, '?'):<22}{cells}")

    for name, s in sets.items():
        k = sum(1 for v in s.values() if v <= RESTRICTED_MAX)
        lo, hi = wilson(k, len(s))
        report[name] = {"n": len(s), "restricted": k,
                        "fraction": k / len(s), "ci": [lo, hi],
                        "distribution": {str(r): sum(1 for v in s.values()
                                                     if v == r)
                                         for r in range(10)}}
        print(f"\n  {name}: {k}/{len(s)} = {100 * k / len(s):.1f}% "
              f"[{100 * lo:.1f}, {100 * hi:.1f}] lineage-restricted "
              f"(<= {RANK_NAME[RESTRICTED_MAX]})")

    for other, why in (("background",
                        "whole proteome — conflates restriction with "
                        "the selection criterion"),
                       ("control",
                        "uncharacterised but not STAT-dependent — "
                        "selection criterion held constant")):
        if "targets" not in report or other not in report:
            continue
        a = report["targets"]["restricted"]
        b = report["targets"]["n"] - a
        c = report[other]["restricted"]
        d = report[other]["n"] - c
        p = fisher(a, b, c, d)
        orr = ((a + .5) * (d + .5)) / ((b + .5) * (c + .5))
        report[f"vs_{other}"] = {"odds_ratio": orr, "fisher_p": p,
                                 "interpretation": why}
        print(f"\n  targets vs {other}: OR {orr:.2f}, p = {p:.3g}")
        print(f"    ({why})")
    if "vs_control" not in report:
        print("\n  No uncharacterised non-STAT control was scored. The "
              "comparison\n  against the whole proteome is partly "
              "circular: these proteins were\n  selected for being "
              "uncharacterised, which is itself a record of\n  homology "
              "search having failed. Report the DESCRIPTIVE figure (how "
              "many\n  of the targets are restricted) rather than the "
              "enrichment.")

    # ---- the detection-bias check --------------------------------------
    # Every stratum above rests on sequence search, which this paper
    # spends its length showing to be unreliable on exactly these
    # proteins. The structural arm asks the one question that can be
    # asked without it: of the targets called lineage-restricted because
    # no sequence homolog was found outside Arthropoda, how many have a
    # STRUCTURAL relative out there after all? Those are over-calls, and
    # the rate belongs beside the headline figure rather than in a log.
    sm8 = out / "targets_vs_far.m8"
    if "targets" in sets and sm8.exists():
        deep = deep_structural(sm8, args.evalue)
        ids = list(sets["targets"].keys())
        resolved = sum(1 for i in ids if i in deep)
        # A namespace mismatch would make every target look structurally
        # unreachable and the over-call rate read as a reassuring zero.
        # Refuse that reading rather than print it.
        if deep and resolved == 0:
            print(f"\n  STRUCTURAL CHECK SKIPPED: none of the "
                  f"{len(deep)} query identifiers in {sm8.name} match "
                  f"the {len(ids)} target identifiers.")
            print(f"    m8 query  : {sorted(deep)[0]}")
            print(f"    target id : {ids[0]}")
            report["structural_check"] = {
                "status": "identifier namespaces differ; not reported"}
        else:
            rest = [i for i, r in sets["targets"].items()
                    if r <= RESTRICTED_MAX]
            far = [i for i, r in sets["targets"].items()
                   if r > RESTRICTED_MAX]
            a = sum(1 for i in rest if i in deep)
            c = sum(1 for i in far if i in deep)
            lo, hi = wilson(a, len(rest)) if rest else (0.0, 0.0)
            print("\n" + "=" * 78)
            print("STRUCTURAL CHECK ON THE SEQUENCE STRATA")
            print("=" * 78)
            print(f"  far panel: ranks >= {args.max_struct_rank}, so any "
                  f"hit is a relative outside Arthropoda")
            print(f"  {len(deep)} of {len(ids)} targets have a structural "
                  f"relative there")
            print(f"\n  called lineage-restricted by sequence: "
                  f"{len(rest)}")
            if rest:
                print(f"    of these, {a} ({100 * a / len(rest):.1f}% "
                      f"[{100 * lo:.1f}, {100 * hi:.1f}]) have a "
                      f"structural relative beyond Arthropoda")
                print(f"    -> sequence phylostratigraphy over-calls "
                      f"restriction for at least these {a}")
            print(f"\n  not called restricted: {len(far)}")
            if far:
                print(f"    of these, {c} ({100 * c / len(far):.1f}%) "
                      f"have one, as expected since they already had a "
                      f"deep sequence hit")
            adj = (len(rest) - a) / len(ids) if ids else 0.0
            print(f"\n  restricted fraction, sequence only : "
                  f"{100 * len(rest) / len(ids):.1f}%")
            print(f"  restricted fraction, after structure: "
                  f"{100 * adj:.1f}%")
            report["structural_check"] = {
                "status": "reported",
                "max_struct_rank": args.max_struct_rank,
                "n_targets": len(ids),
                "with_deep_structural_hit": len(deep),
                "seq_restricted": len(rest),
                "seq_restricted_with_deep_structure": a,
                "overcall_fraction": (a / len(rest)) if rest else None,
                "overcall_ci": [lo, hi] if rest else None,
                "seq_unrestricted": len(far),
                "seq_unrestricted_with_deep_structure": c,
                "restricted_fraction_sequence_only":
                    len(rest) / len(ids) if ids else None,
                "restricted_fraction_after_structure": adj,
                "note": ("An over-call is a protein with no detectable "
                         "sequence homolog outside Arthropoda but a "
                         "structural relative there. The adjusted "
                         "fraction is a lower bound on true restriction: "
                         "structure search has its own detection limit, "
                         "so proteins counted as restricted after it may "
                         "still have undetected relatives."),
            }

    report["caveat"] = (
        "Strata are inferred from detectable sequence homology. Absence "
        "of detection is not absence of homolog, and this study shows "
        "sequence search failing on exactly these proteins. Confirm with "
        "the structural arm before claiming lineage restriction.")
    if (report.get("structural_check") or {}).get("status") == "reported":
        sc = report["structural_check"]
        if sc.get("overcall_fraction") is not None:
            report["caveat"] += (
                f" The structural arm quantifies this: "
                f"{sc['seq_restricted_with_deep_structure']} of the "
                f"{sc['seq_restricted']} targets called restricted have a "
                f"structural relative beyond Arthropoda "
                f"({100 * sc['overcall_fraction']:.1f}%), so the "
                f"sequence-only figure of "
                f"{100 * sc['restricted_fraction_sequence_only']:.1f}% "
                f"falls to at most "
                f"{100 * sc['restricted_fraction_after_structure']:.1f}%.")
    p = res / "lineage.json"
    p.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\n  wrote {p}")
    print("\n  NEXT: these strata rest on sequence search alone, which "
          "this paper\n  shows is unreliable here. The structural "
          "confirmation is the\n  separate stage, worth building only if "
          "the enrichment above is real.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
