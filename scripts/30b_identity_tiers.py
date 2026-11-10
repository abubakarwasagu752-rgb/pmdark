#!/usr/bin/env python3
"""Sequence-identity tiers — the benchmark that actually answers our question.

    python scripts/30b_identity_tiers.py --threads 16

WHY THIS REPLACES THE TAXONOMY TIERS AS THE HEADLINE

  Taxonomy is a proxy, and on this dataset a weak one: Arthropoda is only
  1,087 of 19,304 references, so T3_no_phylum removes 5.6% of the pool and
  barely changes what the model can see. A query whose nearest neighbours
  were non-arthropod anyway scores identically at T0 and T3, and the
  contrast washes out.

  What we actually want to know is not "is there a relative in the same
  phylum" but "is there a SIMILAR SEQUENCE". That is the condition the
  P. monodon dark proteome is in: proteins whose best homolog, anywhere,
  is below the detection threshold of BLAST-style search.

  So we align every query against every reference with DIAMOND and build
  exclusion lists at decreasing identity ceilings:

    I90   nothing above 90% identity may be used   (near-duplicates gone)
    I50   nothing above 50%                        (clear homologs gone)
    I30   nothing above 30%                        (twilight zone)
    I20   nothing above 20%                        (no detectable homology)
    Inone no DIAMOND hit at all may be used        (the true orphan case)

  I30 and below is where homology transfer fails and where a PLM has to
  earn its keep. This is also how the CAFA community benchmarks, so a
  reviewer will recognise it immediately.

  Both tier families are reported. Taxonomy answers "does phylogenetic
  distance matter", identity answers "does sequence divergence matter".
  The second is the one the abstract quotes.
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod  # noqa: E402

# Ceilings, strictest last. A tier excludes every reference whose identity
# to the query is AT OR ABOVE the ceiling.
IDENTITY_TIERS = [
    ("I90", 90.0, "exclude near-duplicates (>=90% id)"),
    ("I50", 50.0, "exclude clear homologs (>=50% id)"),
    ("I30", 30.0, "exclude the twilight zone (>=30% id)"),
    ("I20", 20.0, "exclude anything detectably similar (>=20% id)"),
    ("Inone", 0.0, "exclude every DIAMOND hit — the true orphan case"),
]


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def find_diamond(explicit: str | None = None) -> str | None:
    """Locate the diamond binary without requiring a PATH change.

    Prepending another conda env's bin/ to PATH also shadows `python`,
    which breaks this script's own imports. So we search known env
    locations directly instead.
    """
    import os
    for cand in (explicit, os.environ.get("DIAMOND"), shutil.which("diamond")):
        if cand and Path(cand).exists():
            return cand
    home = Path.home()
    for env in ("busco", "pmbio"):
        p = home / "pmdark_work" / "envs" / env / "bin" / "diamond"
        if p.exists():
            return str(p)
    return None


def run_diamond(query: Path, reference: Path, out: Path, threads: int,
                evalue: float, max_targets: int,
                diamond: str | None = None) -> bool:
    dmnd = find_diamond(diamond)
    if not dmnd:
        log("ERROR: diamond not found.")
        log("  expected at ~/pmdark_work/envs/busco/bin/diamond")
        log("  override with --diamond /path/to/diamond or DIAMOND=/path/...")
        return False
    log(f"diamond: {dmnd}")

    db = out.parent / "reference.dmnd"
    if not db.exists():
        log(f"building DIAMOND db from {reference.name}")
        rc = subprocess.run([dmnd, "makedb", "--in", str(reference),
                             "-d", str(db.with_suffix("")), "--threads", str(threads)
                             ]).returncode
        if rc != 0:
            return False

    log(f"aligning {query.name} vs reference "
        f"(evalue {evalue}, max {max_targets} targets)")
    # --ultra-sensitive matters here: we specifically need the weak hits in
    # the 20-30% range, which default mode misses.
    rc = subprocess.run([
        dmnd, "blastp",
        "-q", str(query), "-d", str(db.with_suffix("")),
        "-o", str(out),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
        "--evalue", str(evalue),
        "--max-target-seqs", str(max_targets),
        "--ultra-sensitive",
        "--threads", str(threads),
    ]).returncode
    return rc == 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--evalue", type=float, default=10.0,
                    help="permissive on purpose — we need weak hits to EXCLUDE them")
    ap.add_argument("--max-targets", type=int, default=2000)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--diamond", default=None,
                    help="explicit path to the diamond binary")
    args = ap.parse_args()

    cfg = cfgmod.load()
    bench = cfg.p("proteins", "benchmark")
    query = bench / "query.faa"
    reference = bench / "reference.faa"
    for p in (query, reference):
        if not p.exists():
            raise SystemExit(f"missing {p} — run scripts/30_build_benchmark.py --build")

    hits_path = bench / "diamond_hits.tsv"
    if args.force or not hits_path.exists() or hits_path.stat().st_size == 0:
        if not run_diamond(query, reference, hits_path,
                           args.threads, args.evalue, args.max_targets,
                           args.diamond):
            return 1
    else:
        log(f"reusing {hits_path.name} ({hits_path.stat().st_size/1e6:.1f} MB)")

    # ---- best identity per (query, reference) pair ---------------------
    log("parsing hits")
    best: dict[str, dict[str, float]] = defaultdict(dict)
    n_lines = 0
    with hits_path.open() as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            q, s, pid = p[0], p[1], float(p[2])
            if q == s:
                continue                     # self-hit, always excluded anyway
            if pid > best[q].get(s, -1.0):
                best[q][s] = pid
            n_lines += 1
    log(f"  {n_lines:,} alignments, {len(best):,} queries with at least one hit")

    # ---- exclusion lists per tier --------------------------------------
    exclusions: dict[str, dict[str, list[str]]] = {}
    for name, ceiling, desc in IDENTITY_TIERS:
        per_query: dict[str, list[str]] = {}
        sizes = []
        for q, hits in best.items():
            if name == "Inone":
                drop = list(hits.keys())
            else:
                drop = [s for s, pid in hits.items() if pid >= ceiling]
            per_query[q] = drop
            sizes.append(len(drop))
        exclusions[name] = per_query
        sizes.sort()
        med = sizes[len(sizes) // 2] if sizes else 0
        mx = sizes[-1] if sizes else 0
        log(f"  {name:<6} {desc:<46} median {med:>5} max {mx:>6} excluded")

    out = bench / "identity_tiers.json.gz"
    with gzip.open(out, "wt") as fh:
        json.dump({
            "tiers": [{"name": n, "ceiling": c, "desc": d}
                      for n, c, d in IDENTITY_TIERS],
            "evalue": args.evalue,
            "exclusions": exclusions,
        }, fh)
    log(f"wrote {out.name} ({out.stat().st_size/1e6:.1f} MB)")

    # ---- how hard is this benchmark, really ----------------------------
    # If most queries already have no homolog above 30%, the dataset is
    # intrinsically hard and the identity tiers will show a clean gradient.
    # If most have a >90% match, the easy tiers are nearly meaningless.
    all_q = [l[1:].split()[0] for l in query.read_text().splitlines()
             if l.startswith(">")]
    print()
    print("  best-hit identity distribution across all queries:")
    buckets = {"no hit": 0, "<20": 0, "20-30": 0, "30-50": 0, "50-90": 0, ">=90": 0}
    for q in all_q:
        hits = best.get(q)
        if not hits:
            buckets["no hit"] += 1
            continue
        top = max(hits.values())
        if top < 20:    buckets["<20"] += 1
        elif top < 30:  buckets["20-30"] += 1
        elif top < 50:  buckets["30-50"] += 1
        elif top < 90:  buckets["50-90"] += 1
        else:           buckets[">=90"] += 1
    for k, v in buckets.items():
        bar = "#" * int(50 * v / max(1, len(all_q)))
        print(f"    {k:<8} {v:>5} ({100*v/max(1,len(all_q)):4.1f}%) {bar}")

    print("\nNext:")
    print("  python scripts/31_calibrate.py --model esm2 \\")
    print("      --tiers T0_full,I90,I50,I30,I20,Inone")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
