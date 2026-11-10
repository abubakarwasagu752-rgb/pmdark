#!/usr/bin/env python3
"""Evidence-based transcript filtering that does not delete short proteins.

    python scripts/23c_filter_transcripts.py --min-libraries 3

THE PROBLEM WITH -m 500

  The merge sweep showed that StringTie's minimum-length flag does nearly
  all the noise removal:

      -m 200  -> 100,781 loci, 27.2% multi-exon
      -m 300  ->  70,376 loci, 32.4%
      -m 500  ->  29,625 loci, 52.9%

  Setting -m 500 lands beautifully on the reference gene count (1.01x) but
  achieves it by deleting every transcript under 500 bp — roughly 40,000 of
  them between m=300 and m=500.

  That is unacceptable for THIS project. Penaeid antimicrobial peptides
  (penaeidins, crustins, ALFs) mature at 5-12 kDa. A 100 aa precursor is
  300 bp of coding sequence; with short UTRs the transcript may never reach
  500 bp. Aim 3 is specifically hunting lineage-restricted small secreted
  effectors, so a blanket length cutoff would remove the target class and
  then report its absence as a finding.

  Length is a proxy for spuriousness, and a bad one. Independent support is
  the real signal.

THE RULE

  Merge permissively (-m 200), then keep a transcript if ANY of:

    1. length >= --long           (500 bp; long transcripts are kept as-is)
    2. multi-exon                 (spliced structure is strong evidence)
    3. supported in >= --min-libraries independent libraries

  Rule 3 is what rescues short single-exon AMP candidates: a 320 bp
  transcript seen in 3 of 12 libraries is real; one seen in a single
  low-mapping library is not. Support is counted by exon-level overlap
  against the per-library GTFs, which were produced independently.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod  # noqa: E402

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "mg", Path(__file__).resolve().parent / "23_merge_assembly.py")
_mg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mg)
gtf_stats = _mg.gtf_stats


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def parse_attr(attrs: str, key: str) -> str | None:
    tag = f'{key} "'
    i = attrs.find(tag)
    if i < 0:
        return None
    j = attrs.find('"', i + len(tag))
    return attrs[i + len(tag):j]


def read_transcripts(gtf: Path) -> dict:
    """transcript_id -> {chrom, strand, gene, exons:[(s,e)], length, n_exons}"""
    tx: dict[str, dict] = {}
    with gtf.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "exon":
                continue
            tid = parse_attr(p[8], "transcript_id")
            if not tid:
                continue
            rec = tx.setdefault(tid, {"chrom": p[0], "strand": p[6],
                                      "gene": parse_attr(p[8], "gene_id"),
                                      "exons": []})
            rec["exons"].append((int(p[3]), int(p[4])))
    for rec in tx.values():
        rec["exons"].sort()
        rec["length"] = sum(e - s + 1 for s, e in rec["exons"])
        rec["n_exons"] = len(rec["exons"])
    return tx


def build_index(gtf: Path, binsize: int = 10000) -> dict:
    """Binned exon index per (chrom, strand) for fast overlap queries."""
    idx: dict[tuple, dict[int, list]] = defaultdict(lambda: defaultdict(list))
    with gtf.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "exon":
                continue
            s, e = int(p[3]), int(p[4])
            key = (p[0], p[6])
            for b in range(s // binsize, e // binsize + 1):
                idx[key][b].append((s, e))
    return idx


def has_overlap(idx: dict, chrom: str, strand: str, exons: list,
                binsize: int = 10000, min_frac: float = 0.5) -> bool:
    """Does at least min_frac of this transcript's bases overlap the index?"""
    key = (chrom, strand)
    if key not in idx:
        return False
    bins = idx[key]
    covered = total = 0
    for s, e in exons:
        total += e - s + 1
        seen = set()
        for b in range(s // binsize, e // binsize + 1):
            for (qs, qe) in bins.get(b, ()):
                lo, hi = max(s, qs), min(e, qe)
                if lo <= hi:
                    seen.update(range(lo, hi + 1))
        covered += len(seen)
    return total > 0 and covered / total >= min_frac


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--merged", type=Path, default=None)
    ap.add_argument("--long", type=int, default=500,
                    help="length at or above which a transcript is kept outright")
    ap.add_argument("--min-libraries", type=int, default=3,
                    help="independent libraries needed to rescue a short "
                         "single-exon transcript")
    ap.add_argument("--min-frac", type=float, default=0.5,
                    help="fraction of a transcript that must overlap a "
                         "library's exons to count as support")
    ap.add_argument("--remerge", action="store_true",
                    help="rebuild the permissive -m 200 merge first")
    args = ap.parse_args()

    cfg = cfgmod.load()
    out = cfg.p("proteins")
    merged = args.merged or (out / "merged_permissive.gtf")

    if args.remerge or not merged.exists():
        import shutil, subprocess
        if not shutil.which("stringtie"):
            raise SystemExit("stringtie not on PATH")
        lst = out / "gtf_list.txt"
        log("building permissive merge (-T 1.0 -m 200 -f 0.01)")
        subprocess.run(["stringtie", "--merge", "-p", "16",
                        "-T", "1.0", "-m", "200", "-f", "0.01",
                        "-o", str(merged), str(lst)], check=False)

    log(f"reading {merged.name}")
    tx = read_transcripts(merged)
    log(f"  {len(tx):,} transcripts")

    lib_gtfs = sorted(p for p in cfg.p("align").glob("*/[SE]RR*.gtf"))
    log(f"indexing {len(lib_gtfs)} per-library GTFs for support counting")
    indexes = []
    for g in lib_gtfs:
        indexes.append((g.parent.name, build_index(g)))
        log(f"  {g.parent.name}")

    log("counting independent support per transcript")
    support: dict[str, int] = {}
    for n, (tid, rec) in enumerate(tx.items(), 1):
        c = 0
        for _, idx in indexes:
            if has_overlap(idx, rec["chrom"], rec["strand"], rec["exons"],
                           min_frac=args.min_frac):
                c += 1
        support[tid] = c
        if n % 20000 == 0:
            log(f"  {n:,}/{len(tx):,}")

    keep: set[str] = set()
    reasons = {"long": 0, "multi_exon": 0, "replicated": 0}
    for tid, rec in tx.items():
        if rec["length"] >= args.long:
            keep.add(tid); reasons["long"] += 1
        elif rec["n_exons"] > 1:
            keep.add(tid); reasons["multi_exon"] += 1
        elif support[tid] >= args.min_libraries:
            keep.add(tid); reasons["replicated"] += 1

    # The number that matters: short single-exon transcripts rescued by
    # replication. These are the AMP-candidate space a -m 500 cutoff would
    # have silently deleted.
    rescued = reasons["replicated"]

    print("\n" + "=" * 74)
    print("EVIDENCE-BASED FILTER")
    print("=" * 74)
    print(f"  input transcripts        {len(tx):,}")
    print(f"  kept, length >= {args.long}      {reasons['long']:,}")
    print(f"  kept, multi-exon         {reasons['multi_exon']:,}")
    print(f"  kept, >= {args.min_libraries} libraries     {rescued:,}   "
          "<-- short single-exon RESCUED")
    print(f"  total kept               {len(keep):,} "
          f"({100*len(keep)/max(1,len(tx)):.1f}%)")
    print(f"  discarded                {len(tx)-len(keep):,}")

    short_kept = sum(1 for t in keep if tx[t]["length"] < args.long)
    print(f"\n  transcripts under {args.long} bp retained: {short_kept:,}")
    print("  (a blanket -m 500 would have removed all of these, including")
    print("   any penaeidin / crustin / ALF precursors)")

    filt = out / "merged_filtered.gtf"
    written = 0
    with merged.open() as fh, filt.open("w") as o:
        for line in fh:
            if line.startswith("#"):
                o.write(line); continue
            p = line.split("\t")
            if len(p) < 9:
                continue
            tid = parse_attr(p[8], "transcript_id")
            if tid in keep:
                o.write(line); written += 1
    log(f"wrote {filt.name} ({written:,} feature lines)")

    s = gtf_stats(filt)
    print("\n" + "=" * 74)
    print("FILTERED ASSEMBLY")
    print("=" * 74)
    for k, v in s.items():
        print(f"  {k:<20} {v:,}" if isinstance(v, int) else f"  {k:<20} {v}")

    gff = next(iter(cfg.p("genome").glob("*_genomic.gff")), None)
    if gff:
        ref_genes = sum(1 for line in gff.open()
                        if not line.startswith("#")
                        and len(line.split("\t")) > 2
                        and line.split("\t")[2] == "gene")
        print(f"\n  reference genes {ref_genes:,}  -> "
              f"{s['n_loci']/max(1,ref_genes):.2f}x")

    rp = cfg.p("results") / "transcript_filter.json"
    rp.write_text(json.dumps({
        "params": {"long": args.long, "min_libraries": args.min_libraries,
                   "min_frac": args.min_frac},
        "input_transcripts": len(tx),
        "kept_by_reason": reasons,
        "total_kept": len(keep),
        "short_transcripts_retained": short_kept,
        "filtered_stats": s,
    }, indent=2))
    print(f"\nwrote {rp}")
    print("\nNext: regenerate the FASTA, then TransDecoder at 70 aa")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
