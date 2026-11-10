#!/usr/bin/env python3
"""Merge the 12 per-library GTFs and test the assembly against its gates.

    python scripts/23_merge_assembly.py --threads 16

WHAT THIS DECIDES

  docs/ASSEMBLY_STRATEGY.md commits to beating the upstream de novo
  assembly on specific, pre-stated numbers. This is where that claim is
  settled:

    gate              upstream        target          why
    ----------------  --------------  --------------  ----------------------
    transcript N50    417 bp          > 1,200 bp      fragmentation resolved
    loci              586,928         25,000-45,000   collapse to real genes
    mean length       409 bp          > 900 bp        full-length transcripts

  586,928 "unigenes" is 15-20x the plausible gene count for a penaeid.
  If the merged genome-guided assembly lands in the 25-45k range, that
  number alone demonstrates the de novo assembly was shattered — which is
  the methodological finding the paper leads with.

  ALL 12 LIBRARIES ARE MERGED, including the failed SRR28572179 (830k
  reads). Shallow data still contributes evidence for transcript structure;
  it is only excluded from differential expression, where its depth would
  distort dispersion estimates.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod  # noqa: E402

GATES = {
    "transcript_n50":   {"upstream": 417,     "target": 1200,  "cmp": "gt"},
    "n_loci":           {"upstream": 586928,  "target": (25000, 45000), "cmp": "range"},
    "mean_transcript":  {"upstream": 409,     "target": 900,   "cmp": "gt"},
}


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def run(cmd: list[str]) -> int:
    log("$ " + " ".join(str(c) for c in cmd)[:220])
    return subprocess.run([str(c) for c in cmd]).returncode


def gtf_stats(gtf: Path) -> dict:
    """Transcript lengths (sum of exons) and locus counts, straight from GTF."""
    tx_len: dict[str, int] = defaultdict(int)
    tx_gene: dict[str, str] = {}
    tx_exons: dict[str, int] = defaultdict(int)

    with gtf.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "exon":
                continue
            attrs = p[8]
            tid = gid = None
            for field in attrs.split(";"):
                field = field.strip()
                if field.startswith('transcript_id "'):
                    tid = field[15:-1]
                elif field.startswith('gene_id "'):
                    gid = field[9:-1]
            if not tid:
                continue
            tx_len[tid] += int(p[4]) - int(p[3]) + 1
            tx_exons[tid] += 1
            if gid:
                tx_gene[tid] = gid

    lens = sorted(tx_len.values(), reverse=True)
    total = sum(lens)
    n50 = 0
    acc = 0
    for L in lens:
        acc += L
        if acc >= total / 2:
            n50 = L
            break

    multi = sum(1 for t, n in tx_exons.items() if n > 1)
    return {
        "n_transcripts": len(lens),
        "n_loci": len(set(tx_gene.values())),
        "total_bp": total,
        "transcript_n50": n50,
        "mean_transcript": round(total / max(1, len(lens)), 1),
        "median_transcript": lens[len(lens) // 2] if lens else 0,
        "max_transcript": lens[0] if lens else 0,
        "min_transcript": lens[-1] if lens else 0,
        "multi_exon": multi,
        "multi_exon_pct": round(100 * multi / max(1, len(lens)), 1),
    }


def check_gates(stats: dict) -> list[tuple[str, bool, str]]:
    out = []
    for key, g in GATES.items():
        val = stats.get(key, 0)
        if g["cmp"] == "gt":
            ok = val > g["target"]
            msg = f"{val:,} vs target >{g['target']:,} (upstream {g['upstream']:,})"
        else:
            lo, hi = g["target"]
            ok = lo <= val <= hi
            msg = f"{val:,} vs target {lo:,}-{hi:,} (upstream {g['upstream']:,})"
        out.append((key, ok, msg))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--min-tpm", type=float, default=1.0,
                    help="stringtie --merge -T; drop transcripts below this TPM")
    ap.add_argument("--min-len", type=int, default=200,
                    help="stringtie --merge -m")
    ap.add_argument("--with-reference", action="store_true",
                    help="supply the NCBI GFF as -G. OFF by default: novel "
                         "unannotated loci are the point of this project, and "
                         "a reference guide biases toward known genes.")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    cfg = cfgmod.load()
    align = cfg.p("align")
    out = cfg.p("proteins")
    out.mkdir(parents=True, exist_ok=True)

    gtfs = sorted(align.glob("*/[SE]RR*.gtf"))
    if not gtfs:
        gtfs = sorted(p for p in align.glob("*/*.gtf") if p.parent.name != "_unmapped")
    if len(gtfs) < 2:
        raise SystemExit(f"found only {len(gtfs)} GTFs under {align}")

    log(f"{len(gtfs)} per-library GTFs")
    for g in gtfs:
        n = sum(1 for line in g.open() if "\ttranscript\t" in line)
        log(f"  {g.parent.name:<14} {n:>7,} transcripts")

    if not shutil.which("stringtie"):
        raise SystemExit("stringtie not on PATH — conda activate "
                         "~/pmdark_work/envs/pmbio")

    lst = out / "gtf_list.txt"
    lst.write_text("\n".join(str(g) for g in gtfs) + "\n")

    merged = out / "merged.gtf"
    if merged.exists() and not args.force:
        log(f"reusing {merged.name} (--force to rebuild)")
    else:
        cmd = ["stringtie", "--merge", "-p", args.threads,
               "-m", args.min_len, "-T", args.min_tpm,
               "-o", merged, lst]
        if args.with_reference:
            gff = next(iter(cfg.p("genome").glob("*_genomic.gff")), None)
            if gff:
                cmd[2:2] = ["-G", str(gff)]
                log(f"guiding with {gff.name}")
        if run(cmd) != 0:
            return 1

    log("computing assembly statistics")
    stats = gtf_stats(merged)

    print("\n" + "=" * 74)
    print("MERGED ASSEMBLY")
    print("=" * 74)
    for k, v in stats.items():
        print(f"  {k:<20} {v:,}" if isinstance(v, int) else f"  {k:<20} {v}")

    print("\n" + "=" * 74)
    print("VALIDATION GATES  (docs/ASSEMBLY_STRATEGY.md)")
    print("=" * 74)
    results = check_gates(stats)
    for key, ok, msg in results:
        print(f"  [{'PASS' if ok else 'FAIL'}] {key:<18} {msg}")
    print("-" * 74)

    n_pass = sum(1 for _, ok, _ in results if ok)
    if n_pass == len(results):
        print(f"  All {len(results)} gates passed. The genome-guided assembly")
        print("  beats the upstream de novo assembly on every pre-stated")
        print("  metric. This comparison is Table 1 of the paper.")
    else:
        print(f"  {n_pass}/{len(results)} gates passed. Do NOT proceed to ORF")
        print("  calling until the failures are understood — a bad assembly")
        print("  produces a bad dark proteome, which is the one thing this")
        print("  project cannot afford.")

    # ---- transcript FASTA for TransDecoder ------------------------------
    fa = out / "merged_transcripts.fa"
    if shutil.which("gffread"):
        genome = next(iter(cfg.p("genome").glob("*_genomic.fna")), None)
        if genome and (not fa.exists() or args.force):
            log("extracting transcript sequences with gffread")
            run(["gffread", "-w", fa, "-g", genome, merged])
        if fa.exists():
            log(f"  {fa.name}: {fa.stat().st_size/1e6:.1f} MB")
    else:
        log("gffread not found; skipping FASTA extraction")

    report = {
        "n_input_gtfs": len(gtfs),
        "merge_params": {"min_tpm": args.min_tpm, "min_len": args.min_len,
                         "reference_guided": args.with_reference},
        "stats": stats,
        "gates": {k: {"pass": ok, "detail": m} for k, ok, m in results},
        "upstream_comparison": {
            "upstream_transcripts": 778680,
            "upstream_unigenes": 586928,
            "upstream_n50": 417,
            "upstream_mean_length": 409,
        },
    }
    rp = cfg.p("results") / "assembly_merge.json"
    rp.parent.mkdir(parents=True, exist_ok=True)
    rp.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {rp}")

    if n_pass == len(results):
        print("\nNext: TransDecoder at 70 aa (and 100 aa for comparison)")
        print(f"  TransDecoder -t {fa} -m 70 --genome ... --gtf {merged}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
