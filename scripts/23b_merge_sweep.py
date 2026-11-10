#!/usr/bin/env python3
"""Tune the StringTie merge against the reference gene count.

    python scripts/23b_merge_sweep.py --threads 16

WHY

  The first merge (-T 1.0 -m 200) produced 100,781 loci and 111,884
  transcripts, of which only 27.2% are multi-exon and the median length is
  396 bp. That is six times better than the upstream de novo assembly
  (586,928 unigenes, 417 bp N50) but still roughly 3x too many loci for a
  penaeid, and the single-exon fraction gives the reason: low-coverage
  pileups in the worst-mapping libraries became "transcripts".

  Per-library transcript counts track mapping rate inversely:

      SRR28572183  84,640 transcripts   56.5% mapped
      SRR28572172  56,772               67.6%
      SRR28572180  48,744               60.7%
      SRR28572173  19,955               78.2%

  So this sweeps the merge thresholds and reports, for each setting, how
  the result compares with the NCBI annotation of the SAME genome — which
  is the only defensible target. A guessed range is not.

  Merging takes ~8 s, so the whole sweep is a couple of minutes.

IMPORTANT

  Thresholds are chosen to match the reference gene count and to raise the
  multi-exon fraction, NOT to maximise the number of novel loci. Choosing
  the setting that yields the most unannotated transcripts would
  manufacture the paper's own result.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod  # noqa: E402

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "mg", Path(__file__).resolve().parent / "23_merge_assembly.py")
_mg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mg)
gtf_stats = _mg.gtf_stats

GRID_T = (1.0, 2.0, 5.0, 10.0)     # stringtie --merge -T, min TPM
GRID_M = (200, 300, 500)           # -m, min transcript length
GRID_F = (0.01, 0.05, 0.10)        # -f, min isoform fraction


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def reference_counts(gff: Path) -> dict:
    """Gene and mRNA counts from the NCBI annotation of this same genome."""
    genes = mrnas = 0
    with gff.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 3:
                continue
            if p[2] == "gene":
                genes += 1
            elif p[2] in ("mRNA", "transcript"):
                mrnas += 1
    return {"ref_genes": genes, "ref_mrnas": mrnas}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--keep", action="store_true",
                    help="keep every intermediate GTF (default: only the best)")
    args = ap.parse_args()

    cfg = cfgmod.load()
    out = cfg.p("proteins")
    lst = out / "gtf_list.txt"
    if not lst.exists():
        raise SystemExit("run scripts/23_merge_assembly.py first")
    if not shutil.which("stringtie"):
        raise SystemExit("stringtie not on PATH")

    gff = next(iter(cfg.p("genome").glob("*_genomic.gff")), None)
    ref = reference_counts(gff) if gff else {}
    if ref:
        log(f"NCBI annotation of the same genome: "
            f"{ref['ref_genes']:,} genes, {ref['ref_mrnas']:,} mRNAs")
        log("  (the target: our locus count should land near ref_genes, and")
        log("   exceed it only by however many novel loci genuinely exist)")
    else:
        log("no reference GFF found; comparing on structure only")

    sweep_dir = out / "merge_sweep"
    sweep_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    total = len(GRID_T) * len(GRID_M) * len(GRID_F)
    n = 0
    for T in GRID_T:
        for m in GRID_M:
            for f in GRID_F:
                n += 1
                tag = f"T{T}_m{m}_f{f}"
                gtf = sweep_dir / f"merged_{tag}.gtf"
                if not gtf.exists():
                    subprocess.run([
                        "stringtie", "--merge", "-p", str(args.threads),
                        "-T", str(T), "-m", str(m), "-f", str(f),
                        "-o", str(gtf), str(lst),
                    ], capture_output=True)
                s = gtf_stats(gtf)
                s.update({"T": T, "m": m, "f": f})
                if ref:
                    s["loci_vs_ref"] = round(s["n_loci"] / max(1, ref["ref_genes"]), 2)
                rows.append(s)
                log(f"  {n}/{total}  {tag:<22} loci {s['n_loci']:>7,}  "
                    f"N50 {s['transcript_n50']:>5,}  "
                    f"multi-exon {s['multi_exon_pct']:>5.1f}%")

    print("\n" + "=" * 92)
    print("MERGE SWEEP")
    if ref:
        print(f"reference: {ref['ref_genes']:,} genes / {ref['ref_mrnas']:,} mRNAs "
              "(NCBI Annotation Release 100, same assembly)")
    print("=" * 92)
    print(f"{'T':>6}{'m':>5}{'f':>6}{'loci':>9}{'xRef':>7}{'tx':>9}"
          f"{'N50':>7}{'mean':>8}{'median':>8}{'multiEx%':>10}")
    print("-" * 92)
    for s in sorted(rows, key=lambda r: -r["multi_exon_pct"]):
        xr = f"{s.get('loci_vs_ref', 0):.2f}" if ref else "-"
        print(f"{s['T']:>6.1f}{s['m']:>5}{s['f']:>6.2f}{s['n_loci']:>9,}"
              f"{xr:>7}{s['n_transcripts']:>9,}{s['transcript_n50']:>7,}"
              f"{s['mean_transcript']:>8.0f}{s['median_transcript']:>8,}"
              f"{s['multi_exon_pct']:>10.1f}")
    print("-" * 92)

    # Selection: closest to the reference gene count, breaking ties by the
    # multi-exon fraction. Deliberately NOT "most loci".
    if ref:
        # Bucket the locus distance to 2% of the reference count before
        # breaking ties on quality. The first version compared exact
        # distances and so preferred 28,612 loci over 28,618 — a 6-locus
        # difference — even though the latter had a higher N50 (2,180 vs
        # 2,042) and more multi-exon transcripts (50.6% vs 47.5%). At that
        # resolution the locus count is noise and the structure is signal.
        tol = max(1, int(0.02 * ref["ref_genes"]))
        best = min(rows, key=lambda r: (abs(r["n_loci"] - ref["ref_genes"]) // tol,
                                        -r["multi_exon_pct"],
                                        -r["transcript_n50"]))
    else:
        best = max(rows, key=lambda r: r["multi_exon_pct"])

    print(f"\nSELECTED  -T {best['T']}  -m {best['m']}  -f {best['f']}")
    print(f"  loci {best['n_loci']:,}"
          + (f"  ({best['loci_vs_ref']}x the reference gene count)" if ref else ""))
    print(f"  transcripts {best['n_transcripts']:,}   N50 {best['transcript_n50']:,} bp"
          f"   multi-exon {best['multi_exon_pct']}%")

    chosen = sweep_dir / f"merged_T{best['T']}_m{best['m']}_f{best['f']}.gtf"
    final = out / "merged.gtf"
    shutil.copy2(chosen, final)
    print(f"\ncopied to {final}")

    if not args.keep:
        for g in sweep_dir.glob("merged_*.gtf"):
            if g != chosen:
                g.unlink()
        print("removed other sweep GTFs (--keep to retain)")

    rp = cfg.p("results") / "merge_sweep.json"
    rp.write_text(json.dumps({"reference": ref, "selected": best,
                              "all": rows}, indent=2))
    print(f"wrote {rp}")

    print("\nRegenerate the transcript FASTA from the selected GTF:")
    print("  python scripts/23_merge_assembly.py --threads 16   # reuses merged.gtf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
