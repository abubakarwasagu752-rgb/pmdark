#!/usr/bin/env python3
"""Detection among genes that are actually expressed. The real denominator.

    python scripts/30_expressed_denominator.py
    python scripts/30_expressed_denominator.py --min-tpm 1 --min-libraries 3

THE FLAW IN EVERY RECOVERY NUMBER SO FAR

  36.7% detected, 21.0% complete, 32.4% before that, 39.6% before that —
  every one of them divides by 24,011, which is every protein-coding gene
  NCBI annotates in this genome. That denominator includes opsins, cuticle
  proteins, vitellogenin and chemoreceptors. 28_recovery_by_family.py
  already showed those come back at 0-11%, correctly, because a hemocyte
  does not transcribe them. Dividing by them measures tissue specificity
  and calls it pipeline performance.

  The denominator we need is the set of genes with READS in these
  libraries. We had never measured it: every figure to date came from our
  own assembly, and an assembly cannot distinguish a gene with ten reads
  from a gene with none. Both yield no transcript.

  21_align_assemble.py now runs a second StringTie pass with -e -G, which
  quantifies the reference annotation and only the reference annotation.
  This reads those tables.

WHY THE LENGTH PATTERN MADE THIS UNAVOIDABLE

  29_fragmentation.py found detection RISING with protein length while
  coverage fell:

      0-150 aa     3,579 genes   detected 19.8%   median coverage 100%
      150-300      7,252         detected 28.5%   median coverage  99%
      300-600      7,953         detected 41.1%   median coverage  89%
      600-1200     3,891         detected 52.7%   median coverage  75%
      1200+        1,336         detected 54.2%   median coverage  56%

  Fragmentation predicts the opposite of the first column, and where we do
  detect a gene we recover 88% of it (median) with 93% of its length. So
  the assembly is not badly broken. Short genes are simply not being seen,
  and the two candidate reasons — they are not expressed in hemocytes, or
  they are expressed and we missed them — are exactly what a read-level
  denominator separates. The answer also differs per length band, which is
  why this reports the breakdown rather than one number.

  PRE-STATED READING, fixed before the run:

    detection among expressed genes >= 85%  -> the pipeline is sound and
        the low overall figures are tissue specificity. Report detection
        against the expressed denominator, and state both numbers.
    detection among expressed 70-85%        -> acceptable with the caveat
        quantified per length band.
    detection among expressed < 70%         -> we are losing expressed
        genes. Fix that before the dark proteome is defined.

  A gene counts as expressed at >= --min-tpm in >= --min-libraries of the
  11 usable libraries. SRR28572179 is excluded: 830k reads, flagged
  use=0 in config/samples.tsv. Both thresholds are reported as a
  sensitivity grid rather than a single choice, because the headline must
  not depend on a cutoff picked after seeing the result.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import statistics
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod  # noqa: E402


def _load(name: str):
    sp = importlib.util.spec_from_file_location(
        "m" + name[:2], Path(__file__).resolve().parent / name)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


_rx = _load("27_recovery_exact.py")
_qc = _load("25_protein_qc.py")
_fam = _load("28_recovery_by_family.py")


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def pct(n: int, d: int) -> str:
    return f"{100*n/max(1,d):.1f}%"


def read_abund(p: Path) -> dict[str, float]:
    """StringTie -A gene table: Gene ID, Gene Name, Reference, Strand,
    Start, End, Coverage, FPKM, TPM."""
    out: dict[str, float] = {}
    with p.open() as fh:
        head = fh.readline().rstrip("\n").split("\t")
        try:
            gi = head.index("Gene ID")
            ti = head.index("TPM")
        except ValueError:
            gi, ti = 0, len(head) - 1
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(gi, ti):
                continue
            try:
                v = float(f[ti])
            except ValueError:
                continue
            g = f[gi]
            if v > out.get(g, -1.0):
                out[g] = v
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-tpm", type=float, default=1.0)
    ap.add_argument("--min-libraries", type=int, default=3)
    ap.add_argument("--min-bp", type=int, default=30)
    ap.add_argument("--min-identity", type=float, default=70.0)
    ap.add_argument("--detect-qcov", type=float, default=10.0)
    ap.add_argument("--complete-qcov", type=float, default=80.0)
    args = ap.parse_args()

    cfg = cfgmod.load()
    out, res = cfg.p("proteins"), cfg.p("results")

    tabs = sorted(cfg.p("align").glob("*/*.ref_gene_abund.tsv"))
    if not tabs:
        print("No *.ref_gene_abund.tsv found under", cfg.p("align"))
        print()
        print("These come from the -e -G pass added to 21_align_assemble.py.")
        print("The BAMs were cleaned up, so the libraries must be re-run:")
        print()
        print("    sbatch slurm/21b_realign.sbatch")
        print()
        print("That also switches the assembly to -G (guided, novel loci")
        print("still reported), which is the other correction from 29.")
        return 2

    log(f"{len(tabs)} reference-quantification tables")
    per_lib: dict[str, dict[str, float]] = {}
    for t in tabs:
        acc = t.name.split(".")[0]
        per_lib[acc] = read_abund(t)
        log(f"  {acc:<14} {len(per_lib[acc]):>7,} genes")

    # The failed library contributes transcript structure, never statistics.
    # Parse by HEADER NAME. An earlier version took the last column, which
    # in this sheet is `note`, not `use` — so SRR28572179's note of
    # FAILED_LIBRARY_16x_shallow was read as a use flag and the row was
    # kept. Only a hardcoded fallback was excluding it.
    bad: set[str] = set()
    st = Path(__file__).resolve().parents[1] / "config" / "samples.tsv"
    if st.exists():
        rows = [ln for ln in st.read_text().splitlines()
                if ln.strip() and not ln.startswith("#")]
        if rows:
            head = rows[0].split("\t")
            try:
                acc_i, use_i = head.index("run_accession"), head.index("use")
            except ValueError:
                acc_i, use_i = 0, -1
            if use_i >= 0:
                for ln in rows[1:]:
                    f = ln.split("\t")
                    if len(f) > max(acc_i, use_i) and \
                            f[use_i].strip() in ("0", "FALSE", "false", "no"):
                        bad.add(f[acc_i].strip())
    if not bad:
        log("WARNING no use=0 rows parsed from samples.tsv; "
            "falling back to the known failed library")
        bad = {"SRR28572179"}
    use = {k: v for k, v in per_lib.items() if k not in bad}
    log(f"using {len(use)} libraries; excluded {sorted(bad & set(per_lib))}")

    log("building the annotation model")
    ref_gff = next(iter(cfg.p("genome").glob("*_genomic.gff")), None)
    A = _rx.annotation(ref_gff)
    coding = A["coding"]
    prod = _fam.gene_products(ref_gff)
    refp = next(iter(cfg.p("genome").glob("*_protein.faa")), None)
    ref_len = _qc.read_fasta_lengths(refp) if refp else {}

    # StringTie -A reports the GFF's gene name/id; map onto our gene keys.
    alias: dict[str, str] = {}
    for g in coding:
        alias[g] = g
        if g.startswith("gene-"):
            alias[g[5:]] = g

    count: dict[str, int] = defaultdict(int)
    maxtpm: dict[str, float] = defaultdict(float)
    for acc, tab in use.items():
        for gid, tpm in tab.items():
            g = alias.get(gid) or alias.get("gene-" + gid)
            if g is None:
                continue
            if tpm >= args.min_tpm:
                count[g] += 1
            maxtpm[g] = max(maxtpm[g], tpm)
    expressed = {g for g, n in count.items() if n >= args.min_libraries}
    log(f"  expressed: {len(expressed):,} of {len(coding):,} coding genes "
        f"(>={args.min_tpm} TPM in >={args.min_libraries} libraries)")

    # ---- our detection sets, locus-aware -------------------------------
    qc = json.loads((res / "protein_qc.json").read_text())
    cds = _rx.predicted_cds_blocks(Path(qc["genome_gff3"]))
    our_len = _qc.read_fasta_lengths(Path(qc["pep"]))
    prot_genes = {}
    for pid in our_len:
        r = cds.get(pid)
        prot_genes[pid] = (_rx.exonic_overlap(A["exon_idx"], r[0], r[1],
                                              r[2], args.min_bp) & coding
                           if r else set())

    tsv = out / "ref_vs_ours_multi.tsv"
    if not tsv.exists():
        raise SystemExit(f"{tsv} missing — run scripts/29_fragmentation.py")
    best: dict[str, tuple[float, float]] = {}
    with tsv.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            g = A["prot2gene"].get(f[0])
            if g is None or g not in prot_genes.get(f[1], ()):
                continue
            try:
                ident, qcov = float(f[2]), float(f[3])
            except ValueError:
                continue
            cur = best.get(g)
            if cur is None or (qcov, ident) > (cur[1], cur[0]):
                best[g] = (ident, qcov)
    detected = {g for g, (i, q) in best.items()
                if i >= args.min_identity and q >= args.detect_qcov}
    complete = {g for g, (i, q) in best.items()
                if i >= args.min_identity and q >= args.complete_qcov}

    tot, nexp = len(coding), len(expressed)
    print("\n" + "=" * 78)
    print("DETECTION AGAINST TWO DENOMINATORS")
    print("=" * 78)
    print(f"  {'':<34}{'all annotated':>16}{'expressed only':>18}")
    print(f"  {'genes':<34}{tot:>16,}{nexp:>18,}")
    print("-" * 78)
    for lab, s in (("detected at its own locus", detected),
                   ("complete (>=80% of protein)", complete)):
        print(f"  {lab:<34}{pct(len(s), tot):>16}"
              f"{pct(len(s & expressed), nexp):>18}")
    print("-" * 78)
    miss = expressed - detected
    print(f"  expressed but NOT detected        {len(miss):>16,}")
    print("  that is the only number describing genes we actually lost")

    # ---- sensitivity to the cutoffs ------------------------------------
    print("\n" + "=" * 78)
    print("SENSITIVITY — the headline must not depend on a chosen cutoff")
    print("=" * 78)
    print(f"  {'min TPM':>9}{'min libs':>10}{'expressed':>12}"
          f"{'detected':>11}{'complete':>11}")
    print("-" * 78)
    grid = []
    for tpm in (0.5, 1.0, 2.0, 5.0):
        for ml in (1, 3, 6):
            c = defaultdict(int)
            for tab in use.values():
                for gid, v in tab.items():
                    g = alias.get(gid) or alias.get("gene-" + gid)
                    if g is not None and v >= tpm:
                        c[g] += 1
            ex = {g for g, n in c.items() if n >= ml}
            row = {"min_tpm": tpm, "min_libs": ml, "expressed": len(ex),
                   "detected_pct": round(100 * len(detected & ex) / max(1, len(ex)), 1),
                   "complete_pct": round(100 * len(complete & ex) / max(1, len(ex)), 1)}
            grid.append(row)
            print(f"  {tpm:>9.1f}{ml:>10}{len(ex):>12,}"
                  f"{row['detected_pct']:>10.1f}%{row['complete_pct']:>10.1f}%")

    # ---- per length band, where the question actually lives ------------
    glen: dict[str, int] = {}
    for p, g in A["prot2gene"].items():
        if p in ref_len:
            glen[g] = max(glen.get(g, 0), ref_len[p])
    print("\n" + "=" * 78)
    print("BY PROTEIN LENGTH, EXPRESSED GENES ONLY")
    print("=" * 78)
    print(f"  {'length':<14}{'annotated':>11}{'expressed':>11}"
          f"{'% expressed':>13}{'detected':>11}{'complete':>11}")
    print("-" * 78)
    bands = [(0, 150), (150, 300), (300, 600), (600, 1200), (1200, 10 ** 9)]
    by_band = []
    for lo, hi in bands:
        gs = {g for g, L in glen.items() if lo <= L < hi} & coding
        if not gs:
            continue
        ex = gs & expressed
        tag = f"{lo}-{hi} aa" if hi < 10 ** 9 else f"{lo}+ aa"
        print(f"  {tag:<14}{len(gs):>11,}{len(ex):>11,}"
              f"{pct(len(ex), len(gs)):>13}"
              f"{pct(len(detected & ex), len(ex)):>11}"
              f"{pct(len(complete & ex), len(ex)):>11}")
        by_band.append({"band": tag, "annotated": len(gs), "expressed": len(ex),
                        "detected": len(detected & ex),
                        "complete": len(complete & ex)})
    print("-" * 78)
    print("  the '% expressed' column is the test: if short genes are simply")
    print("  not expressed, detection within the band will be normal and the")
    print("  19.8% in 29's table was tissue biology, not a pipeline failure")

    overall = 100 * len(detected & expressed) / max(1, nexp)
    print("\n" + "=" * 78)
    print("VERDICT  (rule fixed in this file's docstring before the run)")
    print("=" * 78)
    print(f"  detection among expressed genes  {overall:.1f}%")
    if overall >= 85:
        v = ("PASS — the pipeline is sound. The low figures against all "
             "24,011\n  genes are tissue specificity. Report both "
             "denominators.")
    elif overall >= 70:
        v = ("ACCEPTABLE — carry the caveat, quantified per length band.")
    else:
        v = ("FAIL — expressed genes are being lost. Fix that before the\n"
             "  dark proteome is defined on top of it.")
    print(f"  {v}")

    # ---- what are the expressed-but-missed genes? ---------------------
    if miss:
        hi = sorted(miss, key=lambda g: -maxtpm.get(g, 0))[:15]
        print("\n" + "=" * 78)
        print("MOST HIGHLY EXPRESSED GENES WE FAILED TO DETECT")
        print("=" * 78)
        for g in hi:
            print(f"  {maxtpm.get(g, 0):>9.1f} TPM  "
                  f"{(prod.get(g, '') or g)[:56]}")
        print("\n  a highly expressed gene we missed is a real bug; a list of")
        print("  1-2 TPM genes is the detection limit doing its job")

    report = {
        "coding_genes": tot,
        "thresholds": {"min_tpm": args.min_tpm,
                       "min_libraries": args.min_libraries,
                       "identity": args.min_identity,
                       "detect_qcov": args.detect_qcov,
                       "complete_qcov": args.complete_qcov},
        "libraries_used": sorted(use), "libraries_excluded": sorted(bad),
        "expressed": nexp,
        "detected_all": len(detected), "complete_all": len(complete),
        "detected_expressed": len(detected & expressed),
        "complete_expressed": len(complete & expressed),
        "detected_pct_expressed": round(overall, 1),
        "expressed_not_detected": len(miss),
        "sensitivity_grid": grid, "by_length_band": by_band,
        "median_max_tpm_of_missed": round(
            statistics.median([maxtpm.get(g, 0) for g in miss]), 2) if miss else None,
    }
    rp = res / "expressed_denominator.json"
    rp.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {rp}")

    f = res / "genes_expressed_not_detected.txt"
    f.write_text("\n".join(
        f"{g}\t{maxtpm.get(g, 0):.2f}\t{prod.get(g, '')}"
        for g in sorted(miss, key=lambda g: -maxtpm.get(g, 0))) + "\n")
    print(f"wrote {f.name} ({len(miss):,} genes, by descending TPM)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
