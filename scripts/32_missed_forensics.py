#!/usr/bin/env python3
"""Why did we miss a gene at 25,496 TPM? Classify every expressed miss.

    python scripts/32_missed_forensics.py
    python scripts/32_missed_forensics.py --top 40

THE ANOMALY

  30_expressed_denominator.py raised detection from 36.7% against all
  annotated genes to 78.1% against the 10,959 that are expressed, and the
  per-length table showed why: short genes are expressed in hemocytes
  roughly half as often as long ones (33.2% vs 64.6%). Most of the
  gradient I previously called a pipeline failure is tissue biology.

  But the list of expressed-yet-undetected genes is led by:

      25,496 TPM  uncharacterized LOC119584357
      14,979 TPM  small proline-rich protein 2D-like
      14,961 TPM  uncharacterized LOC119584358
      14,711 TPM  chelonianin-like
      12,192 TPM  uncharacterized LOC119576431

  A gene at 25,496 TPM is ~2.5% of a library. It is not possible to fail
  to ASSEMBLE that. So "not detected" here means our detection TEST failed,
  not our assembly — and until we know which step, 78.1% is not a number
  to report. Six of the top fifteen are "uncharacterized LOC...", which is
  precisely the dark-proteome target set the project now rests on, so this
  cannot be left as a footnote.

SEVEN MECHANICAL EXPLANATIONS, AND WHICH ONE IT IS

  Detection required three things at once: an ORF of ours on the gene's
  annotated exons, on the SAME strand, matching one of the gene's proteins
  at >= 70% identity over >= 10% of its length. Any of these can fail
  independently:

    NO_TRANSCRIPT        nothing we assembled touches its exons at all
    TRANSCRIPT_NO_ORF    a transcript is there, no ORF was called on it
    ORF_WRONG_STRAND     our ORF is on its exons but ANTISENSE. 49% of our
                         transcripts are unstranded (StringTie cannot
                         orient a single-exon transcript), so TransDecoder
                         picks a strand and can pick the wrong one. This is
                         my leading suspicion for short, highly expressed,
                         single-exon genes.
    ORF_UNPLACED         our protein matches it by sequence but carries no
                         genomic coordinates — 567 proteins are in that
                         state and can never pass a locus-aware test
    SEQ_MATCH_ELSEWHERE  matched, placed, but at a different locus
    LOW_IDENTITY         right place, right strand, identity below 70%
    LOW_COVERAGE         right place, right strand, identity fine, but the
                         alignment covers under 10% of the protein

  The first two are real misses. The next three are bookkeeping, and if
  the top-expressed genes land there, our detection figure is understated
  and the fix costs nothing but code.

NOTHING HERE CHANGES ANY OUTPUT. It only measures.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
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

FLIP = {"+": "-", "-": "+", ".": "."}

ORDER = ["NO_TRANSCRIPT", "TRANSCRIPT_NO_ORF", "MITO_GENETIC_CODE",
         "ORF_WRONG_STRAND", "ORF_UNPLACED", "SEQ_MATCH_ELSEWHERE",
         "LOW_IDENTITY", "LOW_COVERAGE", "UNEXPLAINED"]

REAL_MISS = {"NO_TRANSCRIPT", "TRANSCRIPT_NO_ORF"}

# Mitochondrially encoded genes are called with the WRONG GENETIC CODE.
#
# TransDecoder ran with -G universal. Invertebrate mitochondria use
# translation table 5, where TGA is tryptophan and AGA/AGG are serine, so
# under the universal code a mito ORF hits premature stops and no ORF is
# called. Every one of the loudest "real misses" was mitochondrial:
# COX1/2/3, cytochrome b, ATP synthase F0-6, NADH dehydrogenase 1/4/5 —
# all single-exon, all among the most abundant transcripts in any library.
#
# There are ~13 protein-coding mito genes, so the count is small, but they
# sit at the top of every expression ranking and were inflating the REAL
# MISS tally with something that is not an assembly failure at all. They
# are separated out here. The reference annotation already provides these
# proteins, so re-calling them is optional; if we do, it needs a second
# TransDecoder pass restricted to the mitochondrial contig with the
# invertebrate mitochondrial code.
MITO_CONTIGS = {"NC_002184.1"}


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top", type=int, default=30,
                    help="how many of the highest-expressed misses to list")
    ap.add_argument("--min-bp", type=int, default=30)
    ap.add_argument("--min-identity", type=float, default=70.0)
    ap.add_argument("--detect-qcov", type=float, default=10.0)
    args = ap.parse_args()

    cfg = cfgmod.load()
    out, res = cfg.p("proteins"), cfg.p("results")
    missed_file = res / "genes_expressed_not_detected.txt"
    if not missed_file.exists():
        raise SystemExit(f"{missed_file} missing — run "
                         "scripts/30_expressed_denominator.py")

    missed: dict[str, tuple[float, str]] = {}
    for line in missed_file.read_text().splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        try:
            missed[f[0]] = (float(f[1]), f[2] if len(f) > 2 else "")
        except (IndexError, ValueError):
            continue
    log(f"{len(missed):,} expressed-but-undetected genes to explain")

    ref_gff = next(iter(cfg.p("genome").glob("*_genomic.gff")), None)
    qc = json.loads((res / "protein_qc.json").read_text())
    pep, gff3 = Path(qc["pep"]), Path(qc["genome_gff3"])

    log("building the annotation model")
    A = _rx.annotation(ref_gff)
    gene_prots: dict[str, set[str]] = defaultdict(set)
    for p, g in A["prot2gene"].items():
        gene_prots[g].add(p)

    log("placing our ORFs on both strands")
    cds = _rx.predicted_cds_blocks(gff3)
    our_len = _qc.read_fasta_lengths(pep)
    unplaced = {p for p in our_len if p not in cds}
    log(f"  {len(unplaced):,} of our proteins have no genomic coordinates")

    same: dict[str, set[str]] = defaultdict(set)
    opp: dict[str, set[str]] = defaultdict(set)
    where: dict[str, set[str]] = defaultdict(set)
    for pid, (chrom, strand, blocks) in cds.items():
        s = _rx.exonic_overlap(A["exon_idx"], chrom, strand, blocks,
                               args.min_bp)
        o = _rx.exonic_overlap(A["exon_idx"], chrom, FLIP[strand], blocks,
                               args.min_bp)
        for g in s:
            same[g].add(pid)
        for g in o - s:
            opp[g].add(pid)
        where[pid] = s

    log("overlapping assembled transcripts, either strand")
    tx = _rx.gtf_exon_blocks(out / "merged_filtered.gtf")
    has_tx: set[str] = set()
    for chrom, strand, blocks in tx.values():
        for st in (("+", "-") if strand == "." else (strand, FLIP[strand])):
            has_tx |= _rx.exonic_overlap(A["exon_idx"], chrom, st, blocks,
                                         args.min_bp)

    log("reading the diamond table (any locus)")
    tsv = out / "ref_vs_ours_multi.tsv"
    if not tsv.exists():
        raise SystemExit(f"{tsv} missing — run scripts/29_fragmentation.py")
    best_any: dict[str, tuple[str, float, float]] = {}
    with tsv.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            g = A["prot2gene"].get(f[0])
            if g is None:
                continue
            try:
                ident, qcov = float(f[2]), float(f[3])
            except ValueError:
                continue
            cur = best_any.get(g)
            if cur is None or (qcov, ident) > (cur[2], cur[1]):
                best_any[g] = (f[1], ident, qcov)

    # Which missed genes live on the mitochondrial genome?
    mito_genes = {g for g, sp in A["span"].items()
                  if sp[0] in MITO_CONTIGS}
    log(f"  {len(mito_genes)} annotated genes on "
        f"{sorted(MITO_CONTIGS)} (mitochondrion)")

    log("classifying")
    cls: dict[str, str] = {}
    for g in missed:
        hit = best_any.get(g)
        on_same = same.get(g) or set()
        on_opp = opp.get(g) or set()

        if g in mito_genes and not on_same:
            cls[g] = "MITO_GENETIC_CODE"
        elif hit and hit[0] in unplaced:
            cls[g] = "ORF_UNPLACED"
        elif on_same:
            # placed on its exons, right strand — so a threshold rejected it
            if hit is None:
                cls[g] = "LOW_IDENTITY"
            elif hit[1] < args.min_identity:
                cls[g] = "LOW_IDENTITY"
            elif hit[2] < args.detect_qcov:
                cls[g] = "LOW_COVERAGE"
            else:
                # matched well but the matching protein is not one of the
                # ones sitting on this gene
                cls[g] = ("SEQ_MATCH_ELSEWHERE" if hit[0] not in on_same
                          else "UNEXPLAINED")
        elif on_opp:
            cls[g] = "ORF_WRONG_STRAND"
        elif hit is not None:
            cls[g] = "SEQ_MATCH_ELSEWHERE"
        elif g in has_tx:
            cls[g] = "TRANSCRIPT_NO_ORF"
        else:
            cls[g] = "NO_TRANSCRIPT"

    counts = defaultdict(int)
    for v in cls.values():
        counts[v] += 1
    n = len(missed)

    print("\n" + "=" * 78)
    print(f"WHY {n:,} EXPRESSED GENES WERE NOT DETECTED")
    print("=" * 78)
    note = {
        "NO_TRANSCRIPT": "REAL MISS — nothing assembled there",
        "TRANSCRIPT_NO_ORF": "REAL MISS — transcript but no ORF",
        "ORF_WRONG_STRAND": "bookkeeping — antisense ORF call",
        "ORF_UNPLACED": "bookkeeping — our protein has no coordinates",
        "SEQ_MATCH_ELSEWHERE": "bookkeeping — matched at another locus",
        "LOW_IDENTITY": "threshold — under 70% identity",
        "LOW_COVERAGE": "threshold — under 10% coverage",
        "MITO_GENETIC_CODE": "wrong code — table 5, not universal",
        "UNEXPLAINED": "needs a closer look",
    }
    for k in ORDER:
        if counts[k]:
            print(f"  {k:<22}{counts[k]:>7,}  {100*counts[k]/n:>5.1f}%  "
                  f"{note[k]}")
    print("-" * 78)
    real = sum(counts[k] for k in REAL_MISS)
    book = n - real
    print(f"  {'REAL MISSES':<22}{real:>7,}  {100*real/n:>5.1f}%")
    print(f"  {'recoverable':<22}{book:>7,}  {100*book/n:>5.1f}%")

    # What detection would be if the recoverable classes were fixed.
    expd = json.loads((res / "expressed_denominator.json").read_text())
    nexp = expd["expressed"]
    det = expd["detected_expressed"]
    print(f"\n  detection now                 {100*det/nexp:>5.1f}%  "
          f"({det:,}/{nexp:,})")
    print(f"  if the recoverable are fixed  {100*(det+book)/nexp:>5.1f}%  "
          f"({det+book:,}/{nexp:,})")
    print("  the second number is an upper bound, not a result: it assumes")
    print("  every recoverable case really is our gene. Each class needs its")
    print("  own fix before any of it can be claimed.")

    # ---- weight by expression: the loud misses matter most -------------
    print("\n" + "=" * 78)
    print(f"THE {args.top} HIGHEST-EXPRESSED MISSES, CLASSIFIED")
    print("=" * 78)
    refp = next(iter(cfg.p("genome").glob("*_protein.faa")), None)
    ref_len = _qc.read_fasta_lengths(refp) if refp else {}
    glen: dict[str, int] = {}
    for g, ps in gene_prots.items():
        L = [ref_len[p] for p in ps if p in ref_len]
        if L:
            glen[g] = max(L)

    top = sorted(missed, key=lambda g: -missed[g][0])[:args.top]
    print(f"  {'TPM':>10}  {'class':<21}{'aa':>6}{'ex':>4}  product")
    print("-" * 78)
    for g in top:
        tpm, prodname = missed[g]
        aa = glen.get(g, 0)
        nex = len(A["exons"].get(g, []))
        print(f"  {tpm:>10.1f}  {cls[g]:<21}{aa:>6}{nex:>4}  {prodname[:32]}")

    hi = [g for g in missed if missed[g][0] >= 1000]
    if hi:
        hc = defaultdict(int)
        for g in hi:
            hc[cls[g]] += 1
        print("\n" + "=" * 78)
        print(f"GENES ABOVE 1000 TPM ({len(hi):,} of them) — these cannot be "
              "assembly failures")
        print("=" * 78)
        for k in ORDER:
            if hc[k]:
                print(f"  {k:<22}{hc[k]:>7,}  {100*hc[k]/len(hi):>5.1f}%")
        hreal = sum(hc[k] for k in REAL_MISS)
        print("-" * 78)
        print(f"  real misses among them {hreal:,} ({100*hreal/len(hi):.1f}%)")
        if hreal / max(1, len(hi)) < 0.25:
            print("  -> the loud misses are overwhelmingly bookkeeping, which")
            print("     is the expected answer and the one we can fix")
        else:
            print("  -> a quarter or more are genuine. That is a pipeline")
            print("     problem at high expression and must be fixed first")

    report = {
        "n_missed": n,
        "classes": dict(counts),
        "real_misses": real,
        "recoverable": book,
        "detection_now_pct": round(100 * det / nexp, 1),
        "detection_upper_bound_pct": round(100 * (det + book) / nexp, 1),
        "above_1000_tpm": {"n": len(hi),
                           "classes": {k: hc[k] for k in ORDER if hc[k]}}
        if hi else None,
        "unplaced_proteins": len(unplaced),
    }
    rp = res / "missed_forensics.json"
    rp.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {rp}")

    f = res / "missed_classified.tsv"
    f.write_text("gene\tmax_tpm\tclass\tn_exons\tproduct\n" + "\n".join(
        f"{g}\t{missed[g][0]:.2f}\t{cls[g]}\t{len(A['exons'].get(g, []))}\t"
        f"{missed[g][1]}"
        for g in sorted(missed, key=lambda g: -missed[g][0])) + "\n")
    print(f"wrote {f.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
