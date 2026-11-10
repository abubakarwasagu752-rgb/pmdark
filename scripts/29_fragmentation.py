#!/usr/bin/env python3
"""Locus-aware recovery, and how fragmented the assembly actually is.

    python scripts/29_fragmentation.py --threads 16

TWO THINGS TO FIX, ONE OF THEM A REAL BUG IN MY METRIC

  1. PARALOG CROSS-MATCHING INVALIDATED THE SEQUENCE CRITERION.
     28_recovery_by_family.py reported, for actin: 34 genes transcribed,
     32 with an ORF on their exons, and 99 "recovered by sequence". That
     is impossible if the criteria nest — and they did not nest, because I
     built the sequence set from the diamond table alone, with no reference
     to WHERE our matching protein sits in the genome. P. monodon has 109
     annotated actin genes. One assembled actin protein matches most of
     them at over 70% identity, so one real recovery was credited as 99.
     Same for histones (160 genes, 91 "recovered" from 58 ORFs), cuticle
     proteins (632 genes, 35 from 7), tubulins, chaperones, caspases.

     The fix is to require LOCUS AGREEMENT: a gene counts as recovered
     only if one of ITS OWN annotated proteins matches one of our proteins
     AND that protein of ours sits on that gene's exons. Paralogous
     matches elsewhere in the genome no longer count. This needs more than
     one hit per query, so the search is re-run with --max-target-seqs 25;
     with 1, a reference protein whose best hit was a paralog's locus had
     no chance to find its own.

  2. DETECTION AND COMPLETENESS WERE CONFLATED.
     "Recovered" required >= 50% query coverage, which answers "did we
     assemble this gene fully", not "did we see it". Those are different
     questions and the data separates them cleanly:

         1,682 of 1,976 failures: high identity, LOW coverage
         median identity 100%, median query coverage 36%

     100% identity means the sequence is right. 36% coverage means we have
     a third of each protein. Nothing is lost — it is fragmented. And the
     positive control is monotonic in protein length, which is the
     signature of exactly that:

         ribosomal protein   ~150 aa   94.4%
         proteasome          ~250 aa   89.5%
         ubiquitin           ~250 aa   69.8%
         spliceosome         ~400 aa   60.6%
         translation factor  ~850 aa   50.3%

     So this reports DETECTED (locus-aware, >= --min-identity, any
     meaningful coverage) and COMPLETE (>= --complete-qcov) separately,
     with the coverage distribution between them.

WHY THE ASSEMBLY IS FRAGMENTED, AND MY PART IN IT

  21_align_assemble.py runs StringTie per library with no -G, and
  23_merge_assembly.py defaults --with-reference off. I wrote in that file
  that a reference guide "biases toward known genes", and used that to
  justify omitting it. That is wrong, and it is worth being precise about
  why: in StringTie it is -e that restricts output to the reference. -G on
  its own supplies the annotation as a GUIDE while still assembling and
  reporting novel transcripts and novel loci. So -G without -e would have
  bridged low-coverage gaps in known genes — exactly the 36%-coverage
  problem — at no cost to novel discovery. The caution I applied belonged
  to a flag I was not using.

  This script quantifies what that cost us, and checks whether the BAMs
  still exist, because that decides whether the fix is a StringTie re-run
  or a full re-alignment of all 12 libraries.

NOTHING HERE CHANGES ANY OUTPUT. It only measures.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import statistics
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod  # noqa: E402


def _load(name: str):
    sp = importlib.util.spec_from_file_location(
        name[:2], Path(__file__).resolve().parent / name)
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


def hist(vals: list[float], edges: list[float], label: str) -> None:
    print(f"\n  {label}")
    n = len(vals)
    for lo, hi in zip(edges, edges[1:]):
        c = sum(1 for v in vals if lo <= v < hi)
        bar = "#" * int(round(40 * c / max(1, n)))
        print(f"    {lo:>3.0f}-{hi:<3.0f}% {c:>7,} {pct(c, n):>7}  {bar}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--min-identity", type=float, default=70.0)
    ap.add_argument("--detect-qcov", type=float, default=10.0,
                    help="coverage above which a gene counts as DETECTED")
    ap.add_argument("--complete-qcov", type=float, default=80.0,
                    help="coverage above which a gene counts as COMPLETE")
    ap.add_argument("--min-bp", type=int, default=30)
    ap.add_argument("--max-targets", type=int, default=25)
    ap.add_argument("--force-search", action="store_true")
    args = ap.parse_args()

    cfg = cfgmod.load()
    out, res = cfg.p("proteins"), cfg.p("results")
    qc = json.loads((res / "protein_qc.json").read_text())
    pep, gff3 = Path(qc["pep"]), Path(qc["genome_gff3"])
    ref_gff = next(iter(cfg.p("genome").glob("*_genomic.gff")), None)
    refp = next(iter(cfg.p("genome").glob("*_protein.faa")), None)
    if ref_gff is None or refp is None:
        raise SystemExit("reference GFF / proteome not staged — run 27 first")

    log("building the annotation model")
    A = _rx.annotation(ref_gff)
    coding = A["coding"]
    prod = _fam.gene_products(ref_gff)

    gene_prots: dict[str, set[str]] = defaultdict(set)
    for p, g in A["prot2gene"].items():
        gene_prots[g].add(p)

    log("reading protein lengths")
    ref_len = _qc.read_fasta_lengths(refp)
    our_len = _qc.read_fasta_lengths(pep)
    log(f"  reference {len(ref_len):,}   ours {len(our_len):,}")

    # ---- where is each of our proteins? -------------------------------
    log("placing our proteins on annotated exons")
    cds = _rx.predicted_cds_blocks(gff3)
    prot_genes: dict[str, set[str]] = {}
    for pid in our_len:
        r = cds.get(pid)
        prot_genes[pid] = (_rx.exonic_overlap(A["exon_idx"], r[0], r[1],
                                              r[2], args.min_bp) & coding
                           if r else set())
    g_orf = set().union(*prot_genes.values()) if prot_genes else set()
    log(f"  {len(g_orf):,} coding genes carry one of our ORFs on their exons")

    # ---- the search, with room for paralogs ---------------------------
    tsv = out / "ref_vs_ours_multi.tsv"
    dm = _rx.find_diamond()
    if dm is None:
        raise SystemExit("diamond not found")
    if args.force_search or not tsv.exists():
        db = out / "ourprots"
        if not Path(str(db) + ".dmnd").exists():
            log("diamond makedb")
            subprocess.run([dm, "makedb", "--in", str(pep), "-d", str(db),
                            "--threads", str(args.threads)],
                           capture_output=True)
        log(f"diamond blastp, --max-target-seqs {args.max_targets} "
            "(one hit per query cannot see past a paralog)")
        rc = subprocess.run([
            dm, "blastp", "--query", str(refp), "--db", str(db),
            "--out", str(tsv), "--outfmt", "6", "qseqid", "sseqid",
            "pident", "qcovhsp",
            "--max-target-seqs", str(args.max_targets),
            "--evalue", "1e-5", "--threads", str(args.threads), "--quiet",
        ]).returncode
        if rc != 0:
            log(f"diamond exited {rc}")
            return 1

    # ---- locus-aware join ---------------------------------------------
    log("joining hits to genomic placement (locus agreement required)")
    best_local: dict[str, tuple[float, float, str, str]] = {}
    best_any: dict[str, tuple[float, float]] = {}
    n_rows = n_local = 0
    with tsv.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            try:
                ident, qcov = float(f[2]), float(f[3])
            except ValueError:
                continue
            n_rows += 1
            rp, op = f[0], f[1]
            g = A["prot2gene"].get(rp)
            if g is None:
                continue
            cur = best_any.get(g)
            if cur is None or (qcov, ident) > (cur[1], cur[0]):
                best_any[g] = (ident, qcov)
            if g in prot_genes.get(op, ()):            # <-- the locus test
                n_local += 1
                cur = best_local.get(g)
                if cur is None or (qcov, ident) > (cur[1], cur[0]):
                    best_local[g] = (ident, qcov, rp, op)
    log(f"  {n_rows:,} hit rows, {n_local:,} with locus agreement")

    detected = {g for g, (i, q, _, _) in best_local.items()
                if i >= args.min_identity and q >= args.detect_qcov}
    complete = {g for g, (i, q, _, _) in best_local.items()
                if i >= args.min_identity and q >= args.complete_qcov}
    paralog_only = {g for g, (i, q) in best_any.items()
                    if i >= args.min_identity and q >= 50.0} - set(best_local)

    tot = len(coding)
    print("\n" + "=" * 78)
    print("RECOVERY, LOCUS-AWARE")
    print("=" * 78)
    print(f"  annotated coding genes                      {tot:>7,}")
    print(f"  exons carry one of our ORFs                 {len(g_orf):>7,}"
          f"  {pct(len(g_orf), tot):>7}")
    print(f"  DETECTED  (>={args.min_identity:.0f}% id at its own locus)      "
          f"{len(detected):>7,}  {pct(len(detected), tot):>7}")
    print(f"  COMPLETE  (>={args.complete_qcov:.0f}% of the protein)          "
          f"{len(complete):>7,}  {pct(len(complete), tot):>7}")
    print("-" * 78)
    print(f"  genes matched ONLY by a paralog elsewhere   {len(paralog_only):>7,}")
    print("  those were counted as recovered by the old metric and are the")
    print("  reason actin read 99 recovered from 32 ORFs")

    # ---- the fragmentation measurement --------------------------------
    qcovs = [q for (_, q, _, _) in best_local.values()]
    ratios = []
    for g, (i, q, rp, op) in best_local.items():
        rl, ol = ref_len.get(rp), our_len.get(op)
        if rl and ol:
            ratios.append(100.0 * ol / rl)
    print("\n" + "=" * 78)
    print("HOW FRAGMENTED, EXACTLY")
    print("=" * 78)
    if qcovs:
        print(f"  genes with a locus-agreeing hit   {len(qcovs):,}")
        print(f"  median query coverage             "
              f"{statistics.median(qcovs):.0f}%")
        print(f"  mean query coverage               "
              f"{statistics.fmean(qcovs):.0f}%")
        hist(qcovs, [0, 10, 25, 50, 75, 90, 101],
             "coverage of the annotated protein by our best local match")
    if ratios:
        print(f"\n  median (our length / annotated length) "
              f"{statistics.median(ratios):.0f}%")

    # ---- is it length-dependent? --------------------------------------
    bands = [(0, 150), (150, 300), (300, 600), (600, 1200), (1200, 10 ** 9)]
    print("\n" + "=" * 78)
    print("COMPLETENESS BY ANNOTATED PROTEIN LENGTH")
    print("=" * 78)
    print(f"  {'protein length':<18}{'genes':>8}{'detected':>11}"
          f"{'complete':>11}{'med qcov':>10}")
    print("-" * 78)
    glen: dict[str, int] = {}
    for g, ps in gene_prots.items():
        L = [ref_len[p] for p in ps if p in ref_len]
        if L:
            glen[g] = max(L)
    by_band = []
    for lo, hi in bands:
        gs = {g for g, L in glen.items() if lo <= L < hi} & coding
        if not gs:
            continue
        d = len(gs & detected)
        c = len(gs & complete)
        qs = [best_local[g][1] for g in gs if g in best_local]
        mq = statistics.median(qs) if qs else 0
        tag = f"{lo}-{hi} aa" if hi < 10 ** 9 else f"{lo}+ aa"
        print(f"  {tag:<18}{len(gs):>8,}{pct(d, len(gs)):>11}"
              f"{pct(c, len(gs)):>11}{mq:>9.0f}%")
        by_band.append({"band": tag, "genes": len(gs), "detected": d,
                        "complete": c, "median_qcov": round(mq, 1)})
    print("-" * 78)
    print("  if completeness falls monotonically with length, the assembly is")
    print("  fragmented and -G would bridge it; if it is flat, it is depth")

    # ---- how many pieces per gene? ------------------------------------
    pieces: dict[str, int] = defaultdict(int)
    for pid, gs in prot_genes.items():
        for g in gs:
            pieces[g] += 1
    if pieces:
        multi = sum(1 for n in pieces.values() if n > 1)
        print(f"\n  our ORFs per annotated gene: median "
              f"{statistics.median(pieces.values()):.0f}, "
              f"{multi:,} genes carry more than one "
              f"({pct(multi, len(pieces))})")
        print("  several partial ORFs on one gene is fragmentation showing")
        print("  up in the protein set rather than the transcript set")

    # ---- can we fix it cheaply? ---------------------------------------
    print("\n" + "=" * 78)
    print("COST OF THE FIX  (-G without -e: guided, novel loci still reported)")
    print("=" * 78)
    bams = sorted(cfg.p("align").glob("**/*.bam"))
    tot_gb = sum(b.stat().st_size for b in bams) / 1e9 if bams else 0.0
    if bams:
        print(f"  {len(bams)} BAM(s) still on disk, {tot_gb:.1f} GB")
        for b in bams[:14]:
            print(f"    {b.parent.name}/{b.name}  "
                  f"{b.stat().st_size/1e9:.1f} GB")
        print("  -> StringTie can be re-run with -G directly. Minutes, not hours.")
    else:
        print("  no BAMs on disk — 21_align_assemble.py cleaned them up.")
        print("  -> the fix needs re-alignment of all 12 libraries:")
        print("     fastp + HISAT2 + StringTie -G, roughly the original cost.")
        print("     Reads come from ENA again; nothing in the archive is touched.")

    report = {
        "coding_genes": tot,
        "orf_on_exons": len(g_orf),
        "detected": len(detected),
        "complete": len(complete),
        "paralog_only_old_metric": len(paralog_only),
        "thresholds": {"identity": args.min_identity,
                       "detect_qcov": args.detect_qcov,
                       "complete_qcov": args.complete_qcov},
        "median_qcov": round(statistics.median(qcovs), 1) if qcovs else None,
        "median_length_ratio": round(statistics.median(ratios), 1) if ratios else None,
        "by_length_band": by_band,
        "genes_with_multiple_orfs": sum(1 for n in pieces.values() if n > 1),
        "bams_present": len(bams), "bam_gb": round(tot_gb, 1),
    }
    rp = res / "fragmentation.json"
    rp.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {rp}")

    # Genes we detected but did not complete: the recoverable set.
    frag = sorted(detected - complete)
    f = res / "genes_detected_not_complete.txt"
    f.write_text("\n".join(f"{g}\t{prod.get(g, '')}" for g in frag) + "\n")
    print(f"wrote {f.name} ({len(frag):,} genes we saw but did not finish)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
