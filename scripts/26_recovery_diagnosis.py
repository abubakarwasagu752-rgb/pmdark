#!/usr/bin/env python3
"""Why only 39.6% of annotated genes? Locate the loss before building on it.

    python scripts/26_recovery_diagnosis.py --threads 16
    python scripts/26_recovery_diagnosis.py --no-diamond    # ~3 min, no search

THE PROBLEM

  scripts/25_protein_qc.py found that 18,892 predicted proteins sit on only
  9,507 distinct annotated protein-coding genes — 1.99 proteins per gene,
  and 39.6% of the 24,011 annotated coding genes. The pre-stated rule said
  60-75% would be unremarkable for a single tissue. 39.6% is not, so the
  dark proteome does not get defined until the missing 14,504 genes are
  accounted for. BUSCO would have answered the completeness half of this,
  but the busco install on this cluster is broken ("BUSCO must be installed
  before it is run"), so completeness is measured here a better way anyway:
  against the conspecific reference proteome rather than a cross-phylum
  marker set.

THREE CANDIDATE EXPLANATIONS, AND THE MEASUREMENT THAT SEPARATES THEM

  Every missing gene fell out at exactly one of three steps, and each step
  leaves a different fingerprint:

    1. NOT SEQUENCED / NOT EXPRESSED. Hemocytes are one tissue, and if a
       few highly expressed loci soaked up the reads, effective depth for
       everything else is low. Fingerprint: no assembled transcript
       overlaps the gene at all. Measured by locus-level overlap, plus the
       share of each library's TPM held by its top 10 / 100 / 1000
       transcripts — a concentrated library is a shallow library in
       disguise.

    2. ASSEMBLED BUT NO ORF CALLED. Fingerprint: a transcript overlaps the
       gene, but no predicted protein does. This is the dangerous category,
       because it is our own pipeline losing real genes, and it would
       silently deflate every denominator downstream.

    3. ASSEMBLED AND CALLED BUT DIVERGENT. Fingerprint: a protein overlaps
       the gene's locus, but the sequence does not match the annotated
       protein. Either mis-assembly or a genuinely different ORF.

  Reported as a cross-tabulation of (transcript? x protein? x sequence
  match?) over all 24,011 annotated coding genes, so each gene lands in
  exactly one cell and the cells sum to the total.

WHY THE SEQUENCE CHECK RUNS THE DIRECTION IT DOES

  diamond blastp with the ANNOTATED proteins as query and ours as database.
  The question is "what did we miss", so the reference has to be the query:
  every annotated protein gets the chance to find itself in our set, and a
  gene counts as recovered only if one of its proteins is found at
  >= --min-identity over >= --min-qcov of its length. Running it the other
  way would answer "is what we built real", which is a different question
  and the one 25_protein_qc.py already settled at 88.1%.

NOTHING HERE CHANGES ANY OUTPUT. It only measures.
"""

from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "qc", Path(__file__).resolve().parent / "25_protein_qc.py")
_qc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_qc)

BIN = _qc.BIN


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def _attr(s: str, key: str) -> str | None:
    m = re.search(rf"(?:^|;){key}=([^;]+)", s)
    return m.group(1) if m else None


def _gtf_attr(s: str, key: str) -> str | None:
    m = re.search(rf'{key} "([^"]+)"', s)
    return m.group(1) if m else None


def find_diamond(explicit: str | None = None) -> str | None:
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


# ------------------------------------------------------------- reference
def reference_detail(gff: Path) -> dict:
    """Coding-gene index, exon counts, and the protein_id -> gene map."""
    idx: dict[str, dict[int, list]] = defaultdict(lambda: defaultdict(list))
    coding: set[str] = set()
    span: dict[str, tuple] = {}
    rna2gene: dict[str, str] = {}
    exons: dict[str, int] = defaultdict(int)
    prot2gene: dict[str, str] = {}

    with gff.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            kind, a = p[2], p[8]
            if kind == "gene":
                gid = _attr(a, "ID") or f"{p[0]}:{p[3]}"
                if (_attr(a, "gene_biotype") or "") == "protein_coding":
                    coding.add(gid)
                    s, e = int(p[3]), int(p[4])
                    span[gid] = (p[0], s, e, p[6])
                    for b in range(s // BIN, e // BIN + 1):
                        idx[p[0]][b].append((s, e, p[6], gid))
            elif kind in ("mRNA", "transcript"):
                rid, par = _attr(a, "ID"), _attr(a, "Parent")
                if rid and par:
                    rna2gene[rid] = par
            elif kind == "exon":
                par = _attr(a, "Parent")
                if par:
                    exons[par] += 1
            elif kind == "CDS":
                pid, par = _attr(a, "protein_id"), _attr(a, "Parent")
                if pid and par and par in rna2gene:
                    prot2gene[pid] = rna2gene[par]

    gene_exons: dict[str, int] = defaultdict(int)
    for rid, n in exons.items():
        g = rna2gene.get(rid)
        if g:
            gene_exons[g] = max(gene_exons[g], n)

    return {"index": idx, "coding": coding, "span": span,
            "gene_exons": gene_exons, "prot2gene": prot2gene}


def sense_hits(idx: dict, chrom: str, start: int, end: int,
               strand: str) -> set[str]:
    """Coding genes overlapped. A '.' strand matches either orientation —
    StringTie leaves single-exon transcripts unstranded, and treating that
    as a mismatch would discard them all."""
    out: set[str] = set()
    bins = idx.get(chrom)
    if not bins:
        return out
    for b in range(start // BIN, end // BIN + 1):
        for (s, e, st, gid) in bins.get(b, ()):
            if e < start or s > end:
                continue
            if strand == "." or st == "." or st == strand:
                out.add(gid)
    return out


def genes_from_gtf(gtf: Path, idx: dict) -> tuple[set[str], int, int]:
    """Coding genes overlapped by assembled transcripts."""
    hit: set[str] = set()
    n = nostrand = 0
    with gtf.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "transcript":
                continue
            n += 1
            if p[6] == ".":
                nostrand += 1
            hit |= sense_hits(idx, p[0], int(p[3]), int(p[4]), p[6])
    return hit, n, nostrand


# ------------------------------------------------------- library skew
def library_concentration(gtfs: list[Path]) -> list[dict]:
    rows = []
    for g in gtfs:
        tpms = []
        with g.open() as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 9 or p[2] != "transcript":
                    continue
                v = _gtf_attr(p[8], "TPM")
                if v:
                    try:
                        tpms.append(float(v))
                    except ValueError:
                        pass
        if not tpms:
            continue
        tpms.sort(reverse=True)
        tot = sum(tpms) or 1.0
        rows.append({
            "library": g.parent.name,
            "n_transcripts": len(tpms),
            "top10_pct": round(100 * sum(tpms[:10]) / tot, 1),
            "top100_pct": round(100 * sum(tpms[:100]) / tot, 1),
            "top1000_pct": round(100 * sum(tpms[:1000]) / tot, 1),
        })
    return rows


# ------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--no-diamond", action="store_true")
    ap.add_argument("--diamond", default=None)
    ap.add_argument("--min-identity", type=float, default=70.0)
    ap.add_argument("--min-qcov", type=float, default=50.0)
    ap.add_argument("--evalue", type=float, default=1e-5)
    args = ap.parse_args()

    cfg = cfgmod.load()
    out, res = cfg.p("proteins"), cfg.p("results")
    qc = json.loads((res / "protein_qc.json").read_text())
    pep = Path(qc["pep"])
    gff3 = Path(qc["genome_gff3"])
    ref_gff = next(iter(cfg.p("genome").glob("*_genomic.gff")), None)
    if ref_gff is None:
        raise SystemExit("reference GFF not found")

    log("indexing the reference annotation (genes, exons, protein ids)")
    R = reference_detail(ref_gff)
    coding = R["coding"]
    log(f"  {len(coding):,} protein-coding genes, "
        f"{len(R['prot2gene']):,} annotated proteins mapped to genes")

    # ---- step 1: which coding genes did the ASSEMBLY touch at all? -----
    filt = out / "merged_filtered.gtf"
    log(f"overlapping {filt.name} against coding genes")
    g_tx, n_tx, n_nostrand = genes_from_gtf(filt, R["index"])
    g_tx &= coding
    log(f"  {n_tx:,} transcripts ({n_nostrand:,} unstranded) "
        f"touch {len(g_tx):,} coding genes")

    perm = out / "merged_permissive.gtf"
    g_perm = set()
    if perm.exists():
        log(f"overlapping {perm.name} (pre-filter) for comparison")
        g_perm, n_perm, _ = genes_from_gtf(perm, R["index"])
        g_perm &= coding
        log(f"  {n_perm:,} transcripts touch {len(g_perm):,} coding genes")

    # ---- step 2: which coding genes carry a predicted PROTEIN? ---------
    log("placing predicted proteins")
    lens = _qc.read_fasta_lengths(pep)
    spans = _qc.predicted_spans(gff3)
    g_orf: set[str] = set()
    for pid in lens:
        sp = spans.get(pid)
        if sp:
            g_orf |= sense_hits(R["index"], *sp)
    g_orf &= coding
    log(f"  {len(lens):,} proteins touch {len(g_orf):,} coding genes")

    # ---- step 3: which annotated proteins are RECOVERED by sequence? ---
    g_seq: set[str] = set()
    dmnd_ran = False
    if not args.no_diamond:
        dm = find_diamond(args.diamond)
        refp = next((p for p in cfg.p("genome").glob("*_protein.faa*")), None)
        if dm is None:
            log("diamond not found — skipping the sequence check")
        elif refp is None:
            log("reference proteome (*_protein.faa) not in the work tier — "
                "skipping the sequence check")
        else:
            qfa = out / "ref_proteome.faa"
            if not qfa.exists():
                log(f"staging {refp.name}")
                op = gzip.open if refp.suffix == ".gz" else open
                with op(refp, "rt") as i, qfa.open("w") as o:
                    shutil.copyfileobj(i, o)
            db = out / "ourprots"
            tsv = out / "ref_vs_ours.tsv"
            if not tsv.exists():
                log("diamond makedb on our predicted proteins")
                subprocess.run([dm, "makedb", "--in", str(pep), "-d", str(db),
                                "--threads", str(args.threads)],
                               capture_output=True)
                log("diamond blastp: annotated proteome -> our proteins")
                rc = subprocess.run([
                    dm, "blastp", "--query", str(qfa), "--db", str(db),
                    "--out", str(tsv), "--outfmt", "6", "qseqid", "sseqid",
                    "pident", "qcovhsp", "evalue",
                    "--max-target-seqs", "1", "--evalue", str(args.evalue),
                    "--threads", str(args.threads), "--quiet",
                ]).returncode
                if rc != 0:
                    log(f"diamond blastp exited {rc}")
            if tsv.exists():
                dmnd_ran = True
                n = 0
                with tsv.open() as fh:
                    for line in fh:
                        f = line.rstrip("\n").split("\t")
                        if len(f) < 4:
                            continue
                        try:
                            pid_, ident, qcov = f[0], float(f[2]), float(f[3])
                        except ValueError:
                            continue
                        n += 1
                        if ident >= args.min_identity and qcov >= args.min_qcov:
                            g = R["prot2gene"].get(pid_)
                            if g:
                                g_seq.add(g)
                g_seq &= coding
                log(f"  {n:,} hits; {len(g_seq):,} coding genes recovered "
                    f"at >={args.min_identity}% id / >={args.min_qcov}% qcov")

    # ---- cross-tabulate ------------------------------------------------
    cells: dict[str, list[str]] = defaultdict(list)
    for g in coding:
        t, o = g in g_tx, g in g_orf
        s = g in g_seq
        if not t and not o:
            cells["no_transcript"].append(g)
        elif t and not o:
            cells["transcript_no_orf"].append(g)
        elif o and (s or not dmnd_ran):
            cells["recovered"].append(g)
        else:
            cells["orf_no_sequence_match"].append(g)

    tot = len(coding)
    print("\n" + "=" * 78)
    print("WHERE THE 24,011 ANNOTATED CODING GENES WENT")
    print("=" * 78)
    rows = [
        ("recovered", "recovered (protein called"
                      + (" + sequence match)" if dmnd_ran else ")"),
         "this is the real recovery figure"),
        ("orf_no_sequence_match", "protein at the locus, no sequence match",
         "mis-assembly or a divergent ORF"),
        ("transcript_no_orf", "transcript assembled, NO protein called",
         "OUR pipeline losing real genes"),
        ("no_transcript", "no transcript at all",
         "not expressed in hemocytes, or too shallow"),
    ]
    for key, lab, why in rows:
        n = len(cells[key])
        print(f"  {lab:<44} {n:>7,}  ({100*n/tot:>4.1f}%)")
        print(f"    {why}")
    print("-" * 78)
    print(f"  {'total':<44} {sum(len(v) for v in cells.values()):>7,}")

    if g_perm:
        lost = len(g_perm - g_tx)
        print(f"\n  coding genes present before the evidence filter but not "
              f"after: {lost:,}")
        print("  (if this is large, 23c_filter_transcripts.py is too strict)")

    # ---- what do the missed genes look like? ---------------------------
    def profile(gs: list[str]) -> str:
        if not gs:
            return "n/a"
        L = sorted((R["span"][g][2] - R["span"][g][1]) for g in gs
                   if g in R["span"])
        ex = [R["gene_exons"].get(g, 0) for g in gs]
        single = sum(1 for e in ex if e == 1)
        med = L[len(L) // 2] if L else 0
        return (f"median locus {med:,} bp, "
                f"single-exon {100*single/max(1,len(ex)):.0f}%")

    print("\n" + "=" * 78)
    print("PROFILE OF THE MISSED GENES  (loss mechanism leaves a fingerprint)")
    print("=" * 78)
    for key, lab, _ in rows:
        print(f"  {lab:<44} {profile(cells[key])}")

    # ---- library skew --------------------------------------------------
    lib_gtfs = sorted(p for p in cfg.p("align").glob("*/[SE]RR*.gtf"))
    skew = library_concentration(lib_gtfs)
    if skew:
        print("\n" + "=" * 78)
        print("LIBRARY CONCENTRATION  (share of TPM held by the top transcripts)")
        print("=" * 78)
        print(f"  {'library':<16}{'transcripts':>12}{'top10':>9}"
              f"{'top100':>9}{'top1000':>9}")
        print("-" * 78)
        for r in sorted(skew, key=lambda x: -x["top10_pct"]):
            print(f"  {r['library']:<16}{r['n_transcripts']:>12,}"
                  f"{r['top10_pct']:>8.1f}%{r['top100_pct']:>8.1f}%"
                  f"{r['top1000_pct']:>8.1f}%")
        med10 = sorted(r["top10_pct"] for r in skew)[len(skew) // 2]
        print("-" * 78)
        print(f"  median top-10 share {med10:.1f}% — above ~40% means a few")
        print("  loci (rRNA, mito, haemocyanin) are consuming the depth and")
        print("  effective coverage of everything else is far below nominal")

    report = {
        "reference_coding_genes": tot,
        "genes_touched_by_transcripts": len(g_tx),
        "genes_touched_by_transcripts_prefilter": len(g_perm) or None,
        "genes_with_protein": len(g_orf),
        "genes_recovered_by_sequence": len(g_seq) if dmnd_ran else None,
        "sequence_check_ran": dmnd_ran,
        "thresholds": {"identity": args.min_identity, "qcov": args.min_qcov},
        "crosstab": {k: len(v) for k, v in cells.items()},
        "unstranded_transcripts": n_nostrand,
        "library_concentration": skew,
    }
    rp = res / "recovery_diagnosis.json"
    rp.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {rp}")

    for key in ("transcript_no_orf", "no_transcript"):
        f = res / f"genes_{key}.txt"
        f.write_text("\n".join(sorted(cells[key])) + "\n")
        print(f"wrote {f.name} ({len(cells[key]):,} gene ids)")

    print("\nRead the cross-tabulation first. It names the step to fix, and")
    print("the fix differs completely between the four rows.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
