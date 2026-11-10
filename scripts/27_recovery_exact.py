#!/usr/bin/env python3
"""Recovery measured exactly: exonic overlap and sequence identity, not spans.

    python scripts/27_recovery_exact.py --threads 16

WHY 26's NUMBERS CANNOT BE USED AS THEY STAND

  Two defects, both mine, both visible in its own output.

  1. THE SEQUENCE CHECK NEVER RAN. The reference proteome was not in the
     work tier, so the "protein at the locus, no sequence match" row read
     0 (0.0%) — not because nothing failed, but because nothing was tested.
     Every gene with an overlapping ORF was credited as recovered by
     default. That row was vacuous and the 9,507 figure is unvalidated.
     This script stages the proteome (from the archive, or NCBI) and
     actually runs it.

  2. OVERLAP WAS MEASURED AGAINST GENE SPANS, WHICH FAVOURS LONG GENES.
     A gene span includes its introns. An ORF sitting inside the intron of
     a 17 kb gene "overlaps" it while having nothing to do with it. The
     profile 26 printed is the signature of exactly that artifact:

         recovered               median locus 10,963 bp
         transcript, no protein  median locus 17,425 bp   <- longest
         no transcript            median locus  3,069 bp   <- shortest

     Recovered loci are 3.6x longer than untranscribed ones, and the class
     that supposedly has a transcript but no protein is longer still. Long
     loci are simply the easiest to hit by chance. Any recovery figure
     built on span overlap is inflated by an amount that grows with gene
     length, so it cannot be reported.

WHAT THIS MEASURES INSTEAD

  Three nested criteria, each strictly inside the last, so the drop at each
  tightening is itself the measurement:

    SPAN      an assembled transcript or predicted ORF overlaps the gene's
              span anywhere, introns included.        (26's criterion)
    EXONIC    the predicted ORF's CDS blocks overlap the gene's ANNOTATED
              EXONS by >= --min-bp bases, same strand. Intronic ORFs are
              excluded by construction.
    SEQUENCE  diamond blastp finds one of the gene's annotated proteins in
              our set at >= --min-identity over >= --min-qcov of its
              length. Reference as query, because the question is what we
              missed.

  SEQUENCE is the number that goes in the paper. The others exist to show
  how much the looser criteria overstate it.

ON THE UNSTRANDED HALF OF THE ASSEMBLY

  29,146 of 59,795 transcripts (49%) carry no strand: StringTie cannot
  orient a single-exon transcript. For transcript-level overlap a '.'
  strand is therefore allowed to match either orientation — rejecting them
  would discard half the assembly. Predicted ORFs always carry a definite
  strand from TransDecoder, so the ORF-level and sequence-level figures are
  unaffected by this.

NOTHING HERE CHANGES ANY OUTPUT. It only measures.
"""

from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
import os
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

BIN = 100_000
NCBI = "https://ftp.ncbi.nlm.nih.gov/genomes/all"


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def _attr(s: str, key: str) -> str | None:
    m = re.search(rf"(?:^|;){key}=([^;]+)", s)
    return m.group(1) if m else None


def find_diamond(explicit: str | None = None) -> str | None:
    for cand in (explicit, os.environ.get("DIAMOND"), shutil.which("diamond")):
        if cand and Path(cand).exists():
            return cand
    home = Path.home()
    for env in ("busco", "pmbio"):
        p = home / "pmdark_work" / "envs" / env / "bin" / "diamond"
        if p.exists():
            return str(p)
    return None


def gzip_ok(p: Path) -> bool:
    """A truncated download is still a file. UniProt taught us this."""
    try:
        with gzip.open(p, "rb") as fh:
            while fh.read(1 << 20):
                pass
        return True
    except (OSError, EOFError):
        return False


def stage_proteome(cfg, asm_prefix: str) -> Path | None:
    """Put <asm>_protein.faa in the work tier: archive first, then NCBI."""
    gdir = cfg.p("genome")
    plain = gdir / f"{asm_prefix}_protein.faa"
    if plain.exists() and plain.stat().st_size > 0:
        log(f"reference proteome already staged: {plain.name}")
        return plain

    # 1. anywhere in the archive
    for root in (cfg.p("raw"), cfg.p("interim"), gdir):
        try:
            cands = list(root.glob("**/*_protein.faa*"))
        except OSError:
            cands = []
        for c in cands:
            if c.suffix == ".gz":
                if not gzip_ok(c):
                    log(f"  {c} is a truncated gzip; ignoring")
                    continue
                log(f"decompressing {c}")
                with gzip.open(c, "rt") as i, plain.open("w") as o:
                    shutil.copyfileobj(i, o)
            else:
                log(f"copying {c}")
                shutil.copy2(c, plain)
            return plain

    # 2. NCBI, built from the accession in the assembly prefix
    m = re.match(r"(GC[AF])_(\d{3})(\d{3})(\d{3})\.(\d+)_", asm_prefix + "_")
    if not m:
        log(f"cannot build an NCBI path from '{asm_prefix}'")
        return None
    url = (f"{NCBI}/{m.group(1)}/{m.group(2)}/{m.group(3)}/{m.group(4)}/"
           f"{asm_prefix}/{asm_prefix}_protein.faa.gz")
    gz = gdir / f"{asm_prefix}_protein.faa.gz"
    log(f"fetching {url}")
    rc = subprocess.run([
        "curl", "-fSL", "--http1.1", "--retry", "5", "--retry-delay", "3",
        "--speed-limit", "1000", "--speed-time", "60",
        "-o", str(gz), url,
    ]).returncode
    if rc != 0 or not gz.exists():
        log(f"download failed (curl {rc})")
        return None
    if not gzip_ok(gz):
        log("downloaded file is a truncated gzip; discarding")
        gz.unlink(missing_ok=True)
        return None
    with gzip.open(gz, "rt") as i, plain.open("w") as o:
        shutil.copyfileobj(i, o)
    log(f"staged {plain.name} ({plain.stat().st_size/1e6:.1f} MB)")
    return plain


# ------------------------------------------------------- annotation model
def annotation(gff: Path) -> dict:
    """Coding genes with exon blocks, spans, and protein->gene mapping."""
    span: dict[str, tuple] = {}
    coding: set[str] = set()
    desc: dict[str, str] = {}
    rna2gene: dict[str, str] = {}
    gene_exons: dict[str, list] = defaultdict(list)
    prot2gene: dict[str, str] = {}
    pending_exons: list[tuple[str, int, int]] = []

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
                    span[gid] = (p[0], int(p[3]), int(p[4]), p[6])
                    d = _attr(a, "description") or _attr(a, "Name") or ""
                    desc[gid] = d.replace("%2C", ",").replace("%3B", ";")
            elif kind in ("mRNA", "transcript"):
                rid, par = _attr(a, "ID"), _attr(a, "Parent")
                if rid and par:
                    rna2gene[rid] = par
            elif kind == "exon":
                par = _attr(a, "Parent")
                if par:
                    pending_exons.append((par, int(p[3]), int(p[4])))
            elif kind == "CDS":
                pid, par = _attr(a, "protein_id"), _attr(a, "Parent")
                if pid and par:
                    g = rna2gene.get(par)
                    if g:
                        prot2gene[pid] = g

    for par, s, e in pending_exons:
        g = rna2gene.get(par)
        if g in coding:
            gene_exons[g].append((s, e))

    # Merge each gene's exons so overlap is counted once across isoforms.
    merged: dict[str, list] = {}
    exon_idx: dict[str, dict[int, list]] = defaultdict(lambda: defaultdict(list))
    for g, blocks in gene_exons.items():
        blocks.sort()
        out = [list(blocks[0])]
        for s, e in blocks[1:]:
            if s <= out[-1][1] + 1:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        merged[g] = [tuple(b) for b in out]
        chrom, _, _, strand = span[g]
        for s, e in merged[g]:
            for b in range(s // BIN, e // BIN + 1):
                exon_idx[chrom][b].append((s, e, strand, g))

    span_idx: dict[str, dict[int, list]] = defaultdict(lambda: defaultdict(list))
    for g in coding:
        chrom, s, e, strand = span[g]
        for b in range(s // BIN, e // BIN + 1):
            span_idx[chrom][b].append((s, e, strand, g))

    return {"coding": coding, "span": span, "desc": desc,
            "exons": merged, "exon_idx": exon_idx, "span_idx": span_idx,
            "prot2gene": prot2gene,
            "n_exons": {g: len(v) for g, v in merged.items()}}


def span_hits(idx, chrom, start, end, strand) -> set[str]:
    out = set()
    bins = idx.get(chrom)
    if not bins:
        return out
    for b in range(start // BIN, end // BIN + 1):
        for (s, e, st, g) in bins.get(b, ()):
            if e < start or s > end:
                continue
            if strand == "." or st == "." or st == strand:
                out.add(g)
    return out


def exonic_overlap(idx, chrom, strand, blocks, min_bp: int) -> set[str]:
    """Genes whose annotated exons share >= min_bp bases with these blocks."""
    acc: dict[str, int] = defaultdict(int)
    bins = idx.get(chrom)
    if not bins:
        return set()
    for (qs, qe) in blocks:
        for b in range(qs // BIN, qe // BIN + 1):
            for (s, e, st, g) in bins.get(b, ()):
                if st != strand:
                    continue
                lo, hi = max(qs, s), min(qe, e)
                if lo <= hi:
                    acc[g] += hi - lo + 1
    return {g for g, n in acc.items() if n >= min_bp}


def predicted_cds_blocks(gff3: Path) -> dict[str, tuple]:
    """ORF id -> (chrom, strand, [(s,e), ...]) from TransDecoder's GFF3."""
    rec: dict[str, list] = {}
    with gff3.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS":
                continue
            ident = _attr(p[8], "Parent")
            if not ident:
                continue
            if "~~" in ident:
                ident = ident.split("~~")[-1]
            ident = ident.split(",")[0]
            r = rec.setdefault(ident, [p[0], p[6], []])
            r[2].append((int(p[3]), int(p[4])))
    return {k: (v[0], v[1], sorted(v[2])) for k, v in rec.items()}


def gtf_exon_blocks(gtf: Path) -> dict[str, tuple]:
    """transcript_id -> (chrom, strand, [(s,e), ...])"""
    rec: dict[str, list] = {}
    with gtf.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "exon":
                continue
            m = re.search(r'transcript_id "([^"]+)"', p[8])
            if not m:
                continue
            r = rec.setdefault(m.group(1), [p[0], p[6], []])
            r[2].append((int(p[3]), int(p[4])))
    return {k: (v[0], v[1], sorted(v[2])) for k, v in rec.items()}


# ------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--min-bp", type=int, default=30,
                    help="exonic bases an ORF must share with a gene")
    ap.add_argument("--min-identity", type=float, default=70.0)
    ap.add_argument("--min-qcov", type=float, default=50.0)
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--diamond", default=None)
    ap.add_argument("--force-search", action="store_true")
    args = ap.parse_args()

    cfg = cfgmod.load()
    out, res = cfg.p("proteins"), cfg.p("results")
    qc = json.loads((res / "protein_qc.json").read_text())
    pep, gff3 = Path(qc["pep"]), Path(qc["genome_gff3"])

    genome = next(iter(p for p in cfg.p("genome").glob("*_genomic.fna")
                       if "_from_genomic" not in p.name), None)
    ref_gff = next(iter(cfg.p("genome").glob("*_genomic.gff")), None)
    if ref_gff is None or genome is None:
        raise SystemExit("reference GFF / genome FASTA not found")
    asm = genome.name.replace("_genomic.fna", "")
    log(f"assembly prefix: {asm}")

    log("building the annotation model (exon blocks per coding gene)")
    A = annotation(ref_gff)
    coding = A["coding"]
    log(f"  {len(coding):,} coding genes, "
        f"{sum(len(v) for v in A['exons'].values()):,} merged exons, "
        f"{len(A['prot2gene']):,} proteins mapped")

    # ---- transcript level ---------------------------------------------
    filt = out / "merged_filtered.gtf"
    log(f"exonic overlap: {filt.name}")
    tx = gtf_exon_blocks(filt)
    g_tx_span: set[str] = set()
    g_tx_exon: set[str] = set()
    for chrom, strand, blocks in tx.values():
        lo = blocks[0][0]
        hi = max(e for _, e in blocks)
        g_tx_span |= span_hits(A["span_idx"], chrom, lo, hi, strand)
        if strand == ".":
            for st in ("+", "-"):
                g_tx_exon |= exonic_overlap(A["exon_idx"], chrom, st,
                                            blocks, args.min_bp)
        else:
            g_tx_exon |= exonic_overlap(A["exon_idx"], chrom, strand,
                                        blocks, args.min_bp)
    g_tx_span &= coding
    g_tx_exon &= coding
    log(f"  span {len(g_tx_span):,}   exonic {len(g_tx_exon):,}")

    # ---- ORF level -----------------------------------------------------
    log("exonic overlap: predicted ORFs")
    cds = predicted_cds_blocks(gff3)
    lens = _qc.read_fasta_lengths(pep)
    g_orf_span: set[str] = set()
    g_orf_exon: set[str] = set()
    novel_exon: set[str] = set()
    for pid in lens:
        r = cds.get(pid)
        if not r:
            continue
        chrom, strand, blocks = r
        lo = blocks[0][0]
        hi = max(e for _, e in blocks)
        s_hits = span_hits(A["span_idx"], chrom, lo, hi, strand) & coding
        e_hits = exonic_overlap(A["exon_idx"], chrom, strand, blocks,
                                args.min_bp) & coding
        g_orf_span |= s_hits
        g_orf_exon |= e_hits
        if not e_hits:
            novel_exon.add(pid)
    log(f"  span {len(g_orf_span):,}   exonic {len(g_orf_exon):,}")
    log(f"  proteins with NO exonic overlap of any coding gene: "
        f"{len(novel_exon):,}")

    # ---- sequence level ------------------------------------------------
    g_seq: set[str] = set()
    ran = False
    refp = stage_proteome(cfg, asm)
    dm = find_diamond(args.diamond)
    if refp is None:
        log("no reference proteome — the sequence check cannot run")
    elif dm is None:
        log("diamond not found — the sequence check cannot run")
    else:
        db = out / "ourprots"
        tsv = out / "ref_vs_ours.tsv"
        if args.force_search or not tsv.exists():
            log("diamond makedb on our predicted proteins")
            subprocess.run([dm, "makedb", "--in", str(pep), "-d", str(db),
                            "--threads", str(args.threads)],
                           capture_output=True)
            log("diamond blastp: annotated proteome -> ours "
                "(reference as query: the question is what we missed)")
            rc = subprocess.run([
                dm, "blastp", "--query", str(refp), "--db", str(db),
                "--out", str(tsv), "--outfmt", "6", "qseqid", "sseqid",
                "pident", "qcovhsp", "evalue",
                "--max-target-seqs", "1", "--evalue", str(args.evalue),
                "--threads", str(args.threads), "--quiet",
            ]).returncode
            if rc != 0:
                log(f"diamond blastp exited {rc}")
        if tsv.exists():
            ran = True
            nhit = npass = 0
            with tsv.open() as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 4:
                        continue
                    try:
                        ident, qcov = float(f[2]), float(f[3])
                    except ValueError:
                        continue
                    nhit += 1
                    if ident >= args.min_identity and qcov >= args.min_qcov:
                        npass += 1
                        g = A["prot2gene"].get(f[0])
                        if g:
                            g_seq.add(g)
            g_seq &= coding
            log(f"  {nhit:,} hits, {npass:,} passing, "
                f"{len(g_seq):,} genes recovered by sequence")

    # ---- report --------------------------------------------------------
    tot = len(coding)
    print("\n" + "=" * 78)
    print("RECOVERY OF ANNOTATED CODING GENES, BY CRITERION")
    print("=" * 78)
    print(f"  {'criterion':<46}{'genes':>8}{'of':>9}")
    print("-" * 78)
    rows = [("transcript overlaps the gene SPAN", len(g_tx_span)),
            ("transcript overlaps annotated EXONS", len(g_tx_exon)),
            ("predicted ORF overlaps the gene SPAN", len(g_orf_span)),
            ("predicted ORF overlaps annotated EXONS", len(g_orf_exon))]
    if ran:
        rows.append(("annotated protein found in our set (SEQUENCE)",
                     len(g_seq)))
    for lab, n in rows:
        print(f"  {lab:<46}{n:>8,}{100*n/tot:>8.1f}%")
    print("-" * 78)
    infl = len(g_orf_span) - len(g_orf_exon)
    print(f"  span overstates ORF-level recovery by {infl:,} genes "
          f"({100*infl/max(1,len(g_orf_span)):.1f}%)")
    print("  those are ORFs inside the introns of long genes")
    if ran:
        print(f"\n  REPORTABLE RECOVERY: {100*len(g_seq)/tot:.1f}% "
              f"({len(g_seq):,}/{tot:,})")

    # ---- where the rest went, on the exonic criterion ------------------
    cells: dict[str, list[str]] = defaultdict(list)
    for g in coding:
        t, o = g in g_tx_exon, g in g_orf_exon
        s = g in g_seq
        if ran and s:
            cells["recovered"].append(g)
        elif ran and o:
            cells["orf_no_sequence"].append(g)
        elif not ran and o:
            cells["recovered_unverified"].append(g)
        elif t:
            cells["transcript_no_orf"].append(g)
        else:
            cells["no_transcript"].append(g)

    print("\n" + "=" * 78)
    print("WHERE THE GENES WENT  (exonic criterion, mutually exclusive)")
    print("=" * 78)
    for k, lab in (("recovered", "recovered, sequence-confirmed"),
                   ("recovered_unverified", "ORF on the exons (unverified)"),
                   ("orf_no_sequence", "ORF on the exons, sequence FAILED"),
                   ("transcript_no_orf", "exons transcribed, no ORF called"),
                   ("no_transcript", "exons never transcribed")):
        n = len(cells[k])
        if n:
            print(f"  {lab:<46}{n:>8,}{100*n/tot:>8.1f}%")
    print("-" * 78)
    print(f"  {'total':<46}{sum(len(v) for v in cells.values()):>8,}")

    def prof(gs: list[str]) -> str:
        if not gs:
            return "n/a"
        L = sorted(A["span"][g][2] - A["span"][g][1] for g in gs)
        ex = [A["n_exons"].get(g, 0) for g in gs]
        return (f"median locus {L[len(L)//2]:,} bp, "
                f"median exons {sorted(ex)[len(ex)//2]}, "
                f"single-exon {100*sum(1 for e in ex if e==1)/len(ex):.0f}%")

    print("\n" + "=" * 78)
    print("PROFILE")
    print("=" * 78)
    for k, lab in (("recovered", "recovered"),
                   ("recovered_unverified", "ORF on exons (unverified)"),
                   ("orf_no_sequence", "sequence failed"),
                   ("transcript_no_orf", "no ORF called"),
                   ("no_transcript", "never transcribed")):
        if cells[k]:
            print(f"  {lab:<30}{prof(cells[k])}")

    # What ARE the never-transcribed genes? Tally their descriptions.
    miss = cells["no_transcript"]
    if miss:
        words: dict[str, int] = defaultdict(int)
        for g in miss:
            for w in re.findall(r"[A-Za-z]{4,}", A["desc"].get(g, "").lower()):
                if w not in ("like", "protein", "containing", "domain",
                             "family", "putative", "type", "subunit"):
                    words[w] += 1
        print("\n" + "=" * 78)
        print(f"THE {len(miss):,} NEVER-TRANSCRIBED GENES — commonest "
              "description words")
        print("=" * 78)
        top = sorted(words.items(), key=lambda kv: -kv[1])[:18]
        for i in range(0, len(top), 3):
            print("   " + "".join(f"{w:<16}{n:>6,}   " for w, n in top[i:i+3]))
        print("\n  if these are sensory, cuticular, reproductive or")
        print("  digestive families, their absence from hemocytes is biology,")
        print("  not a pipeline failure")

    report = {
        "coding_genes": tot,
        "min_bp_exonic": args.min_bp,
        "thresholds": {"identity": args.min_identity, "qcov": args.min_qcov},
        "transcript_span": len(g_tx_span),
        "transcript_exonic": len(g_tx_exon),
        "orf_span": len(g_orf_span),
        "orf_exonic": len(g_orf_exon),
        "sequence": len(g_seq) if ran else None,
        "sequence_check_ran": ran,
        "span_inflation_genes": infl,
        "crosstab": {k: len(v) for k, v in cells.items()},
        "proteins_no_exonic_overlap": len(novel_exon),
    }
    rp = res / "recovery_exact.json"
    rp.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {rp}")

    # The dark pool, now defined by exonic overlap rather than span.
    nf = out / "proteins_novel_exonic.faa"
    w = False
    with pep.open() as fh, nf.open("w") as o:
        for line in fh:
            if line.startswith(">"):
                w = line[1:].split()[0] in novel_exon
            if w:
                o.write(line)
    print(f"wrote {nf.name} ({len(novel_exon):,} sequences)")
    print("\n  NOTE this novel set is LARGER than 25_protein_qc's 1,648,")
    print("  because an ORF in the intron of an annotated gene now counts as")
    print("  novel rather than as that gene. Which definition the paper uses")
    print("  must be stated explicitly; intronic ORFs are a real class.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
