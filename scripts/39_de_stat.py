#!/usr/bin/env python3
"""STAT-dependence of the reference transcriptome, from the one clean contrast.

    python scripts/39_de_stat.py
    python scripts/39_de_stat.py --read-length 100
    python scripts/39_de_stat.py --inventory      # report inputs, compute nothing

WHY THIS CONTRAST AND NO OTHER

  docs/DESIGN_LIMITS.md settles it: there is no naive arm and no
  non-targeting dsRNA arm, so every contrast except one confounds infection
  with knockdown. The only single-variable axis in the dataset is

      WSSV (n=5)   vs   PmSTAT_dsRNA_WSSV (n=4)

  which isolates STAT knockdown UNDER infection. That is also the axis the
  paper's title claims, so the design limitation and the central claim
  happen to coincide.

  The PmSTAT_dsRNA-only arm is n=2 after SRR28572179 is excluded as a
  failed library (830k reads against 10-16M elsewhere). It is not used for
  the primary contrast.

DIRECTION, STATED ONCE

  Positive log2FoldChange = higher in WSSV than in knockdown+WSSV
                          = suppressed when STAT is knocked down
                          = STAT-DEPENDENT

  pydeseq2 fixes the reference level alphabetically and ignores ref_level,
  so the coefficient is always group[T.WSSV]. This script DISCOVERS the
  coefficient from the design matrix and derives the contrast direction
  from its name rather than hardcoding either, then asserts that the
  shrunk fold changes correlate POSITIVELY with the unshrunk ones.
  Without that guard a mismatched coeff silently inverts every fold change
  and the STAT-dependent set comes back as its own complement.

THE POSITIVE CONTROL THAT DECIDES EVERYTHING

  PmSTAT's own transcript must be down in the knockdown arm. If it is not,
  the knockdown failed and nothing downstream means anything. This is
  checked and printed before any other result, and the script says plainly
  when it cannot find the gene.

WHERE THE COUNTS COME FROM, AND WHAT THAT COSTS

  The BAMs are gone -- 21_align_assemble.py only keeps them with
  --keep-bam, which is off by default -- so counts cannot come from
  featureCounts without re-aligning all 12 libraries.

  Instead they are derived from the per-library {acc}.ref_expr.gtf written
  by the -e -G quantification pass, as StringTie's own prepDE.py does:

      count = coverage * transcript_length / read_length

  summed over a gene's transcripts. These are ESTIMATED counts, not
  observed read counts, and that is a real limitation to state in the
  methods. --read-length enters as a global constant, identical for every
  gene and library, so it cancels in size-factor normalisation; it shifts
  absolute count magnitude only.

  Only reference genes are quantified here, because -e -G restricts output
  to the annotation. The 2,108 novel-locus proteins are NOT in this matrix
  and need their own quantification against the merged assembly. The
  "dark" fraction tested here is the uncharacterised (LOC-only) share of
  the annotated transcriptome.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod                              # noqa: E402

# NOTE: an earlier version classified genes by whether Name= matched
# ^LOC\d+$. That is dead and must not come back: in this annotation
# every protein-coding gene is LOC-named, so the test was vacuous.
# Characterisation comes from product=; see gff_annotation().

# Two patterns on purpose. The first version used r"\bstat\b", which fails
# on every real STAT symbol -- "stat5b", "STAT1", "PmSTAT" -- because \b
# needs a non-word character after "stat" and a digit is a word character.
# That would have reported the knockdown as unverifiable while the gene sat
# in the table. A false negative here silently invalidates the whole
# analysis, so STRICT is permissive about prefixes and suffixes, and LOOSE
# prints anything else name-matching for a human to eyeball.
#
# STRICT anchors on the END of the name, which is what lets "PmSTAT" and
# "stat5b" through while keeping "stathmin", "statherin" and "prostatin"
# out. "STATH" (statherin's symbol) does slip through; it is printed with
# its name and fold change, so a human sees immediately that it is not the
# control.
STAT_STRICT = re.compile(
    r"stat\d*[a-z]?$|signal transducer and activator", re.I)
UNCHAR = re.compile(r"uncharacterized|uncharacterised|hypothetical", re.I)
STAT_LOOSE = re.compile(r"stat|signal transducer and activator", re.I)
ATTR = re.compile(r'(\S+) "([^"]*)"')


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


# ---------------------------------------------------------------------------
# sample sheet
# ---------------------------------------------------------------------------

def read_samples(path: Path) -> list[dict]:
    """Parse samples.tsv BY HEADER NAME.

    scripts/30 had a bug here: it read `use` positionally as the last
    column when `note` is last, so every library's inclusion flag came
    from the wrong field. Columns are addressed by name only.
    """
    rows: list[dict] = []
    header: list[str] | None = None
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if header is None:
            header = [h.strip() for h in f]
            for need in ("run_accession", "group", "use"):
                if need not in header:
                    raise SystemExit(
                        f"{path.name}: no '{need}' column; found {header}")
            continue
        rows.append(dict(zip(header, f)))
    if not rows:
        raise SystemExit(f"{path.name}: no data rows")
    return rows


# ---------------------------------------------------------------------------
# annotation: the characterisation lives in product=, not Name=
# ---------------------------------------------------------------------------

GFF_ATTR = re.compile(r"([^=;]+)=([^;]*)")


def is_uncharacterised(product: str) -> bool:
    return bool(UNCHAR.search(product or ""))


def gff_annotation(gff: Path) -> dict[str, dict]:
    """gene_id -> {name, product, biotype} from the reference GFF.

    WHY NOT Name=, WHICH IS WHAT THIS SCRIPT USED FIRST

      In GCF_015228065.2 only 1,929 of 28,261 genes carry a symbol, and
      those turn out to be tRNAs (Trnar-acg, Trnam-cau, ...). Every
      protein-coding gene is named LOC<digits> whether it is well
      understood or completely unknown -- there is no actb, no tubb. So
      Name= cannot separate characterised from uncharacterised, and the
      first version of this script produced a 13-gene contingency table
      of tRNAs and reported a p-value off it.

      The characterisation is in product=, on the mRNA features: NCBI
      gives a LOC gene either "uncharacterized protein LOC119587821" or a
      real description like "signal transducer and activator of
      transcription 5B-like". That split is 8,780 / 24,023 here.

      A gene with ANY described isoform counts as characterised, so a
      described product always wins over an uncharacterised one.
    """
    info: dict[str, dict] = {}
    with gff.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            feat, attrs = f[2], f[8]
            if feat in ("gene", "pseudogene"):
                a = dict(GFF_ATTR.findall(attrs))
                gid = a.get("ID")
                if gid:
                    d = info.setdefault(gid, {})
                    d["name"] = a.get("Name", "")
                    d["biotype"] = a.get("gene_biotype", "")
            elif feat in ("mRNA", "transcript", "ncRNA", "tRNA", "rRNA"):
                a = dict(GFF_ATTR.findall(attrs))
                parent, prod = a.get("Parent", ""), a.get("product", "")
                # mRNA parents are gene-*; CDS parents are rna-*, skip those
                if not prod or not parent.startswith("gene-"):
                    continue
                d = info.setdefault(parent, {})
                cur = d.get("product", "")
                if not cur or (is_uncharacterised(cur)
                               and not is_uncharacterised(prod)):
                    d["product"] = prod
    return info


# ---------------------------------------------------------------------------
# counts from the -e -G quantification
# ---------------------------------------------------------------------------

def gene_counts(gtf: Path, read_length: float) -> tuple[dict, dict]:
    """(gene -> estimated count, gene -> name) from one ref_expr.gtf.

    prepDE.py's conversion: a transcript's coverage times its exonic
    length, divided by the read length, estimates how many reads produced
    it. Gene counts sum over transcripts.
    """
    tx_cov: dict[str, float] = {}
    tx_gene: dict[str, str] = {}
    tx_len: dict[str, int] = defaultdict(int)
    names: dict[str, str] = {}

    with gtf.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            feat, start, end, attrs = f[2], f[3], f[4], f[8]
            a = dict(ATTR.findall(attrs))
            tid = a.get("transcript_id")
            if not tid:
                continue
            gid = a.get("ref_gene_id") or a.get("gene_id") or tid
            tx_gene[tid] = gid
            nm = a.get("ref_gene_name") or a.get("gene_name")
            if nm:
                names.setdefault(gid, nm)
            if feat == "transcript":
                try:
                    tx_cov[tid] = float(a.get("cov", "0"))
                except ValueError:
                    tx_cov[tid] = 0.0
            elif feat == "exon":
                try:
                    tx_len[tid] += int(end) - int(start) + 1
                except ValueError:
                    pass

    counts: dict[str, float] = defaultdict(float)
    for tid, gid in tx_gene.items():
        cov, length = tx_cov.get(tid, 0.0), tx_len.get(tid, 0)
        if cov > 0 and length > 0:
            counts[gid] += cov * length / read_length
    return dict(counts), names


def build_matrix(align_dir: Path, samples: list[dict], read_length: float):
    """(libraries, genes, integer count matrix, gene names)."""
    per_lib: dict[str, dict] = {}
    names: dict[str, str] = {}
    missing: list[str] = []
    for s in samples:
        acc = s["run_accession"]
        gtf = align_dir / acc / f"{acc}.ref_expr.gtf"
        if not gtf.exists():
            alt = sorted(align_dir.glob(f"**/{acc}.ref_expr.gtf"))
            gtf = alt[0] if alt else gtf
        if not gtf.exists():
            missing.append(acc)
            continue
        c, nm = gene_counts(gtf, read_length)
        per_lib[acc] = c
        names.update({k: v for k, v in nm.items() if k not in names})
        log(f"  {acc}: {len(c):,} genes, {sum(c.values()):,.0f} est. reads")
    if missing:
        print("\nMissing the -e -G quantification for: "
              + ", ".join(missing))
        print("Regenerate with, per library:")
        print("  python scripts/21_align_assemble.py --only <ACC> "
              "--ref-quant --force")
        if not per_lib:
            raise SystemExit("no usable libraries")
    libs = [s["run_accession"] for s in samples if s["run_accession"] in per_lib]
    genes = sorted({g for c in per_lib.values() for g in c})
    mat = np.zeros((len(libs), len(genes)), dtype=np.int64)
    idx = {g: j for j, g in enumerate(genes)}
    for i, acc in enumerate(libs):
        for g, v in per_lib[acc].items():
            mat[i, idx[g]] = int(round(v))
    return libs, genes, mat, names


# ---------------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------------

def qc(libs, groups, mat, size_factors=None) -> dict:
    """Library sizes, a PCA, and the mandated depth-outlier check."""
    tot = mat.sum(axis=1)
    print("\n" + "-" * 78)
    print("LIBRARY QC")
    print("-" * 78)
    med = float(np.median(tot))
    for lib, g, t in zip(libs, groups, tot):
        ratio = t / med if med else 0.0
        flag = "  <- DEPTH OUTLIER" if ratio > 2 or ratio < 0.5 else ""
        print(f"  {lib:<14} {g:<20} {t:>12,} est. reads  "
              f"{ratio:>4.2f}x median{flag}")

    # log-CPM PCA on the most variable genes, via SVD (no sklearn needed)
    cpm = mat / np.maximum(tot[:, None], 1) * 1e6
    lg = np.log2(cpm + 1)
    keep = np.argsort(-lg.var(axis=0))[:2000]
    X = lg[:, keep]
    X = X - X.mean(axis=0, keepdims=True)
    U, S, _ = np.linalg.svd(X, full_matrices=False)
    pc = U[:, :2] * S[:2]
    var = (S ** 2 / (S ** 2).sum())[:2] * 100
    print(f"\n  PCA on 2,000 most variable genes "
          f"(PC1 {var[0]:.1f}%, PC2 {var[1]:.1f}%)")
    for lib, g, p in zip(libs, groups, pc):
        print(f"    {lib:<14} {g:<20} PC1 {p[0]:+8.2f}  PC2 {p[1]:+8.2f}")

    # does PC1 separate the groups we are contrasting?
    uniq = sorted(set(groups))
    if len(uniq) >= 2:
        a = pc[[i for i, g in enumerate(groups) if g == uniq[0]], 0]
        b = pc[[i for i, g in enumerate(groups) if g == uniq[1]], 0]
        sep = abs(a.mean() - b.mean()) / max(1e-9, pc[:, 0].std())
        print(f"\n  PC1 group separation: {sep:.2f} SD "
              f"({uniq[0]} vs {uniq[1]})")
        if sep < 0.5:
            print("  WARNING: the groups do not separate on PC1. Treat every")
            print("  downstream result with suspicion.")
    return {"library_sizes": {l: int(t) for l, t in zip(libs, tot)},
            "pc1": {l: float(p[0]) for l, p in zip(libs, pc)},
            "pc2": {l: float(p[1]) for l, p in zip(libs, pc)},
            "pct_var": [float(v) for v in var]}


# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--contrast", default="WSSV:PmSTAT_dsRNA_WSSV",
                    help="'A:B'. The default is the only clean "
                         "single-variable axis in the dataset.")
    ap.add_argument("--read-length", type=float, default=150.0,
                    help="used only to convert coverage to counts; a global "
                         "constant that cancels in normalisation")
    ap.add_argument("--min-count", type=int, default=10,
                    help="drop genes below this summed count")
    ap.add_argument("--min-libraries", type=int, default=3,
                    help="...in at least this many libraries")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--lfc-report", type=float, default=2.0,
                    help="second, stringent threshold reported only for "
                         "concordance with the published DEG counts")
    ap.add_argument("--stat-gene", default=None,
                    help="gene id of the STAT ortholog, e.g. "
                         "gene-LOC119587821, to use as the positive control "
                         "instead of searching product descriptions")
    ap.add_argument("--exclude", default="",
                    help="comma-separated run accessions to drop, for "
                         "leave-one-out checks")
    ap.add_argument("--inventory", action="store_true",
                    help="report what inputs exist and stop")
    args = ap.parse_args()

    cfg = cfgmod.load()
    root = Path(__file__).resolve().parents[1]
    sheet = root / "config" / "samples.tsv"
    if not sheet.exists():
        raise SystemExit(f"missing {sheet}")
    samples = read_samples(sheet)
    align = cfg.p("align")

    want_a, _, want_b = args.contrast.partition(":")
    drop = {x.strip() for x in args.exclude.split(",") if x.strip()}
    use = [s for s in samples if s.get("use", "1").strip() == "1"
           and s["group"] in (want_a, want_b)
           and s["run_accession"] not in drop]
    if drop:
        print(f"  dropped by --exclude: {', '.join(sorted(drop))}")
    excluded = [s for s in samples if s.get("use", "1").strip() != "1"]

    print("=" * 78)
    print(f"CONTRAST  {want_a}  vs  {want_b}")
    print("=" * 78)
    print(f"  positive log2FoldChange = higher in {want_a}")
    print("                          = suppressed when STAT is knocked down")
    print("                          = STAT-DEPENDENT")
    counts_by_group = defaultdict(int)
    for s in use:
        counts_by_group[s["group"]] += 1
    for g, n in sorted(counts_by_group.items()):
        print(f"  {g:<22} n = {n}")
    for s in excluded:
        print(f"  excluded: {s['run_accession']} ({s.get('note', '-')})")
    if len(counts_by_group) < 2:
        raise SystemExit("the contrast needs two groups with use=1 libraries")

    gff = next(iter(cfg.p("genome").glob("*_genomic.gff")), None)
    if gff is None:
        print("\nNo *_genomic.gff in the genome directory. Gene products")
        print("cannot be read, so the characterised/uncharacterised split")
        print("and the STAT positive control are both unavailable.")
        ann: dict[str, dict] = {}
    else:
        log(f"reading annotation from {gff.name}")
        ann = gff_annotation(gff)
        n_unc = sum(1 for v in ann.values()
                    if is_uncharacterised(v.get("product", "")))
        n_desc = sum(1 for v in ann.values()
                     if v.get("product") and
                     not is_uncharacterised(v["product"]))
        log(f"  {len(ann):,} genes; {n_desc:,} with a described product, "
            f"{n_unc:,} uncharacterised")

    log(f"reading -e -G quantification from {align}")
    libs, genes, mat, names = build_matrix(align, use, args.read_length)
    groups = [next(s["group"] for s in use if s["run_accession"] == l)
              for l in libs]
    print(f"\n  {len(libs)} libraries x {len(genes):,} reference genes")
    if args.inventory:
        qc(libs, groups, mat)
        print("\n(--inventory: stopping before DE)")
        return 0

    qc_report = qc(libs, groups, mat)

    # ---- filter -------------------------------------------------------
    expressed = ((mat >= args.min_count).sum(axis=0) >= args.min_libraries)
    keep = np.flatnonzero(expressed)
    log(f"{len(keep):,} of {len(genes):,} genes pass "
        f">={args.min_count} counts in >={args.min_libraries} libraries")
    mat_f = mat[:, keep]
    genes_f = [genes[j] for j in keep]

    # ---- DE -----------------------------------------------------------
    try:
        import pandas as pd
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as e:
        print(f"\nmissing {e.name}. Install into the project env:")
        print("    pip install pydeseq2")
        return 2

    counts_df = pd.DataFrame(mat_f, index=libs, columns=genes_f)
    meta = pd.DataFrame({"group": groups}, index=libs)
    log("fitting the negative-binomial model")
    dds = DeseqDataSet(counts=counts_df, metadata=meta, design="~group",
                       refit_cooks=True, quiet=True)
    dds.deseq2()

    # Discover the coefficient and take the contrast direction FROM it.
    coeffs = [c for c in dds.varm["LFC"].columns if c != "Intercept"]
    if len(coeffs) != 1:
        raise SystemExit(f"expected one group coefficient, got {coeffs}")
    coeff = coeffs[0]
    target = coeff.split("[T.")[-1].rstrip("]")
    ref = want_a if target == want_b else want_b
    log(f"coefficient {coeff} -> contrast [{target} vs {ref}]")
    flipped = (target != want_a)
    if flipped:
        print(f"  NOTE: pydeseq2 fixes the reference level alphabetically, "
              f"so the\n  fitted direction is {target} vs {ref}. Signs are "
              f"flipped below to\n  keep 'positive = STAT-dependent'.")

    st = DeseqStats(dds, contrast=["group", target, ref],
                    alpha=args.alpha, quiet=True)
    st.summary()
    raw = st.results_df.copy()

    st.lfc_shrink(coeff=coeff)
    shrunk = st.results_df.copy()

    ok = np.isfinite(raw.log2FoldChange) & np.isfinite(shrunk.log2FoldChange)
    r = float(np.corrcoef(raw.log2FoldChange[ok],
                          shrunk.log2FoldChange[ok])[0, 1])
    if r <= 0.5:
        raise SystemExit(
            f"shrunk and unshrunk fold changes correlate {r:+.3f}. The "
            "coefficient does not match the contrast and every sign would "
            "be inverted. Refusing to report.")
    log(f"sign guard: corr(raw, shrunk) = {r:+.3f}")

    res = shrunk
    if flipped:
        res = res.copy()
        res["log2FoldChange"] = -res["log2FoldChange"]
    res["gene_name"] = [ann.get(g, {}).get("name", "") or names.get(g, "")
                        for g in res.index]
    res["product"] = [ann.get(g, {}).get("product", "") for g in res.index]
    res["biotype"] = [ann.get(g, {}).get("biotype", "") for g in res.index]
    res["uncharacterised"] = [is_uncharacterised(p) if p else None
                              for p in res["product"]]

    # ---- the positive control -----------------------------------------
    print("\n" + "=" * 78)
    print("POSITIVE CONTROL — DID THE KNOCKDOWN WORK?")
    print("=" * 78)
    if args.stat_gene:
        strict_hits = [g for g in res.index if g == args.stat_gene]
        if not strict_hits:
            print(f"  --stat-gene {args.stat_gene} is not in the tested set.")
    else:
        strict_hits = [g for g in res.index
                       if STAT_STRICT.search(
                           res.at[g, "product"] or res.at[g, "gene_name"] or g)]
    loose_hits = [g for g in res.index
                  if g not in set(strict_hits)
                  and STAT_LOOSE.search(
                      res.at[g, "product"] or res.at[g, "gene_name"] or g)]
    ctrl = None
    if not strict_hits:
        print("  No gene whose NAME looks like a STAT symbol is in the")
        print("  tested set. The knockdown is therefore UNVERIFIED here.")
        print("  Two innocent explanations before concluding it failed:")
        print("    - the P. monodon STAT ortholog may be annotated only as")
        print("      LOC<digits>, in which case no name search can find it")
        print("    - it may have been filtered out by --min-count")
        print("  Identify it by homology instead (the dsRNA target sequence")
        print("  against the reference proteome) and re-check that gene by")
        print("  ID. Until then, every result below assumes a knockdown")
        print("  this script could not confirm.")
    if loose_hits:
        print("  Other name matches, for inspection (not used as the")
        print("  control):")
        for g in loose_hits[:10]:
            row = res.loc[g]
            print(f"    {names.get(g, g):<28} log2FC "
                  f"{row.log2FoldChange:+6.2f}  padj {row.padj:.3g}")
    norm = mat / np.maximum(mat.sum(axis=1, keepdims=True), 1) * 1e6
    gidx = {g: j for j, g in enumerate(genes)}
    for g in strict_hits:
        row = res.loc[g]
        fin = np.isfinite(row.padj)
        mark = "" if not fin else (
            "  <- significant" if row.padj < args.alpha else "  <- NOT significant")
        print(f"  {g}")
        print(f"    product: {res.at[g, 'product'][:68]}")
        print(f"    log2FC {row.log2FoldChange:+6.2f}   "
              f"padj {row.padj if fin else float('nan'):.3g}{mark}")
        # Per-library CPM, so the knockdown is visible rather than inferred
        if g in gidx:
            print("    CPM per library:")
            for lib, grp, v in zip(libs, groups, norm[:, gidx[g]]):
                print(f"      {lib:<14} {grp:<20} {v:9.1f}")
        if ctrl is None:
            ctrl = {"gene": g, "product": res.at[g, "product"],
                    "log2FoldChange": float(row.log2FoldChange),
                    "padj": float(row.padj) if fin else None,
                    "cpm": {lib: float(v) for lib, v in
                            zip(libs, norm[:, gidx[g]])} if g in gidx else {}}
    if ctrl and ctrl["padj"] is not None:
        if ctrl["log2FoldChange"] > 0 and ctrl["padj"] < args.alpha:
            print("\n  STAT itself is suppressed in the knockdown arm, in the")
            print("  expected direction. The contrast behaves as designed.")
        else:
            print("\n  WARNING: STAT is not significantly suppressed in the")
            print("  expected direction. If the knockdown did not take, the")
            print("  STAT-dependent set below is not interpretable. Resolve")
            print("  this before building anything on it.")

    # ---- results ------------------------------------------------------
    sig = res[(res.padj < args.alpha)]
    up = sig[sig.log2FoldChange > 0]      # STAT-dependent
    dn = sig[sig.log2FoldChange < 0]
    strict = res[(res.padj < args.alpha)
                 & (res.log2FoldChange.abs() > args.lfc_report)]
    s_up = strict[strict.log2FoldChange > 0]
    s_dn = strict[strict.log2FoldChange < 0]

    print("\n" + "=" * 78)
    print("DIFFERENTIAL EXPRESSION")
    print("=" * 78)
    print(f"  pre-registered primary (apeglm shrinkage, no LFC threshold, "
          f"BH < {args.alpha}):")
    print(f"    {len(sig):,} genes    {len(up):,} STAT-dependent "
          f"(down on knockdown)    {len(dn):,} STAT-repressed")
    print(f"\n  stringent, for concordance with the published counts "
          f"(|LFC| > {args.lfc_report:g}):")
    print(f"    {len(strict):,} genes    {len(s_up):,} up in {want_a}    "
          f"{len(s_dn):,} up in {want_b}")
    print("    published: 675 up in WSSV, 138 up in dsRNA+WSSV")
    print("    These will not match exactly: the published analysis used")
    print("    observed counts and no shrinkage, this one uses")
    print("    coverage-derived counts and apeglm. Concordance of direction")
    print("    and rough magnitude is what matters.")

    # ---- is the dark fraction enriched? -------------------------------
    print("\n" + "-" * 78)
    print("IS THE UNCHARACTERISED FRACTION STAT-DEPENDENT?")
    print("-" * 78)
    # Protein-coding only, and only genes with a product. tRNAs carry the
    # only gene symbols in this annotation and would otherwise dominate.
    pc = res[(res.biotype == "protein_coding") & res["product"].astype(bool)]
    enrich = None
    if len(pc) < 50 or pc.uncharacterised.isna().all():
        print(f"  Only {len(pc)} protein-coding genes with a product in the")
        print("  tested set. Too few to test; check the annotation parse.")
    else:
        dep = (pc.padj < args.alpha) & (pc.log2FoldChange > 0)
        unc = pc.uncharacterised.fillna(False).astype(bool)
        a = int((dep & unc).sum())        # uncharacterised, STAT-dependent
        b = int((~dep & unc).sum())       # uncharacterised, not
        c = int((dep & ~unc).sum())       # described, STAT-dependent
        d = int((~dep & ~unc).sum())      # described, not
        print(f"  protein-coding genes tested                 : {len(pc):>6}")
        print(f"  uncharacterised, STAT-dependent             : {a:>6}")
        print(f"  uncharacterised, not STAT-dependent         : {b:>6}")
        print(f"  described, STAT-dependent                   : {c:>6}")
        print(f"  described, not STAT-dependent               : {d:>6}")
        p_unc = a / max(1, a + b)
        p_dsc = c / max(1, c + d)
        print(f"\n  {100 * p_unc:.1f}% of uncharacterised are STAT-dependent")
        print(f"  {100 * p_dsc:.1f}% of described are STAT-dependent")
        enrich = {"table": {"unc_dep": a, "unc_not": b,
                            "desc_dep": c, "desc_not": d},
                  "frac_uncharacterised": p_unc, "frac_described": p_dsc,
                  "n_protein_coding": int(len(pc))}

        # A zero or near-zero cell makes the odds ratio meaningless. The
        # first version of this script printed "OR 0.00, p = 2.72e-08" off
        # a table whose fourth cell was 0 and called it a finding.
        cells = [a, b, c, d]
        if min(cells) < 5:
            print(f"\n  NOT TESTED: smallest cell is {min(cells)}. An odds")
            print("  ratio from a table with an empty or near-empty cell is")
            print("  not interpretable, so none is reported.")
            enrich["refused"] = f"smallest cell {min(cells)}"
        else:
            try:
                from scipy.stats import fisher_exact
                odds, pv = fisher_exact([[a, b], [c, d]])
                enrich.update({"odds_ratio": float(odds),
                               "p_value": float(pv)})
                direction = ("ENRICHED among uncharacterised" if odds > 1
                             else "DEPLETED among uncharacterised")
                print(f"  Fisher exact: OR {odds:.2f}, p = {pv:.3g}")
                if pv < 0.05:
                    print(f"  STAT-dependence is {direction}.")
                    if odds < 1:
                        print("  That is the OPPOSITE of the title's framing:")
                        print("  the dark fraction is LESS STAT-dependent than")
                        print("  the characterised one. The claim must change.")
                else:
                    print("  Not significant. The uncharacterised fraction is")
                    print("  not preferentially STAT-dependent, and the")
                    print("  title's framing would need to change.")
            except ImportError:
                print("  (scipy unavailable; odds ratio not computed)")

        print("\n  Caveat: 'uncharacterised' is NCBI's product annotation,")
        print("  which correlates with expression level, gene length and")
        print("  conservation. An expression- and length-matched background")
        print("  is the next refinement before this is a strong claim.")

    # ---- outputs ------------------------------------------------------
    out_dir = cfg.p("results")
    out_dir.mkdir(parents=True, exist_ok=True)
    table = out_dir / "de_stat_dependence.csv"
    res.sort_values("padj").to_csv(table)
    summary = {
        "contrast": {"a": want_a, "b": want_b,
                     "fitted_coeff": coeff, "sign_flipped": bool(flipped),
                     "sign_guard_corr": r},
        "libraries": {l: g for l, g in zip(libs, groups)},
        "excluded": [{"run": s["run_accession"], "note": s.get("note", "")}
                     for s in excluded],
        "read_length": args.read_length,
        "n_genes_quantified": len(genes),
        "n_genes_tested": len(genes_f),
        "alpha": args.alpha,
        "primary": {"n_sig": len(sig), "n_stat_dependent": len(up),
                    "n_stat_repressed": len(dn)},
        "stringent": {"lfc": args.lfc_report, "n_sig": len(strict),
                      "n_up_a": len(s_up), "n_up_b": len(s_dn),
                      "published_up_a": 675, "published_up_b": 138},
        "positive_control": ctrl,
        "uncharacterised_enrichment": enrich,
        "qc": qc_report,
    }
    (out_dir / "de_stat_summary.json").write_text(json.dumps(summary, indent=2))
    log(f"wrote {table}")
    log(f"wrote {out_dir / 'de_stat_summary.json'}")

    print("\n  Counts are coverage-derived (prepDE-style), not observed read")
    print("  counts, because the BAMs were not retained. State that in the")
    print("  methods. Re-aligning with --keep-bam and running featureCounts")
    print("  would remove the approximation if a reviewer asks.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
