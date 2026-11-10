#!/usr/bin/env python3
"""Is 32.4% recovery biology or breakage? Ask the genes that cannot be absent.

    python scripts/28_recovery_by_family.py --threads 16

THE QUESTION

  27_recovery_exact.py put sequence-confirmed recovery at 7,770 / 24,011 =
  32.4%, with this attrition:

      coding genes                24,011   100.0%
      exons assembled             11,000    45.8%     -13,011
      ORF called on those exons     9,161    38.2%      -1,839
      annotated protein recovered   7,770    32.4%      -1,391

  A single percentage cannot say whether that is a hemocyte transcriptome
  behaving normally or a pipeline destroying real genes. But some genes
  cannot be absent from a living cell. Ribosomal proteins, the proteasome,
  actin, tubulin, elongation factors, histones, ATP synthase, the glycolytic
  enzymes — every nucleated cell transcribes these constitutively and at
  high abundance. They are the internal positive control this dataset
  already contains, and a conspecific one, which is strictly better here
  than BUSCO's cross-phylum markers.

  The mirror image is just as informative. Odorant and gustatory receptors,
  cuticle proteins, vitellogenin, opsins, sperm and testis-specific genes:
  a hemocyte has no business transcribing any of them, and their absence is
  the expected result, not a loss.

  PRE-STATED DECISION RULE, fixed before the run:

    housekeeping recovery >= 90%  -> the pipeline is sound. 32.4% overall
                                     is single-tissue expression plus
                                     depth, and is reported as such.
    housekeeping recovery 75-90%  -> partial loss; quantify it and carry
                                     the caveat into every denominator.
    housekeeping recovery < 75%   -> the pipeline is losing real genes.
                                     The dark proteome does not get defined
                                     until that is fixed.

  The same table also answers a question the paper needs on its own terms:
  how much of the penaeid immune repertoire did we actually capture. Toll,
  Dscam, C-type lectins, crustins, penaeidins, ALFs, lysozyme,
  prophenoloxidase, the JAK/STAT components. If those are recovered well,
  Aim 4 has a substrate.

TWO FIXES TO 27's OUTPUT

  * Its description tally printed nothing. NCBI gene lines here carry
    Name=LOC115123456, and a 4+ letter word regex finds nothing in "LOC".
    Gene function lives in product= on the mRNA lines, which is what this
    reads.

  * Its 1,972 "ORF on the exons but sequence FAILED" genes were never
    broken down. Partial assembly and real divergence look completely
    different in the two numbers diamond already reported: fragmentation
    gives high identity with low query coverage, divergence the opposite.
    This splits them.

NOTHING HERE CHANGES ANY OUTPUT. It only measures.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod  # noqa: E402

def _load(name: str):
    """Same helper the other scripts use; 28 predates it and I referenced
    it here without defining it, which crashed at the expression block."""
    sp = importlib.util.spec_from_file_location(
        "m" + name[:2], Path(__file__).resolve().parent / name)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


_spec = importlib.util.spec_from_file_location(
    "rx", Path(__file__).resolve().parent / "27_recovery_exact.py")
_rx = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_rx)

_qspec = importlib.util.spec_from_file_location(
    "qc", Path(__file__).resolve().parent / "25_protein_qc.py")
_qc = importlib.util.module_from_spec(_qspec)
_qspec.loader.exec_module(_qc)

BIN = _rx.BIN

# Genes no living cell can be without. The positive control.
HOUSEKEEPING = {
    "ribosomal protein": r"\bribosomal protein\b",
    "proteasome": r"\bproteasome\b",
    "ubiquitin": r"\bubiquitin\b",
    "actin": r"\bactin\b",
    "tubulin": r"\btubulin\b",
    "translation factor": r"\b(elongation factor|translation initiation factor)\b",
    "RNA polymerase": r"\bRNA polymerase\b",
    "histone": r"\bhistone\b",
    "ATP synthase": r"\bATP synthase\b",
    "cytochrome c oxidase": r"\bcytochrome c oxidase\b",
    "glycolytic": r"\b(glyceraldehyde-3-phosphate dehydrogenase|enolase|"
                  r"phosphoglycerate kinase|pyruvate kinase)\b",
    "chaperone": r"\b(heat shock protein|chaperonin)\b",
    "spliceosome": r"\b(small nuclear ribonucleoprotein|splicing factor)\b",
}

# A hemocyte should not be transcribing these. The negative control.
NOT_IN_HEMOCYTES = {
    "chemoreceptor": r"\b(odorant receptor|gustatory receptor|"
                     r"chemosensory|olfactory)\b",
    "cuticle": r"\bcuticle|cuticular\b",
    "vitellogenin / egg": r"\b(vitellogenin|vitelline|chorion)\b",
    "opsin / visual": r"\b(opsin|rhodopsin)\b",
    "sperm / testis": r"\b(sperm|testis-specific|spermatogenesis)\b",
    "muscle myosin": r"\bmyosin heavy chain, muscle\b",
}

# What the paper is actually about.
IMMUNE = {
    "C-type lectin": r"\b(c-type lectin|lectin)\b",
    "Toll / TLR": r"\btoll\b",
    "Dscam": r"\bdscam|down syndrome cell adhesion\b",
    "crustin": r"\bcrustin\b",
    "penaeidin": r"\bpenaeidin\b",
    "anti-LPS factor": r"\b(anti-lipopolysaccharide|alf\b)",
    "lysozyme": r"\blysozyme\b",
    "prophenoloxidase": r"\b(prophenoloxidase|phenoloxidase)\b",
    "hemocyanin": r"\b(hemocyanin|haemocyanin)\b",
    "JAK / STAT": r"\b(janus kinase|signal transducer and activator)\b",
    "NF-kB / Relish / IKK": r"\b(relish|nf-kappa|inhibitor of nuclear factor)\b",
    "caspase": r"\bcaspase\b",
    "serpin": r"\bserpin\b",
    "clip-domain protease": r"\bclip\b",
    "transglutaminase": r"\btransglutaminase\b",
    "thioester / TEP": r"\bthioester\b",
}

STOP = {"like", "protein", "containing", "domain", "family", "putative",
        "type", "subunit", "isoform", "uncharacterized", "partial",
        "homolog", "related", "precursor", "chain", "factor", "with"}


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def _attr(s: str, key: str) -> str | None:
    m = re.search(rf"(?:^|;){key}=([^;]+)", s)
    return m.group(1) if m else None


def unquote(s: str) -> str:
    for a, b in (("%2C", ","), ("%3B", ";"), ("%3D", "="), ("%25", "%"),
                 ("%28", "("), ("%29", ")"), ("%2F", "/")):
        s = s.replace(a, b)
    return s


def gene_products(gff: Path) -> dict[str, str]:
    """gene id -> product string, taken from its mRNA/CDS product= fields.

    This is where NCBI puts function. The gene line only has Name=LOC...,
    which is why 27's word tally came back empty.
    """
    rna2gene: dict[str, str] = {}
    prod: dict[str, str] = {}
    pend: list[tuple[str, str]] = []
    with gff.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            a = p[8]
            if p[2] in ("mRNA", "transcript"):
                rid, par = _attr(a, "ID"), _attr(a, "Parent")
                if rid and par:
                    rna2gene[rid] = par
                pr = _attr(a, "product")
                if pr and par:
                    prod.setdefault(par, unquote(pr))
            elif p[2] == "CDS":
                pr, par = _attr(a, "product"), _attr(a, "Parent")
                if pr and par:
                    pend.append((par, unquote(pr)))
    for par, pr in pend:
        g = rna2gene.get(par)
        if g:
            prod.setdefault(g, pr)
    return prod


def match_family(prod: dict[str, str], patterns: dict[str, str],
                 universe: set[str]) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for name, pat in patterns.items():
        rx = re.compile(pat, re.I)
        out[name] = {g for g in universe if rx.search(prod.get(g, ""))}
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--min-identity", type=float, default=70.0)
    ap.add_argument("--min-qcov", type=float, default=50.0)
    ap.add_argument("--min-bp", type=int, default=30)
    ap.add_argument("--min-tpm", type=float, default=1.0,
                    help="expression threshold for the per-family denominator")
    ap.add_argument("--min-libraries", type=int, default=3,
                    help="libraries a gene must reach --min-tpm in")
    args = ap.parse_args()

    cfg = cfgmod.load()
    out, res = cfg.p("proteins"), cfg.p("results")
    qc = json.loads((res / "protein_qc.json").read_text())
    pep, gff3 = Path(qc["pep"]), Path(qc["genome_gff3"])
    ref_gff = next(iter(cfg.p("genome").glob("*_genomic.gff")), None)
    tsv = out / "ref_vs_ours.tsv"
    if ref_gff is None:
        raise SystemExit("reference GFF not found")
    if not tsv.exists():
        raise SystemExit(f"{tsv} missing — run scripts/27_recovery_exact.py")

    log("building the annotation model")
    A = _rx.annotation(ref_gff)
    coding = A["coding"]
    log("reading product= names from the mRNA lines")
    prod = gene_products(ref_gff)
    named = sum(1 for g in coding if prod.get(g))
    log(f"  {named:,}/{len(coding):,} coding genes carry a product name")

    # ---- recovery sets ------------------------------------------------
    log("recomputing exonic transcript and ORF overlap")
    tx = _rx.gtf_exon_blocks(out / "merged_filtered.gtf")
    g_tx: set[str] = set()
    for chrom, strand, blocks in tx.values():
        strands = ("+", "-") if strand == "." else (strand,)
        for st in strands:
            g_tx |= _rx.exonic_overlap(A["exon_idx"], chrom, st, blocks,
                                       args.min_bp)
    g_tx &= coding

    cds = _rx.predicted_cds_blocks(gff3)
    lens = _qc.read_fasta_lengths(pep)
    prot_genes: dict[str, set[str]] = {}
    g_orf: set[str] = set()
    for pid in lens:
        r = cds.get(pid)
        hits = (_rx.exonic_overlap(A["exon_idx"], r[0], r[1], r[2],
                                   args.min_bp) & coding) if r else set()
        prot_genes[pid] = hits
        g_orf |= hits

    # LOCUS AGREEMENT IS REQUIRED HERE TOO.
    #
    # Without it this table credited 9 real opsin ORFs to 65 annotated
    # opsin paralogs (74.7%) and 83 cuticle ORFs to 410 cuticle genes
    # (64.9%), which destroyed the negative control: housekeeping 79.3%
    # against hemocyte-irrelevant 65.7% does not separate, and the script
    # correctly warned about it. The ORF column, which IS locus-aware,
    # separates cleanly in the same run: ribosomal 97.5% vs opsin 10.3%.
    # 29_fragmentation.py fixed this for the headline figure; this is the
    # back-port into the control panel.
    multi = out / "ref_vs_ours_multi.tsv"
    if multi.exists():
        tsv = multi          # one hit per query cannot see past a paralog
    log(f"reading {tsv.name} with locus agreement required")
    best: dict[str, tuple[float, float]] = {}
    for line in tsv.open():
        f = line.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        g = A["prot2gene"].get(f[0])
        if not g or g not in prot_genes.get(f[1], ()):
            continue
        try:
            ident, qcov = float(f[2]), float(f[3])
        except ValueError:
            continue
        # A gene may have several annotated proteins; keep its best hit,
        # ranked on coverage first because coverage is what distinguishes
        # a fragmented assembly from a divergent sequence.
        cur = best.get(g)
        if cur is None or (qcov, ident) > (cur[1], cur[0]):
            best[g] = (ident, qcov)
    g_seq = {g for g, (i, q) in best.items()
             if i >= args.min_identity and q >= args.min_qcov} & coding
    log(f"  transcript {len(g_tx):,}  ORF {len(g_orf):,}  "
        f"sequence {len(g_seq):,}")

    tot = len(coding)

    # ---- EXPRESSED denominator, per family ----------------------------
    #
    # The family table used to divide by every annotated member, which is
    # the same denominator error 30_expressed_denominator.py fixed at the
    # global level. Actin has 109 annotated paralogs and only 40 are
    # transcribed in hemocytes: dividing by 109 scored 33.9%, dividing by
    # what is actually present scored 92.5%.
    #
    # The two questions have to be read off different columns:
    #   % expressed            -> BIOLOGY. Housekeeping high, irrelevant low.
    #   recovery among expressed -> PIPELINE. Should be high for BOTH.
    # Reading the verdict off recovery-over-all-members conflated them.
    expressed: set[str] = set()
    tabs = sorted(cfg.p("align").glob("*/*.ref_gene_abund.tsv"))
    if tabs:
        _ed = _load("30_expressed_denominator.py")
        bad: set[str] = set()
        st = Path(__file__).resolve().parents[1] / "config" / "samples.tsv"
        if st.exists():
            rws = [ln for ln in st.read_text().splitlines()
                   if ln.strip() and not ln.startswith("#")]
            if rws:
                hd = rws[0].split("\t")
                try:
                    ai, ui = hd.index("run_accession"), hd.index("use")
                    for ln in rws[1:]:
                        fl = ln.split("\t")
                        if len(fl) > max(ai, ui) and fl[ui].strip() in (
                                "0", "FALSE", "false", "no"):
                            bad.add(fl[ai].strip())
                except ValueError:
                    pass
        alias = {g: g for g in coding}
        alias.update({g[5:]: g for g in coding if g.startswith("gene-")})
        cnt: dict[str, int] = defaultdict(int)
        used = 0
        for t in tabs:
            acc = t.name.split(".")[0]
            if acc in bad:
                continue
            used += 1
            for gid, tpm in _ed.read_abund(t).items():
                g = alias.get(gid) or alias.get("gene-" + gid)
                if g is not None and tpm >= args.min_tpm:
                    cnt[g] += 1
        expressed = {g for g, n in cnt.items() if n >= args.min_libraries}
        log(f"  expressed: {len(expressed):,} genes "
            f"(>={args.min_tpm} TPM in >={args.min_libraries} of {used} libs)")
    else:
        log("  no *.ref_gene_abund.tsv — falling back to all-member "
            "denominators, which understate recovery for paralogous families")

    # ---- the control panel --------------------------------------------
    def table(title: str, patterns: dict, note: str) -> dict:
        fams = match_family(prod, patterns, coding)
        print("\n" + "=" * 80)
        print(title)
        print("=" * 80)
        print(f"  {'family':<24}{'genes':>6}{'expr':>6}{'%expr':>7}"
              f"{'seq':>6}{'rec|all':>9}{'rec|expr':>10}")
        print("-" * 80)
        rows = {}
        agg_n = agg_s = agg_e = agg_se = 0
        for name, _ in sorted(patterns.items()):
            gset = fams[name]
            n = len(gset)
            if n == 0:
                continue
            e = len(gset & expressed) if expressed else 0
            s = len(gset & g_seq)
            se = len(gset & g_seq & expressed) if expressed else 0
            agg_n += n
            agg_s += s
            agg_e += e
            agg_se += se
            rx = f"{100*se/e:>9.1f}%" if e else f"{'n/a':>10}"
            print(f"  {name:<24}{n:>6,}{e:>6,}{100*e/n:>6.1f}%"
                  f"{s:>6,}{100*s/n:>8.1f}%{rx}")
            rows[name] = {"n": n, "expressed": e, "sequence": s,
                          "sequence_expressed": se,
                          "transcribed": len(gset & g_tx),
                          "orf": len(gset & g_orf)}
        print("-" * 80)
        pct = 100 * agg_s / max(1, agg_n)
        pe = 100 * agg_se / max(1, agg_e) if agg_e else 0.0
        pxp = 100 * agg_e / max(1, agg_n)
        print(f"  {'ALL OF THE ABOVE':<24}{agg_n:>6,}{agg_e:>6,}"
              f"{pxp:>6.1f}%{agg_s:>6,}{pct:>8.1f}%{pe:>9.1f}%")
        print(f"  {note}")
        return {"rows": rows, "n": agg_n, "sequence": agg_s,
                "pct": round(pct, 1), "expressed": agg_e,
                "sequence_expressed": agg_se,
                "pct_expressed_of_family": round(pxp, 1),
                "pct_recovery_of_expressed": round(pe, 1)}

    hk = table("POSITIVE CONTROL — genes no living cell can be without",
               HOUSEKEEPING,
               "every cell transcribes these constitutively; "
               "this is the pipeline test")
    neg = table("NEGATIVE CONTROL — genes a hemocyte should NOT transcribe",
                NOT_IN_HEMOCYTES,
                "low recovery here is the expected result, not a loss")
    imm = table("THE IMMUNE REPERTOIRE — what the paper is about",
                IMMUNE,
                "this is whether Aim 4 has a substrate")

    print("\n" + "=" * 80)
    print("VERDICT  (rule fixed in this file's docstring before the run)")
    print("=" * 80)
    print(f"  overall recovery          {100*len(g_seq)/tot:>6.1f}%  "
          f"({len(g_seq):,}/{tot:,})")
    print("\n  BIOLOGY — what fraction of each family is expressed at all")
    print(f"    housekeeping            {hk['pct_expressed_of_family']:>6.1f}%")
    print(f"    hemocyte-irrelevant     {neg['pct_expressed_of_family']:>6.1f}%"
          "   <- should be LOW; that is the control working")
    print(f"    immune repertoire       {imm['pct_expressed_of_family']:>6.1f}%")
    sep = hk["pct_expressed_of_family"] - neg["pct_expressed_of_family"]
    print(f"    separation              {sep:>6.1f} points")
    print("\n  PIPELINE — recovery AMONG EXPRESSED genes")
    print(f"    housekeeping            {hk['pct_recovery_of_expressed']:>6.1f}%"
          f"  ({hk['sequence_expressed']:,}/{hk['expressed']:,})  <- decides it")
    print(f"    hemocyte-irrelevant     {neg['pct_recovery_of_expressed']:>6.1f}%"
          "   should ALSO be high: if a hemocyte does")
    print("                                    transcribe an opsin, we should "
          "find it")
    print(f"    immune repertoire       {imm['pct_recovery_of_expressed']:>6.1f}%")
    print(f"    given exons assembled   "
          f"{100*len(g_seq & g_tx)/max(1,len(g_tx)):>6.1f}%")
    print("-" * 80)
    hk_pipeline = (hk["pct_recovery_of_expressed"] if expressed
                   else hk["pct"])
    if not expressed:
        print("  NOTE no expression tables; verdict falls back to the "
              "all-member\n  denominator, which understates paralogous "
              "families.")
    if hk_pipeline >= 90:
        v = ("PASS — the pipeline is sound. Overall 32.4% reflects "
             "single-tissue\n  expression and depth, and is reported as such.")
    elif hk_pipeline >= 75:
        v = ("MARGINAL — partial loss of real genes. Quantify it and carry "
             "the\n  caveat into every downstream denominator.")
    else:
        v = ("FAIL — real genes are being lost. Do not define the dark "
             "proteome\n  on top of this; fix assembly or ORF calling first.")
    print(f"  {v}")
    if expressed and sep < 20:
        print("\n  WARNING the controls barely separate on EXPRESSION. That "
              "pattern\n  means the measurement is not discriminating, "
              "whatever the\n  absolute values say.")

    # ---- why did 1,972 genes fail the sequence check? -----------------
    failed = (g_orf - g_seq)
    nohit = {g for g in failed if g not in best}
    withhit = failed - nohit
    lowq = {g for g in withhit
            if best[g][0] >= args.min_identity and best[g][1] < args.min_qcov}
    lowi = {g for g in withhit
            if best[g][0] < args.min_identity and best[g][1] >= args.min_qcov}
    both = withhit - lowq - lowi
    print("\n" + "=" * 80)
    print(f"THE {len(failed):,} GENES WITH AN ORF ON THEIR EXONS BUT NO "
          "SEQUENCE MATCH")
    print("=" * 80)
    print(f"  no diamond hit at all                      {len(nohit):>7,}")
    print(f"  high identity, LOW coverage (fragmented)   {len(lowq):>7,}"
          "   <- our assembly is partial")
    print(f"  good coverage, LOW identity (divergent)    {len(lowi):>7,}"
          "   <- different sequence")
    print(f"  both below threshold                       {len(both):>7,}")
    if withhit:
        qs = sorted(best[g][1] for g in withhit)
        ids = sorted(best[g][0] for g in withhit)
        print(f"\n  median query coverage {qs[len(qs)//2]:.0f}%   "
              f"median identity {ids[len(ids)//2]:.0f}%")
        print("  fragmentation and divergence are distinguishable here, and")
        print("  only fragmentation is something we could fix")

    # ---- what ARE the never-transcribed genes? ------------------------
    never = coding - g_tx
    words: dict[str, int] = defaultdict(int)
    for g in never:
        for w in re.findall(r"[A-Za-z]{4,}", prod.get(g, "").lower()):
            if w not in STOP:
                words[w] += 1
    print("\n" + "=" * 80)
    print(f"THE {len(never):,} GENES WHOSE EXONS WERE NEVER TRANSCRIBED")
    print("=" * 80)
    top = sorted(words.items(), key=lambda kv: -kv[1])[:24]
    if top:
        for i in range(0, len(top), 3):
            print("   " + "".join(f"{w:<18}{n:>6,}   " for w, n in top[i:i+3]))
    else:
        print("   (no product names found — check the GFF attributes)")
    unch = sum(1 for g in never
               if "uncharacterized" in prod.get(g, "").lower())
    print(f"\n  'uncharacterized' among them: {unch:,} "
          f"({100*unch/max(1,len(never)):.0f}%)")
    unch_all = sum(1 for g in coding
                   if "uncharacterized" in prod.get(g, "").lower())
    print(f"  'uncharacterized' across all coding genes: {unch_all:,} "
          f"({100*unch_all/tot:.0f}%)")
    print("  if the never-transcribed set is enriched for uncharacterized")
    print("  genes, many are weakly-supported annotations rather than genes")
    print("  we failed to detect")

    report = {
        "coding_genes": tot,
        "recovery": {"transcript_exonic": len(g_tx), "orf_exonic": len(g_orf),
                     "sequence": len(g_seq),
                     "conditional_pct": round(
                         100 * len(g_seq & g_tx) / max(1, len(g_tx)), 1)},
        "housekeeping": hk, "not_in_hemocytes": neg, "immune": imm,
        "sequence_failure_modes": {
            "no_hit": len(nohit), "fragmented_low_qcov": len(lowq),
            "divergent_low_identity": len(lowi), "both": len(both)},
        "never_transcribed": {
            "n": len(never), "uncharacterized": unch,
            "uncharacterized_all_coding": unch_all,
            "top_words": dict(top)},
        "verdict_housekeeping_pct": hk["pct"],
        "housekeeping_recovery_of_expressed": hk["pct_recovery_of_expressed"],
        "control_separation_points": round(sep, 1),
        "expressed_denominator_used": bool(expressed),
    }
    rp = res / "recovery_by_family.json"
    rp.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {rp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
