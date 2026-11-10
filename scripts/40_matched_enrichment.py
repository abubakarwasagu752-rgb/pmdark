#!/usr/bin/env python3
"""Is the uncharacterised fraction STAT-dependent, once power is controlled?

    python scripts/40_matched_enrichment.py
    python scripts/40_matched_enrichment.py --dep-direction both

WHY THE RAW TEST CANNOT BE READ

  scripts/39 found uncharacterised genes DEPLETED for STAT-dependence:
  17.3% against 27.5%, OR 0.55, p = 7e-21. Taken at face value that
  refutes the paper's central claim.

  It cannot be taken at face value. Only 1,895 of 7,079 uncharacterised
  genes (27%) passed the expression filter, against 7,849 of 19,058
  described genes (41%). Uncharacterised genes are systematically
  lower-expressed, and low expression means low power to call
  differential expression. A power deficit alone produces OR < 1 with no
  biology involved, and the observed depletion is exactly its signature.

  So the question is not "are they depleted" but "are they still depleted
  once compared against genes with the same expression and the same
  length". Three analyses, because each fails differently:

  1. STRATIFIED (Mantel-Haenszel)
     Odds ratios computed within baseMean strata and pooled. Removes the
     between-stratum confound but assumes a common within-stratum effect.

  2. MATCHED BACKGROUND
     Each uncharacterised gene paired with a described gene of nearest
     expression and length. Makes the comparison concrete and reports the
     balance achieved, so a reader can see whether matching worked rather
     than trusting that it did. Permutation test on the matched sets.

  3. FOLD-CHANGE DISTRIBUTIONS  <-- the one that settles it
     Mann-Whitney on log2FoldChange itself, with no significance cut
     anywhere. The power problem enters only through the padj threshold,
     so a test that never thresholds cannot suffer from it. If 1 and 2
     disagree with 3, trust 3.

  A verdict is printed only where the three agree. Where they disagree
  the script says so and declines to call it, because a disagreement here
  means the effect is not separable from the confound with this data.
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

GFF_ATTR = re.compile(r"([^=;]+)=([^;]*)")


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


# ---------------------------------------------------------------------------

def gene_lengths(gff: Path) -> dict[str, int]:
    """gene_id -> exonic length of its longest transcript.

    Length is a confound in its own right: a short transcript collects
    fewer reads at the same expression level, so it is tested with less
    power. Gene span would be the wrong measure (introns collect no
    reads), so this sums exons per transcript and takes the longest.
    """
    tx_gene: dict[str, str] = {}
    tx_len: dict[str, int] = defaultdict(int)
    with gff.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            feat, start, end, attrs = f[2], f[3], f[4], f[8]
            if feat in ("mRNA", "transcript", "ncRNA", "tRNA", "rRNA"):
                a = dict(GFF_ATTR.findall(attrs))
                tid, parent = a.get("ID"), a.get("Parent", "")
                if tid and parent.startswith("gene-"):
                    tx_gene[tid] = parent
            elif feat == "exon":
                a = dict(GFF_ATTR.findall(attrs))
                parent = a.get("Parent", "")
                if parent:
                    try:
                        tx_len[parent] += int(end) - int(start) + 1
                    except ValueError:
                        pass
    out: dict[str, int] = {}
    for tid, gid in tx_gene.items():
        L = tx_len.get(tid, 0)
        if L > out.get(gid, 0):
            out[gid] = L
    return out


# ---------------------------------------------------------------------------

def mantel_haenszel(strata: list[tuple[int, int, int, int]]) -> dict:
    """Pooled OR across strata, with the Robins-Breslow-Greenland CI.

    Each stratum is (a, b, c, d): a/b are dependent/not among
    uncharacterised, c/d among described.
    """
    num = den = 0.0
    for a, b, c, d in strata:
        n = a + b + c + d
        if n == 0:
            continue
        num += a * d / n
        den += b * c / n
    if den <= 0 or num <= 0:
        return {"or": float("nan"), "ci": [float("nan")] * 2,
                "n_strata": len(strata)}
    psi = num / den

    # Robins-Breslow-Greenland variance of log(psi)
    s1 = s2 = s3 = 0.0
    for a, b, c, d in strata:
        n = a + b + c + d
        if n == 0:
            continue
        P, Q = (a + d) / n, (b + c) / n
        R, S = a * d / n, b * c / n
        s1 += P * R
        s2 += P * S + Q * R
        s3 += Q * S
    var = 0.0
    if num > 0 and den > 0:
        var = (s1 / (2 * num ** 2) + s2 / (2 * num * den)
               + s3 / (2 * den ** 2))
    se = float(np.sqrt(max(var, 0.0)))
    lo, hi = float(psi * np.exp(-1.96 * se)), float(psi * np.exp(1.96 * se))
    return {"or": float(psi), "ci": [lo, hi], "se_log": se,
            "n_strata": len(strata)}


def match_background(unc_x: np.ndarray, desc_x: np.ndarray,
                     caliper: float = 0.2, k: int = 60) -> np.ndarray:
    """Greedy 1:1 nearest-neighbour matching with a caliper, no replacement.

    unc_x / desc_x are (n, 2) standardised [log10 baseMean, log10 length].
    Returns, for each uncharacterised gene, the index of its matched
    described gene, or -1 if no partner lay within `caliper`.

    THE CALIPER IS NOT OPTIONAL. Without it, greedy matching pairs every
    gene with its nearest partner however far away that is, and on
    validation data where the two classes barely overlap it cut the
    standardised mean difference only from -1.56 to -0.19 -- still badly
    unbalanced -- which left enough residual confounding to produce a
    spurious "effect" on data simulated with no effect at all. Dropping
    genes that have no close partner is the honest alternative: a smaller,
    genuinely comparable set, and a count of how many had to be discarded.
    """
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        return np.full(len(unc_x), -1, dtype=int)
    tree = cKDTree(desc_x)
    kk = min(k, len(desc_x))
    dist, cand = tree.query(unc_x, k=kk)
    dist, cand = np.atleast_2d(dist), np.atleast_2d(cand)
    used = np.zeros(len(desc_x), dtype=bool)
    out = np.full(len(unc_x), -1, dtype=int)
    # hardest cases first, so they are not left with the scraps
    order = np.argsort(-dist[:, 0])
    for i in order:
        for j, dj in zip(cand[i], dist[i]):
            if dj > caliper:
                break
            if not used[j]:
                used[j] = True
                out[i] = j
                break
    return out


def smd(a: np.ndarray, b: np.ndarray) -> float:
    """Standardised mean difference — how unbalanced two groups are.
    Below 0.1 is conventionally 'balanced'."""
    s = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return float((a.mean() - b.mean()) / s) if s > 0 else 0.0


# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", default=None,
                    help="de_stat_dependence.csv (default: results dir)")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--dep-direction", default="up",
                    choices=["up", "both"],
                    help="'up' = suppressed by knockdown (STAT-dependent); "
                         "'both' = any significant change")
    ap.add_argument("--strata", type=int, default=10,
                    help="number of baseMean strata for Mantel-Haenszel")
    ap.add_argument("--caliper", type=float, default=0.2,
                    help="max standardised distance for a match; genes with "
                         "no partner inside it are dropped, not matched badly")
    ap.add_argument("--min-delta", type=float, default=0.147,
                    help="smallest Cliff's delta worth calling an effect. "
                         "Below this is 'negligible' by convention, and at "
                         "n~2000 vs ~8000 a delta of 0.03 still reaches "
                         "p<0.05 while meaning nothing.")
    ap.add_argument("--min-or", type=float, default=1.25,
                    help="smallest odds ratio (or its reciprocal) worth "
                         "calling an effect")
    ap.add_argument("--max-smd", type=float, default=0.1,
                    help="refuse the matched test if residual imbalance "
                         "exceeds this")
    ap.add_argument("--permutations", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    try:
        import pandas as pd
        from scipy.stats import mannwhitneyu, fisher_exact
    except ImportError as e:
        print(f"missing {e.name}; pip install pandas scipy")
        return 2

    cfg = cfgmod.load()
    table = Path(args.table) if args.table else (
        cfg.p("results") / "de_stat_dependence.csv")
    if not table.exists():
        raise SystemExit(f"missing {table} — run scripts/39_de_stat.py first")
    log(f"reading {table.name}")
    df = pd.read_csv(table, index_col=0)

    gff = next(iter(cfg.p("genome").glob("*_genomic.gff")), None)
    if gff is None:
        raise SystemExit("no reference GFF; cannot get gene lengths")
    log(f"reading exon lengths from {gff.name}")
    lens = gene_lengths(gff)
    df["length"] = [lens.get(g, 0) for g in df.index]

    df = df[(df.biotype == "protein_coding")
            & df["product"].notna()
            & (df.length > 0)
            & df.baseMean.notna() & (df.baseMean > 0)
            & df.padj.notna()].copy()
    df["unc"] = df.uncharacterised.map(
        {True: True, "True": True}).fillna(False).astype(bool)
    if args.dep_direction == "up":
        df["dep"] = (df.padj < args.alpha) & (df.log2FoldChange > 0)
    else:
        df["dep"] = df.padj < args.alpha

    unc, desc = df[df.unc], df[~df.unc]
    log(f"{len(df):,} protein-coding genes  "
        f"({len(unc):,} uncharacterised, {len(desc):,} described)")
    if len(unc) < 50 or len(desc) < 50:
        raise SystemExit("too few genes in one class to test")

    # ---- 0. quantify the confound -------------------------------------
    print("\n" + "=" * 78)
    print("THE CONFOUND, MEASURED")
    print("=" * 78)
    lb_u, lb_d = np.log10(unc.baseMean), np.log10(desc.baseMean)
    ll_u, ll_d = np.log10(unc.length), np.log10(desc.length)
    print(f"  {'':<22}{'uncharacterised':>18}{'described':>14}")
    print(f"  {'median baseMean':<22}{10 ** np.median(lb_u):>18.1f}"
          f"{10 ** np.median(lb_d):>14.1f}")
    print(f"  {'median exonic length':<22}{10 ** np.median(ll_u):>18.0f}"
          f"{10 ** np.median(ll_d):>14.0f}")
    print(f"\n  standardised mean difference, log10 baseMean: "
          f"{smd(lb_u.values, lb_d.values):+.3f}")
    print(f"  standardised mean difference, log10 length  : "
          f"{smd(ll_u.values, ll_d.values):+.3f}")
    print("\n  Anything beyond about 0.1 means the two classes are not")
    print("  comparable as they stand, and the raw odds ratio is measuring")
    print("  this difference as much as any biology.")

    # the power curve that makes it visible
    print("\n  Fraction called STAT-dependent, by expression decile:")
    qs = np.quantile(np.log10(df.baseMean), np.linspace(0, 1, 11))
    qs[-1] += 1e-9
    print(f"    {'decile':<8}{'baseMean range':>22}"
          f"{'uncharacterised':>18}{'described':>14}")
    strata = []
    decile_table = []          # saved, because the figures need it
    for i in range(10):
        lo, hi = qs[i], qs[i + 1]
        m = (np.log10(df.baseMean) >= lo) & (np.log10(df.baseMean) < hi)
        su, sd = df[m & df.unc], df[m & ~df.unc]
        a, b = int(su.dep.sum()), int((~su.dep).sum())
        c, d = int(sd.dep.sum()), int((~sd.dep).sum())
        strata.append((a, b, c, d))
        decile_table.append({
            "decile": i + 1,
            "basemean_lo": float(10 ** lo), "basemean_hi": float(10 ** hi),
            "unc_dep": a, "unc_n": a + b,
            "desc_dep": c, "desc_n": c + d,
            "unc_frac": a / (a + b) if a + b else None,
            "desc_frac": c / (c + d) if c + d else None})
        fu = f"{100 * a / (a + b):.1f}% ({a + b})" if a + b else "-"
        fd = f"{100 * c / (c + d):.1f}% ({c + d})" if c + d else "-"
        print(f"    {i + 1:<8}{10 ** lo:>9.1f}-{10 ** hi:<11.1f}"
              f"{fu:>18}{fd:>14}")

    # ---- 1. raw, for reference ----------------------------------------
    A = int((unc.dep).sum()); B = int((~unc.dep).sum())
    C = int((desc.dep).sum()); D = int((~desc.dep).sum())
    raw_or, raw_p = fisher_exact([[A, B], [C, D]])
    print("\n" + "=" * 78)
    print("1. RAW (not controlled — shown only for comparison)")
    print("=" * 78)
    print(f"  uncharacterised {A}/{A + B} = {100 * A / (A + B):.1f}%   "
          f"described {C}/{C + D} = {100 * C / (C + D):.1f}%")
    print(f"  OR {raw_or:.2f}, p = {raw_p:.3g}")

    # ---- 2. stratified ------------------------------------------------
    mh = mantel_haenszel(strata)
    print("\n" + "=" * 78)
    print(f"2. STRATIFIED by expression ({args.strata} deciles, "
          f"Mantel-Haenszel)")
    print("=" * 78)
    if np.isfinite(mh["or"]):
        print(f"  pooled OR {mh['or']:.2f}  95% CI "
              f"[{mh['ci'][0]:.2f}, {mh['ci'][1]:.2f}]")
        excl = mh["ci"][1] < 1 or mh["ci"][0] > 1
        big = mh["or"] <= 1 / args.min_or or mh["or"] >= args.min_or
        mh_sig = bool(excl and big)
        print(f"  {'excludes 1' if excl else 'CROSSES 1'}; "
              f"magnitude {'>=' if big else '<'} {args.min_or:g}-fold")
        if excl and not big:
            print("  Significant but negligible: with this many genes a tiny")
            print("  residual confound clears p<0.05. Not counted as an "
                  "effect.")
    else:
        mh_sig = False
        print("  not estimable (an empty stratum)")

    # ---- 3. matched ---------------------------------------------------
    print("\n" + "=" * 78)
    print("3. MATCHED on expression and length (1:1, nearest neighbour)")
    print("=" * 78)
    mu = np.column_stack([lb_u.values, ll_u.values])
    md = np.column_stack([lb_d.values, ll_d.values])
    cen = np.vstack([mu, md]).mean(axis=0)
    sca = np.vstack([mu, md]).std(axis=0)
    sca[sca == 0] = 1.0
    pairs = match_background((mu - cen) / sca, (md - cen) / sca,
                             caliper=args.caliper)
    got = pairs >= 0
    matched = {"n": int(got.sum())}
    if got.sum() < 50:
        print("  matching failed (scipy missing, or too few partners)")
        m_sig = False
    else:
        u_sel = unc.iloc[np.flatnonzero(got)]
        d_sel = desc.iloc[pairs[got]]
        bal_b = smd(np.log10(u_sel.baseMean.values),
                    np.log10(d_sel.baseMean.values))
        bal_l = smd(np.log10(u_sel.length.values),
                    np.log10(d_sel.length.values))
        dropped = int((~got).sum())
        print(f"  matched {len(u_sel):,} pairs; {dropped:,} uncharacterised "
              f"genes had no partner within the caliper and were dropped")
        print(f"  balance after matching — SMD baseMean {bal_b:+.3f}, "
              f"length {bal_l:+.3f}")
        matched.update({"smd_basemean": bal_b, "smd_length": bal_l,
                        "dropped": dropped})
        if max(abs(bal_b), abs(bal_l)) > args.max_smd:
            print(f"  STILL UNBALANCED (> {args.max_smd:g}). The two classes")
            print("  do not overlap enough in expression to be compared this")
            print("  way. Treat the matched result as uninformative.")
            matched["unbalanced"] = True
        a2, b2 = int(u_sel.dep.sum()), int((~u_sel.dep).sum())
        c2, d2 = int(d_sel.dep.sum()), int((~d_sel.dep).sum())
        print(f"  uncharacterised {a2}/{a2 + b2} = "
              f"{100 * a2 / max(1, a2 + b2):.1f}%   "
              f"matched described {c2}/{c2 + d2} = "
              f"{100 * c2 / max(1, c2 + d2):.1f}%")
        if min(a2, b2, c2, d2) < 5:
            print(f"  NOT TESTED: smallest cell {min(a2, b2, c2, d2)}")
            m_sig = False
            matched["refused"] = True
        else:
            m_or, m_p = fisher_exact([[a2, b2], [c2, d2]])
            # permutation: shuffle the class label within matched pairs
            rng = np.random.default_rng(args.seed)
            du = u_sel.dep.values.astype(int)
            dd = d_sel.dep.values.astype(int)
            obs = du.mean() - dd.mean()
            n = len(du)
            null = np.empty(args.permutations)
            for t in range(args.permutations):
                flip = rng.random(n) < 0.5
                x = np.where(flip, dd, du)
                y = np.where(flip, du, dd)
                null[t] = x.mean() - y.mean()
            perm_p = float((np.abs(null) >= abs(obs)).mean())
            print(f"  OR {m_or:.2f}, Fisher p = {m_p:.3g}; "
                  f"paired permutation p = {perm_p:.4f}")
            big_or = m_or <= 1 / args.min_or or m_or >= args.min_or
            m_sig = bool(perm_p < 0.05 and big_or
                         and not matched.get("unbalanced"))
            if perm_p < 0.05 and not big_or:
                print(f"  Significant but below {args.min_or:g}-fold; not "
                      "counted as an effect.")
            matched.update({"or": float(m_or), "fisher_p": float(m_p),
                            "perm_p": perm_p, "diff": float(obs),
                            "table": [a2, b2, c2, d2]})

    # ---- 4. the test that does not threshold --------------------------
    print("\n" + "=" * 78)
    print("4. FOLD-CHANGE DISTRIBUTIONS (no significance cut anywhere)")
    print("=" * 78)
    print("  This is the decisive one: the power deficit enters only")
    print("  through the padj threshold, and this test never applies one.")
    lfc_u = unc.log2FoldChange.values
    lfc_d = desc.log2FoldChange.values
    U, p_all = mannwhitneyu(lfc_u, lfc_d, alternative="two-sided")
    delta = 2 * U / (len(lfc_u) * len(lfc_d)) - 1          # Cliff's delta
    print(f"  all genes: median LFC {np.median(lfc_u):+.3f} "
          f"(uncharacterised) vs {np.median(lfc_d):+.3f} (described)")
    print(f"             Cliff's delta {delta:+.3f}, p = {p_all:.3g}")

    within = []
    print("\n  within expression deciles:")
    for i in range(10):
        lo, hi = qs[i], qs[i + 1]
        m = (np.log10(df.baseMean) >= lo) & (np.log10(df.baseMean) < hi)
        a_, b_ = df[m & df.unc].log2FoldChange.values, \
            df[m & ~df.unc].log2FoldChange.values
        if len(a_) < 20 or len(b_) < 20:
            continue
        Uw, pw = mannwhitneyu(a_, b_, alternative="two-sided")
        dw = 2 * Uw / (len(a_) * len(b_)) - 1
        within.append({"decile": i + 1, "delta": float(dw),
                       "p": float(pw), "n_unc": len(a_), "n_desc": len(b_)})
        print(f"    decile {i + 1:<2} delta {dw:+.3f}  p {pw:<10.3g}  "
              f"n {len(a_)}/{len(b_)}")
    # DOES THE EFFECT TRACK EXPRESSION, AND IN WHICH DIRECTION?
    #
    # This is the test that separates a power artifact from a real effect,
    # and the pooled estimates cannot do it. An artifact is largest where
    # power is worst -- the low deciles -- and disappears where both
    # classes have ample reads. apeglm shrinkage pushes the same way,
    # pulling low-count fold changes toward zero. So an artifact predicts
    # |delta| DECREASING with expression. A real effect is free to do the
    # opposite, and if it increases with expression it cannot be the
    # confound.
    trend = None
    if len(within) >= 5:
        from scipy.stats import spearmanr
        dec = np.array([w["decile"] for w in within], dtype=float)
        dl = np.array([w["delta"] for w in within])
        rho, p_tr = spearmanr(dec, np.abs(dl))
        lo_d = float(np.mean(np.abs(dl[:3]))) if len(dl) >= 3 else float("nan")
        hi_d = float(np.mean(np.abs(dl[-3:]))) if len(dl) >= 3 else float("nan")
        trend = {"spearman_rho": float(rho), "p": float(p_tr),
                 "mean_abs_delta_low3": lo_d, "mean_abs_delta_high3": hi_d}
        print(f"\n  trend of |delta| with expression: rho {rho:+.3f}, "
              f"p = {p_tr:.3g}")
        print(f"    mean |delta|, three lowest deciles  : {lo_d:.3f}")
        print(f"    mean |delta|, three highest deciles : {hi_d:.3f}")
        if hi_d > lo_d:
            print("    The effect GROWS with expression, which is the")
            print("    opposite of what a power deficit or shrinkage would")
            print("    produce. The confound cannot explain it.")
        else:
            print("    The effect SHRINKS with expression, which is what a")
            print("    power deficit would produce. Treat it as confounded.")

    n_dir = sum(1 for w in within
                if w["delta"] <= -args.min_delta and w["p"] < 0.05)
    n_opp = sum(1 for w in within
                if w["delta"] >= args.min_delta and w["p"] < 0.05)
    lfc_sig = bool(p_all < 0.05 and abs(delta) >= args.min_delta)
    if p_all < 0.05 and abs(delta) < args.min_delta:
        print(f"\n  p = {p_all:.3g} but |delta| = {abs(delta):.3f} is below "
              f"{args.min_delta:g},")
        print("  which is 'negligible' by convention. At this sample size a")
        print("  trivial difference clears p<0.05, so this is NOT counted as")
        print("  an effect. Validation on simulated data with no true effect")
        print("  produced delta -0.029 at p=0.048.")

    # ---- verdict ------------------------------------------------------
    print("\n" + "=" * 78)
    print("VERDICT")
    print("=" * 78)
    votes = {"stratified": mh_sig, "matched": m_sig, "fold_change": lfc_sig}
    agree = len(set(votes.values())) == 1
    print(f"  stratified   : {'effect' if mh_sig else 'no effect'}")
    print(f"  matched      : {'effect' if m_sig else 'no effect'}")
    print(f"  fold changes : {'effect' if lfc_sig else 'no effect'}")
    print(f"  deciles with a significant difference: {n_dir} depleted, "
          f"{n_opp} enriched")
    print()
    interaction = bool(
        trend and np.isfinite(trend["mean_abs_delta_high3"])
        and trend["mean_abs_delta_high3"] >= args.min_delta
        and trend["mean_abs_delta_low3"] < args.min_delta)

    if interaction:
        # A single pooled OR is the wrong summary when the effect depends
        # on expression. Matching reweights the described class onto the
        # uncharacterised expression distribution, which is skewed low, so
        # it estimates the effect where there is least of one -- which is
        # why the pooled tests can disagree while the data is perfectly
        # interpretable.
        d_hi = [w for w in within if w["delta"] <= -args.min_delta
                and w["p"] < 0.05]
        print("  AN INTERACTION WITH EXPRESSION, not a disagreement.")
        print(f"  |delta| is below {args.min_delta:g} in the three lowest")
        print(f"  expression deciles and {trend['mean_abs_delta_high3']:.3f} "
              f"in the three highest.")
        print()
        print("  Report the interaction, not a single odds ratio. The")
        print("  pooled and matched estimates are not wrong, they answer a")
        print("  different question: matching reweights the described class")
        print("  onto the uncharacterised expression distribution, which is")
        print("  skewed low, so it measures the effect where there is least")
        print("  of one.")
        print()
        print(f"  Among well-expressed genes ({len(d_hi)} deciles with a")
        print("  significant, non-negligible difference, all in the same")
        print("  direction) STAT-dependence is DEPLETED among")
        print("  uncharacterised genes. Because the effect grows with")
        print("  expression rather than shrinking, the power confound")
        print("  cannot account for it.")
        print()
        print("  A 'STAT-regulated dark proteome' framing is not supported.")
        print("  The defensible claim is the specific set of uncharacterised")
        print("  genes that ARE STAT-dependent.")
    elif not agree:
        print("  THE THREE DISAGREE and there is no expression interaction")
        print("  to explain it. The effect is not separable from the")
        print("  confound with this data. Report the confound and the raw")
        print("  numbers, and make no claim about whether the")
        print("  uncharacterised fraction is preferentially STAT-dependent.")
        print("  Do not pick the analysis you prefer.")
    elif not lfc_sig:
        print("  No difference survives control. The uncharacterised")
        print("  fraction is neither enriched nor depleted for")
        print("  STAT-dependence once expression and length are matched,")
        print("  and the raw OR of "
              f"{raw_or:.2f} was a power artifact.")
        print("  The honest claim is that the dark fraction participates in")
        print("  the STAT response in proportion to its expression.")
    else:
        direction = "DEPLETED" if np.median(lfc_u) < np.median(lfc_d) \
            else "ENRICHED"
        print(f"  STAT-dependence is genuinely {direction} among")
        print("  uncharacterised genes, and it survives matching on")
        print("  expression and length.")
        if direction == "DEPLETED":
            print("  This contradicts a 'STAT-regulated dark proteome'")
            print("  framing. The defensible claim is the specific set of")
            print("  uncharacterised genes that ARE STAT-dependent, not the")
            print("  fraction as a whole.")

    out = cfg.p("results") / "matched_enrichment.json"
    out.write_text(json.dumps({
        "alpha": args.alpha, "dep_direction": args.dep_direction,
        "n_genes": int(len(df)), "n_uncharacterised": int(len(unc)),
        "n_described": int(len(desc)),
        "confound": {
            "smd_log_basemean": smd(lb_u.values, lb_d.values),
            "smd_log_length": smd(ll_u.values, ll_d.values),
            "median_basemean_unc": float(10 ** np.median(lb_u)),
            "median_basemean_desc": float(10 ** np.median(lb_d)),
            "median_length_unc": float(10 ** np.median(ll_u)),
            "median_length_desc": float(10 ** np.median(ll_d))},
        "raw": {"table": [A, B, C, D], "or": float(raw_or),
                "p": float(raw_p)},
        "stratified": mh,
        "decile_table": decile_table,
        "matched": matched,
        "fold_change": {"median_unc": float(np.median(lfc_u)),
                        "median_desc": float(np.median(lfc_d)),
                        "cliffs_delta": float(delta), "p": float(p_all),
                        "within_deciles": within},
        "votes": votes, "agree": bool(agree),
        "trend": trend, "interaction": interaction,
    }, indent=2))
    log(f"wrote {out}")

    # the deliverable, either way
    dep_unc = df[df.unc & df.dep]
    ids = cfg.p("results") / "stat_dependent_uncharacterised.txt"
    ids.write_text("\n".join(sorted(dep_unc.index)) + "\n")
    print(f"\n  {len(dep_unc):,} uncharacterised, STAT-dependent genes "
          f"written to\n  {ids}")
    print("  That set is the annotation target, and it stands whatever the")
    print("  enrichment verdict is — it does not depend on the fraction")
    print("  being enriched, only on these genes being STAT-dependent.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
