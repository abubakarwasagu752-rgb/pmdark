#!/usr/bin/env python3
"""Collect every result into one document, with nothing retyped.

    python scripts/42_consolidate.py
    python scripts/42_consolidate.py --out /path/to/RESULTS.md

WHY

  The project's numbers now live in about a dozen JSON and CSV files.
  Writing a paper by reading them off a screen is how a figure ends up
  disagreeing with its own caption, and how a number that was superseded
  three runs ago survives into a submission.

  So this file contains NO numbers. Every value in the output is read
  from a results file at run time, and every section names the file it
  came from and when that file was last written. If a result is missing
  the document says so in place rather than quietly leaving a gap, and if
  two files disagree about the same quantity it says that too.

  Re-run it after any analysis changes. If the document changes, the
  paper needs changing.

WHAT IT ALSO DOES

  Runs a set of consistency checks that encode decisions already made --
  that the hardened orphan calibration is the one applied, that both
  homology baselines used the same transfer threshold, that the knockdown
  positive control is significant in the right direction -- and reports
  PASS/FAIL for each. A FAIL means something regressed, not that the
  document is wrong.

  It also lists, from the data rather than from memory, the claims that
  are NOT supported: every paired difference whose interval crosses zero,
  every tier that is undefined by construction, every test that was
  refused. That list is as much a part of the results as the table.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod                              # noqa: E402


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


class Source:
    """A results file, loaded once, with provenance recorded."""

    def __init__(self, root: Path):
        self.root = root
        self.seen: dict[str, dict] = {}
        self.missing: list[str] = []

    def json(self, name: str) -> dict | None:
        if name in self.seen:
            return self.seen[name].get("data")
        p = self.root / name
        if not p.exists():
            self.missing.append(name)
            return None
        try:
            data = json.loads(p.read_text())
        except ValueError as e:
            self.missing.append(f"{name} (unreadable: {e})")
            return None
        self.seen[name] = {
            "data": data,
            "mtime": datetime.fromtimestamp(
                p.stat().st_mtime, tz=timezone.utc).isoformat(timespec="minutes"),
            "bytes": p.stat().st_size,
        }
        return data

    def glob_json(self, pattern: str) -> dict[str, dict]:
        out = {}
        for p in sorted(self.root.glob(pattern)):
            d = self.json(p.name)
            if d is not None:
                out[p.name] = d
        return out

    def lines(self, name: str) -> list[str] | None:
        p = self.root / name
        if not p.exists():
            self.missing.append(name)
            return None
        self.seen[name] = {
            "data": None,
            "mtime": datetime.fromtimestamp(
                p.stat().st_mtime, tz=timezone.utc).isoformat(timespec="minutes"),
            "bytes": p.stat().st_size,
        }
        return [ln for ln in p.read_text().splitlines() if ln.strip()]


def f(x, nd: int = 3) -> str:
    """Format a number, or say plainly that it is absent."""
    if x is None:
        return "—"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if v != v:
        return "nan"
    if v in (float("inf"), float("-inf")):
        return "inf"
    return f"{v:.{nd}f}"


def ci(pair, nd: int = 3) -> str:
    if not pair or len(pair) != 2:
        return "—"
    return f"[{f(pair[0], nd)}, {f(pair[1], nd)}]"


# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    cfg = cfgmod.load()
    res = cfg.p("results")
    src = Source(res)
    out_path = Path(args.out) if args.out else (res / "RESULTS.md")

    L: list[str] = []
    checks: list[tuple[str, str, str]] = []     # (status, name, detail)

    def h(level: int, text: str) -> None:
        L.append(("#" * level) + " " + text)
        L.append("")

    def p(text: str = "") -> None:
        L.append(text)

    L.append("# PmDark — consolidated results")
    L.append("")
    p(f"Generated {datetime.now(tz=timezone.utc).isoformat(timespec='minutes')} "
      f"by `scripts/42_consolidate.py`. Every number below is read from a "
      f"results file; none is written into this script. Re-run after any "
      f"analysis change.")
    p()

    # ---- load everything ----------------------------------------------
    boot = src.json("bootstrap_ci.json")
    hardening = src.json("orphan_hardening.json")
    base_d = src.json("homology_baseline_DIAMOND.json")
    base_j = src.json("homology_baseline_jackhmmer.json")
    cmp_d = src.json("plm_vs_homology_DIAMOND.json")
    cmp_j = src.json("plm_vs_homology_jackhmmer.json")
    de = src.json("de_stat_summary.json")
    match = src.json("matched_enrichment.json")
    annot = src.json("annotation_summary.json")
    calibs = src.glob_json("calibration_*.json")
    targets = src.lines("stat_dependent_uncharacterised.txt")

    # =====================================================================
    h(2, "1. Method — anisotropy correction")
    if not boot:
        p("**MISSING** `bootstrap_ci.json` — run `scripts/36_bootstrap_baseline.py`.")
        p()
    else:
        p(f"Paired bootstrap, {boot.get('replicates')} replicates, "
          f"seed {boot.get('seed')}, alpha {boot.get('alpha')}. "
          f"Source: `bootstrap_ci.json`.")
        p()
        p("Per-cell Fmax with 95% CI:")
        p()
        p("| model / transform / tier | n | Fmax | 95% CI | coverage |")
        p("|---|---|---|---|---|")
        for k, v in sorted(boot.get("cells", {}).items()):
            p(f"| {k} | {v.get('n')} | {f(v.get('fmax'))} | "
              f"{ci(v.get('fmax_ci'))} | {f(v.get('coverage'), 2)} |")
        p()
        p("Paired differences (same resampled queries in both arms):")
        p()
        p("| comparison | n | Δ Fmax | 95% CI | p | verdict |")
        p("|---|---|---|---|---|---|")
        for k, v in boot.get("deltas", {}).items():
            verdict = "**excludes 0**" if v.get("excludes_zero") else "crosses 0"
            p(f"| {k} | {v.get('n')} | {f(v.get('delta'))} | "
              f"{ci(v.get('ci'))} | {v.get('p_display', f(v.get('p_two_sided')))} "
              f"| {verdict} |")
        p()

    # =====================================================================
    h(2, "2. Calibration grid")
    if not calibs:
        p("**MISSING** `calibration_*.json`.")
        p()
    else:
        p(f"Sources: {', '.join('`' + k + '`' for k in sorted(calibs))}.")
        p()
        tiers: list[str] = []
        for d in calibs.values():
            for t in d.get("results", {}):
                if t not in tiers:
                    tiers.append(t)
        p("| model | transform | " + " | ".join(tiers) + " |")
        p("|---" * (2 + len(tiers)) + "|")
        for name, d in sorted(calibs.items()):
            row = [d.get("model", "?"),
                   d.get("params", {}).get("transform", "?")]
            for t in tiers:
                blk = d.get("results", {}).get(t)
                row.append(f(blk["metrics"]["overall"]["fmax"]) if blk else "—")
            p("| " + " | ".join(row) + " |")
        p()
        prm = next(iter(calibs.values())).get("params", {})
        p(f"Frozen hyperparameters: k={prm.get('k')}, "
          f"sim_floor={prm.get('sim_floor')}, power={prm.get('power')}. "
          f"`sim_floor` is an absolute cosine cut, so it must be identical "
          f"across transforms or the comparison is a units error.")
        p()

    # =====================================================================
    h(2, "3. Orphan set, hardened")
    if not hardening:
        p("**MISSING** `orphan_hardening.json` — run `scripts/37_harden_orphans.py`.")
        p()
    else:
        p(f"Orphan defined as no DIAMOND hit at e-value "
          f"{hardening.get('diamond_evalue')} under `--ultra-sensitive`, then "
          f"hardened against jackhmmer "
          f"({hardening.get('jackhmmer', {}).get('iterations')} iterations). "
          f"Source: `orphan_hardening.json`.")
        p()
        p("| jackhmmer E-value cut | with a hit | still orphan |")
        p("|---|---|---|")
        for cut, v in hardening.get("survival", {}).items():
            p(f"| {cut} | {v.get('with_hit')} | {v.get('still_orphan')} |")
        p()
        p(f"Hardening threshold applied: E ≤ {hardening.get('harden_at')}, "
          f"leaving {hardening.get('n_survivors')} of "
          f"{hardening.get('n_orphan_diamond')}.")
        p()
        p("| stratum | n | Fmax | 95% CI |")
        p("|---|---|---|---|")
        for lab, v in hardening.get("strata", {}).items():
            p(f"| {lab} | {v.get('n')} | {f(v.get('fmax'))} | "
              f"{ci(v.get('ci'))} |")
        p()

    # =====================================================================
    h(2, "4. Homology baselines")
    for tag, b in (("DIAMOND", base_d), ("jackhmmer", base_j)):
        if not b:
            p(f"**MISSING** `homology_baseline_{tag}.json`.")
            p()
            continue
        p(f"**{tag}** — transfer threshold "
          f"`{b.get('transfer_evalue', 'unfiltered')}`, k={b.get('k')}. "
          f"Source: `homology_baseline_{tag}.json`.")
        p()
        p("| tier | best-hit Fmax | top-k Fmax | answerable |")
        p("|---|---|---|---|")
        for tier, row in b.get("tiers", {}).items():
            bb, tt = row.get("best", {}), row.get("topk", {})
            ans = (f"{bb.get('answerable')}/{bb.get('n_scored')} "
                   f"({100 * bb['answerable_frac']:.1f}%)"
                   if "answerable_frac" in bb else "—")
            p(f"| {tier} | {f(bb.get('fmax'))} | {f(tt.get('fmax'))} | {ans} |")
        und = b.get("undefined_by_construction") or []
        if und:
            p()
            p(f"Undefined by construction (the baseline cannot answer, so "
              f"Fmax there is not an accuracy measurement): "
              f"{', '.join('`' + x + '`' for x in und)}.")
        p()

    # consistency: both baselines must share a transfer threshold
    if base_d and base_j:
        a = base_d.get("transfer_evalue", "unfiltered")
        bq = base_j.get("transfer_evalue", "unfiltered")
        checks.append(("PASS" if a == bq else "FAIL",
                       "baselines share a transfer threshold",
                       f"DIAMOND {a}, jackhmmer {bq}"))

    # =====================================================================
    h(2, "5. PLM vs homology, decomposed")
    p("`all` mixes coverage and accuracy and must not be quoted as an "
      "accuracy difference. `hit` is the accuracy comparison, conditional "
      "on the baseline finding a hit. `nohit` is the PLM's standalone Fmax "
      "where the baseline found nothing — the coverage claim.")
    p()
    for tag, c in (("DIAMOND", cmp_d), ("jackhmmer", cmp_j)):
        if not c:
            p(f"**MISSING** `plm_vs_homology_{tag}.json`.")
            p()
            continue
        p(f"**vs {tag}** (transform {c.get('transform')}, "
          f"baseline mode {c.get('baseline_mode')}). "
          f"Source: `plm_vs_homology_{tag}.json`.")
        p()
        p("| tier | model | all Δ | all CI | hit Δ | hit CI | hit verdict | "
          "nohit n | nohit PLM Fmax |")
        p("|---|---|---|---|---|---|---|---|---|")
        for tier, row in c.get("tiers", {}).items():
            for model, m in (row.get("models") or {}).items():
                al, ac = m.get("all_queries") or {}, m.get("diamond_answerable")
                ms = m.get("diamond_unanswerable")
                hv = ("**excludes 0**" if ac and ac.get("excludes_zero")
                      else "crosses 0" if ac else "—")
                p(f"| {tier} | {model} | {f(al.get('delta'))} | "
                  f"{ci(al.get('ci'))} | "
                  f"{f(ac.get('delta')) if ac else '—'} | "
                  f"{ci(ac.get('ci')) if ac else '—'} | {hv} | "
                  f"{ms.get('n') if ms else '—'} | "
                  f"{f(ms.get('fmax')) if ms else '—'} |")
        sk = c.get("skipped") or []
        if sk:
            p()
            p("Excluded from the accuracy comparison: " + "; ".join(sk) + ".")
        p()

    # =====================================================================
    h(2, "6. Differential expression — STAT dependence")
    if not de:
        p("**MISSING** `de_stat_summary.json` — run `scripts/39_de_stat.py`.")
        p()
    else:
        con = de.get("contrast", {})
        p(f"Contrast **{con.get('a')} vs {con.get('b')}** — the only clean "
          f"single-variable axis in the dataset. Fitted coefficient "
          f"`{con.get('fitted_coeff')}`"
          + (", signs flipped to keep positive = STAT-dependent"
             if con.get("sign_flipped") else "")
          + f". Sign guard corr(raw, shrunk) = "
            f"{f(con.get('sign_guard_corr'))}. Source: `de_stat_summary.json`.")
        p()
        libs = de.get("libraries", {})
        excl = ", ".join(f"{x['run']} ({x['note']})"
                         for x in de.get("excluded", []))
        p(f"Libraries: {len(libs)} used"
          + (f"; excluded {excl}" if excl else ""))
        p()
        ctrl = de.get("positive_control")
        if ctrl:
            p(f"**Positive control** — `{ctrl.get('gene')}`, "
              f"{ctrl.get('product', '')}: log2FC {f(ctrl.get('log2FoldChange'), 2)}, "
              f"padj {ctrl.get('padj'):.3g}."
              if ctrl.get("padj") is not None else
              f"**Positive control** — `{ctrl.get('gene')}`: "
              f"log2FC {f(ctrl.get('log2FoldChange'), 2)}, padj —.")
            ok = (ctrl.get("padj") is not None and ctrl["padj"] < 0.05
                  and (ctrl.get("log2FoldChange") or 0) > 0)
            checks.append(("PASS" if ok else "FAIL",
                           "knockdown confirmed (STAT suppressed, significant)",
                           f"log2FC {f(ctrl.get('log2FoldChange'), 2)}, "
                           f"padj {ctrl.get('padj')}"))
        else:
            p("**Positive control: ABSENT.** The knockdown is unverified and "
              "every result below assumes it worked.")
            checks.append(("FAIL", "knockdown confirmed", "no control recorded"))
        p()
        pr, st = de.get("primary", {}), de.get("stringent", {})
        p(f"Pre-registered primary (apeglm shrinkage, no LFC threshold, "
          f"BH < {de.get('alpha')}): **{pr.get('n_sig')}** genes — "
          f"{pr.get('n_stat_dependent')} STAT-dependent, "
          f"{pr.get('n_stat_repressed')} STAT-repressed.")
        p()
        p(f"Stringent (|LFC| > {st.get('lfc')}), for concordance with the "
          f"published analysis: **{st.get('n_sig')}** genes — "
          f"{st.get('n_up_a')} up in {con.get('a')}, "
          f"{st.get('n_up_b')} up in {con.get('b')}. "
          f"Published: {st.get('published_up_a')} and "
          f"{st.get('published_up_b')}.")
        p()
        p(f"Counts are coverage-derived (prepDE-style, read length "
          f"{de.get('read_length')}) because the BAMs were not retained — "
          f"state this in the methods. "
          f"{de.get('n_genes_tested'):,} of "
          f"{de.get('n_genes_quantified'):,} genes tested.")
        p()

    # =====================================================================
    h(2, "7. Is the uncharacterised fraction STAT-dependent?")
    if not match:
        p("**MISSING** `matched_enrichment.json` — run "
          "`scripts/40_matched_enrichment.py`.")
        p()
    else:
        cf = match.get("confound", {})
        p(f"Source: `matched_enrichment.json`. "
          f"{match.get('n_genes'):,} protein-coding genes "
          f"({match.get('n_uncharacterised'):,} uncharacterised, "
          f"{match.get('n_described'):,} described).")
        p()
        p(f"Confound: standardised mean difference "
          f"{f(cf.get('smd_log_basemean'))} in log10 baseMean and "
          f"{f(cf.get('smd_log_length'))} in log10 length. Median baseMean "
          f"{f(cf.get('median_basemean_unc'), 1)} vs "
          f"{f(cf.get('median_basemean_desc'), 1)}.")
        p()
        raw, mh, mt = (match.get("raw", {}), match.get("stratified", {}),
                       match.get("matched", {}))
        p("| analysis | estimate | interval / p |")
        p("|---|---|---|")
        p(f"| raw (uncontrolled) | OR {f(raw.get('or'), 2)} | "
          f"p = {raw.get('p'):.3g} |" if raw.get("p") is not None else
          f"| raw | OR {f(raw.get('or'), 2)} | — |")
        p(f"| stratified (Mantel–Haenszel) | OR {f(mh.get('or'), 2)} | "
          f"95% CI {ci(mh.get('ci'), 2)} |")
        if mt.get("or") is not None:
            p(f"| matched (caliper, {mt.get('n')} pairs) | "
              f"OR {f(mt.get('or'), 2)} | permutation p = "
              f"{mt.get('perm_p'):.4f} |")
            p(f"| matched balance | SMD {f(mt.get('smd_basemean'))} | "
              f"{mt.get('dropped')} dropped |")
        fc = match.get("fold_change", {})
        p(f"| fold-change distributions | Cliff's δ "
          f"{f(fc.get('cliffs_delta'))} | p = {fc.get('p'):.3g} |"
          if fc.get("p") is not None else "")
        p()
        tr = match.get("trend")
        if tr:
            p(f"Trend of |δ| with expression: Spearman ρ "
              f"{f(tr.get('spearman_rho'), 2)}, p = {tr.get('p'):.3g}. "
              f"Mean |δ| in the three lowest deciles "
              f"{f(tr.get('mean_abs_delta_low3'))}, three highest "
              f"{f(tr.get('mean_abs_delta_high3'))}.")
            p()
            grows = (tr.get("mean_abs_delta_high3") or 0) > \
                (tr.get("mean_abs_delta_low3") or 0)
            p("The effect **grows** with expression. A power deficit or "
              "apeglm shrinkage would do the opposite, so the confound "
              "cannot account for it."
              if grows else
              "The effect **shrinks** with expression, which is what a "
              "power deficit produces. Treat it as confounded.")
            p()
            checks.append(("PASS" if grows else "WARN",
                           "depletion is not a power artifact",
                           f"low3 {f(tr.get('mean_abs_delta_low3'))}, "
                           f"high3 {f(tr.get('mean_abs_delta_high3'))}"))
        if match.get("interaction"):
            p("**Interaction with expression.** No single pooled odds ratio "
              "is the right summary: matching reweights the described class "
              "onto the uncharacterised expression distribution, which is "
              "skewed low, so it measures the effect where there is least "
              "of one. Report the decile table.")
            p()
        p("Per-decile fold-change comparison:")
        p()
        p("| decile | Cliff's δ | p | n unc / desc |")
        p("|---|---|---|---|")
        for w in fc.get("within_deciles", []):
            p(f"| {w.get('decile')} | {f(w.get('delta'))} | "
              f"{w.get('p'):.3g} | {w.get('n_unc')} / {w.get('n_desc')} |")
        p()

    # =====================================================================
    h(2, "8. Annotation of the target set")
    if not annot:
        p("**MISSING** `annotation_summary.json` — run "
          "`scripts/41_annotate_targets.py`.")
        p()
    else:
        p(f"Source: `annotation_summary.json`. Model "
          f"{annot.get('model')} + {annot.get('transform')}, "
          f"k={annot.get('params', {}).get('k')}, "
          f"sim_floor={annot.get('params', {}).get('sim_floor')}, "
          f"power={annot.get('params', {}).get('power')}.")
        p()
        p(f"{annot.get('n_targets')} targets "
          f"({annot.get('n_targets_requested')} requested, "
          f"{annot.get('n_without_protein')} without a coding isoform). "
          f"{annot.get('n_orphan')} have no reachable homolog "
          f"(DIAMOND {annot.get('diamond_hits')} hits, "
          f"jackhmmer {annot.get('jackhmmer_hits')}).")
        p()
        y = annot.get("yield", {})
        if y:
            n = annot.get("n_targets") or 1
            p("**Useful yield** — quote this, not the 'any term' figure. "
              "A gene with only cellular-component terms has been "
              "localised, not functionally annotated.")
            p()
            p("| assignment | n | % of targets |")
            p("|---|---|---|")
            for lab, key in (("molecular function", "molecular_function"),
                             ("biological process", "biological_process"),
                             ("cellular component only",
                              "cellular_component_only")):
                v = y.get(key)
                p(f"| {lab} | {v} | {100 * v / n:.1f}% |"
                  if v is not None else f"| {lab} | — | — |")
            p()
        p("Calibration applied per protein:")
        p()
        p("| stratum | source | n (benchmark) | Fmax | threshold | "
          "precision | recall | BP | MF | CC |")
        p("|---|---|---|---|---|---|---|---|---|---|")
        for lab, c in (annot.get("calibration") or {}).items():
            p(f"| {lab} | {c.get('source', 'calibration file')} | "
              f"{c.get('n_benchmark')} | {f(c.get('fmax'))} | "
              f"{f(c.get('threshold'), 2)} | {f(c.get('precision'))} | "
              f"{f(c.get('recall'))} | {f(c.get('BP'))} | {f(c.get('MF'))} | "
              f"{f(c.get('CC'))} |")
        p()
        tf = annot.get("transferability", {})
        if tf.get("calibration_median_top_sim") is not None:
            shift = (tf["target_median_top_sim"]
                     - tf["calibration_median_top_sim"])
            p(f"Transferability: targets' median top similarity "
              f"{f(tf['target_median_top_sim'])} against the "
              f"{tf.get('calibration_stratum')} calibration set's "
              f"{f(tf['calibration_median_top_sim'])} "
              f"(shift {shift:+.3f}). "
              + ("The calibration is conservative for these proteins."
                 if shift > 0.05 else
                 "The calibration is optimistic for these proteins; treat "
                 "the numbers as an upper bound." if shift < -0.05 else
                 "Comparable; the calibration transfers."))
            p()
            checks.append(("PASS" if shift > -0.05 else "WARN",
                           "calibration transfers to the targets",
                           f"shift {shift:+.3f}"))
        # the orphan calibration must be the hardened one
        orph = (annot.get("calibration") or {}).get("orphan", {})
        checks.append(("PASS" if "hardened" in str(orph.get("source", ""))
                       else "FAIL",
                       "hardened orphan calibration applied",
                       str(orph.get("source", "un-hardened"))))

    # target-count consistency
    if annot and targets is not None:
        checks.append(("PASS" if annot.get("n_targets_requested") == len(targets)
                       else "WARN",
                       "annotation covers the whole target list",
                       f"{annot.get('n_targets_requested')} requested vs "
                       f"{len(targets)} in the list"))

    # orphan n consistency between calibration and hardening
    if hardening and calibs:
        for name, d in calibs.items():
            blk = d.get("results", {}).get("ORPHAN_nohit")
            if blk and d.get("model") == "esm2":
                n_cal = blk.get("n")
                n_h = hardening.get("n_orphan_diamond")
                checks.append(("PASS" if n_cal == n_h else "WARN",
                               "orphan set size agrees across files",
                               f"calibration {n_cal}, hardening {n_h}"))
                break

    # =====================================================================
    h(2, "9. Claims that are NOT supported")
    p("Read from the results files, not from memory. Each of these was "
      "tested and failed to reach a conclusion, and none may be written "
      "as an effect.")
    p()
    any_neg = False
    if boot:
        crossing = [k for k, v in boot.get("deltas", {}).items()
                    if not v.get("excludes_zero")]
        if crossing:
            any_neg = True
            p("**Paired differences whose interval crosses zero:**")
            p()
            for k in crossing:
                v = boot["deltas"][k]
                p(f"- `{k}` — Δ {f(v.get('delta'))} {ci(v.get('ci'))}")
            p()
    for tag, c in (("DIAMOND", cmp_d), ("jackhmmer", cmp_j)):
        if not c:
            continue
        rows = []
        for tier, row in c.get("tiers", {}).items():
            for model, m in (row.get("models") or {}).items():
                ac = m.get("diamond_answerable")
                if ac and not ac.get("excludes_zero"):
                    rows.append(f"- `{tier}` / {model} — accuracy Δ vs {tag} "
                                f"{f(ac.get('delta'))} {ci(ac.get('ci'))}")
        if rows:
            any_neg = True
            p(f"**No accuracy advantage over {tag} at these tiers:**")
            p()
            L.extend(rows)
            p()
    for tag, b in (("DIAMOND", base_d), ("jackhmmer", base_j)):
        und = (b or {}).get("undefined_by_construction") or []
        if und:
            any_neg = True
            p(f"**Undefined by construction for {tag}** (coverage claims "
              f"only, never accuracy): "
              f"{', '.join('`' + x + '`' for x in und)}.")
            p()
    if match and match.get("interaction"):
        any_neg = True
        p("**No single enrichment estimate** for the uncharacterised "
          "fraction: the effect interacts with expression, so the pooled "
          "and matched odds ratios answer different questions. Report the "
          "decile table.")
        p()
    if hardening:
        any_neg = True
        p(f"**Structural homology unchecked.** The orphan set is hardened "
          f"against sequence search only (DIAMOND e-value "
          f"{hardening.get('diamond_evalue')}, jackhmmer "
          f"{hardening.get('jackhmmer', {}).get('iterations')} iterations). "
          f"Foldseek needs predicted structures that do not exist. "
          f"Limitations, not implied.")
        p()
    if not any_neg:
        p("_(nothing recorded — suspicious; check that the results files "
          "were read)_")
        p()

    # =====================================================================
    h(2, "10. Consistency checks")
    if not checks:
        p("_(no checks ran — results files missing)_")
        p()
    else:
        p("| check | status | detail |")
        p("|---|---|---|")
        for status, name, detail in checks:
            mark = {"PASS": "PASS", "FAIL": "**FAIL**",
                    "WARN": "_WARN_"}.get(status, status)
            p(f"| {name} | {mark} | {detail} |")
        p()

    # =====================================================================
    h(2, "Provenance")
    p("| file | last written (UTC) | bytes |")
    p("|---|---|---|")
    for name, meta in sorted(src.seen.items()):
        p(f"| `{name}` | {meta['mtime']} | {meta['bytes']:,} |")
    p()
    if src.missing:
        p("**Missing or unreadable:** "
          + ", ".join("`" + m + "`" for m in src.missing))
        p()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(L) + "\n")
    log(f"wrote {out_path} ({out_path.stat().st_size / 1000:.1f} kB)")

    n_fail = sum(1 for s, _, _ in checks if s == "FAIL")
    n_warn = sum(1 for s, _, _ in checks if s == "WARN")
    print(f"\n  {len(checks)} consistency checks: "
          f"{len(checks) - n_fail - n_warn} pass, {n_warn} warn, "
          f"{n_fail} FAIL")
    for s, name, detail in checks:
        if s != "PASS":
            print(f"    [{s}] {name} — {detail}")
    if src.missing:
        print(f"  {len(src.missing)} results file(s) missing; the document "
              f"says so in place:")
        for m in src.missing:
            print(f"    {m}")
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
