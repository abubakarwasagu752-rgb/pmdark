# PmDark — consolidated results

Generated 2026-11-09T03:40+00:00 by `scripts/42_consolidate.py`. Every number below is read from a results file; none is written into this script. Re-run after any analysis change.

## 1. Method — anisotropy correction

Paired bootstrap, 2000 replicates, seed 0, alpha 0.05. Source: `bootstrap_ci.json`.

Per-cell Fmax with 95% CI:

| model / transform / tier | n | Fmax | 95% CI | coverage |
|---|---|---|---|---|
| esm2/center/I30 | 1087 | 0.656 | [0.639, 0.674] | 0.98 |
| esm2/center/I50 | 1087 | 0.690 | [0.674, 0.709] | 0.98 |
| esm2/center/ORPHAN_nohit | 94 | 0.589 | [0.519, 0.665] | 0.93 |
| esm2/center/T0_full | 1087 | 0.727 | [0.711, 0.745] | 1.00 |
| esm2/none/I30 | 1087 | 0.648 | [0.629, 0.667] | 0.98 |
| esm2/none/I50 | 1087 | 0.681 | [0.663, 0.699] | 0.98 |
| esm2/none/ORPHAN_nohit | 94 | 0.580 | [0.515, 0.657] | 0.94 |
| esm2/none/T0_full | 1087 | 0.708 | [0.690, 0.726] | 0.99 |
| esm2/whiten/I30 | 1087 | 0.721 | [0.705, 0.737] | 0.99 |
| esm2/whiten/I50 | 1087 | 0.775 | [0.762, 0.791] | 0.99 |
| esm2/whiten/ORPHAN_nohit | 94 | 0.664 | [0.599, 0.732] | 0.93 |
| esm2/whiten/T0_full | 1087 | 0.815 | [0.801, 0.829] | 0.99 |
| prostt5/center/I30 | 1087 | 0.663 | [0.646, 0.681] | 0.98 |
| prostt5/center/I50 | 1087 | 0.695 | [0.678, 0.713] | 0.97 |
| prostt5/center/ORPHAN_nohit | 94 | 0.595 | [0.530, 0.675] | 0.93 |
| prostt5/center/T0_full | 1087 | 0.729 | [0.713, 0.747] | 0.98 |
| prostt5/none/I30 | 1087 | 0.662 | [0.644, 0.680] | 0.99 |
| prostt5/none/I50 | 1087 | 0.689 | [0.671, 0.708] | 0.97 |
| prostt5/none/ORPHAN_nohit | 94 | 0.609 | [0.529, 0.684] | 0.86 |
| prostt5/none/T0_full | 1087 | 0.719 | [0.701, 0.738] | 0.97 |
| prostt5/whiten/I30 | 1087 | 0.705 | [0.689, 0.722] | 1.00 |
| prostt5/whiten/I50 | 1087 | 0.763 | [0.748, 0.779] | 0.98 |
| prostt5/whiten/ORPHAN_nohit | 94 | 0.634 | [0.572, 0.705] | 0.88 |
| prostt5/whiten/T0_full | 1087 | 0.806 | [0.792, 0.821] | 0.98 |
| prott5/center/I30 | 1087 | 0.677 | [0.661, 0.695] | 0.99 |
| prott5/center/I50 | 1087 | 0.714 | [0.699, 0.731] | 1.00 |
| prott5/center/ORPHAN_nohit | 94 | 0.611 | [0.553, 0.686] | 1.00 |
| prott5/center/T0_full | 1087 | 0.751 | [0.735, 0.767] | 1.00 |
| prott5/none/I30 | 1087 | 0.668 | [0.653, 0.686] | 0.99 |
| prott5/none/I50 | 1087 | 0.700 | [0.685, 0.718] | 1.00 |
| prott5/none/ORPHAN_nohit | 94 | 0.617 | [0.552, 0.686] | 0.97 |
| prott5/none/T0_full | 1087 | 0.733 | [0.718, 0.750] | 0.99 |
| prott5/whiten/I30 | 1087 | 0.710 | [0.694, 0.727] | 0.99 |
| prott5/whiten/I50 | 1087 | 0.770 | [0.756, 0.786] | 0.99 |
| prott5/whiten/ORPHAN_nohit | 94 | 0.679 | [0.615, 0.744] | 0.97 |
| prott5/whiten/T0_full | 1087 | 0.810 | [0.797, 0.825] | 0.99 |

Paired differences (same resampled queries in both arms):

| comparison | n | Δ Fmax | 95% CI | p | verdict |
|---|---|---|---|---|---|
| esm2/ORPHAN_nohit: center - none | 94 | 0.009 | [-0.011, 0.022] | 0.385 | crosses 0 |
| esm2/ORPHAN_nohit: whiten - none | 94 | 0.084 | [0.038, 0.128] | <0.0010 | **excludes 0** |
| prott5/ORPHAN_nohit: center - none | 94 | -0.006 | [-0.013, 0.014] | 0.913 | crosses 0 |
| prott5/ORPHAN_nohit: whiten - none | 94 | 0.063 | [0.029, 0.098] | <0.0010 | **excludes 0** |
| prostt5/ORPHAN_nohit: center - none | 94 | -0.013 | [-0.022, 0.011] | 0.289 | crosses 0 |
| prostt5/ORPHAN_nohit: whiten - none | 94 | 0.026 | [-0.016, 0.078] | 0.196 | crosses 0 |
| esm2/T0_full: center - none | 1087 | 0.020 | [0.014, 0.029] | <0.0010 | **excludes 0** |
| esm2/T0_full: whiten - none | 1087 | 0.108 | [0.096, 0.120] | <0.0010 | **excludes 0** |
| prott5/T0_full: center - none | 1087 | 0.018 | [0.012, 0.022] | <0.0010 | **excludes 0** |
| prott5/T0_full: whiten - none | 1087 | 0.078 | [0.067, 0.086] | <0.0010 | **excludes 0** |
| prostt5/T0_full: center - none | 1087 | 0.010 | [0.005, 0.018] | <0.0010 | **excludes 0** |
| prostt5/T0_full: whiten - none | 1087 | 0.087 | [0.076, 0.099] | <0.0010 | **excludes 0** |
| esm2/I50: center - none | 1087 | 0.009 | [0.004, 0.017] | <0.0010 | **excludes 0** |
| esm2/I50: whiten - none | 1087 | 0.094 | [0.084, 0.107] | <0.0010 | **excludes 0** |
| prott5/I50: center - none | 1087 | 0.014 | [0.009, 0.019] | <0.0010 | **excludes 0** |
| prott5/I50: whiten - none | 1087 | 0.070 | [0.059, 0.080] | <0.0010 | **excludes 0** |
| prostt5/I50: center - none | 1087 | 0.006 | [0.001, 0.013] | 0.029 | **excludes 0** |
| prostt5/I50: whiten - none | 1087 | 0.074 | [0.063, 0.087] | <0.0010 | **excludes 0** |
| esm2/I30: center - none | 1087 | 0.008 | [0.002, 0.015] | 0.011 | **excludes 0** |
| esm2/I30: whiten - none | 1087 | 0.073 | [0.061, 0.085] | <0.0010 | **excludes 0** |
| prott5/I30: center - none | 1087 | 0.009 | [0.005, 0.013] | <0.0010 | **excludes 0** |
| prott5/I30: whiten - none | 1087 | 0.042 | [0.030, 0.051] | <0.0010 | **excludes 0** |
| prostt5/I30: center - none | 1087 | 0.001 | [-0.005, 0.007] | 0.717 | crosses 0 |
| prostt5/I30: whiten - none | 1087 | 0.043 | [0.032, 0.055] | <0.0010 | **excludes 0** |
| none/ORPHAN_nohit: prott5 - esm2 | 94 | 0.037 | [0.001, 0.073] | 0.045 | **excludes 0** |
| center/ORPHAN_nohit: prott5 - esm2 | 94 | 0.022 | [-0.008, 0.069] | 0.146 | crosses 0 |
| whiten/ORPHAN_nohit: prott5 - esm2 | 94 | 0.015 | [-0.007, 0.041] | 0.177 | crosses 0 |
| none/T0_full: prott5 - esm2 | 1087 | 0.025 | [0.019, 0.036] | <0.0010 | **excludes 0** |
| center/T0_full: prott5 - esm2 | 1087 | 0.024 | [0.015, 0.031] | <0.0010 | **excludes 0** |
| whiten/T0_full: prott5 - esm2 | 1087 | -0.005 | [-0.009, 0.001] | 0.100 | crosses 0 |
| none/I50: prott5 - esm2 | 1087 | 0.019 | [0.012, 0.029] | <0.0010 | **excludes 0** |
| center/I50: prott5 - esm2 | 1087 | 0.024 | [0.015, 0.032] | <0.0010 | **excludes 0** |
| whiten/I50: prott5 - esm2 | 1087 | -0.005 | [-0.013, 0.000] | 0.075 | crosses 0 |
| none/I30: prott5 - esm2 | 1087 | 0.020 | [0.012, 0.031] | <0.0010 | **excludes 0** |
| center/I30: prott5 - esm2 | 1087 | 0.021 | [0.012, 0.030] | <0.0010 | **excludes 0** |
| whiten/I30: prott5 - esm2 | 1087 | -0.011 | [-0.018, -0.002] | 0.009 | **excludes 0** |

## 2. Calibration grid

Sources: `calibration_esm2.json`, `calibration_esm2_center.json`, `calibration_esm2_none.json`, `calibration_esm2_whiten.json`, `calibration_prostt5_center.json`, `calibration_prostt5_none.json`, `calibration_prostt5_whiten.json`, `calibration_prott5_center.json`, `calibration_prott5_none.json`, `calibration_prott5_whiten.json`.

| model | transform | T0_full | I90 | I50 | I30 | Inone | ORPHAN_nohit |
|---|---|---|---|---|---|---|---|
| esm2 | ? | 0.708 | 0.702 | 0.681 | 0.648 | 0.632 | 0.580 |
| esm2 | center | 0.727 | 0.720 | 0.690 | 0.656 | 0.640 | 0.589 |
| esm2 | none | 0.708 | 0.702 | 0.681 | 0.648 | 0.632 | 0.580 |
| esm2 | whiten | 0.815 | 0.812 | 0.775 | 0.721 | 0.697 | 0.664 |
| prostt5 | center | 0.729 | 0.725 | 0.695 | 0.663 | 0.653 | 0.595 |
| prostt5 | none | 0.719 | 0.712 | 0.689 | 0.662 | 0.649 | 0.609 |
| prostt5 | whiten | 0.806 | 0.803 | 0.763 | 0.705 | 0.679 | 0.634 |
| prott5 | center | 0.751 | 0.745 | 0.714 | 0.677 | 0.657 | 0.611 |
| prott5 | none | 0.733 | 0.727 | 0.700 | 0.668 | 0.654 | 0.617 |
| prott5 | whiten | 0.810 | 0.808 | 0.770 | 0.710 | 0.688 | 0.679 |

Frozen hyperparameters: k=25, sim_floor=0.5, power=3.0. `sim_floor` is an absolute cosine cut, so it must be identical across transforms or the comparison is a units error.

## 3. Orphan set, hardened

Orphan defined as no DIAMOND hit at e-value 10.0 under `--ultra-sensitive`, then hardened against jackhmmer (5 iterations). Source: `orphan_hardening.json`.

| jackhmmer E-value cut | with a hit | still orphan |
|---|---|---|
| 1e-20 | 18 | 76 |
| 1e-10 | 22 | 72 |
| 1e-05 | 22 | 72 |
| 0.001 | 24 | 70 |
| 0.01 | 30 | 64 |
| 0.1 | 43 | 51 |
| 1 | 50 | 44 |
| 10 | 51 | 43 |

Hardening threshold applied: E ≤ 0.001, leaving 70 of 94.

| stratum | n | Fmax | 95% CI |
|---|---|---|---|
| DIAMOND-orphan | 94 | 0.664 | [0.599, 0.732] |
| HARDENED | 70 | 0.625 | [0.542, 0.706] |
| dropped (had a profile hit) | 24 | 0.734 | [0.647, 0.863] |

## 4. Homology baselines

**DIAMOND** — transfer threshold `0.001`, k=25. Source: `homology_baseline_DIAMOND.json`.

| tier | best-hit Fmax | top-k Fmax | answerable |
|---|---|---|---|
| T0_full | 0.768 | 0.768 | 927/1087 (85.3%) |
| I90 | 0.757 | 0.761 | 916/1087 (84.3%) |
| I50 | 0.667 | 0.661 | 795/1087 (73.1%) |
| I30 | 0.316 | 0.323 | 378/1087 (34.8%) |
| I20 | 0.014 | 0.015 | 25/1087 (2.3%) |
| Inone | 0.000 | 0.000 | 0/1087 (0.0%) |
| ORPHAN_nohit | 0.000 | 0.000 | 0/94 (0.0%) |

Undefined by construction (the baseline cannot answer, so Fmax there is not an accuracy measurement): `I20`, `Inone`, `ORPHAN_nohit`.

**jackhmmer** — transfer threshold `0.001`, k=25. Source: `homology_baseline_jackhmmer.json`.

| tier | best-hit Fmax | top-k Fmax | answerable |
|---|---|---|---|
| T0_full | 0.764 | 0.760 | 974/1087 (89.6%) |
| I90 | 0.757 | 0.748 | 969/1087 (89.1%) |
| I50 | 0.692 | 0.692 | 902/1087 (83.0%) |
| I30 | 0.598 | 0.615 | 794/1087 (73.0%) |
| I20 | 0.548 | 0.563 | 695/1087 (63.9%) |
| Inone | 0.547 | 0.562 | 691/1087 (63.6%) |
| ORPHAN_nohit | 0.251 | 0.271 | 24/94 (25.5%) |

## 5. PLM vs homology, decomposed

`all` mixes coverage and accuracy and must not be quoted as an accuracy difference. `hit` is the accuracy comparison, conditional on the baseline finding a hit. `nohit` is the PLM's standalone Fmax where the baseline found nothing — the coverage claim.

**vs DIAMOND** (transform whiten, baseline mode topk). Source: `plm_vs_homology_DIAMOND.json`.

| tier | model | all Δ | all CI | hit Δ | hit CI | hit verdict | nohit n | nohit PLM Fmax |
|---|---|---|---|---|---|---|---|---|
| T0_full | esm2 | 0.047 | [0.036, 0.058] | 0.005 | [-0.000, 0.011] | crosses 0 | 160 | 0.647 |
| T0_full | prott5 | 0.042 | [0.032, 0.054] | 0.001 | [-0.006, 0.007] | crosses 0 | 160 | 0.648 |
| T0_full | prostt5 | 0.038 | [0.028, 0.049] | -0.001 | [-0.008, 0.004] | crosses 0 | 160 | 0.621 |
| I90 | esm2 | 0.051 | [0.040, 0.063] | 0.002 | [-0.003, 0.008] | crosses 0 | 171 | 0.672 |
| I90 | prott5 | 0.047 | [0.036, 0.059] | -0.001 | [-0.008, 0.005] | crosses 0 | 171 | 0.669 |
| I90 | prostt5 | 0.042 | [0.031, 0.054] | -0.004 | [-0.011, 0.002] | crosses 0 | 171 | 0.645 |
| I50 | esm2 | 0.114 | [0.098, 0.132] | -0.003 | [-0.009, 0.007] | crosses 0 | 292 | 0.738 |
| I50 | prott5 | 0.110 | [0.092, 0.126] | -0.010 | [-0.019, -0.002] | **excludes 0** | 292 | 0.736 |
| I50 | prostt5 | 0.102 | [0.085, 0.120] | -0.013 | [-0.020, -0.003] | **excludes 0** | 292 | 0.719 |
| I30 | esm2 | 0.398 | [0.372, 0.424] | 0.005 | [-0.008, 0.019] | crosses 0 | 709 | 0.763 |
| I30 | prott5 | 0.387 | [0.361, 0.413] | -0.005 | [-0.020, 0.010] | crosses 0 | 709 | 0.750 |
| I30 | prostt5 | 0.381 | [0.356, 0.409] | -0.005 | [-0.015, 0.010] | crosses 0 | 709 | 0.742 |

Excluded from the accuracy comparison: I20 (DIAMOND answers 2.3% — undefined by construction).

**vs jackhmmer** (transform whiten, baseline mode topk). Source: `plm_vs_homology_jackhmmer.json`.

| tier | model | all Δ | all CI | hit Δ | hit CI | hit verdict | nohit n | nohit PLM Fmax |
|---|---|---|---|---|---|---|---|---|
| T0_full | esm2 | 0.056 | [0.044, 0.066] | 0.034 | [0.025, 0.041] | **excludes 0** | 113 | 0.577 |
| T0_full | prott5 | 0.051 | [0.040, 0.061] | 0.029 | [0.020, 0.037] | **excludes 0** | 113 | 0.585 |
| T0_full | prostt5 | 0.046 | [0.036, 0.057] | 0.027 | [0.018, 0.034] | **excludes 0** | 113 | 0.559 |
| I90 | esm2 | 0.064 | [0.052, 0.074] | 0.038 | [0.028, 0.046] | **excludes 0** | 118 | 0.598 |
| I90 | prott5 | 0.060 | [0.048, 0.070] | 0.034 | [0.024, 0.043] | **excludes 0** | 118 | 0.603 |
| I90 | prostt5 | 0.054 | [0.043, 0.064] | 0.031 | [0.021, 0.039] | **excludes 0** | 118 | 0.571 |
| I50 | esm2 | 0.084 | [0.072, 0.098] | 0.030 | [0.022, 0.040] | **excludes 0** | 185 | 0.664 |
| I50 | prott5 | 0.079 | [0.066, 0.092] | 0.025 | [0.015, 0.034] | **excludes 0** | 185 | 0.658 |
| I50 | prostt5 | 0.071 | [0.060, 0.085] | 0.022 | [0.013, 0.032] | **excludes 0** | 185 | 0.627 |
| I30 | esm2 | 0.106 | [0.090, 0.122] | 0.022 | [0.009, 0.032] | **excludes 0** | 293 | 0.642 |
| I30 | prott5 | 0.095 | [0.080, 0.112] | 0.008 | [-0.004, 0.020] | crosses 0 | 293 | 0.639 |
| I30 | prostt5 | 0.089 | [0.075, 0.107] | 0.010 | [-0.002, 0.022] | crosses 0 | 293 | 0.610 |
| I20 | esm2 | 0.135 | [0.117, 0.152] | 0.024 | [0.009, 0.037] | **excludes 0** | 392 | 0.608 |
| I20 | prott5 | 0.125 | [0.108, 0.144] | 0.013 | [-0.002, 0.026] | crosses 0 | 392 | 0.609 |
| I20 | prostt5 | 0.118 | [0.100, 0.136] | 0.008 | [-0.006, 0.022] | crosses 0 | 392 | 0.588 |

## 6. Differential expression — STAT dependence

Contrast **WSSV vs PmSTAT_dsRNA_WSSV** — the only clean single-variable axis in the dataset. Fitted coefficient `group[T.WSSV]`. Sign guard corr(raw, shrunk) = 0.958. Source: `de_stat_summary.json`.

Libraries: 9 used; excluded SRR28572179 (FAILED_LIBRARY_16x_shallow)

**Positive control** — `gene-LOC119587821`, signal transducer and activator of transcription 5B-like: log2FC 1.12, padj 3.77e-08.

Pre-registered primary (apeglm shrinkage, no LFC threshold, BH < 0.05): **5152** genes — 2683 STAT-dependent, 2469 STAT-repressed.

Stringent (|LFC| > 2.0), for concordance with the published analysis: **813** genes — 631 up in WSSV, 182 up in PmSTAT_dsRNA_WSSV. Published: 675 and 138.

Counts are coverage-derived (prepDE-style, read length 150.0) because the BAMs were not retained — state this in the methods. 10,231 of 25,226 genes tested.

## 7. Is the uncharacterised fraction STAT-dependent?

Source: `matched_enrichment.json`. 9,741 protein-coding genes (1,894 uncharacterised, 7,847 described).

Confound: standardised mean difference -0.372 in log10 baseMean and -0.008 in log10 length. Median baseMean 51.5 vs 103.1.

| analysis | estimate | interval / p |
|---|---|---|
| raw (uncontrolled) | OR 0.55 | p = 6.93e-21 |
| stratified (Mantel–Haenszel) | OR 0.72 | 95% CI [0.63, 0.83] |
| matched (caliper, 1888 pairs) | OR 0.80 | permutation p = 0.0040 |
| matched balance | SMD -0.001 | 6 dropped |
| fold-change distributions | Cliff's δ -0.154 | p = 2e-25 |

Trend of |δ| with expression: Spearman ρ 0.84, p = 0.00222. Mean |δ| in the three lowest deciles 0.051, three highest 0.238.

The effect **grows** with expression. A power deficit or apeglm shrinkage would do the opposite, so the confound cannot account for it.

**Interaction with expression.** No single pooled odds ratio is the right summary: matching reweights the described class onto the uncharacterised expression distribution, which is skewed low, so it measures the effect where there is least of one. Report the decile table.

Per-decile fold-change comparison:

| decile | Cliff's δ | p | n unc / desc |
|---|---|---|---|
| 1 | -0.045 | 0.263 | 297 / 677 |
| 2 | -0.040 | 0.334 | 264 / 710 |
| 3 | -0.069 | 0.112 | 232 / 742 |
| 4 | -0.136 | 0.00249 | 211 / 763 |
| 5 | -0.041 | 0.401 | 175 / 799 |
| 6 | -0.105 | 0.032 | 167 / 808 |
| 7 | -0.127 | 0.00983 | 165 / 808 |
| 8 | -0.186 | 0.000336 | 146 / 828 |
| 9 | -0.214 | 5.76e-05 | 137 / 837 |
| 10 | -0.313 | 2.85e-07 | 100 / 875 |

## 8. Annotation of the target set

Source: `annotation_summary.json`. Model esm2 + whiten, k=25, sim_floor=0.0, power=3.0.

328 targets (328 requested, 0 without a coding isoform). 190 have no reachable homolog (DIAMOND 135 hits, jackhmmer 64).

**Useful yield** — quote this, not the 'any term' figure. A gene with only cellular-component terms has been localised, not functionally annotated.

| assignment | n | % of targets |
|---|---|---|
| molecular function | 136 | 41.5% |
| biological process | 127 | 38.7% |
| cellular component only | 149 | 45.4% |

Calibration applied per protein:

| stratum | source | n (benchmark) | Fmax | threshold | precision | recall | BP | MF | CC |
|---|---|---|---|---|---|---|---|---|---|
| orphan | hardened (DIAMOND + jackhmmer) | 70 | 0.625 | 0.30 | 0.862 | 0.491 | 0.324 | 0.382 | 0.767 |
| homolog | calibration file | 1087 | 0.815 | 0.38 | 0.862 | 0.773 | 0.656 | 0.837 | 0.882 |

Transferability: targets' median top similarity 0.375 against the ORPHAN_nohit calibration set's 0.351 (shift +0.023). Comparable; the calibration transfers.

## 9. Claims that are NOT supported

Read from the results files, not from memory. Each of these was tested and failed to reach a conclusion, and none may be written as an effect.

**Paired differences whose interval crosses zero:**

- `esm2/ORPHAN_nohit: center - none` — Δ 0.009 [-0.011, 0.022]
- `prott5/ORPHAN_nohit: center - none` — Δ -0.006 [-0.013, 0.014]
- `prostt5/ORPHAN_nohit: center - none` — Δ -0.013 [-0.022, 0.011]
- `prostt5/ORPHAN_nohit: whiten - none` — Δ 0.026 [-0.016, 0.078]
- `prostt5/I30: center - none` — Δ 0.001 [-0.005, 0.007]
- `center/ORPHAN_nohit: prott5 - esm2` — Δ 0.022 [-0.008, 0.069]
- `whiten/ORPHAN_nohit: prott5 - esm2` — Δ 0.015 [-0.007, 0.041]
- `whiten/T0_full: prott5 - esm2` — Δ -0.005 [-0.009, 0.001]
- `whiten/I50: prott5 - esm2` — Δ -0.005 [-0.013, 0.000]

**No accuracy advantage over DIAMOND at these tiers:**

- `T0_full` / esm2 — accuracy Δ vs DIAMOND 0.005 [-0.000, 0.011]
- `T0_full` / prott5 — accuracy Δ vs DIAMOND 0.001 [-0.006, 0.007]
- `T0_full` / prostt5 — accuracy Δ vs DIAMOND -0.001 [-0.008, 0.004]
- `I90` / esm2 — accuracy Δ vs DIAMOND 0.002 [-0.003, 0.008]
- `I90` / prott5 — accuracy Δ vs DIAMOND -0.001 [-0.008, 0.005]
- `I90` / prostt5 — accuracy Δ vs DIAMOND -0.004 [-0.011, 0.002]
- `I50` / esm2 — accuracy Δ vs DIAMOND -0.003 [-0.009, 0.007]
- `I30` / esm2 — accuracy Δ vs DIAMOND 0.005 [-0.008, 0.019]
- `I30` / prott5 — accuracy Δ vs DIAMOND -0.005 [-0.020, 0.010]
- `I30` / prostt5 — accuracy Δ vs DIAMOND -0.005 [-0.015, 0.010]

**No accuracy advantage over jackhmmer at these tiers:**

- `I30` / prott5 — accuracy Δ vs jackhmmer 0.008 [-0.004, 0.020]
- `I30` / prostt5 — accuracy Δ vs jackhmmer 0.010 [-0.002, 0.022]
- `I20` / prott5 — accuracy Δ vs jackhmmer 0.013 [-0.002, 0.026]
- `I20` / prostt5 — accuracy Δ vs jackhmmer 0.008 [-0.006, 0.022]

**Undefined by construction for DIAMOND** (coverage claims only, never accuracy): `I20`, `Inone`, `ORPHAN_nohit`.

**No single enrichment estimate** for the uncharacterised fraction: the effect interacts with expression, so the pooled and matched odds ratios answer different questions. Report the decile table.

**Structural homology unchecked.** The orphan set is hardened against sequence search only (DIAMOND e-value 10.0, jackhmmer 5 iterations). Foldseek needs predicted structures that do not exist. Limitations, not implied.

## 10. Consistency checks

| check | status | detail |
|---|---|---|
| baselines share a transfer threshold | PASS | DIAMOND 0.001, jackhmmer 0.001 |
| knockdown confirmed (STAT suppressed, significant) | PASS | log2FC 1.12, padj 3.7714681133180874e-08 |
| depletion is not a power artifact | PASS | low3 0.051, high3 0.238 |
| calibration transfers to the targets | PASS | shift +0.023 |
| hardened orphan calibration applied | PASS | hardened (DIAMOND + jackhmmer) |
| annotation covers the whole target list | PASS | 328 requested vs 328 in the list |
| orphan set size agrees across files | PASS | calibration 94, hardening 94 |

## Provenance

| file | last written (UTC) | bytes |
|---|---|---|
| `annotation_summary.json` | 2026-11-09T03:27+00:00 | 1,398 |
| `bootstrap_ci.json` | 2026-11-08T06:37+00:00 | 24,981 |
| `calibration_esm2.json` | 2026-11-07T11:23+00:00 | 8,866 |
| `calibration_esm2_center.json` | 2026-11-08T05:44+00:00 | 8,885 |
| `calibration_esm2_none.json` | 2026-11-08T05:43+00:00 | 8,907 |
| `calibration_esm2_whiten.json` | 2026-11-08T05:45+00:00 | 9,038 |
| `calibration_prostt5_center.json` | 2026-11-08T05:50+00:00 | 8,953 |
| `calibration_prostt5_none.json` | 2026-11-08T05:49+00:00 | 8,901 |
| `calibration_prostt5_whiten.json` | 2026-11-08T05:51+00:00 | 9,043 |
| `calibration_prott5_center.json` | 2026-11-08T05:47+00:00 | 8,887 |
| `calibration_prott5_none.json` | 2026-11-08T05:46+00:00 | 8,938 |
| `calibration_prott5_whiten.json` | 2026-11-08T05:48+00:00 | 8,940 |
| `de_stat_summary.json` | 2026-11-08T13:05+00:00 | 3,045 |
| `homology_baseline_DIAMOND.json` | 2026-11-08T08:09+00:00 | 6,456 |
| `homology_baseline_jackhmmer.json` | 2026-11-08T08:11+00:00 | 6,844 |
| `matched_enrichment.json` | 2026-11-09T03:04+00:00 | 3,137 |
| `orphan_hardening.json` | 2026-11-09T03:27+00:00 | 2,432 |
| `plm_vs_homology_DIAMOND.json` | 2026-11-08T08:10+00:00 | 16,583 |
| `plm_vs_homology_jackhmmer.json` | 2026-11-08T08:12+00:00 | 20,544 |
| `stat_dependent_uncharacterised.txt` | 2026-11-09T03:04+00:00 | 5,904 |

