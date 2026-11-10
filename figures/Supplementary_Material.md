# Supplementary Material

All values are generated from the analysis outputs by `scripts/43_figures.py`; none is transcribed.

## Table S1. Search sensitivity and tier leakage

Queries with any permitted hit under each search, at a common transfer threshold. Conditions where a search can answer almost nothing are undefined for that search by construction and carry no accuracy interpretation.

| condition | DIAMOND | jackhmmer | DIAMOND Fmax | jackhmmer Fmax |
|---|---|---|---|---|
| T0 | 927/1087 (85.3%) | 974/1087 (89.6%) | 0.768 | 0.764 |
| I90 | 916/1087 (84.3%) | 969/1087 (89.1%) | 0.757 | 0.757 |
| I50 | 795/1087 (73.1%) | 902/1087 (83.0%) | 0.667 | 0.692 |
| I30 | 378/1087 (34.8%) | 794/1087 (73.0%) | 0.316 | 0.598 |
| I20 | 25/1087 (2.3%) | 695/1087 (63.9%) | 0.014 | 0.548 |
| Inone | 0/1087 (0.0%) | 691/1087 (63.6%) | 0.000 | 0.547 |
| orphan | 0/94 (0.0%) | 24/94 (25.5%) | 0.000 | 0.251 |

Undefined by construction: I20, Inone, orphan.

## Table S2. Differential expression summary

| quantity | value |
|---|---|
| contrast | WSSV vs PmSTAT_dsRNA_WSSV |
| libraries used | 9 |
| libraries excluded | SRR28572179 |
| genes quantified | 25226 |
| genes tested | 10231 |
| knockdown control gene | gene-LOC119587821 |
| control log2 fold change | 1.12 |
| control adjusted P | 3.77e-08 |
| significant, no fold-change threshold | 5152 |
| of which STAT-dependent | 2683 |
| significant at |log2FC| > 2.0 | 813 |
| previously reported, same threshold | 675 and 138 |

## Table S3. STAT dependence by expression decile

| decile | baseMean range | uncharacterized | described | Cliff's δ | P |
|---|---|---|---|---|---|
| 1 | 2.3–11.9 | 10/297 (3.4%) | 34/677 (5.0%) | -0.045 | 0.263 |
| 2 | 11.9–21.1 | 16/264 (6.1%) | 33/710 (4.6%) | -0.040 | 0.334 |
| 3 | 21.1–35.5 | 35/232 (15.1%) | 84/742 (11.3%) | -0.069 | 0.112 |
| 4 | 35.5–58.9 | 28/211 (13.3%) | 149/763 (19.5%) | -0.136 | 0.00249 |
| 5 | 58.9–91.5 | 42/175 (24.0%) | 201/799 (25.2%) | -0.041 | 0.401 |
| 6 | 91.5–142.5 | 39/167 (23.4%) | 221/808 (27.4%) | -0.105 | 0.032 |
| 7 | 142.5–233.2 | 45/165 (27.3%) | 280/808 (34.7%) | -0.127 | 0.00983 |
| 8 | 233.2–405.9 | 41/146 (28.1%) | 330/828 (39.9%) | -0.186 | 0.000336 |
| 9 | 405.9–956.5 | 44/137 (32.1%) | 366/837 (43.7%) | -0.214 | 5.76e-05 |
| 10 | 956.5–83908.2 | 28/100 (28.0%) | 458/875 (52.3%) | -0.313 | 2.85e-07 |

Confounding: standardized mean difference -0.372 in log10 expression and -0.008 in log10 transcript length. Absolute effect size against expression decile, Spearman ρ = 0.84 (P = 0.00222).

## Table S4. Calibration applied to the target proteins

| subset | source | n (benchmark) | Fmax | score threshold | precision | recall | BP | MF | CC |
|---|---|---|---|---|---|---|---|---|---|
| orphan | hardened (DIAMOND + jackhmmer) | 70 | 0.625 | 0.30 | 0.862 | 0.491 | 0.324 | 0.382 | 0.767 |
| homolog | benchmark | 1087 | 0.815 | 0.38 | 0.862 | 0.773 | 0.656 | 0.837 | 0.882 |

Target proteins' median nearest-neighbour similarity 0.375 against 0.351 for the calibration subset.

Per-protein predictions are provided as `annotation_targets.csv` (one row per protein) and `annotation_targets_long.csv` (one row per protein and term).

