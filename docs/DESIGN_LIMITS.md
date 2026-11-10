# Design limits of the upstream dataset

Read this before writing a single claim. These constraints come from the
source Data Descriptor (Jatuyosporn et al., *Sci Data* 13:829, 2026) and they
bound what this project can honestly conclude. Three of them differ from what
we assumed when sketching the aims.

---

## 1. There is no uninfected control arm

The three RNA-seq conditions are:

| Arm | n | dsRNA | WSSV |
|---|---|---|---|
| WSSV | 5 | PBS | yes |
| PmSTAT dsRNA | 3 | PmSTAT | PBS |
| PmSTAT dsRNA + WSSV | 4 | PmSTAT | yes |

There is **no naive arm and no GFP-dsRNA (non-targeting) arm**. Consequences:

- "WSSV-responsive" cannot be defined from this dataset alone. Any contrast
  involving the WSSV arm and the dsRNA-only arm confounds infection with
  knockdown, because the dsRNA-only arm also lacks the injection-matched
  viral control.
- The **only clean single-variable contrast is WSSV vs PmSTAT dsRNA + WSSV**,
  which isolates the effect of STAT knockdown *under infection*.

This is not fatal — it is in fact lucky, because the one clean axis is exactly
the novel axis the paper is built on. But the aims must be rewritten so that
the STAT-dependence contrast is the primary and the only causal claim, with
"infection-associated" treated as descriptive and clearly caveated.

**Action:** rewrite Aim 4. Drop the 2x2 language. The partition becomes
*STAT-dependent under infection* vs *STAT-independent under infection*, not a
four-cell classification.

---

## 2. The published DEG set is small and the threshold is unusually stringent

Reported with |log2FC| > 2 and BH-FDR < 0.05:

- PmSTAT dsRNA vs PmSTAT dsRNA + WSSV: **zero** significant DEGs
- WSSV vs PmSTAT dsRNA + WSSV: 675 up in WSSV, 138 up in dsRNA+WSSV

And after protein-coding filtering, the usable set collapses further: of the
138 upregulated in dsRNA+WSSV only **47 encode proteins**; of the 349
resampling-robust WSSV-up genes only **88 encode proteins**. That is ~135
protein-coding DEGs in total — too thin to build a regulon on.

An LFC threshold of 2 is very aggressive for n = 3-5. A defensible reanalysis
using apeglm shrinkage, `lfcThreshold = 0`, and independent hypothesis
weighting will recover far more signal at controlled FDR.

**Action:** treat DESeq2 reanalysis as a required, creditable contribution,
not a convenience. Report both the original stringent set (for concordance)
and the shrinkage-based set (for discovery). Pre-register the thresholds in
the methods so it does not read as p-hacking.

---

## 3. RESOLVED — we rebuild the protein set independently

Data Records deposits:

- Raw reads: **SRA SRP499992**
- Pipelines, scripts, and annotation tables: the upstream GitHub repository
  (`LuiguiGallardo/p_monodon_transcriptome`, `Additional_Tables/`)

There is **no TSA accession and no deposited assembly or protein FASTA**, and
the Data Availability statement defers additional data to the corresponding
author on request.

**Decision: we do not request it. We rebuild from SRA.**

This removes the dependency, removes the scooping exposure, and — because a
chromosome-level *P. monodon* reference genome exists — lets us produce a
*better* protein set than the one we would have been given. The upstream
assembly has a 417 bp N50 and 586,928 unigenes, which is a shattered
assembly; genome-guided reconstruction fixes the single largest technical
risk to the project.

Full rationale, pipeline and compute budget: **`docs/ASSEMBLY_STRATEGY.md`**.

Two consequences worth stating up front:

- The rebuild becomes a **stated contribution**, not hidden plumbing. "De
  novo assembly inflates the apparent dark proteome" is itself a publishable,
  generalisable finding.
- We control the ORF cutoff. Upstream used >= 100 aa, which by construction
  excludes most penaeid AMP precursors — exactly the class Aim 3 hunts. We
  use 70 aa and retain the >= 100 aa subset as a comparable stratum.

**Residual risk:** we cannot reproduce their exact 28,423-protein set, so
direct protein-ID-level comparison to their annotation tables is not
possible. Comparison is at the level of assembly statistics, locus counts and
annotation rates (see the validation gate table in `ASSEMBLY_STRATEGY.md`).
Their `Additional_Tables/` remain useful for cross-checking annotation rates
and DEG concordance, and the repo can be cloned without contacting anyone.

---

## 4. "Dark" has three defensible definitions — report all three

The 15,966 "annotated" figure is a **union across evidence types**, and the
per-type coverage is much thinner:

| Evidence | Proteins | % |
|---|---:|---:|
| Swiss-Prot hit | 5,138 | 18.1 |
| GO term | 7,372 | 25.9 |
| eggNOG / COG | 9,313 | 32.8 |
| InterProScan | 13,361 | 47.0 |
| KEGG KO | 6,898 | 24.3 |
| **Union (any)** | **15,966** | **56.2** |

So:

- `dark_union` = no evidence of any kind = **12,457** proteins
- `dark_named` = no Swiss-Prot hit = **23,285** proteins
- `dark_strict` = no Swiss-Prot *and* no InterPro = computed in stage 01

Most of the proteome has no *named* homolog. That is a much stronger framing
than "44% unannotated" and it is the honest one. Lead with `dark_union` as the
conservative headline and report `dark_named` as the scope of the problem.

---

## 5. Assembly fragmentation is the dominant technical risk

Contig N50 = 417 bp; median contig length = 325 bp. TransDecoder retained
ORFs >= 100 aa. A short-ORF-heavy, highly fragmented assembly means a
substantial fraction of "orphans" will be fragments, UTR-frame artifacts, and
chimeras rather than real proteins.

BUSCO on the 1-isoform set was C = 85.2% (S = 84.6%, D = 0.6%), F = 6.9%,
M = 7.9% for Arthropoda — so conserved gene content is well captured, but
that says nothing about the non-conserved fraction, which is precisely our
target.

**Action:** the proteinness filter (stage 02) is not optional, and the
proteogenomic MS evidence (stage 03) is the only hard answer. Run stage 03
early; it determines whether this is a discovery paper or a benchmark paper.

---

## 6. Scooping risk is real and should be managed, not ignored

The source paper explicitly states the resource supports "future studies of
JAK/STAT-dependent antiviral programs and host-WSSV interactions" and
"training and evaluation of machine-learning models". The authoring groups
(Chulalongkorn, UNAM-IBT) have both the data and the motive.

Because the protein FASTA likely requires contacting them anyway, the
pragmatic move is to open that conversation framed as collaboration rather
than silent reuse. Worst case they decline and we re-assemble from SRA; best
case we gain co-authors with wet-lab capacity and remove the biggest risk to
the project.
