# Assembly strategy: why we rebuild, and why genome-guided

The upstream protein FASTA is not deposited anywhere public, and we are not
requesting it. We therefore rebuild the protein set from SRA SRP499992.

This is not a workaround. It is the single biggest upgrade available to the
project, for four reasons set out below.

---

## 1. The upstream assembly is severely fragmented, and that is our main risk

Reported upstream statistics:

| Metric | Value |
|---|---|
| Transcripts | 778,680 |
| Unigenes (longest isoform/gene) | 586,928 |
| Contig N50 | **417 bp** |
| Median contig length | **325 bp** |
| Mean isoform length | 409 bp |
| Predicted proteins (TransDecoder, >=100 aa) | 28,423 |

A 417 bp N50 on a ~2.4 Gb-genome organism means most "transcripts" are
fragments. 586,928 unigenes is roughly 15-20x the expected gene count for a
penaeid, which is the signature of a shattered de novo assembly, not of
biological complexity.

The consequence for us is direct: **an unknown but large fraction of the
12,457 "orphans" are assembly artifacts** — fragments, UTR-frame ORFs, and
chimeras. If we inherit that assembly, the central claim of the paper
("the dark proteome is functionally structured, not noise") is built on
sand, and the first competent reviewer will say so.

We cannot fix this by filtering alone. We fix it at the source.

---

## 2. A chromosome-level reference genome exists

Uengwetwanit et al. 2021, *Molecular Ecology Resources*,
"A chromosome-level assembly of the black tiger shrimp (*Penaeus monodon*)
genome facilitates the identification of growth-associated genes"
(doi:10.1111/1755-0998.13357). The genome is in NCBI and indexed by KEGG as
organism `pmoo`.

The upstream study assembled de novo anyway. With a chromosome-level
reference available, genome-guided reconstruction is simply the better
method: reads are placed on a scaffold rather than stitched from overlaps,
so transcript boundaries, isoform structure and ORF frames are all far more
reliable.

`scripts/02_fetch_refs.py` resolves the current best assembly accession at
runtime via the NCBI datasets API rather than hard-coding one, so the
pipeline stays correct as the assembly is updated.

---

## 3. Genome-guided assembly gives us two analyses for free

**5' flanking regions for STAT motif scanning.** The original plan had a
separate, awkward step: map de novo transcripts back to the genome to recover
promoters. With a genome-guided assembly every transcript already has
coordinates. GAS-motif scanning upstream of STAT-dependent transcripts
becomes a coordinate lookup, not a mapping project.

**Locus-level deduplication.** 586,928 unigenes collapse to real loci. Dark
proteins that are actually fragments of the same gene merge instead of
inflating family counts — which would otherwise produce spurious "novel
families" in Aim 3.

---

## 4. We control the ORF cutoff, and short ORFs are where the biology is

TransDecoder upstream retained ORFs **>= 100 aa**. Penaeid antimicrobial
peptides are routinely shorter: penaeidins, crustins and ALFs mature at
~5-12 kDa, and many precursors sit near or below that cutoff.

So the upstream protein set is, by construction, depleted of exactly the
class of lineage-restricted small secreted effectors that Aim 3 is hunting.
Rebuilding at a 70 aa cutoff (with the >=100 aa subset retained as a
comparable stratum) recovers that space.

This also gives a clean, quotable result: *"lowering the ORF threshold from
100 to 70 aa recovers N candidate small secreted proteins, of which M carry
peptide-level MS evidence."*

---

## Pipeline

```
SRA SRP499992 (12 runs, 179.4 M pairs, 2x75 bp)
        |
        v  fastp  (QC; matches upstream intent, modern tool)
        |
        +--> HISAT2 --> P. monodon reference genome
        |        |
        |        +--> mapped reads   --> StringTie (per sample)
        |        |                         --> StringTie --merge
        |        |                         --> genome-guided transcripts
        |        |
        |        +--> unmapped reads --> rnaSPAdes
        |                                  --> de novo transcripts
        |                                      (captures loci absent from
        |                                       the reference assembly)
        v
   merged transcript set
        |
        v  TransDecoder (--min_protein_length 70)
        |
        v  protein set  ->  three dark-set definitions  ->  stage 02
```

**Why hybrid rather than genome-guided alone.** A reference genome is never
complete, and fast-evolving immune loci are disproportionately likely to be
missing or mis-assembled in it. Assembling the unmapped fraction de novo
keeps those recoverable. The fraction of final dark proteins coming from the
unmapped arm is itself a reportable number.

**Why not Trinity.** Trinity on 179 M read pairs needs roughly 1 GB RAM per
million reads and commonly runs for days. The `cpu` partition caps at 24 h
per job. HISAT2 + StringTie is dramatically cheaper and, with a reference
available, more accurate. rnaSPAdes on the unmapped residue only is a small
job by comparison.

---

## Compute budget (estimates — confirm against actual node memory)

| Step | Partition | CPU | Est. RAM | Est. walltime | Notes |
|---|---|---|---|---|---|
| SRA prefetch + fasterq-dump | `cpu` | 8 | 8 GB | 4-8 h | network-bound; resumable per run |
| fastp QC, 12 libraries | `cpu` | 16 | 8 GB | 1 h | trivially parallel |
| hisat2-build (2.4 Gb genome) | `cpu` | 32 | **60-100 GB** | 4-10 h | **memory risk — see below** |
| HISAT2 align, 12 libraries | `cpu` | 32 | 12 GB | 3-6 h | per-sample, resumable |
| StringTie + merge | `cpu` | 16 | 16 GB | 2 h | |
| rnaSPAdes on unmapped | `cpu` | 32 | 60-120 GB | 8-20 h | **memory risk** |
| TransDecoder | `cpu` | 16 | 16 GB | 2-4 h | |

**Two steps can exceed node memory and must be checked first.** Run:

```bash
sinfo -o "%P %m %c %G %l" | column -t
```

If no node offers >= 100 GB:

- `hisat2-build`: pass `--noauto -p <cpus>` with explicit `--bmax` and
  `--dcv` to trade time for memory; or build per-chromosome and merge.
  `scripts/02_fetch_refs.py` prints the exact fallback invocation.
- `rnaSPAdes`: run with `--memory <GB>` set below the node limit, and
  normalise the unmapped reads first (`bbnorm` target=60) to cut the graph
  size. Digital normalisation of the unmapped fraction alone is safe —
  these reads are not used for quantification.

Quantification for differential expression uses the **full read set** against
the final transcript set (salmon or RSEM), never the normalised reads.

---

## Validation gates

Report all of these against the upstream assembly, as Table 1 of the paper:

| Gate | Upstream | Target | Meaning |
|---|---|---|---|
| Contig N50 | 417 bp | > 1,200 bp | fragmentation resolved |
| Loci / unigenes | 586,928 | 25,000-45,000 | collapse to real genes |
| BUSCO Arthropoda C | 85.2% | >= 90% | completeness not lost |
| BUSCO duplicated | 0.6% | < 3% | no isoform inflation |
| Predicted proteins | 28,423 | report both cutoffs | 70 aa and 100 aa strata |
| Reads mapped to final set | 76-92% | >= 90% | representation improved |

If the rebuilt assembly does not beat the upstream one on N50 and locus
count, something is wrong with the pipeline — stop and diagnose rather than
proceeding.

---

## What this does to the paper

The rebuild becomes a stated contribution rather than a hidden dependency.
It supports an additional, defensible result that strengthens Aim 1:

> *De novo assembly inflates the apparent dark proteome. Of N orphan ORFs
> called from the fragmented de novo assembly, M resolve to fragments of
> annotated loci under genome-guided reconstruction, while K represent
> genuine lineage-restricted proteins supported by peptide evidence.*

That is a cautionary, methodological finding of general interest to anyone
annotating a non-model transcriptome — and it is exactly the kind of result
that lands in *Briefings in Bioinformatics*.
