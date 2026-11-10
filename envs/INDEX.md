# Conda environments — captured 2026-11-09T08:46:28+00:00

Active at capture time: `pmbio`

This supersedes `ENVIRONMENT.txt` in the same directory, which recorded only the environment that happened to be active. That file is left in place rather than overwritten, because it is a true record of what the freeze captured.

## Environments

| env | conda pkgs | pip pkgs | prefix |
|---|---:|---:|---|
| `busco` | 253 | 27 | `/storage/slurm/home/68090800301@cpe.kmutt.ac.th/pmdark_work/envs/busco` |
| `esmprobe` | 96 | 66 | `/storage/slurm/home/68090800301@cpe.kmutt.ac.th/miniconda3/envs/esmprobe` |
| `gmx_gpu` | 17 | 0 | `/storage/slurm/home/68090800301@cpe.kmutt.ac.th/miniconda3/envs/gmx_gpu` |
| `miniconda3` | 129 | 71 | `/storage/slurm/home/68090800301@cpe.kmutt.ac.th/miniconda3` |
| `pmbio` | 231 | 46 | `/storage/slurm/home/68090800301@cpe.kmutt.ac.th/pmdark_work/envs/pmbio` |
| `serprot` | 128 | 84 | `/storage/slurm/home/68090800301@cpe.kmutt.ac.th/miniconda3/envs/serprot` |

## Where each tool lives

The question someone rebuilding this will actually ask.

| tool | environment(s) |
|---|---|
| `diamond` | `busco` |
| `jackhmmer` | `busco`, `pmbio` |
| `hmmsearch` | `busco`, `pmbio` |
| `hmmbuild` | `busco`, `pmbio` |
| `hisat2` | `pmbio` |
| `stringtie` | `pmbio` |
| `samtools` | `busco`, `pmbio` |
| `foldseek` | — |
| `mmseqs` | `serprot` |
| `prodigal` | `busco`, `pmbio` |
| `TransDecoder.LongOrfs` | — |
| `busco` | `busco` |
| `fasterq-dump` | — |
| `zstd` | `busco`, `miniconda3`, `pmbio` |

## Where each package lives

| package | environment(s) |
|---|---|
| `pytorch` | — |
| `torch` | `esmprobe`, `serprot` |
| `transformers` | `esmprobe`, `serprot` |
| `sentencepiece` | `serprot` |
| `pydeseq2` | `pmbio` |
| `biopython` | `busco`, `serprot` |
| `numpy` | `busco`, `esmprobe`, `pmbio`, `serprot` |
| `scipy` | `esmprobe`, `pmbio`, `serprot` |
| `pandas` | `busco`, `esmprobe`, `pmbio`, `serprot` |
| `matplotlib` | `esmprobe`, `pmbio`, `serprot` |
| `scikit-learn` | `esmprobe`, `pmbio`, `serprot` |
| `statsmodels` | `serprot` |
| `pyhmmer` | — |
| `python-docx` | `pmbio`, `serprot` |
| `openpyxl` | `esmprobe`, `serprot` |
| `pyyaml` | `busco`, `esmprobe`, `pmbio`, `serprot` |
| `goatools` | — |
| `anndata` | `pmbio` |

## Capture problems

- `gmx_gpu`: no python in this prefix

## Rebuilding

```bash
conda create --name <env> --file <env>.explicit.txt
# then, if that env has a pip file:
pip install -r <env>.pip.txt
```

