#!/bin/bash
# ===========================================================================
# PmDark — PLM embedding, checkpointed job array
#
# The gpul140 partition (4x L40 48GB) caps at 12 h and allows 30 CPU/job.
# This script runs up to the cap, exits cleanly at 11h48m, and resubmits
# itself if work remains. Running it twice is always safe.
#
# Submit:
#   sbatch --export=ALL,MODEL=esm2,FASTA=$SHARE/01_interim/all_proteins.faa \
#          slurm/10_embed.sbatch
# ===========================================================================
#SBATCH --job-name=pmd-embed
#SBATCH --partition=gpul140
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=11:55:00
#SBATCH --array=0-3
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

PROJECT="${PROJECT:-$HOME/pmdark}"
MODEL="${MODEL:?set MODEL=esm2|prott5|prostt5|saprot|esm2_3b}"
FASTA="${FASTA:?set FASTA=/path/to/proteins.faa}"
WALLTIME_MIN="${WALLTIME_MIN:-715}"
MAX_RESUBMIT="${MAX_RESUBMIT:-6}"
RESUBMIT_COUNT="${RESUBMIT_COUNT:-0}"

cd "$PROJECT"

# --- environment ----------------------------------------------------------
# Keep every cache off /home — it has ~69 GB free and HF caches are large.
SHARE_DIR="$(python -c 'import sys;sys.path.insert(0,"src");from pmdark import config;print(config.load().p("share"))')"
export HF_HOME="$SHARE_DIR/02_models/hf"
export TRANSFORMERS_CACHE="$HF_HOME"
export TORCH_HOME="$SHARE_DIR/02_models/torch"
export TOKENIZERS_PARALLELISM=false
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

OUT="$SHARE_DIR/04_embeddings/$MODEL"
mkdir -p "$OUT" "$HF_HOME" "$TORCH_HOME"

if [[ -n "${CONDA_ENV:-}" ]]; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$CONDA_ENV"
elif [[ -f "$PROJECT/.venv/bin/activate" ]]; then
  source "$PROJECT/.venv/bin/activate"
fi

echo "=== $(date -Is) | task ${SLURM_ARRAY_TASK_ID:-0} | resubmit $RESUBMIT_COUNT ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# --- run ------------------------------------------------------------------
set +e
python scripts/10_embed.py \
  --fasta "$FASTA" \
  --model "$MODEL" \
  --out   "$OUT" \
  --walltime-min "$WALLTIME_MIN"
RC=$?
set -e

# --- resubmit if incomplete ----------------------------------------------
if [[ $RC -eq 99 ]]; then
  if (( RESUBMIT_COUNT >= MAX_RESUBMIT )); then
    echo "Hit MAX_RESUBMIT=$MAX_RESUBMIT; stopping. Inspect $OUT and resubmit manually."
    exit 1
  fi
  # Only array task 0 resubmits, so we don't fork N chains per wave
  if [[ "${SLURM_ARRAY_TASK_ID:-0}" == "0" ]]; then
    echo "Work remains — resubmitting (wave $((RESUBMIT_COUNT + 1)))"
    sbatch --export=ALL,MODEL="$MODEL",FASTA="$FASTA",\
RESUBMIT_COUNT=$((RESUBMIT_COUNT + 1)),PROJECT="$PROJECT" \
      "$PROJECT/slurm/10_embed.sbatch"
  fi
  exit 0
fi

exit $RC
