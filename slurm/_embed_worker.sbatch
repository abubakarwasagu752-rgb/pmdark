#!/bin/bash
# ===========================================================================
# Internal worker, launched by slurm/embed_fanout.sh — not run directly.
#
# All resource flags (partition, gres, mem, cpus, array, time) come from the
# fan-out submitter's sbatch command line, so this file carries none of its
# own. The same script runs unchanged on an L40 and on an RTX 3080.
#
# TWO THINGS THIS GETS RIGHT, both learned the hard way:
#
#  1. Output goes to the WORK tier, never the archive. /shared/a2-public is
#     not mounted on compute nodes. The previous version wrote to
#     $SHARE/04_embeddings and every worker died instantly, leaving all ten
#     shards unclaimed.
#
#  2. The interpreter is pinned by absolute path. `conda activate` inside a
#     non-interactive batch script is unreliable — the job inherits the
#     submitting shell's env, so it may keep the wrong python and fail on
#     `import yaml`.
# ===========================================================================
set -uo pipefail

PROJECT="${PROJECT:-$HOME/pmdark}"
WORK="${PMDARK_WORK:-$HOME/pmdark_work}"
MODEL="${MODEL:?}"
FASTA="${FASTA:?}"
BATCH="${BATCH:-16}"

cd "$PROJECT"

# serprot carries torch + transformers; pmbio does not.
PY="${PMDARK_PY:-$HOME/miniconda3/envs/serprot/bin/python}"
[ -x "$PY" ] || PY="$(command -v python)"

# Caches off the small HOME quota root, but still on compute-visible storage.
export HF_HOME="$WORK/models/hf"
export HF_HUB_CACHE="$HF_HOME/hub"
export TORCH_HOME="$WORK/models/torch"
export TOKENIZERS_PARALLELISM=false
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export PYTHONPATH="$PROJECT/src:${PYTHONPATH:-}"
mkdir -p "$HF_HOME" "$TORCH_HOME" "$WORK/logs"

# Resolve the output directory through the config so it always tracks the
# WORK tier, and give each FASTA its own subdirectory (query/ vs reference/)
# because 31_calibrate.py loads them separately.
STEM="$(basename "$FASTA")"; STEM="${STEM%.*}"
OUT="$WORK/embeddings/$MODEL/$STEM"
mkdir -p "$OUT"

echo "=========================================================="
echo "host      : $(hostname)"
echo "partition : ${SLURM_JOB_PARTITION:-?}"
echo "job/task  : ${SLURM_JOB_ID:-?} / ${SLURM_ARRAY_TASK_ID:-0}"
echo "python    : $PY"
echo "model     : $MODEL   batch $BATCH"
echo "fasta     : $FASTA"
echo "out       : $OUT"
echo "started   : $(date -Is)"
echo "=========================================================="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# ---------------------------------------------------------------------------
# GPU HEALTH GATE. Exit BEFORE claiming a shard if this node has no usable
# device.
#
# Observed on this cluster:
#   Unable to determine the device handle for GPU0: 0000:00:10.0: Unknown Error
#   CUDA initialization: CUDA unknown error ... Setting the available
#   devices to be zero
#
# The card is present to Slurm but broken for CUDA. Without this gate the
# worker starts, claims a shard, fails, and the shard then sits claimed
# until `stale_after` (an hour by default) expires — so one sick node can
# stall the whole fan-out. Exiting 0 keeps Slurm from flagging a failure
# and leaves the shard unclaimed for a healthy worker to pick up.
# ---------------------------------------------------------------------------
if ! "$PY" -c "import sys, torch; sys.exit(0 if torch.cuda.is_available() else 1)" 2>/dev/null; then
  echo "NO USABLE GPU on $(hostname) — torch.cuda.is_available() is False."
  echo "Exiting without claiming a shard so a healthy node can take it."
  nvidia-smi || true
  exit 0
fi
echo "gpu check : ok"

# ---------------------------------------------------------------------------
# UNBUFFERED, AND TEED TO OUR OWN FILE.
#
# Slurm buffers a job's stdout/stderr and loses the buffer when the task
# dies, so every failed worker left behind exactly one line —
# "Traceback (most recent call last):" — with the actual exception gone.
# Four debugging rounds were spent on logs that had already been thrown
# away.
#
#   python -u  : no Python-side buffering, so the traceback reaches the fd
#   tee        : writes through to a file we control as it arrives, which
#                survives the task being killed
#   PIPESTATUS : keep the python exit code, not tee's
# ---------------------------------------------------------------------------
# $LOGS is defined in embed_fanout.sh, NOT here, and this script runs
# `set -u`. Referencing it killed every worker two seconds in with
# "LOGS: unbound variable" — the logging fix meant to help us debug was
# itself the failure. Derive it from WORK, which this script does define.
LOGS="$WORK/logs"
mkdir -p "$LOGS"
WLOG="$LOGS/worker_${MODEL}_${STEM}_${SLURM_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID:-0}.log"
echo "worker log: $WLOG"

"$PY" -u scripts/10_embed.py \
  --fasta "$FASTA" \
  --model "$MODEL" \
  --out   "$OUT" \
  --batch "$BATCH" \
  --walltime-min 715 2>&1 | tee "$WLOG"
RC=${PIPESTATUS[0]}

echo "finished  : $(date -Is)  rc=$RC"

# rc 99 means shards remain. With fan-out other workers may still be live,
# so do not self-resubmit; re-run embed_fanout.sh instead (it is idempotent
# and only picks up genuinely unclaimed shards).
exit $(( RC == 99 ? 0 : RC ))
