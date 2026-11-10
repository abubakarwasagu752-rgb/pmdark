#!/bin/bash
# ===========================================================================
# Fan embedding out across EVERY partition that has a free GPU.
#
#   bash slurm/embed_fanout.sh esm2 $SHARE/01_interim/all_proteins.faa
#   bash slurm/embed_fanout.sh esm2 proteins.faa --dry-run
#
# Submits one independent array job per partition. They never talk to each
# other — coordination is entirely through atomic shard claims on the shared
# filesystem (src/pmdark/claim.py). Submitting twice is harmless; extra
# workers find everything claimed and exit.
#
# Per-GPU settings are tuned to card memory. ESM-2 3B only fits on the L40s.
# ===========================================================================
set -uo pipefail
# NOTE: deliberately NOT `set -e`. The GPU-counting loop below uses
# `grep -oP` which exits 1 when a node has no GRES allocated — a normal,
# expected condition. Combined with pipefail inside a `var=$(...)`
# assignment, `set -e` terminated the script silently right after printing
# the table header. Errors are checked explicitly instead.

MODEL="${1:?usage: embed_fanout.sh MODEL FASTA [--dry-run]}"
FASTA="${2:?usage: embed_fanout.sh MODEL FASTA [--dry-run]}"
DRY=0
SKIP_CHECK=0
for a in "${@:3}"; do
  case "$a" in
    --dry-run)    DRY=1 ;;
    --skip-check) SKIP_CHECK=1 ;;
  esac
done

PROJECT="${PROJECT:-$HOME/pmdark}"
LOGS="${PMDARK_WORK:-$HOME/pmdark_work}/logs"
mkdir -p "$LOGS"
cd "$PROJECT"

# ---------------------------------------------------------------------------
# PREFLIGHT. Refuse to submit if the tokenizer cannot even load.
#
# This gate exists because 20 array jobs were once submitted across five
# partitions for ProtT5 and ProstT5 while the tokenizer check was failing:
# serprot had no `sentencepiece`, T5Tokenizer fell back to a TikToken
# extractor and raised. Every one of those jobs was doomed before it
# started, and they occupied GPUs that other work could have used.
#
# A few seconds of CPU here is worth more than a queue of certain failures.
# --skip-check overrides, for when you know better than the check.
# ---------------------------------------------------------------------------
PLM_PY="${PMDARK_PY:-$HOME/miniconda3/envs/serprot/bin/python}"
[ -x "$PLM_PY" ] || PLM_PY="$(command -v python)"
if [ "$SKIP_CHECK" = 0 ] && [ -f "$PROJECT/scripts/33b_verify_tokenization.py" ]; then
  echo "preflight: tokenizer check for $MODEL"
  if ! "$PLM_PY" scripts/33b_verify_tokenization.py --models "$MODEL" \
        > "$LOGS/preflight_$MODEL.log" 2>&1; then
    echo
    echo "ABORTED — the tokenizer for '$MODEL' does not load."
    echo "  log: $LOGS/preflight_$MODEL.log"
    echo
    tail -n 12 "$LOGS/preflight_$MODEL.log" | sed 's/^/    /'
    echo
    echo "  Nothing was submitted. Fix the environment, not the queue."
    echo "  For the T5 models this is usually a missing sentencepiece:"
    echo "    $PLM_PY -m pip install sentencepiece protobuf"
    echo "  Re-run this script afterwards, or pass --skip-check to override."
    exit 1
  fi
  echo "preflight: ok"
  echo
fi

# partition | gpu_type | vram_gb | cpus | mem_mb | max_array | batch_650m
PARTITIONS=(
  "gpul40|l40|48|8|32000|4|64"
  "gpu4090|rtx4090|24|6|20000|2|24"
  "gpua4000|a4000|16|8|24000|1|16"
  "gpu3060|rtx3060|12|8|24000|2|10"
  "gpu3080|rtx3080|10|8|24000|1|8"
  "gpu4070ti|rtx4070ti|12|8|20000|1|10"
)

# Models that need more than ~16 GB of VRAM; small cards are skipped for these.
case "$MODEL" in
  esm2_3b) MIN_VRAM=40 ;;
  prott5|prostt5|saprot) MIN_VRAM=12 ;;
  *) MIN_VRAM=8 ;;
esac

echo "=============================================================="
echo "fan-out embedding"
echo "  model      : $MODEL  (needs >= ${MIN_VRAM} GB VRAM)"
echo "  fasta      : $FASTA"
echo "=============================================================="

[ -f "$FASTA" ] || { echo "FASTA not found: $FASTA"; exit 1; }

printf "\n%-12s %-10s %-6s %-9s %-7s %s\n" \
       PARTITION GPU VRAM "free/tot" BATCH DECISION
echo "--------------------------------------------------------------"

SUBMITTED=0
TOTAL_WORKERS=0

for spec in "${PARTITIONS[@]}"; do
  IFS='|' read -r part gputype vram cpus mem maxarr batch <<< "$spec"

  # Does the partition exist and is it up?
  if ! sinfo -h -p "$part" -o "%P" 2>/dev/null | grep -q .; then
    printf "%-12s %-10s %-6s %-9s %-7s %s\n" "$part" "$gputype" "${vram}G" "-" "-" "absent"
    continue
  fi

  # Count free GPUs on this partition: total GRES minus allocated GRES.
  n_tot=0; n_use=0
  while IFS='|' read -r g b; do
    t=$(echo "$g" | grep -oP 'gpu:[^:]+:\K[0-9]+' 2>/dev/null | head -1 || true)
    n_tot=$(( n_tot + ${t:-0} ))
    u=$(echo "$b" | grep -oP 'gpu:[^:]+:\K[0-9]+' 2>/dev/null | head -1 || true)
    n_use=$(( n_use + ${u:-0} ))
  done < <(sinfo -h -N -p "$part" -o "%G|%b" 2>/dev/null)
  n_free=$(( n_tot - n_use )); (( n_free < 0 )) && n_free=0

  if (( vram < MIN_VRAM )); then
    printf "%-12s %-10s %-6s %-9s %-7s %s\n" \
      "$part" "$gputype" "${vram}G" "${n_free}/${n_tot}" "-" "skip: VRAM too small for $MODEL"
    continue
  fi

  # Queue even when currently busy — Slurm will start us when a GPU frees.
  # Cap the array at the partition's GPU count so we never over-subscribe.
  n_req=$maxarr
  (( n_req > n_tot )) && n_req=$n_tot
  (( n_req < 1 )) && n_req=1

  note="submit array 0-$(( n_req - 1 ))"
  (( n_free == 0 )) && note="$note (queued; no free GPU now)"

  printf "%-12s %-10s %-6s %-9s %-7s %s\n" \
    "$part" "$gputype" "${vram}G" "${n_free}/${n_tot}" "$batch" "$note"

  if (( DRY == 1 )); then
    continue
  fi

  JID=$(sbatch --parsable \
    --job-name="pmd-${MODEL}-${part}" \
    --partition="$part" \
    --gres=gpu:1 \
    --cpus-per-task="$cpus" \
    --mem="${mem}M" \
    --time=11:55:00 \
    --array="0-$(( n_req - 1 ))" \
    --output="${LOGS}/%x_%A_%a.out" \
    --error="${LOGS}/%x_%A_%a.err" \
    --export=ALL,MODEL="$MODEL",FASTA="$FASTA",BATCH="$batch",PROJECT="$PROJECT" \
    "$PROJECT/slurm/_embed_worker.sbatch" 2>&1) || {
      echo "    submit FAILED on $part"; continue; }

  echo "    -> job $JID"
  SUBMITTED=$(( SUBMITTED + 1 ))
  TOTAL_WORKERS=$(( TOTAL_WORKERS + n_req ))
done

echo "--------------------------------------------------------------"
if (( DRY == 1 )); then
  echo "dry run — nothing submitted"
else
  echo "submitted $SUBMITTED array jobs, up to $TOTAL_WORKERS concurrent GPU workers"
  echo
  echo "monitor:"
  echo "  squeue -u \$USER"
  echo "  python scripts/embed_status.py --model $MODEL"
fi
