#!/bin/bash
# ===========================================================================
# RE-RUN of alignment + assembly, with the two corrections from 29.
#
#   sbatch --array=0-11%4 slurm/21b_realign.sbatch   # 4 at a time
#
# WHY RE-RUN
#   1. StringTie now gets -G (guide, NOT -e). The original run omitted it on
#      my mistaken reasoning that a guide suppresses novel loci; -e does
#      that, -G does not. Guided assembly completes transcripts at
#      annotated loci that sparse coverage left as fragments.
#   2. A second StringTie pass with -e -G quantifies the reference
#      annotation, giving the expression denominator that every recovery
#      figure so far has lacked.
#   BAMs are deliberately NOT kept. Both StringTie passes run from the same
#   scratch BAM inside the same job, so retention buys nothing, and 12 BAMs
#   at 3-5 GB each would be 36-60 GB against a ~100 GB quota shared with
#   five other projects. Add --keep-bam only if you want them for IGV.
#
#   --force is passed because each library carries a DONE marker from the
#   un-guided run. Existing outputs are overwritten, nothing else is
#   touched, and /shared is never written by a compute node.
#
# Each task is self-contained: it downloads its own reads from ENA into
# node-local /tmp, aligns there, and writes only small outputs back to the
# quota'd HOME tier. Nothing touches /shared, which compute nodes cannot see.
#
# Runs on the cpu partition (24 h, 32 cores). The GPU partitions would work
# too but there is no GPU step here and they cap at 12 h.
# ===========================================================================
#SBATCH --job-name=pmd-realign
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -uo pipefail

PROJECT="${PROJECT:-$HOME/pmdark}"
WORK="${PMDARK_WORK:-$HOME/pmdark_work}"

cd "$PROJECT"

# Use the pmbio interpreter by ABSOLUTE PATH, not via `conda activate`.
#
# `conda activate` inside a non-interactive batch script is unreliable: the
# job inherits the submitting shell's environment (--export=ALL is the
# default), so if you submit from `serprot` the job may keep serprot's
# python even after sourcing activate_bio.sh. That is what killed job
# 24712 with ModuleNotFoundError: No module named 'yaml' — pmbio HAS
# pyyaml, the job just was not using pmbio's python.
#
# Activation is still sourced, because hisat2/stringtie/samtools must be on
# PATH; but the interpreter is pinned explicitly.
[ -f "$WORK/activate_bio.sh" ] && source "$WORK/activate_bio.sh"

PY="$WORK/envs/pmbio/bin/python"
[ -x "$PY" ] || PY="$(command -v python)"

export TMPDIR="/tmp/pmdark_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$TMPDIR"

echo "=========================================================="
echo "node    : $(hostname)"
echo "task    : ${SLURM_ARRAY_TASK_ID:-0}"
echo "cpus    : ${SLURM_CPUS_PER_TASK:-16}"
echo "tmpdir  : $TMPDIR  ($(df -h /tmp | awk 'NR==2{print $4}') free)"
echo "python  : $PY"
echo "started : $(date -Is)"
echo "=========================================================="

"$PY" scripts/21_align_assemble.py --threads "${SLURM_CPUS_PER_TASK:-16}" \
      --force --guided --ref-quant
RC=$?

echo "cleaning $TMPDIR"
rm -rf "$TMPDIR"

echo "finished: $(date -Is)  rc=$RC"
exit $RC
