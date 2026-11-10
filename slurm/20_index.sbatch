#!/bin/bash
# ===========================================================================
# Build the HISAT2 index for the P. monodon genome (GCF_015228065.2).
#
#   sbatch slurm/20_index.sbatch
#
# --mem=100G forces placement on cpu01, the only node with 128 GB. The other
# cpu nodes have 63 GB. The cpu partition allows 24 h, double the 12 h GPU
# cap, which matters because hisat2-build cannot resume if interrupted.
#
# The genome lives on the ARCHIVE tier, which compute nodes cannot see, so
# this job stages it through HOME first. That copy is done by the submitting
# user on the login node — see the guard below.
# ===========================================================================
#SBATCH --job-name=pmd-index
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=23:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
# NOTE: relative paths here resolve against the submit dir. Submit from
# $HOME/pmdark, never /shared — compute nodes cannot write there and the
# job dies before producing a log (see job 24697).

set -euo pipefail

PROJECT="${PROJECT:-$HOME/pmdark}"
WORK="${PMDARK_WORK:-$HOME/pmdark_work}"
GENOME_DIR="$WORK/genome"
PREFIX="$GENOME_DIR/pmon"

cd "$PROJECT"
[ -f "$WORK/activate_bio.sh" ] && source "$WORK/activate_bio.sh"

echo "=========================================================="
echo "node     : $(hostname)"
echo "memory   : $(free -g | awk 'NR==2{print $2" GB"}')"
echo "cpus     : ${SLURM_CPUS_PER_TASK:-?}"
echo "genome   : $GENOME_DIR"
echo "started  : $(date -Is)"
echo "=========================================================="

# CAREFUL: the NCBI download contains BOTH
#     *_genomic.fna.gz            (the 2.4 Gb genome)
#     *_cds_from_genomic.fna.gz   (13 MB of coding sequence)
# and the second matches the same glob. Alphabetically "cds_from" sorts
# FIRST, so a naive `head -1` silently indexes the CDS file and finishes in
# two minutes with a 358 MB index instead of several GB. That is exactly
# what job 24698 did. Exclude the derived files, then take the largest.
pick() { ls -S "$GENOME_DIR"/$1 2>/dev/null \
           | grep -vE '_(cds|rna)_from_genomic' | head -1 || true; }

FA=$(pick '*_genomic.fna')
if [ -z "$FA" ]; then
  GZ=$(pick '*_genomic.fna.gz')
  if [ -z "$GZ" ]; then
    echo
    echo "ERROR: no genome FASTA in $GENOME_DIR"
    echo
    echo "Compute nodes cannot see /shared, so stage it from the LOGIN node:"
    echo "  mkdir -p $GENOME_DIR"
    echo "  cp /shared/a2-public/abk/pmdark/00_raw/genome/*_genomic.fna.gz \\"
    echo "     /shared/a2-public/abk/pmdark/00_raw/genome/*_genomic.gff.gz \\"
    echo "     $GENOME_DIR/"
    echo "  gunzip $GENOME_DIR/*.gz"
    exit 1
  fi
  echo "decompressing $(basename "$GZ")"
  if command -v pigz >/dev/null; then
    pigz -d -p "${SLURM_CPUS_PER_TASK:-8}" -k "$GZ"
  else
    gunzip -k "$GZ"
  fi
  FA="${GZ%.gz}"
fi

FA_BYTES=$(stat -c %s "$FA")
echo "fasta    : $FA  ($(du -h "$FA" | cut -f1))"

# Sanity gate. The P. monodon genome is ~2.4 Gb; anything under 1 GB means
# we picked a derived file again and the whole pipeline downstream would be
# quietly wrong. Fail loudly instead.
if [ "$FA_BYTES" -lt 1000000000 ]; then
  echo
  echo "ERROR: $FA is only $((FA_BYTES/1000000)) MB — expected ~2400 MB."
  echo "This is a CDS or RNA file, not the genome. Files present:"
  ls -laS "$GENOME_DIR"/*.fna 2>/dev/null
  exit 1
fi

if [ -f "${PREFIX}.1.ht2" ] || [ -f "${PREFIX}.1.ht2l" ]; then
  echo "index already present at $PREFIX — nothing to do"
  ls -la "$GENOME_DIR"/pmon.*.ht2* | head
  exit 0
fi

echo
echo "building HISAT2 index (expect 2-8 h for a 2.4 Gb genome)"
echo "NOTE: hisat2-build cannot resume. If this is killed, it restarts."
/usr/bin/time -v hisat2-build -p "${SLURM_CPUS_PER_TASK:-32}" "$FA" "$PREFIX" \
  2> >(tee "$GENOME_DIR/hisat2-build.time" >&2)

echo
echo "--- index files ---"
ls -la "$GENOME_DIR"/pmon.*.ht2*
du -sh "$GENOME_DIR"
echo
echo "--- quota after build ---"
df -h "$HOME" | tail -1

echo "finished : $(date -Is)"
echo "Next: sbatch --array=0-11 slurm/21_align.sbatch"
