#!/bin/bash
# ===========================================================================
# Move the bulk of the project off the quota'd work tier onto the archive,
# leaving symlinks behind.
#
#   bash scripts/35_retier.sh              # dry run, shows the plan
#   bash scripts/35_retier.sh --apply      # do it
#   bash scripts/35_retier.sh --check      # just test archive write access
#
# WHY THIS EXISTS, AND WHAT I GOT WRONG
#
#   The whole three-tier design rests on a claim in config/config.yaml:
#
#     "The decisive fact: /shared/a2-public is NOT MOUNTED ON COMPUTE NODES.
#        slurmstepd-gpu05: couldn't chdir to `/shared/a2-public/...'
#        mkdir: cannot create directory '/shared': Permission denied
#      It is login-node only. No batch job may read or write it."
#
#   That conclusion was wrong. Tested on gpul40:
#
#     --- read:   /shared/a2-public/abk/pmdark
#     --- write:  WRITE OK
#     --- df:     8.8T  680G  8.1T   8%  /shared/a2-public
#
#   The second error was `mkdir /shared` failing on the filesystem ROOT,
#   which is a permission error about creating a top-level directory, not
#   evidence of a missing mount. I read it as the latter and generalised
#   from a single node, and everything followed from there: a 92 GB
#   allowance instead of 8.1 TB, aggressive mid-pipeline cleanup, deleting
#   BAMs that then forced a full re-alignment of 12 libraries, and finally
#   a quota deadlock that stopped the embeddings outright.
#
# WHY SYMLINKS RATHER THAN EDITING THE CONFIG
#
#   Absolute paths are already baked into results on disk —
#   orf_calling.json records pep_file and genome_gff3, protein_qc.json
#   points at both, the shard stores carry their own paths. Repointing the
#   config would invalidate all of it and require rewriting those files.
#   Moving the data and leaving a symlink at the old location keeps every
#   recorded path valid, costs nothing against the quota, and is trivially
#   reversible.
#
# WHAT IS NOT MOVED
#
#   envs/   conda environments. Interpreter paths are pinned by absolute
#           path in every sbatch; moving them breaks all of it for no gain.
#   logs/   small, and written constantly by running jobs.
#   scratch/ should be deleted, not moved. See the note below.
# ===========================================================================
set -uo pipefail

WORK="${PMDARK_WORK:-$HOME/pmdark_work}"
SHARE="${PMDARK_ARCHIVE:-/shared/a2-public/abk/pmdark}"
DEST="$SHARE/work_tier"
MODE="dry"
for a in "$@"; do
  case "$a" in
    --apply) MODE="apply" ;;
    --check) MODE="check" ;;
  esac
done

MOVE=(models genome align proteins embeddings structures)

hr() { printf '%s\n' "--------------------------------------------------------------"; }

echo "work    : $WORK"
echo "archive : $DEST"
hr

# --- can the archive actually take writes? ---------------------------------
mkdir -p "$DEST" 2>/dev/null
probe="$DEST/.write_probe.$$"
if ! touch "$probe" 2>/dev/null; then
  echo "FAIL: cannot write to $DEST from this node."
  echo "  Nothing was moved."
  exit 1
fi
rm -f "$probe"
echo "archive is writable from $(hostname)"

echo
echo "IMPORTANT: this script runs on the login node, so it has only proved"
echo "the LOGIN node can write there. Compute nodes are what matter. Verify"
echo "every GPU partition before relying on this:"
echo
for p in gpul40 gpu4090 gpua4000 gpu3060 gpu3080 gpu4070ti; do
  echo "  srun -p $p --gres=gpu:1 -t 2 --pty bash -c 'touch $DEST/.probe_$p && echo $p OK'"
done
[ "$MODE" = "check" ] && exit 0

echo
echo "QUOTA NOW"
quota -s 2>/dev/null | tail -2 || echo "  (no quota command)"

echo
echo "PLAN"
hr
total=0
for d in "${MOVE[@]}"; do
  src="$WORK/$d"
  [ -e "$src" ] || continue
  if [ -L "$src" ]; then
    echo "  skip     $d  (already a symlink)"
    continue
  fi
  sz=$(du -sm "$src" 2>/dev/null | cut -f1)
  total=$((total + sz))
  echo "  move     $d  (${sz} MB)  ->  $DEST/$d"
  echo "           then symlink $src -> $DEST/$d"
done
hr
echo "  frees roughly ${total} MB of quota"

if [ -d "$WORK/scratch" ]; then
  ssz=$(du -sm "$WORK/scratch" 2>/dev/null | cut -f1)
  echo
  echo "  NOT handled here: $WORK/scratch (${ssz} MB)."
  echo "  That is dead output from finished jobs — config uses /tmp for"
  echo "  scratch and this directory is only the fallback. Delete it:"
  echo "    rm -rf $WORK/scratch/*"
fi

if [ "$MODE" != "apply" ]; then
  echo
  echo "(dry run — rerun with --apply)"
  exit 0
fi

echo
echo "APPLYING"
hr
for d in "${MOVE[@]}"; do
  src="$WORK/$d"
  [ -e "$src" ] || continue
  [ -L "$src" ] && continue
  tgt="$DEST/$d"
  if [ -e "$tgt" ]; then
    echo "  $tgt already exists — merging with rsync"
    if ! rsync -a --remove-source-files "$src/" "$tgt/"; then
      echo "  FAILED on $d; leaving it in place"
      continue
    fi
    find "$src" -type d -empty -delete 2>/dev/null
  else
    echo "  moving $d ..."
    if ! mv "$src" "$tgt"; then
      echo "  FAILED on $d; leaving it in place"
      continue
    fi
  fi
  if [ -e "$src" ]; then
    echo "  $src still present; not symlinking (inspect manually)"
    continue
  fi
  ln -s "$tgt" "$src" && echo "  linked $src -> $tgt"
done

hr
echo "QUOTA AFTER"
quota -s 2>/dev/null | tail -2 || true
echo
echo "Verify a compute node can follow the symlinks before resubmitting:"
echo "  srun -p gpul40 --gres=gpu:1 -t 2 --pty bash -c '"
echo "    ls \$HOME/pmdark_work/embeddings/ && echo READ-THROUGH OK'"
