#!/bin/bash
# ===========================================================================
# Install the project into compute-visible storage.
#
#   bash scripts/90_bootstrap_home.sh
#
# WHY THIS EXISTS
#   /shared/a2-public is login-node only. A Slurm job cannot read the
#   scripts there, and cannot write its --output file there either. Job
#   24697 died in 1 second for exactly this reason, with no log to show for
#   it, because slurmstepd could not create the output file.
#
#   So the CODE moves to $HOME/pmdark (compute-visible), and /shared keeps
#   its job as the data archive. The code is ~60 KB; this costs nothing
#   against the quota.
#
# After this runs, work from $HOME/pmdark. Not from /shared.
# ===========================================================================
set -euo pipefail

SRC="${SRC:-/shared/a2-public/abk/pmdark}"
DST="${DST:-$HOME/pmdark}"
WORK="${PMDARK_WORK:-$HOME/pmdark_work}"

echo "source (archive) : $SRC"
echo "dest   (code)    : $DST"
echo "work             : $WORK"
echo

mkdir -p "$DST" "$WORK"/{logs,models,genome,align,proteins,embeddings,structures,scratch}

# Code only — never the 00_raw / 06_results data trees.
for d in src scripts slurm config docs env; do
  if [ -d "$SRC/$d" ]; then
    rsync -a --delete "$SRC/$d/" "$DST/$d/" 2>/dev/null \
      || { mkdir -p "$DST/$d"; cp -r "$SRC/$d/." "$DST/$d/"; }
    echo "  copied $d"
  fi
done
[ -f "$SRC/README.md" ] && cp "$SRC/README.md" "$DST/" || true

echo
echo "size of code tree: $(du -sh "$DST" | cut -f1)"

# ---------------------------------------------------------------------------
# Python dependencies for the pmbio env.
# pmbio was built for bioinformatics binaries and has no pyyaml, so every
# pmdark script fails there with ModuleNotFoundError. serprot has torch and
# transformers but no bio tools. Both envs need the small pure-python deps.
# ---------------------------------------------------------------------------
echo
echo "--- python deps ---"
if [ -n "${CONDA_PREFIX:-}" ]; then
  echo "active env: $CONDA_PREFIX"
  python -m pip install -q --no-input pyyaml 2>&1 | tail -2 || true
  python -c "import yaml; print('  pyyaml', yaml.__version__)" 2>/dev/null \
    || echo "  pyyaml still missing — install manually"
else
  echo "  no conda env active; activate one first"
fi

# ---------------------------------------------------------------------------
# Point the config at the right places and verify.
# ---------------------------------------------------------------------------
cd "$DST"
echo
echo "--- resolved paths ---"
python -c "
import sys; sys.path.insert(0,'src')
from pmdark import config as c
print(c.describe(c.load()))
" || echo "config load failed — check pyyaml"

cat > "$DST/env.sh" <<EOF
# source this before working:  source ~/pmdark/env.sh
export PMDARK_HOME="$DST"
export PMDARK_WORK="$WORK"
export PMDARK_ARCHIVE="$SRC"
export PYTHONPATH="$DST/src:\${PYTHONPATH:-}"
export HF_HOME="$WORK/models/hf"
export HF_HUB_CACHE="\$HF_HOME/hub"
export TORCH_HOME="$WORK/models/torch"
export TOKENIZERS_PARALLELISM=false
cd "$DST"
EOF

echo
echo "wrote $DST/env.sh"
echo
echo "=========================================================="
echo "From now on, work from $DST — not from /shared."
echo
echo "  source ~/pmdark/env.sh"
echo "  sbatch slurm/20_index.sbatch"
echo "=========================================================="
