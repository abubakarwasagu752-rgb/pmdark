#!/bin/bash
# ===========================================================================
# Verify the shared folder is genuinely writable and fast enough, WITHOUT
# touching anything that is already there.
#
#   bash scripts/00b_verify_share.sh
#
# /shared/a2-public is mode 755 root:root. The `[ -w ]` shell test returns
# true on NFS mounts even when the user cannot actually write, because it
# inspects mount options rather than server-side permissions. The only
# reliable check is a real mkdir. That is what this does.
#
# This script creates exactly one new directory for you and removes its own
# probe files. It never deletes, moves or modifies anything pre-existing.
# ===========================================================================
set -uo pipefail

SHARED_ROOT="${SHARED_ROOT:-/shared/a2-public}"
# Numeric part of the username — the full USER contains '@' which is awkward
# inside paths and breaks some tools.
ME="${PMDARK_USER:-$(echo "$USER" | cut -d@ -f1)}"
WORKSPACE="$SHARED_ROOT/$ME/pmdark"

echo "=========================================================="
echo "PmDark share verification"
echo "=========================================================="
echo "shared root : $SHARED_ROOT"
echo "user        : $USER"
echo "workspace   : $WORKSPACE"
echo

# --- 1. what already exists (read-only inspection) ------------------------
echo "--- existing contents of $SHARED_ROOT (NOT touched) ---"
ls -la "$SHARED_ROOT" 2>/dev/null || { echo "cannot read $SHARED_ROOT"; exit 1; }
echo
echo "--- per-directory sizes (may take a moment) ---"
du -sh "$SHARED_ROOT"/* 2>/dev/null | sort -h | tail -20
echo

# --- 2. real write test ---------------------------------------------------
echo "--- write test ---"
if [ -e "$WORKSPACE" ]; then
  echo "NOTE: $WORKSPACE already exists; reusing it."
else
  if mkdir -p "$WORKSPACE" 2>/dev/null; then
    echo "created $WORKSPACE"
  else
    echo
    echo "CANNOT CREATE $WORKSPACE"
    echo "  $SHARED_ROOT is root-owned mode 755, so a subdirectory may need"
    echo "  to be created for you by an admin, or you may have an existing"
    echo "  allocation under a name you did not try."
    echo
    echo "  Look for a directory you already own:"
    echo "    find $SHARED_ROOT -maxdepth 2 -user \"\$(id -u)\" 2>/dev/null"
    echo
    echo "  If none, ask the cluster admin for a directory under"
    echo "  $SHARED_ROOT owned by $USER."
    exit 1
  fi
fi

PROBE="$WORKSPACE/.pmdark_probe"
mkdir -p "$PROBE" || { echo "cannot create probe dir"; exit 1; }

if echo ok > "$PROBE/w.txt" 2>/dev/null && [ "$(cat "$PROBE/w.txt")" = "ok" ]; then
  echo "write/read : OK"
else
  echo "write/read : FAILED"
  rm -rf "$PROBE"
  exit 1
fi

# --- 3. throughput --------------------------------------------------------
echo
echo "--- sequential write throughput (1 GB) ---"
dd if=/dev/zero of="$PROBE/speed.bin" bs=1M count=1024 conv=fsync 2>&1 | tail -1
echo "--- sequential read throughput ---"
# Drop-cache is not permitted as a normal user, so this read may be partly
# served from page cache. Treat it as an upper bound.
dd if="$PROBE/speed.bin" of=/dev/null bs=1M 2>&1 | tail -1

# --- 4. small-file / metadata latency ------------------------------------
# PLM embedding writes thousands of small shard files. NFS metadata latency
# matters more here than raw bandwidth.
echo
echo "--- metadata latency (500 small files) ---"
START=$(date +%s.%N)
for i in $(seq 1 500); do echo x > "$PROBE/f_$i"; done
END=$(date +%s.%N)
echo "create 500 files: $(echo "$END - $START" | bc) s"
START=$(date +%s.%N)
ls "$PROBE" > /dev/null
END=$(date +%s.%N)
echo "list dir        : $(echo "$END - $START" | bc) s"

# --- 5. cleanup of our own probe only ------------------------------------
rm -rf "$PROBE"
echo "probe removed (nothing else touched)"

# --- 6. capacity ----------------------------------------------------------
echo
echo "--- capacity ---"
df -h "$WORKSPACE"
echo
echo "--- quota on this filesystem (if any) ---"
quota -s -f "$WORKSPACE" 2>/dev/null \
  || lfs quota -h -u "$USER" "$WORKSPACE" 2>/dev/null \
  || echo "no quota reported for $WORKSPACE (likely unquotaed)"

# --- 7. the line to paste ------------------------------------------------
echo
echo "=========================================================="
echo "Paste this into config/config.yaml under paths:"
echo
echo "  share: $WORKSPACE"
echo "=========================================================="
