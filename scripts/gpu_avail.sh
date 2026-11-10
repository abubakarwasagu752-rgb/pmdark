#!/bin/bash
# ===========================================================================
# Live GPU availability across every partition.
#
#   bash scripts/gpu_avail.sh            # snapshot
#   watch -n 20 bash scripts/gpu_avail.sh
#
# Read-only. Shows allocated vs total GPUs per node so the fan-out submitter
# knows where to aim.
# ===========================================================================
set -uo pipefail

echo "======================================================================"
echo "GPU availability   $(date '+%F %T')"
echo "======================================================================"

printf "%-12s %-16s %-6s %-9s %-7s %-6s %s\n" \
       PARTITION NODE STATE "GPU free" "CPUfree" "MEM_MB" GRES
echo "----------------------------------------------------------------------"

TOTAL=0; FREE=0

while IFS='|' read -r part node state gres gres_used cpus_state mem; do
  [ -z "$node" ] && continue

  # gres  = gpu:l40:2          gres_used = gpu:l40:2(IDX:0-1)
  n_tot=$(echo "$gres"      | grep -oP 'gpu:[^:]+:\K[0-9]+' | head -1)
  n_use=$(echo "$gres_used" | grep -oP 'gpu:[^:]+:\K[0-9]+' | head -1)
  n_tot=${n_tot:-0}; n_use=${n_use:-0}
  n_free=$(( n_tot - n_use ))
  (( n_free < 0 )) && n_free=0

  # cpus_state = allocated/idle/other/total
  cpu_idle=$(echo "$cpus_state" | cut -d/ -f2)

  TOTAL=$(( TOTAL + n_tot )); FREE=$(( FREE + n_free ))

  flag=""
  (( n_free > 0 )) && flag=" <-- available"

  printf "%-12s %-16s %-6s %-9s %-7s %-6s %s%s\n" \
         "$part" "$node" "$state" "${n_free}/${n_tot}" \
         "${cpu_idle}" "$mem" "$gres" "$flag"
done < <(sinfo -h -N -o "%P|%N|%t|%G|%b|%C|%m" 2>/dev/null | sort -u)

echo "----------------------------------------------------------------------"
echo "GPUs free: $FREE / $TOTAL"
echo

echo "-- queue pressure by partition --"
squeue -h -o "%P" 2>/dev/null | sort | uniq -c | sort -rn | sed 's/^/  /'
[ -z "$(squeue -h 2>/dev/null)" ] && echo "  queue empty"
echo

echo "-- my jobs --"
squeue -h -u "$USER" -o "  %.10i %.12P %.18j %.2t %.10M %R" 2>/dev/null \
  || echo "  none"
[ -z "$(squeue -h -u "$USER" 2>/dev/null)" ] && echo "  none"
echo

echo "-- memory ceiling (RAM does NOT pool across nodes) --"
sinfo -h -o "%P %m" | sort -u | sort -k2 -n | awk '
  {printf "  %-14s %8s MB\n", $1, $2}
  END {print "  ^ a single memory-bound job can use at most ONE of these"}'
