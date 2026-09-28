#!/bin/sh
# Larger levels of the L2 diagnostics (notes/v6/l2/REPORT.md, section 6).  Not run locally: N = 128 slip needs
# more than the ~6 GB free on the 2-core machine (it was OOM-killed there).  On the 32-core pod, from the repo root:
#
#   sh notes/v6/l2/run_pod.sh > notes/v6/l2/pod.log 2>&1
#
# Each line appends JSON records to notes/v6/l2/l2_results.jsonl.  Expected cost (SuperLU, 1 thread per job):
# N=128 ~3 min / 8 GB, N=192 ~10 min / 20 GB, N=256 ~25 min / 40 GB per slip job (two factorisations).
set -e
cd "$(dirname "$0")"
for N in 128 192 256; do
  python3 -W ignore l2_numerics.py slip "$N"
  python3 -W ignore l2_numerics.py moments "$N" 100
  python3 -W ignore l2_numerics.py moments "$N" 10000
  python3 -W ignore l2_numerics.py leak "$N" 100
  python3 -W ignore l2_numerics.py leak "$N" 1000
done
