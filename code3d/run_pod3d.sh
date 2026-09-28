#!/bin/sh
# 3D SVN convergence study -- pod script (the finer levels that do not fit the shared 7 GB / 2-core box).
# Run from the repo root:
#     sh code3d/run_pod3d.sh            (all stages, sequential)
#     sh code3d/run_pod3d.sh STAGE      (one stage: box4 | pen | k3a | k3b | k3c | pen4 | k4)
# Needs: python3, numpy, scipy, sympy, pypardiso + MKL (PYPARDISO_MKL_RT points to libmkl_rt.so.3).
# Results append to notes/v6/num3d/results.jsonl; tables: python3 code3d/table3d.py
#
# Memory / time (MEASURED locally, k=3, PARDISO LDL^T of the velocity-only IPM matrix, peak RSS):
#   sphere N=4 :  70k vel dofs  -> 1.0 GB, 30-110 s per method incl. IPM (2 shared cores)
#   box    N=6 :  79k vel dofs  -> 1.15 GB
# EXTRAPOLATED (3D direct-solver fill ~ n^{4/3}, n = velocity dofs), k=3 sphere, L=2, m=N/2:
#   N=6 : 233k dofs  ~ 5-6 GB,   N=8 : 550k dofs ~ 17 GB,   N=10 : 1.07M dofs ~ 40 GB,  N=12 : 1.85M ~ 80 GB.
#   k=4 multiplies the dof count by ~ (4/3)^3 = 2.4 at fixed mesh.
# Give each stage a machine with the memory listed; stages are independent.
set -e
cd "$(dirname "$0")/.."
export PYPARDISO_MKL_RT=${PYPARDISO_MKL_RT:-/usr/local/lib/libmkl_rt.so.3}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8} MKL_NUM_THREADS=${MKL_NUM_THREADS:-8}
LOG=notes/v6/num3d/pod_logs; mkdir -p $LOG
S=${1:-all}
run() { echo "== $*"; "$@"; }

# (1) k=4 box verification (manufactured solution, consistent Nitsche on the box; <= ~4 GB)
if [ $S = all -o $S = box4 ]; then
  run python3 code3d/run3d.py box 4 2 3 4 5 > $LOG/box_k4.log 2>&1
fi
# (2) Z_h coercivity threshold mu* on the sphere meshes (confirms mu=400 > mu*; N=6 ~ 6 GB, N=8 ~ 18 GB, ~15 factorizations each)
if [ $S = all -o $S = pen ]; then
  run python3 code3d/penalty3d.py 3 6 8 > $LOG/penalty_k3_6_8.log 2>&1
fi
# (3) k=3 sphere tests, mu = 400 (also writes the nodal-interpolation record): N=6 (~6 GB), N=8 (~18 GB)
if [ $S = all -o $S = k3a ]; then
  for K in P sph rot; do run python3 code3d/run3d.py sph 3 $K 400 6 8 > $LOG/sph_k3_${K}_6_8.log 2>&1; done
fi
# (4) k=3 sphere, N=10 (~40 GB) and N=12 (~80 GB)
if [ $S = all -o $S = k3b ]; then
  for K in P sph rot; do run python3 code3d/run3d.py sph 3 $K 400 10 12 > $LOG/sph_k3_${K}_10_12.log 2>&1; done
fi
# (5) mu-independence of the leak constant (slip ratio) at mu = 1000, N = 2..8
if [ $S = all -o $S = k3c ]; then
  run python3 code3d/run3d.py sph 3 P 1000 2 4 6 8 > $LOG/sph_k3_P_mu1000.log 2>&1
fi
# (6) k=4 sphere tests, mu = 1000 (local: muZ*(k=4,N=2) = 194, 2.5x the k=3 value; check with the pen4 stage),
#     N=2..6 (N=4 ~ 3.5 GB, N=6 ~ 15-20 GB)
if [ $S = all -o $S = pen4 ]; then
  run python3 code3d/penalty3d.py 4 4 6 > $LOG/penalty_k4_4_6.log 2>&1
fi
if [ $S = all -o $S = k4 ]; then
  for K in P sph rot; do run python3 code3d/run3d.py sph 4 $K 1000 2 4 6 > $LOG/sph_k4_${K}.log 2>&1; done
fi
python3 code3d/table3d.py > notes/v6/num3d/tables_after_pod.txt
