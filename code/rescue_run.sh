#!/bin/bash
# Batch driver for code/rescue.py (notes/v4/rescue.tex). Sequential, 1 thread; peak RSS 2.2 GB (N = 64, SuperLU);
# N = 96 uses PARDISO (1.2 GB; SuperLU would need 5.6 GB).
#   bash code/rescue_run.sh [study|musweep|flux|cond|all]
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
what=${1:-all}
LOG=../logs
mkdir -p $LOG

if [[ $what == cond || $what == all ]]; then
  for N in 16 24 32 48; do
    python3 rescue.py cond $N 100,300,1000,3000,10000,100000
  done > $LOG/rescue_cond.log 2>&1
fi

if [[ $what == flux || $what == all ]]; then
  for K in F F2; do for N in 16 24 32 48 64; do
    python3 rescue.py flux $K $N 1e2,1e4,1e6,1e8 && echo "flux $K $N done"
  done; done > $LOG/rescue_flux.log 2>&1
fi

if [[ $what == musweep || $what == all ]]; then
  for K in A B1; do for N in 32 64; do
    python3 rescue.py musweep $K $N 1e2,1e3,1e4,1e5,1e6,1e7,1e8 && echo "musweep $K $N done"
  done; done > $LOG/rescue_musweep.log 2>&1
fi

if [[ $what == study || $what == all ]]; then
  for N in 16 24 32 48 64 96; do
    for K in A P1 B1 B100; do
      for AC in "0 100" "0.5 100" "0.5 1000" "1 100"; do
        set -- $AC
        if [[ $N == 96 ]]; then export RESCUE_SOLVER=pardiso; else export RESCUE_SOLVER=lu; fi
        python3 rescue.py study $K $N $1 $2 && echo "study $K N=$N alpha=$1 c=$2 done"
      done
    done
  done > $LOG/rescue_study.log 2>&1
fi
