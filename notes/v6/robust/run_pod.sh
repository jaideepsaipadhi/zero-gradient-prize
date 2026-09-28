#!/usr/bin/env bash
# Larger hydrostatic runs (N = 96, 128) for REPORT.md Sec. R4. Not run locally (2 shared cores, no MKL).
# Needs: numpy, scipy, sympy, pypardiso + MKL (libmkl_rt), as for the paper's svn_k.py PARDISO runs.
# Expected cost: SuperLU at N=64 was ~15 s/solve on 2 cores; PARDISO at N=128 should be ~1 min/solve, ~8 GB peak.
# Note: the PARDISO path regularises the pressure block by 1e-12*M (as svn_k.solve_pardiso), which puts a
# floor ~1e-11 on quantities that are exactly zero (CNS*-R; cubic phi with CNS*).
set -euo pipefail
cd "$(dirname "$0")"
export PYPARDISO_MKL_RT=${PYPARDISO_MKL_RT:-/usr/local/lib/libmkl_rt.so.3}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
HYDRO_SOLVER=pardiso python hydro.py hydrolite 96,128 4 > hydro_96_128.jsonl
HYDRO_SOLVER=pardiso python hydro.py rates 96,128 4 > rates_96_128.jsonl
