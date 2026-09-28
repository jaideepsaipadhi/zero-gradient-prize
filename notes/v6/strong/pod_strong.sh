#!/bin/bash
# Larger strong-imposition levels for T3 (need > 7 GB RAM: N=256 ~ 7 GB, N=320 ~ 11 GB, N=384 ~ 16 GB RSS, PARDISO).
# Run from the repo root on a pod with MKL:  bash notes/v6/strong/pod_strong.sh
# Prediction (REPORT.md sec. 6): local H^1 rate 1.54 (192->256), 1.53 (256->320), 1.52 (320->384);
#   E/hG^1.5 -> ~0.90-0.94.  One-locked-vertex mesh: locked part L/hG -> ~0.68-0.70, rate -> 1.
set -e
export PYPARDISO_MKL_RT=${PYPARDISO_MKL_RT:-/usr/local/lib/libmkl_rt.so.3}
export SVN_OUT=$(pwd)/notes/v6/strong/runs
cd code
for N in 256 320 384; do python3 strong_bc.py $N A 100 zero pardiso std; done
cd ../notes/v6/strong
for N in 256 320; do python3 one_lock.py $N pardiso; done
python3 fit_rates.py
