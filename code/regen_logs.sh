#!/usr/bin/env bash
# Regenerates every certificate / lemma-test / v3 / v4 log in ../logs (sequential, single process).
# v2 logs: minutes on 2 cores.  v3 logs: about an hour in total; h1_limit fem needs a few GB (N <= 96, SuperLU).
set -e
cd "$(dirname "$0")"
mkdir -p ../logs
# --- v2
python3 lamstar_exact.py    2>&1 | tee ../logs/lamstar_exact_S3_S4singular.log
python3 lamstar_exact_v2.py 2>&1 | tee ../logs/lamstar_exact_S4prime.log
python3 -W ignore::DeprecationWarning lamstar.py 2>&1 | tee ../logs/lamstar_float.log
python3 lemma_tests.py 16,24,32,48,64 2>&1 | tee ../logs/lemma_tests.log
# --- v3: Clough--Tocher certificates (Section 9)
python3 lamstar_exact_ct.py 2>&1 | tee ../logs/lamstar_exact_ct.log
# --- v3: broader penalty necessity (Section 8.1)
python3 lamstar_family.py bubble    2>&1 | tee ../logs/lamstar_family_bubble.log
python3 lamstar_family.py certify   2>&1 | tee ../logs/lamstar_family_certify.log
python3 lamstar_family.py mesh      2>&1 | tee ../logs/lamstar_family_mesh.log
python3 lamstar_family.py limit     2>&1 | tee ../logs/lamstar_family_limit.log
python3 lamstar_family.py limitstar 2>&1 | tee ../logs/lamstar_family_limitstar.log
python3 lamstar_family.py levels    2>&1 | tee ../logs/lamstar_family_levels.log
# --- v3: H^1 leak constant (Section 6.1)
{ python3 h1_limit.py annulus; python3 h1_limit.py mps; } 2>&1 | tee ../logs/h1_limit_series.log
python3 h1_limit.py fem 2>&1 | tee ../logs/h1_limit_fem.log
# --- v3: lower bound (Section 7.1)
python3 lower_bound_tests.py stars  8,12,16,20,24,32,48,64,96,128,192,256,512,1024 2>&1 | tee ../logs/lower_bound_stars.log
python3 lower_bound_tests.py field  16,32,64,128,256,512 2>&1 | tee ../logs/lower_bound_field.log
python3 lower_bound_tests.py global 16,24,32,48,64 2>&1 | tee ../logs/lower_bound_global.log
# --- v4: growing penalty, mu-sweeps, flux-corrected data, conditioning (Sections 7.2 and 14.10 of the paper):
#     logs/rescue_{cond,flux,musweep,study}.log.  The JSONL records go to results/rerun/v4_rescue
#     so the shipped results/v4/rescue is not overwritten (peak 2.2 GB; N = 96 needs pypardiso).
RESCUE_OUT=../results/rerun/v4_rescue bash rescue_run.sh all
