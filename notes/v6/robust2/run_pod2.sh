#!/usr/bin/env bash
# Larger runs for REPORT.tex Sec. 5 (not run locally: 2 shared cores, 7 GB).
# Needs numpy, scipy (SuperLU), sympy; the repo's code/ directory (svn_k.py, svn.py).
# Estimated cost: N=128 SuperLU ~6-10 min/solve, ~4-6 GB; mucrit N=128 ~3 min, N=256 ~20 min / ~10 GB.
set -euo pipefail
cd "$(dirname "$0")"
python leak.py 128 sin 100,1000        > leak_128.jsonl        # leak law, normalised H1/energy at N=128
python slip_split.py 128 100 sin       > slip_split_128.jsonl  # does the edge-mean normal slip saturate at mu=100?
python mucrit.py 128,256               > mucrit_128_256.jsonl  # critical Nitsche penalty mu_c (limit ~80 expected)
python leak.py 128 r1_4 1000           > leak_128_r1_4.jsonl   # pure-normal potential: H1 order k+1/2 -> higher?
python leak.py cert 256 sin            > cert_Lk_256.jsonl
python leak.py certY 128 sin           > certY_128.jsonl
