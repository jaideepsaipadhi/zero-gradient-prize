#!/usr/bin/env bash
# Regenerates the certificate / lemma-test logs in ../logs (sequential; ~minutes to an hour on 2 cores).
set -e
cd "$(dirname "$0")"
mkdir -p ../logs
python3 lamstar_exact.py    2>&1 | tee ../logs/lamstar_exact_S3_S4singular.log
python3 lamstar_exact_v2.py 2>&1 | tee ../logs/lamstar_exact_S4prime.log
python3 lamstar.py          2>&1 | tee ../logs/lamstar_float.log
python3 lemma_tests.py 16,24,32,48,64 2>&1 | tee ../logs/lemma_tests.log
