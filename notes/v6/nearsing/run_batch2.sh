#!/bin/bash
# strictly after run_batch.sh (waits on its PID): the OOM-killed N=96 runs, then the inf-sup sweeps (N=16/24, dense)
cd "$(dirname "$0")"
while kill -0 3625 2>/dev/null; do sleep 5; done
for a in "altsplit c0.3" "altsplit h1.0" "altsplit q1.0" "altwall c0.3" "altwall h1.0"; do
  set -- $a; n=$(grep -c '"N": 96' runs/$1_$2.jsonl); [ "$n" = 0 ] && python3 ns_fe.py $1 $2 96 >> runs/$1_$2.jsonl 2>> runs/$1_$2.err
done
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
python3 ns_infsup.py onewall 16 0.3 0.03 0.003 0.0003 > runs/infsup_onewall16.jsonl 2>&1
python3 ns_infsup.py onesplit 16 0.3 0.1 0.03 0.01 0.003 > runs/infsup_onesplit16.jsonl 2>&1
python3 ns_infsup.py onewall 24 0.3 0.03 0.003 > runs/infsup_onewall24.jsonl 2>&1
