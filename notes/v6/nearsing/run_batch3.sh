#!/bin/bash
# after run_batch2.sh: N=80 for the families whose N=96 runs were OOM-killed (shared memory cgroup), one at a time
cd "$(dirname "$0")"
while pgrep -f run_batch2.sh > /dev/null; do sleep 5; done
for a in "altsplit h1.0" "altwall h1.0" "altsplit q1.0" "altsplit c0.3" "altwall c0.3" "std c0" "alt c0" "altwall q1.0"; do
  set -- $a; python3 ns_fe.py $1 $2 80 >> runs/$1_$2.jsonl 2>> runs/$1_$2.err
done
