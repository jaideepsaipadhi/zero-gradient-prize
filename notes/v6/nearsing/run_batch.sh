#!/bin/bash
# one queue of FE runs (N <= 96: ~60 s, ~3 GB each; two concurrent N=96 runs exceed the 7 GB cgroup)
cd "$(dirname "$0")"
for a in "altsplit c0.3" "altsplit h1.0" "altsplit q1.0" "altwall c0.3" "altwall h1.0" "altwall q1.0" "onesplit q1.0" "onesplit h1.0" "onewall q1.0" "alt c0" "std c0"; do
  set -- $a; python3 ns_fe.py $1 $2 32 48 64 96 > runs/$1_$2.jsonl 2> runs/$1_$2.err
done
