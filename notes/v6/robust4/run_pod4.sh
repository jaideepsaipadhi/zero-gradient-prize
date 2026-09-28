#!/bin/bash
# robust4 pod job (needs ~10-20 GB RAM; N=128 peaks above the 7 GB sandbox).  ~1-2 h.
# H^1 order of W_*(phi) on the uniform production-type mesh for the invisible (B, R), K-degenerate (A) and visible (sin, Q)
# potentials, N = 128, 192, 256 at mu = 1000 (coercive regime).
cd "$(dirname "$0")"
for N in 128 192 256; do
  for p in B A R sin Q; do
    python run_u.py $N $p 1000 >> u_pod.jsonl 2>> u_pod.err
  done
done
