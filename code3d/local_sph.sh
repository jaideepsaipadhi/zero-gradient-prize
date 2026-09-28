#!/bin/sh
# local (shared 7 GB box) sphere runs, k=3, mu=400, N=2,4
cd "$(dirname "$0")"
for K in P sph rot; do python3 run3d.py sph 3 $K 400 2 4 > ../notes/v6/num3d/sph_k3_${K}_mu400_local.log 2>&1; done
