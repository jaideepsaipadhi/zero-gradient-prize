#!/bin/sh
# mu-dependence of the GS leak (P test): mu = 200 (> muZ* = 148 at N=4) and mu = 1000
cd "$(dirname "$0")"
for MU in 200 1000; do python3 run3d.py sph 3 P $MU 2 4 > ../notes/v6/num3d/sph_k3_P_mu${MU}_local.log 2>&1; done
