#!/bin/sh
# second local batch: P rerun with the tightened IPM stopping rule, inf-sup checks, k=4 small levels
cd "$(dirname "$0")"
python3 run3d.py sph 3 P 400 4 > ../notes/v6/num3d/sph_k3_P_mu400_N4_rerun.log 2>&1
python3 infsup3d.py box 3 1 2 3 > ../notes/v6/num3d/infsup_box_k3.log 2>&1
python3 infsup3d.py sph 3 2 > ../notes/v6/num3d/infsup_sph_k3.log 2>&1
python3 penalty3d.py 4 2 > ../notes/v6/num3d/penalty_k4.log 2>&1
python3 run3d.py box 4 2 3 4 > ../notes/v6/num3d/box_k4.log 2>&1
