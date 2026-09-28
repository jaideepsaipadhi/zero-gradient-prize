#!/bin/sh
cd "$(dirname "$0")"
while pgrep -f "run3d.py box 4" > /dev/null; do sleep 10; done
python3 patch_test.py > ../notes/v6/num3d/patch_test.log 2>&1
for K in P sph rot; do python3 run3d.py interp 3 $K 400 2 4 > ../notes/v6/num3d/interp_k3_$K.log 2>&1; done
python3 -c "
import svn3d
for k in (3,4):
    for N in (2,4):
        S=svn3d.build('box',N,k,L=1.0); print('box k=%d N=%d local sufficient mu (V_h, sym flux) = %.2f'%(k,N,svn3d.local_threshold(S,S.fBox,'sym').max()),flush=True)
    for N in (2,4,6,8):
        S=svn3d.build('sph',N,k,L=2.0); print('sph k=%d N=%d local sufficient mu (V_h, grad flux) = %.2f  h/hG=%.3f'%(k,N,svn3d.local_threshold(S,S.fGamma,'grad').max(),S.h/S.hGamma),flush=True)
" > ../notes/v6/num3d/local_threshold.log 2>&1
