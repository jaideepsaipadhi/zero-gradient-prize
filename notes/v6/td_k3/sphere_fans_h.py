"""FE (full Alfeld space, k=3) energy-normalised constant of the ring witness on sphere-inscribed fans:
z = north pole, m ring points on the unit sphere (rational, via inverse stereographic projection) at chord ~h,
apex q=(1+H h) z, all m faces on Gamma.  kappa_hat = sigma_min(tangential total moment on Y^1)/h^{5/2}."""
from patchlib import *
from fractions import Fraction as Fr
import numpy as np
def sph(u,v):
    s=1+u*u+v*v; return (2*u/s,2*v/s,(1-u*u-v*v)/s)
fr=[(1,0,0),(0,1,0)]
for m in (5,6):
  for H in (Fr(1,2),Fr(1)):
    out=[]
    for h in (Fr(1,2),Fr(1,4),Fr(1,8),Fr(1,16)):
        pts=[]
        for j in range(m):
            th=2*np.pi*j/m+0.1*np.sin(3*j)
            u=Fr(float(h)/2*np.cos(th)).limit_denominator(4000); v=Fr(float(h)/2*np.sin(th)).limit_denominator(4000)
            pts.append(sph(u,v))
        z=(0,0,1); q=(0,0,1+H*h)
        P=Patch([[q,z,pts[j],pts[(j+1)%m]] for j in range(m)],[(j,(1,2,3)) for j in range(m)],3)
        s=kappa_float(P,fr); rk=moment_rank_mod(P,fr)[2]
        out.append(f'h={float(h):.4f}: rank {rk}, kappa_hat={s[-1]/float(h)**2.5:.4f}')
    print(f'm={m} H={float(H)}: '+' | '.join(out),flush=True)
