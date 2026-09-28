"""Regular m-fans with the Gamma-faces tilted: ring (cos th_j, sin th_j, -eps), z = 0, apex q = (0,0,H).
eps = 0 is the flat limit; m=4, eps=1, H=1 is the octahedral fan (rank 0).  Exact ring coordinates are
rational points on the unit circle.  Frame = horizontal plane (tangent plane at z of the sphere through
the vertices, which is horizontal by symmetry)."""
from patchlib import *
from fractions import Fraction as Fr
import numpy as np, sys
def ratcirc(m):
    # rational points on unit circle near angles 2 pi j/m (via t=tan(th/2))
    pts=[]
    for j in range(m):
        th=2*np.pi*j/m
        if abs(np.cos(th/2))<1e-12: pts.append((Fr(-1),Fr(0))); continue
        tt=Fr(np.tan(th/2)).limit_denominator(1000)
        pts.append(((1-tt*tt)/(1+tt*tt),2*tt/(1+tt*tt)))
    return pts
fr=[(1,0,0),(0,1,0)]
for m in (4,5,6):
    C=ratcirc(m)
    for H in (Fr(1,2),Fr(1),Fr(2)):
        row=[]
        for eps in [Fr(0),Fr(1,10),Fr(1,4),Fr(1,2),Fr(3,4),Fr(1),Fr(5,4),Fr(3,2)]:
            ring=[(c[0],c[1],-eps) for c in C]
            P=Patch([[(0,0,H),(0,0,0),ring[j],ring[(j+1)%m]] for j in range(m)],[(j,(1,2,3)) for j in range(m)],3)
            s=kappa_float(P,fr); rk=moment_rank_mod(P,fr)[2]
            row.append(f'{float(eps):.2f}:{s[-1]:.4f}(r{rk})')
        print(f'm={m} H={float(H)}:',' '.join(row),flush=True)
