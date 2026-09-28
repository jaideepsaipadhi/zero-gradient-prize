"""Flat rings: singular-edge configurations (faces around [z,q] in two planes), equal weights."""
from patchlib import *
from fractions import Fraction as Fr
from ring1 import ring
fr=[(1,0,0),(0,1,0)]
for h1,h2 in [(Fr(1,2),Fr(1,2)),(0,Fr(1,2)),(Fr(1,3),Fr(2,3))]:
    P=ring((1,0,0),(0,1,0),[(-1,0,h1),(0,-1,h2)])
    print('singular ring (two planes), c-heights',h1,h2, moment_rank_mod(P,fr), kappa_float(P,fr).round(6))
# singular, but with q tilted within... q generic -> planes z,q,a etc.
q=(Fr(1,5),Fr(1,7),1)
# choose c1 on plane(z,q,a) and c2 on plane(z,q,b): c = -s*a + t*q
a=(1,0,0); b=(0,1,0)
c1=tuple(-a[i]+Fr(1,2)*q[i] for i in range(3)); c2=tuple(-b[i]+Fr(1,2)*q[i] for i in range(3))
P=ring(a,b,[c1,c2],q=q)
print('singular ring, tilted q', moment_rank_mod(P,fr), kappa_float(P,fr).round(6))
