"""Ring = star of the interior edge [z,q] (z on Gamma, q the apex of the macro-tet of F=[z,a,b]).
Only F is a Gamma face in the ring (other link vertices c_i are above the plane).
EQUAL weights on F: single-macro fields contribute nothing; do the ring fields?"""
from patchlib import *
from fractions import Fraction as Fr
import numpy as np, sys
def ring(a,b,cs,q=(0,0,1),z=(0,0,0),w=(1,1,1)):
    link=[a,b]+cs
    m=len(link)
    mac=[[q,z,link[j],link[(j+1)%m]] for j in range(m)]
    return Patch(mac,[(0,(1,2,3))],3,weights={0:w} if w else None)
fr=[(1,0,0),(0,1,0)]
tests={
 'r=2 c above':((1,0,0),(0,1,0),[(-1,Fr(1,3),Fr(1,2)),(Fr(1,4),-1,Fr(1,2))]),
 'r=1 c above':((1,0,0),(0,1,0),[(-1,-1,Fr(1,2))]),
 'r=3 c above':((1,0,0),(0,1,0),[(-1,Fr(1,2),Fr(1,3)),(-1,-1,Fr(1,2)),(Fr(1,2),-1,Fr(1,4))]),
 'r=3 generic q':((1,Fr(1,9),0),(Fr(-1,3),1,0),[(-1,Fr(1,2),Fr(1,3)),(-1,-1,Fr(1,2)),(Fr(1,2),-1,Fr(1,4))]),
}
for nm,(a,b,cs) in tests.items():
    for w in ((1,1,1),None):
        P=ring(a,b,cs,w=w)
        s,d1=kappa_float(P,fr,return_all=True)
        print(f'{nm:15s} w={"equal" if w else "l^2  "} dimY1={d1}', moment_rank_mod(P,fr), s.round(6))
