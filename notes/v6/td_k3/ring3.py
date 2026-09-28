"""Joint edge-moment map (R^6) of the single Gamma face F on ring patches."""
from patchlib import *
from fractions import Fraction as Fr
def ring(a,b,cs,q=(0,0,1),z=(0,0,0)):
    link=[a,b]+cs; m=len(link)
    return Patch([[q,z,link[j],link[(j+1)%m]] for j in range(m)],[(0,(1,2,3))],3)
tests={
 'm=3':((1,0,0),(0,1,0),[(-1,-1,Fr(1,2))]),
 'm=4':((1,0,0),(0,1,0),[(-1,Fr(1,3),Fr(1,2)),(Fr(1,4),-1,Fr(1,2))]),
 'm=4 singular':((1,0,0),(0,1,0),[(-1,0,Fr(1,2)),(0,-1,Fr(1,2))]),
 'm=5':((1,0,0),(0,1,0),[(-1,Fr(1,2),Fr(1,3)),(-1,-1,Fr(1,2)),(Fr(1,2),-1,Fr(1,4))]),
 'm=6':((1,0,0),(0,1,0),[(-1,Fr(1,2),Fr(1,3)),(-1,-Fr(1,3),Fr(1,2)),(-Fr(1,4),-1,Fr(1,2)),(Fr(1,2),-1,Fr(1,4))]),
}
for nm,(a,b,cs) in tests.items():
    P=ring(a,b,cs)
    print(nm,'single-face ring: rank of joint edge-moment map on Y^1 =',rank_on_Y1_mod(P,edge_moment_rows(P,0))[1])
