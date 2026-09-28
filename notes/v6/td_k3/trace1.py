from patchlib import *
from fractions import Fraction as Fr
def ring(a,b,cs,q=(0,0,1),z=(0,0,0)):
    link=[a,b]+cs; m=len(link)
    return Patch([[q,z,link[j],link[(j+1)%m]] for j in range(m)],[(0,(1,2,3))],3)
q=(0,0,1); z=(0,0,0); a=(1,0,0); b=(0,1,0)
P=Patch([[q,z,a,b]],[(0,(1,2,3))],3); r,_=trace_rows(P,0); print('single: trace rank',rank_on_Y_mod(P,r))
P=Patch([[q,z,a,b],[q,z,(-1,-1,0),a]],[(0,(1,2,3)),(1,(1,2,3))],3); r,_=trace_rows(P,0); print('pair: trace rank on F',rank_on_Y_mod(P,r))
for cs in ([(-1,-1,Fr(1,2))],[(-1,Fr(1,3),Fr(1,2)),(Fr(1,4),-1,Fr(1,2))],[(-1,0,Fr(1,2)),(0,-1,Fr(1,2))]):
    P=ring(a,b,cs); r,bt=trace_rows(P,0); print('ring',len(cs)+2,': trace rank',rank_on_Y_mod(P,r))
    # g(z): beta=(2,0,0) (z is first face vertex in macro order q,z,a,b -> face order z,a,b)
    rz=[r[i] for i,(bb,c) in enumerate(bt) if bb==(2,0,0)]
    print('    rank of g(z):',rank_on_Y_mod(P,rz))
