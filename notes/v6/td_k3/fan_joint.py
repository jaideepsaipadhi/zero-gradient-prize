"""Fans (ring around [z,q] with ALL link faces on Gamma): joint trace rank, and the trace space on F_0
of fields whose traces vanish on the other Gamma faces."""
from patchlib import *
from fractions import Fraction as Fr
import numpy as np
def fan(ring,q,z=(0,0,0)):
    m=len(ring); return Patch([[q,z,ring[j],ring[(j+1)%m]] for j in range(m)],[(j,(1,2,3)) for j in range(m)],3)
cases={'flat square, q on axis':([(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)],(0,0,1)),
 'flat pentagon-ish':([(1,0,0),(Fr(3,10),1,0),(-Fr(4,5),Fr(3,5),0),(-Fr(4,5),-Fr(3,5),0),(Fr(3,10),-1,0)],(Fr(1,7),Fr(1,9),1)),
 'octa (non-flat), q on axis':([(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)],(0,0,2)),}
for nm,(ring,q) in cases.items():
    z=(0,0,0) if 'octa' not in nm else (0,0,1)
    P=fan(ring,q,z); m=len(ring)
    T=[trace_rows(P,i)[0] for i in range(m)]
    alltr=sum(T,[])
    print(nm,': joint trace rank',rank_on_Y_mod(P,alltr),'(7m =',7*m,')', ' single-face trace ranks',[rank_on_Y_mod(P,T[i]) for i in range(m)])
    # isolated: fields with zero trace on faces 1..m-1 -> rank of trace on face 0
    D=P.dense(P.div_rows,'mod'); O=P.dense(sum(T[1:],[]),'mod'); T0=P.dense(T[0],'mod')
    base=rank_mod(np.vstack([D,O])); print('   isolated trace rank on F0 =',rank_mod(np.vstack([D,O,T0]))-base)
