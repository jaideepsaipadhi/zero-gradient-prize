"""Two macro-tets sharing the interior face [q,z,a1] which contains the boundary edge E=[z,a1]; flat."""
from patchlib import *
from fractions import Fraction as Fr
import numpy as np
def pair(a0,a1,a2,q,z=(0,0,0),w=None):
    mac=[[q,z,a1,a2],[q,z,a0,a1]]; gf=[(0,(1,2,3)),(1,(1,2,3))]
    return Patch(mac,gf,3,weights=w)
fr=[(1,0,0),(0,1,0)]
cases={'sym':((-1,-1,0),(1,0,0),(-1,1,0),(0,0,1)),
       'gen':((-Fr(1,2),-1,0),(1,Fr(1,7),0),(-Fr(3,4),Fr(6,5),0),(Fr(1,5),Fr(1,9),1)),
       'gen2':((Fr(1,3),-1,0),(1,0,0),(Fr(1,4),1,0),(Fr(1,2),0,Fr(3,4)))}
for nm,(a0,a1,a2,q) in cases.items():
    for w in (None,{0:(1,1,1),1:(1,1,1)}):
        P=pair(a0,a1,a2,q,w=w)
        D=P.dense(P.div_rows,'float'); C=P.dense(P.cons_rows(),'float')
        from patchlib import null_float
        Ny=null_float(D); Ny1=Ny@null_float(C@Ny)
        M=P.dense(P.mom_rows(fr),'float')@Ny1
        u,s,vt=np.linalg.svd(M)
        print(nm,'weights','l^2' if w is None else 'equal','dimY',Ny.shape[1],'dimY1',Ny1.shape[1],moment_rank_mod(P,fr),'sv',s.round(6),'left sing vecs',u.round(4).tolist())
