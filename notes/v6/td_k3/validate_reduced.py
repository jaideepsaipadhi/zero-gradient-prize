"""Validate the reduced model against the full finite-element computation: joint trace rank and
total-moment rank must coincide (the report proves J = J_red)."""
from patchlib import *
from reduced import Reduced
from fractions import Fraction as Fr
import time
def fullP(link,q,z,gam):
    m=len(link); return Patch([[q,z,link[j],link[(j+1)%m]] for j in range(m)],[(j,(1,2,3)) for j in gam],3)
cases=[
 ('single ring m=3',[(1,0,0),(0,1,0),(-1,-1,Fr(1,2))],(0,0,1),(0,0,0),[0]),
 ('single ring m=5 generic q',[(1,Fr(1,9),0),(Fr(-1,3),1,0),(-1,Fr(1,2),Fr(1,3)),(-1,-1,Fr(1,2)),(Fr(1,2),-1,Fr(1,4))],(Fr(1,5),Fr(-1,7),1),(0,0,0),[0]),
 ('flat square fan',[(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)],(0,0,1),(0,0,0),[0,1,2,3]),
 ('octa fan',[(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)],(0,0,2),(0,0,1),[0,1,2,3]),
 ('octa fan off-axis',[(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)],(Fr(1,5),Fr(1,7),2),(0,0,1),[0,1,2,3]),
 ('2 adjacent flat',[(1,0,0),(0,1,0),(-1,0,0),(-Fr(1,2),-1,Fr(1,2)),(Fr(1,2),-1,Fr(1,3))],(0,0,1),(0,0,0),[0,1]),
 ('2 non-adjacent',[(1,0,0),(0,1,0),(-1,Fr(1,3),Fr(1,2)),(-1,-1,0),(Fr(1,2),-1,0)],(Fr(1,5),Fr(1,7),1),(0,0,0),[0,3]),
 ('pentagon fan',[(1,0,0),(Fr(3,10),1,0),(-Fr(4,5),Fr(3,5),0),(-Fr(4,5),-Fr(3,5),0),(Fr(3,10),-1,0)],(Fr(1,7),Fr(1,9),1),(0,0,0),[0,1,2,3,4]),
]
fr=[(1,0,0),(0,1,0)]
for nm,link,q,z,gam in cases:
    t0=time.time()
    R=Reduced(link,q,z,gam); jd,mr,_=R.analyse(fr)
    P=fullP(link,q,z,gam)
    T=sum([trace_rows(P,i)[0] for i in range(len(gam))],[])
    fj=rank_on_Y_mod(P,T); fm=moment_rank_mod(P,fr)[2]
    print(f'{nm:28s} reduced: J dim {jd}, moment rank {mr} | FE: joint trace rank {fj}, moment rank {fm}  ({time.time()-t0:.1f}s)')
