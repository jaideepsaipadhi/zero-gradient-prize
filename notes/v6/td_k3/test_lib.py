import time
from patchlib import *
t0=time.time()
T=[(0,0,3),(0,0,0),(5,1,0),(1,4,0)]
for k in (3,4):
    P=Patch([T],[(0,(1,2,3))],k)
    print('single general k=%d'%k, 'nunk',P.nunk, moment_rank_mod(P,[(1,0,0),(0,1,0)]), kappa_float(P,[(1,0,0),(0,1,0)]))
ring=[(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)]
mac=[[(0,0,2),(0,0,1),ring[j],ring[(j+1)%4]] for j in range(4)]
P=Patch(mac,[(j,(1,2,3)) for j in range(4)],3)
print('octa fan', P.nunk, moment_rank_mod(P,[(1,0,0),(0,1,0)]), kappa_float(P,[(1,0,0),(0,1,0)]))
from fractions import Fraction as Fr
q=(Fr(1,5),Fr(1,7),2)
P=Patch([[q,(0,0,1),ring[j],ring[(j+1)%4]] for j in range(4)],[(j,(1,2,3)) for j in range(4)],3)
print('octa offaxis', moment_rank_mod(P,[(1,0,0),(0,1,0)]), kappa_float(P,[(1,0,0),(0,1,0)]))
print(time.time()-t0)
