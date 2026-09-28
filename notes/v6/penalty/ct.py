"""Clough-Tocher split of a patch; float lambda* via patch.lamstar."""
import numpy as np
from patch import lamstar
def ct_split(P, tris, nit):
    P = [tuple(map(float,p)) for p in P]; T2=[]
    for (a,b,c) in tris:
        g = tuple((np.array(P[a])+np.array(P[b])+np.array(P[c]))/3); P.append(g); gi=len(P)-1
        T2 += [(a,b,gi),(b,c,gi),(c,a,gi)]
    return P, T2, nit
def ct_fan(rim,k):
    P=[(0.,0.)]+list(rim); m=len(rim)-1
    tris=[(0,i+1,i+2) for i in range(m)]
    nit={(0,1):(0,-1),(0,m+1):(0,-1)}
    return lamstar(*ct_split(P,tris,nit),k)
if __name__=="__main__":
    for H in [1,2,4,6,8,12,16]:
        print(H,[(k,round(ct_fan([(1,0),(1,H),(0,H),(-1,0)],k)[0],3)) for k in (2,3,4)],
              "single macro tri k=2:", round(lamstar(*ct_split([(0,0),(1,0),(1,H)],[(0,1,2)],{(0,1):(0,-1)}),2)[0],3), flush=True)
