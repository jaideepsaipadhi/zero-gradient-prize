import numpy as np, sys
from patch import lamstar
from ct import ct_split
def strip_mesh(W,J,H,shear=0.0):
    idx=lambda i,j:j*(W+1)+i
    P=[(i+shear*j*H, j*H) for j in range(J+1) for i in range(W+1)]
    tris=[]
    for j in range(J):
        for i in range(W):
            a,b,c,d=idx(i,j),idx(i+1,j),idx(i+1,j+1),idx(i,j+1)
            tris+=[(a,b,c),(a,c,d)]
    nit={(idx(i,0),idx(i+1,0)):(0,-1) for i in range(W)}
    return P,tris,nit
def lam(W,J,H,k,ct=False,shear=0.0):
    P,t,n=strip_mesh(W,J,H,shear)
    if ct: P,t,n=ct_split(P,t,n)
    return lamstar(P,t,n,k)[0]
if __name__=="__main__":
    mode=sys.argv[1]
    if mode=="p4":
        for H in [6,8,12]:
            print("P4 H",H,[ (W,round(lam(W,1,H,4),3)) for W in (2,4,8,16,24)],flush=True)
    if mode=="ct2":
        for H in [2,4,8]:
            print("CT2 H",H,[ (W,round(lam(W,1,H,2,True),3)) for W in (2,4,8,16)],flush=True)
