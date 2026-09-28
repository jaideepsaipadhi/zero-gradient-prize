"""Referee 8: independent FULL-FE check (patchlib FE, not the reduced model) of the NON-FLAT vertex formula
  sum_i M_i(v_delta) = -(1/360) sum_i c_i G_i delta   (mod the span of moments of the remaining freedom)
on sphere-inscribed stars with true sphere weights l^2, rational points on the unit sphere."""
import sys; sys.path.insert(0,'../td_k3')
from patchlib import *
from fractions import Fraction as Fr
import numpy as np
def sph(s,t):
    s,t=Fr(s),Fr(t); D=1+s*s+t*t; return (2*s/D,2*t/D,(1-s*s-t*t)/D)
def pred(link,q,z,gam,d):
    link=[np.array([float(c) for c in x]) for x in link]; q=np.array([float(c) for c in q]); z=np.array([float(c) for c in z])
    d=np.array(d,float); m=len(link); tot=np.zeros(3)
    for i in gam:
        a=link[i]; b=link[(i+1)%m]; N=np.cross(a-z,b-z); area=np.linalg.norm(N)/2; n=N/np.linalg.norm(N)
        P=np.eye(3)-np.outer(n,n); r=(a-z)+(b-z)+P@(q-z); G=P-np.outer(r,n)/((q-z)@n)
        H=abs((q-z)@n); c=area*((a-b)@(a-b))/H
        tot+=-c/360*(G@d)
    return tot
def run(link,q,z,gam):
    m=len(link)
    P=Patch([[q,z,link[j],link[(j+1)%m]] for j in range(m)],[(j,(1,2,3)) for j in gam],3)
    qF,zF=F3(q),F3(z)
    D=P.dense(P.div_rows,'float'); C=P.dense(P.cons_rows(),'float')
    M=P.dense(P.mom_rows([(1,0,0),(0,1,0),(0,0,1)]),'float')
    for d in [(1,0,0),(0,1,0),(0,0,1)]:
        fix=[]
        for j in range(m):
            x=F3(link[j])
            for pt,val in [(tuple((qF[c]+2*zF[c])/3 for c in range(3)),[Fr(dc)/3 for dc in d]),
                           (tuple((2*qF[c]+zF[c])/3 for c in range(3)),[-Fr(dc)/6 for dc in d]),
                           (tuple((qF[c]+zF[c]+x[c])/3 for c in range(3)),[-Fr(dc)/6 for dc in d])]:
                n=P.key[pt]
                for c in range(3):
                    r=np.zeros(P.nunk); r[3*n+c]=1; fix.append((r,float(val[c])))
        F=np.array([r for r,_ in fix]); fb=np.array([v for _,v in fix])
        A=np.vstack([D,C,F]); b=np.concatenate([np.zeros(len(D)+len(C)),fb])
        x0=np.linalg.lstsq(A,b,rcond=None)[0]; res=np.linalg.norm(A@x0-b)
        Nul=null_float(A); MN=M@Nul
        pr=pred(link,q,z,gam,d); diff=pr-M@x0
        # is diff in range(MN)?
        if MN.shape[1] and np.linalg.norm(MN)>1e-13:
            y=np.linalg.lstsq(MN,diff,rcond=None)[0]; out=np.linalg.norm(MN@y-diff)
            rk=np.linalg.matrix_rank(MN,tol=1e-10*max(1,np.linalg.norm(MN)))
        else: out=np.linalg.norm(diff); rk=0
        print(f' d={d} resid {res:.1e}  FE {np.round(M@x0,8)}  pred {np.round(pr,8)}  free-rank {rk}  dist(pred, FE affine image) {out:.1e}',flush=True)
z=(0,0,1)
cases=[
 ('sphere pentagon, all 5 Gamma',[sph(Fr(1,5),0),sph(Fr(1,16),Fr(3,16)),sph(Fr(-3,20),Fr(1,8)),sph(Fr(-1,7),Fr(-1,9)),sph(Fr(1,20),Fr(-1,5))],(Fr(1,30),Fr(-1,40),Fr(5,4)),[0,1,2,3,4]),
 ('sphere quad, faces 0,2 only',[sph(Fr(1,4),Fr(1,30)),sph(0,Fr(1,5)),sph(Fr(-1,5),Fr(1,50)),(Fr(1,20),Fr(-1,4),Fr(9,10))],(Fr(-1,30),Fr(1,25),Fr(13,10)),[0,2]),
 ('coarse tilted hexagon (tau big)',[sph(Fr(1,2),0),sph(Fr(1,4),Fr(1,2)),sph(Fr(-1,4),Fr(1,2)),sph(Fr(-1,2),0),sph(Fr(-1,4),Fr(-1,2)),sph(Fr(1,4),Fr(-1,2))],(Fr(1,10),0,Fr(3,2)),[0,1,2,3,4,5]),
]
for nm,link,q,gam in cases:
    print(nm,flush=True); run(link,q,z,gam)
