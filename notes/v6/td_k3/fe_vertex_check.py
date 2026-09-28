"""Independent FE check (full Alfeld space, no reduced model) of the flat vertex formula.
Impose the interior-face data of the 'vertex move' (Bernstein coefficients on every interior macro face
f_j=[q,z,x_j]:  c[(q+2z)/3]=d/3, c[(2q+z)/3]=dp/3=-d/6, c[(q+z+x_j)/3]=p_j/6=-d/6), zero mean on every Gamma
face, and compute the total moment over the affine solution set (floats, lstsq + nullspace).
Equal weights (1,1,1) (then D_i=0): prediction  -(1/1440) sum_i (4|F_i|/H) d  exactly (no freedom)."""
from patchlib import *
from fractions import Fraction as Fr
import numpy as np
def check(link,q,gam,d):
    z=(0,0,0); m=len(link)
    P=Patch([[q,z,link[j],link[(j+1)%m]] for j in range(m)],[(j,(1,2,3)) for j in gam],3,weights={fi:(1,1,1) for fi in range(len(gam))})
    qF,zF=F3(q),F3(z)
    rows=[];rhs=[]
    D=P.dense(P.div_rows,'float'); rows.append(D); rhs.append(np.zeros(D.shape[0]))
    C=P.dense(P.cons_rows(),'float'); rows.append(C); rhs.append(np.zeros(C.shape[0]))
    fix=[]
    for j in range(m):
        x=F3(link[j])
        for pt,val in [(tuple((qF[c]+2*zF[c])/3 for c in range(3)),[dc/3 for dc in d]),
                       (tuple((2*qF[c]+zF[c])/3 for c in range(3)),[-dc/6 for dc in d]),
                       (tuple((qF[c]+zF[c]+x[c])/3 for c in range(3)),[-dc/6 for dc in d])]:
            n=P.key[pt]
            for c in range(3):
                r=np.zeros(P.nunk); r[3*n+c]=1; fix.append((r,float(val[c])))
    rows.append(np.array([r for r,_ in fix])); rhs.append(np.array([v for _,v in fix]))
    A=np.vstack(rows); b=np.concatenate(rhs)
    x0,res,rk,sv=np.linalg.lstsq(A,b,rcond=None)
    resid=np.linalg.norm(A@x0-b)
    M=P.dense(P.mom_rows([(1,0,0),(0,1,0)]),'float')
    Nul=null_float(A)
    spread=np.linalg.norm(M@Nul) if Nul.shape[1] else 0.0
    H=float(q[2]); area=lambda i: abs(float(link[i][0]*link[(i+1)%m][1]-link[i][1]*link[(i+1)%m][0]))/2
    kap=sum(4*area(i)/H for i in gam)/1440
    return M@x0, -kap*np.array(d[:2],dtype=float), resid, spread
tests=[([(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)],(Fr(1,5),Fr(1,7),1),[0,1,2,3]),
       ([(1,0,0),(Fr(1,2),Fr(7,8),0),(-Fr(1,2),Fr(7,8),0),(-1,0,0),(-Fr(1,2),-Fr(7,8),0),(Fr(1,2),-Fr(7,8),0)],(0,0,1),[0,1,2,3,4,5]),
       ([(1,0,0),(0,1,0),(-1,Fr(1,3),Fr(1,2)),(-1,-1,0),(Fr(1,2),-1,0)],(Fr(1,5),Fr(1,7),Fr(4,5)),[0,3]),
       ([(1,Fr(1,9),0),(Fr(-1,3),1,0),(-1,-1,Fr(1,2))],(Fr(1,5),Fr(-1,7),1),[0])]
for link,q,gam in tests:
    for d in [(1,0,0),(0,1,0)]:
        got,pred,resid,spread=check(link,q,gam,d)
        print('gam',gam,'d',d[:2],'FE moment',got.round(10),'predicted',pred.round(10),'residual %.1e'%resid,'moment spread over solution set %.1e'%spread)
