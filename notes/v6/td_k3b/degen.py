"""EXACT: for a flat edge star R(z,q) (z=0, link in the plane x3=0), Gamma-face subset I, and sagitta weights
w_ij = h(e_ij,e_ij) with h a symmetric 2x2 form (h11,h12,h22), find all real h != 0 for which the tangential
total moment map on Y^1(R) has rank < 2.  Mt(h) = h11 M1 + h12 M2 + h22 M3 (linear in h).  Rank<2 iff
exists lam=(l1,l2)!=0 with lam^T Mt(h) = 0: bilinear; we eliminate h via 3x3 minors in tau=l2/l1."""
import sys
sys.path.insert(0,'../td_k3')
from reduced import Reduced
from fractions import Fraction as Fr
import sympy as sp
tau=sp.symbols('tau')
HB=[((1,0),(0,0)),((0,1),(1,0)),((0,0),(0,1))]
def hq(H,e): return sum(H[a][b]*e[a]*e[b] for a in range(2) for b in range(2))
def weights(link,z,gam,H):
    m=len(link); W={}
    for i in gam:
        a,b=link[i],link[(i+1)%m]
        ea=(a[0]-z[0],a[1]-z[1]); eb=(b[0]-z[0],b[1]-z[1]); eab=(b[0]-a[0],b[1]-a[1])
        W[i]=(Fr(hq(H,ea)),Fr(hq(H,eb)),Fr(hq(H,eab)))
    return W
def moment_mats(link,q,z,gam):
    Ms=[]
    for H in HB:
        R=Reduced(link,q,z,gam,weights=weights(link,z,gam,H))
        _,_,Mt=R.analyse([(1,0,0),(0,1,0)])
        Ms.append(Mt)
    return Ms
def degenerate_h(Ms):
    """return list of (lam, h-basis) real degeneracies"""
    rows=[Ms[k][c,:] for c in range(2) for k in range(3)]
    S=sp.Matrix.vstack(*rows)
    # coordinates in row space
    rs=S.T.columnspace()
    if not rs: return 'allzero'
    Bas=sp.Matrix.hstack(*rs)          # n x r
    P=(Bas.T*Bas).inv()*Bas.T
    co=[[P*Ms[k][c,:].T for k in range(3)] for c in range(2)]
    out=[]
    # lam=(1,tau)
    Bt=sp.Matrix.hstack(*[co[0][k]+tau*co[1][k] for k in range(3)])
    r=Bt.shape[0]
    import itertools
    g=None
    for rows3 in itertools.combinations(range(r),3):
        d=sp.expand(Bt.extract(list(rows3),[0,1,2]).det())
        if d!=0: g=d if g is None else sp.gcd(g,d)
    if g is None: return 'identically degenerate'
    roots=[rt for rt in sp.Poly(g,tau).real_roots()] if sp.Poly(g,tau).degree()>0 else []
    for rt in set(roots):
        Bn=Bt.subs(tau,rt); ns=Bn.nullspace(); out.append(((1,rt),[list(v) for v in ns]))
    B2=sp.Matrix.hstack(*[co[1][k] for k in range(3)])
    ns=B2.nullspace()
    if ns: out.append(((0,1),[list(v) for v in ns]))
    return out
if __name__=="__main__":
    cases=[
     ('flat square, all 4',[(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)],(Fr(1,5),Fr(1,7),1),[0,1,2,3]),
     ('flat square, 1',[(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)],(Fr(1,5),Fr(1,7),1),[0]),
     ('flat square, 2 adj',[(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)],(Fr(1,5),Fr(1,7),1),[0,1]),
     ('flat square, 2 opp',[(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)],(Fr(1,5),Fr(1,7),1),[0,2]),
     ('flat pentagon all',[(1,0,0),(Fr(3,10),1,0),(-Fr(4,5),Fr(3,5),0),(-Fr(4,5),-Fr(3,5),0),(Fr(3,10),-1,0)],(Fr(1,7),Fr(1,9),1),[0,1,2,3,4]),
     ('flat hexagon all',[(2,0,0),(1,2,0),(-1,2,0),(-2,0,0),(-1,-2,0),(1,-2,0)],(Fr(1,3),Fr(-1,4),1),[0,1,2,3,4,5]),
     ('flat m=3 all',[(1,0,0),(-1,1,0),(-Fr(1,2),-1,0)],(Fr(1,9),0,1),[0,1,2]),
    ]
    for nm,link,q,gam in cases:
        Ms=moment_mats(link,q,(0,0,0),gam)
        print(nm,': degenerate (lam, h in basis h11,h12,h22):',degenerate_h(Ms),flush=True)
