"""Annihilator of the joint edge-moment image for single-face rings, in face-intrinsic coordinates.
F=[z,a,b] ordered as the macro's local (1,2,3) = (z,a,b); E_ij for pairs (z,a),(z,b),(a,b);
each E_ij (tangential) written as alpha (a-z) + beta (b-z).  Exact via sympy on a float-free basis."""
from patchlib import *
from fractions import Fraction as Fr
import sympy as sp
def ring(a,b,cs,q=(0,0,1),z=(0,0,0)):
    link=[a,b]+cs; m=len(link)
    return Patch([[q,z,link[j],link[(j+1)%m]] for j in range(m)],[(0,(1,2,3))],3)
def image_basis(P):
    # exact nullspaces with sympy DomainMatrix
    from sympy.polys.matrices import DomainMatrix
    from sympy import QQ
    def dm(rows):
        return DomainMatrix([[QQ(r.get(c,0).numerator,r.get(c,0).denominator) if c in r else QQ(0) for c in range(P.nunk)] for r in rows],(len(rows),P.nunk),QQ)
    D=dm(P.div_rows); N=D.nullspace()      # rows = basis of Y
    C=dm(P.cons_rows()); K=(C*N.transpose()).nullspace()  # coefficient combos -> Y^1
    Y1=K*N
    E=dm(edge_moment_rows(P,0))
    return (E*Y1.transpose()).to_Matrix()   # 9 x dimY1
def intrinsic(P,Img):
    V=P.fgeom[0]  # face verts in macro order: (z,a,b) since macro = [q,z,a,b]
    z,a,b=[sp.Matrix([sp.Rational(c.numerator,c.denominator) for c in p]) for p in V]
    Bm=sp.Matrix.hstack(a-z,b-z)
    rows=[]
    for i in range(3):
        blk=Img[3*i:3*i+3,:]
        co=(Bm.T*Bm).inv()*Bm.T*blk   # 2 x n
        assert (Bm*co-blk).is_zero_matrix
        rows.append(co)
    M=sp.Matrix.vstack(*rows)   # 6 x n  (order: E_za(alpha,beta),E_zb,E_ab)
    ann=M.T.nullspace()
    return M.rank(),[list(v.T/ (v[[i for i in range(6) if v[i]!=0][0]])) for v in ann]
tests={
 'm=3':((1,0,0),(0,1,0),[(-1,-1,Fr(1,2))]),
 'm=3b':((1,Fr(1,5),0),(-Fr(1,3),1,0),[(-1,-1,Fr(1,3))]),
 'm=4':((1,0,0),(0,1,0),[(-1,Fr(1,3),Fr(1,2)),(Fr(1,4),-1,Fr(1,2))]),
 'm=4b':((2,0,0),(0,1,0),[(-1,Fr(1,3),Fr(1,2)),(Fr(1,4),-1,Fr(1,2))]),
}
for nm,(a,b,cs) in tests.items():
    for q in [(0,0,1),(Fr(1,3),Fr(-1,5),Fr(4,5))]:
        P=ring(a,b,cs,q=q)
        r,ann=intrinsic(P,image_basis(P))
        print(nm,'q=',q,'rank',r,'annihilator (E_za a,b | E_zb a,b | E_ab a,b):',ann)
