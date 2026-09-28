import sympy as sp
from math import factorial as F
from fractions import Fraction as Q
y,t=sp.symbols('y t')
for k in range(4,13):
    B=[(a,b) for a in range(k) for b in range(k-a)]
    m=lambda a,b: sp.Rational(F(a)*F(b),F(a+b+2))
    G=sp.Matrix(len(B),len(B),lambda i,j: m(B[i][0]+B[j][0],B[i][1]+B[j][1]))
    rhs=sp.Matrix([sp.Rational(1,a+1) if b==0 else 0 for a,b in B])
    c=G.LUsolve(rhs)
    g=sp.Poly(sp.expand(2*sp.diff(sp.legendre(k,t),t).subs(t,1-2*y)),y)
    target=sp.Matrix([g.coeff_monomial(y**b) if a==0 else 0 for a,b in B])
    ck=sum(c[i]/(B[i][0]+1) for i in range(len(B)) if B[i][1]==0)
    print(k, (c-target).is_zero_matrix, ck, ck==k*(k+1), flush=True)
