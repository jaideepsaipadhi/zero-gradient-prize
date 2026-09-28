# symbolic (all shapes) check of the triangular-system diagonals and of g_C, and of the Piola field traces
import sys; sys.path.insert(0,'../robust5'); sys.path.insert(0,'../robust6')
import sympy as sp, pt_fields as PF, uref_fix
from fractions import Fraction as Fr
PF.UREF=uref_fix.get()
x1,y1,x2,y2=sp.symbols('x1 y1 x2 y2')
class N:
    exact=True
    @staticmethod
    def q(p,r=1):
        if isinstance(p,Fr): return sp.Rational(p.numerator,p.denominator)/r
        return sp.Rational(p,r) if isinstance(p,int) else p/r
G,U,wA,wC=PF.fields(x1,y1,x2,y2,N)
mA,mC=PF.Mt(wA),PF.Mt(wC); dA,dC=PF.Dt(wA,G),PF.Dt(wC,G)
print('mA,mC,dA,dC =',[sp.factor(sp.simplify(q)) for q in (mA,mC,dA,dC)])
for u in U:
    gA=-PF.Mt(u)/mA; gC=-(PF.Dt(u,G)+gA*dA)/dC
    print('gA =',sp.factor(sp.simplify(gA)),' gC =',sp.factor(sp.simplify(gC)))
