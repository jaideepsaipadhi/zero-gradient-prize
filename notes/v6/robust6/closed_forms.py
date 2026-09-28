"""Closed forms (sympy, all shapes) for the CORRECTED (P^T) fields: T1 Gram pieces as rational functions of
(x1, y1), gamma_C, and E2 = ||grad w_C||^2_{T2} = (x1^2+1) F(x2, y2).  Uses robust5/pt_fields.py (whose Gram
routines are correct; only UREF was wrong) with UREF replaced.  Saves closed_forms.pkl."""
import sys; sys.path.insert(0, '../robust5')
import sympy as sp, pickle
from fractions import Fraction as Fr
import pt_fields as PF, uref_fix
PF.UREF = uref_fix.get()
x1, y1 = sp.symbols('x1 y1', positive=True)
x2, y2 = sp.Integer(1), sp.Integer(1)   # T1 pieces and gamma_C do not depend on T2 (structural); F is checked separately
class NumSym:
    exact = True
    @staticmethod
    def q(p, r=1):
        if isinstance(p, Fr): return sp.Rational(p.numerator, p.denominator) / r
        return sp.Rational(p, r) if isinstance(p, int) else p / r
G, V, info = PF.star_fields(x1, y1, x2, y2, NumSym)
_, U, wA, wC = PF.fields(x1, y1, x2, y2, NumSym)
gC = []
for F in V:
    for comp in range(2):
        nz = [al for al in PF.I4 if sp.simplify(wC[2][comp][al]) != 0]
        if nz:
            gC.append(sp.factor(sp.simplify(F[2][comp][nz[0]] / wC[2][comp][nz[0]]))); break
print("gamma_C =", gC)
T1only = [{1: F[1], 2: [{a: 0 for a in PF.I4} for _ in range(2)]} for F in V]
Q1 = {}
for i in range(2):
    for j in range(i, 2):
        gg = sp.factor(sp.cancel(PF.gram_grad(T1only[i], T1only[j], G)))
        ge = sp.factor(sp.cancel(PF.gram_e(V[i], V[j])))
        gd = sp.factor(sp.cancel(PF.gram_dn(V[i], V[j], G)))
        Q1[(i, j)] = (gg, ge, gd)
        print((i, j), "\n  grad_T1:", gg, "\n  e:", ge, "\n  dn:", gd, flush=True)
pickle.dump(dict(gC=gC, Q1=Q1), open("closed_forms.pkl", "wb"))
