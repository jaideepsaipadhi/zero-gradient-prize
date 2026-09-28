"""Independent SYMBOLIC verification (all shapes) of the two-triangle (P^T) fields v_1, v_2 on Q = T1 u T2.
The BB coefficients come from robust5/pt_fields.py (the object under test), with the reference fields either
robust5's UREF ('old') or the corrected ones from uref_fix.py ('new').  Everything else is independent:
vertices and barycentrics are recomputed here, fields are converted to monomials, and every constraint is checked
by direct sympy integration (no Bernstein integration formulas).
usage: python indep_Q.py old|new [x1 y1 x2 y2]   (symbolic if no shape given)"""
import sys; sys.path.insert(0, '../robust5')
import sympy as sp
from math import factorial as f
from fractions import Fraction as Fr
import pt_fields as PF
mode = sys.argv[1] if len(sys.argv) > 1 else "new"
if mode == "new":
    import uref_fix
    PF.UREF = uref_fix.get()
if len(sys.argv) > 2:   # rational shape given: python indep_Q.py new 7/10 1/2 3/5 4/5
    x1, y1, x2, y2 = [sp.Rational(q) for q in sys.argv[2:6]]
else:
    x1, y1, x2, y2 = sp.symbols('x1 y1 x2 y2', positive=True)
X, Y, t, u, w = sp.symbols('X Y t u w')
class NumSym:
    exact = True
    @staticmethod
    def q(p, r=1):
        if isinstance(p, Fr): return sp.Rational(p.numerator, p.denominator) / r
        return sp.Rational(p, r) if isinstance(p, int) else p / r
G, V, info = PF.star_fields(x1, y1, x2, y2, NumSym)
s1, s2 = x1 + y1, x2 + y2
z = sp.Matrix([0, 0]); a = sp.Matrix([1, 0]); c = sp.Matrix([x1, 1]) / s1
b = (x2 * c + sp.Matrix([-c[1], c[0]])) / s2
def bary(P0, P1, P2):
    M = sp.Matrix([[P0[0], P1[0], P2[0]], [P0[1], P1[1], P2[1]], [1, 1, 1]])
    lam = M.inv() * sp.Matrix([X, Y, 1])
    return [sp.cancel(l) for l in lam]
L1 = bary(z, a, c); L2 = bary(z, c, b)
def poly(coef, L):
    return sp.expand(sum(coef[al] * sp.Rational(f(4), f(al[0]) * f(al[1]) * f(al[2])) * L[0]**al[0] * L[1]**al[1] * L[2]**al[2]
                         for al in PF.I4 if coef[al] != 0))
def Z(e): return sp.simplify(sp.cancel(sp.together(e))) == 0
def zero_poly(e, var):
    e = sp.together(sp.expand(e)); num = sp.numer(e)
    return all(Z(cf) for cf in sp.Poly(sp.expand(num), var).coeffs()) if num != 0 else True
def on_seg(p, P, Q): return [sp.expand(q.subs({X: P[0] + t * (Q[0] - P[0]), Y: P[1] + t * (Q[1] - P[1])})) for q in p]
def int_T1(expr):   # T1 = (0,0),(1,0),c ; X = u + c0 w, Y = c1 w, det = c1
    e = sp.expand(expr.subs({X: u + c[0] * w, Y: c[1] * w}, simultaneous=True))
    return sp.integrate(sp.integrate(e, (u, 0, 1 - w)), (w, 0, 1)) * c[1]
Bs = [sp.binomial(4, j) * X**j * (1 - X)**(4 - j) for j in range(5)]
target = [Bs[1] - Bs[2], Bs[1] - Bs[3]]
ok_all = True
for fi, F in enumerate(V):
    v1 = [poly(F[1][comp], L1) for comp in range(2)]
    v2 = [poly(F[2][comp], L2) for comp in range(2)]
    res = {}
    res["cont [z,c]"] = all(zero_poly(p - q, t) for p, q in zip(on_seg(v1, z, c), on_seg(v2, z, c)))
    res["v=0 on [a,c]"] = all(zero_poly(p, t) for p in on_seg(v1, a, c))
    res["v=0 on [c,b]"] = all(zero_poly(p, t) for p in on_seg(v2, c, b))
    res["v=0 on [z,b]"] = all(zero_poly(p, t) for p in on_seg(v2, z, b))
    div2 = sp.expand(sp.diff(v2[0], X) + sp.diff(v2[1], Y))
    res["div v = 0 on T2"] = zero_poly(div2, [X, Y]) if div2 != 0 else True
    div1 = sp.diff(v1[0], X) + sp.diff(v1[1], Y)
    vn_e = sp.expand(-v1[1].subs(Y, 0))           # n = (0,-1)
    pb = []
    for p in range(4):
        for q in range(4 - p):
            r = X**p * Y**q
            lhs = int_T1(div1 * r); rhs = sp.integrate(sp.expand(vn_e * r.subs(Y, 0)), (X, 0, 1))
            pb.append(Z(lhs - rhs))
    res["pressure-blind on T1 (10 moments)"] = all(pb)
    res["pb detail"] = pb
    res["(M) int_e v = 0"] = all(Z(sp.integrate(vc.subs(Y, 0), (X, 0, 1))) for vc in v1)
    res["(D) int_e d_n v = 0"] = all(Z(sp.integrate((-sp.diff(vc, Y)).subs(Y, 0), (X, 0, 1))) for vc in v1)
    res["trace"] = zero_poly(vn_e - sp.expand(target[fi]), X)
    print(f"[{mode}] field v_{fi+1}:", res, flush=True)
    ok_all &= all(v for k, v in res.items() if k != "pb detail")
print(f"[{mode}] ALL CONSTRAINTS HOLD (" + ("symbolic, all shapes" if len(sys.argv) <= 2 else "at shape " + " ".join(sys.argv[2:6])) + "):", ok_all)
print("mA, mC, dA, dC =", [sp.factor(info[k]) for k in ("mA", "mC", "dA", "dC")])
