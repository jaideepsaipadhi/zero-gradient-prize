"""Symbolic (sympy) verification that the explicit fields satisfy every constraint of fan.py-type for ALL shapes:
run the generic pipeline with x1,y1,x2,y2 as symbols and simplify the residuals to 0.
Also: reference Riesz representer zeta of r -> int_e r (P_3 on T^) has constant trace on e."""
import sympy as sp
import pt_fields as PF
from fractions import Fraction as Fr

x1, y1, x2, y2 = sp.symbols('x1 y1 x2 y2', real=True)

class NumSym:
    exact = True
    @staticmethod
    def q(p, r=1):
        if isinstance(p, Fr):
            return sp.Rational(p.numerator, p.denominator) / r
        return sp.Rational(p, r) if isinstance(p, int) else p / r

G, V, info = PF.star_fields(x1, y1, x2, y2, NumSym)
print("mA =", sp.factor(info["mA"]), " mC =", sp.simplify(info["mC"]), " dA =", sp.factor(info["dA"]), " dC =", sp.factor(info["dC"]))

# independent symbolic constraint check: continuity/BC by construction of BB supports; here check
#  pressure-blindness on T1 and T2, (M), (D) vector, using the BB formulas directly.
from math import factorial as fact
def div_bb(F, t):
    g = G["g%d" % t]
    dx0, dy0 = PF.grad_bb(F[t][0], g, 4); dx1, dy1 = PF.grad_bb(F[t][1], g, 4)
    return {b: dx0[b] + dy1[b] for b in PF.I3}
res = []
for F in V:
    # continuity across [z,c]: T1 (i,0,j) == T2 (i,j,0); BCs: T1 alpha_z = 0; T2 alpha_z = 0 or alpha_b... (z,c,b): [z,b] is lam_c = 0
    for comp in range(2):
        for i in range(5):
            res.append(F[1][comp][(i, 0, 4 - i)] - F[2][comp][(i, 4 - i, 0)])
        for al in PF.I4:
            if al[0] == 0: res += [F[1][comp][al], F[2][comp][al]]
            if al[1] == 0: res.append(F[2][comp][al])
    for t in (1, 2):
        dv = div_bb(F, t); area2 = G["det%d" % t]
        for r in PF.I3:
            val = sum(dv[b] * area2 * NumSym.q(PF.int_T(1, b, r, 3, 3)) for b in PF.I3)
            if t == 1:   # - <v.n, r>_e, n = (0,-1)
                val -= sum(-F[1][1][al] * NumSym.q(PF.int_e(al, r, 4, 3)) for al in PF.I4 if al[2] == 0)
            res.append(val)
    for comp in range(2):
        res.append(sum(F[1][comp][al] for al in PF.I4 if al[2] == 0))
        res.append(PF.Dt(F, G, comp))
bad = [r for r in res if sp.simplify(sp.together(r)) != 0]
print("residuals checked:", len(res), " nonzero:", len(bad))

# zeta on reference
x, y = sp.symbols('x y')
mons = [x**i * y**j for i in range(4) for j in range(4 - i)]
intT = lambda f: sp.integrate(sp.integrate(f, (y, 0, 1 - x)), (x, 0, 1))
Gm = sp.Matrix(len(mons), len(mons), lambda i, j: intT(mons[i] * mons[j]))
rhs = sp.Matrix([sp.integrate(m.subs(y, 0), (x, 0, 1)) for m in mons])
cz = Gm.LUsolve(rhs); zeta = sum(cz[i] * mons[i] for i in range(len(mons)))
print("trace of zeta on e:", sp.expand(zeta.subs(y, 0)))
