"""Closed forms (sympy) for the (P^T) certificate, reducing it to the shape (x1, y1) of T1 plus one number F_max:
   Q = Q1(x1,y1) + 13 * E2 * gC gC^T,   E2 = ||grad w_C||^2_{T2} = (x1^2+1) * F(x2,y2),
   F(x2,y2) = ||D^2 (30 mu_z^2 mu_c^2 mu_b)||^2_{T2n} / s2^2  (T2n: T2 normalised to base [0,1]),
and all remaining pieces (T1 grams, gamma's, G_W) as rational functions of (x1, y1)."""
import sympy as sp, pickle
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
# gamma_C of each field: coefficient of w_C; recover from T2 part: v_f|T2 = gC * wC|T2
_, U, wA, wC = PF.fields(x1, y1, x2, y2, NumSym)
gC = []
for F in V:
    al = (2, 1, 1)   # any T2 coefficient where wC is nonzero
    for comp in range(2):
        for a in PF.I4:
            if sp.simplify(wC[2][comp][a]) != 0:
                gC.append(sp.factor(sp.simplify(F[2][comp][a] / wC[2][comp][a]))); break
        else:
            continue
        break
print("gamma_C:", gC)
T1only = [{1: F[1], 2: [{a: 0 for a in PF.I4} for _ in range(2)]} for F in V]
Q1 = {}
for i in range(2):
    for j in range(i, 2):
        gg = sp.factor(sp.simplify(PF.gram_grad(T1only[i], T1only[j], G)))
        ge = sp.factor(sp.simplify(PF.gram_e(V[i], V[j])))
        gd = sp.factor(sp.simplify(PF.gram_dn(V[i], V[j], G)))
        Q1[(i, j)] = (gg, ge, gd)
        print((i, j), "grad_T1:", gg, "\n    e:", ge, "\n    dn:", gd)
# E2 check: ||grad wC||^2_{T2}
wC2 = {1: [{a: 0 for a in PF.I4} for _ in range(2)], 2: wC[2]}
E2 = sp.factor(sp.simplify(PF.gram_grad(wC2, wC2, G)))
print("E2 =", E2)
F = sp.factor(sp.simplify(E2 / (x1**2 + 1)))
print("F(x2,y2) = E2/(x1^2+1) =", F)
pickle.dump(dict(gC=gC, Q1=Q1, E2=E2, F=F), open("sym_Q.pkl", "wb"))
