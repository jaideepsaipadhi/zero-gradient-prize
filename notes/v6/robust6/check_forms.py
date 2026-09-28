"""Independent exact check, at rational shapes, of (i) the corrected closed forms (closed_forms.pkl) for the T1 Gram
pieces, (ii) the T2 part E2 = gC_f gC_g (x1^2+1) F(x2,y2) with robust5's F, (iii) the pairing matrix P, computed
here by an L^2 projection on the PHYSICAL T1 (not the reference), and (iv) the Frobenius-distance Gram closed form.
Fields are converted to monomials and integrated directly (no Bernstein integration formulas)."""
import sys; sys.path.insert(0, '../robust5')
import sympy as sp, pickle, random
from math import factorial as f
from fractions import Fraction as Fr
import pt_fields as PF, uref_fix, cert_PT as C
PF.UREF = uref_fix.get()
CF = pickle.load(open('closed_forms.pkl', 'rb'))
X, Y, u, w = sp.symbols('X Y u w'); x1s, y1s = sp.symbols('x1 y1')
class NumSym:
    exact = True
    @staticmethod
    def q(p, r=1):
        if isinstance(p, Fr): return sp.Rational(p.numerator, p.denominator) / r
        return sp.Rational(p, r) if isinstance(p, int) else p / r
def run(x1, y1, x2, y2):
    G, V, info = PF.star_fields(x1, y1, x2, y2, NumSym)
    s1, s2 = x1 + y1, x2 + y2
    z = sp.Matrix([0, 0]); a = sp.Matrix([1, 0]); c = sp.Matrix([x1, 1]) / s1
    b = (x2 * c + sp.Matrix([-c[1], c[0]])) / s2
    def bary(P0, P1, P2):
        M = sp.Matrix([[P0[0], P1[0], P2[0]], [P0[1], P1[1], P2[1]], [1, 1, 1]]); return list(M.inv() * sp.Matrix([X, Y, 1]))
    L1, L2 = bary(z, a, c), bary(z, c, b)
    def poly(coef, L): return sp.expand(sum(coef[al] * sp.Rational(f(4), f(al[0]) * f(al[1]) * f(al[2])) * L[0]**al[0] * L[1]**al[1] * L[2]**al[2] for al in PF.I4 if coef[al] != 0))
    def intT(e, P0, P1, P2):
        J = sp.Matrix.hstack(P1 - P0, P2 - P0)
        e = sp.expand(e.subs({X: P0[0] + J[0, 0] * u + J[0, 1] * w, Y: P0[1] + J[1, 0] * u + J[1, 1] * w}, simultaneous=True))
        return sp.integrate(sp.integrate(e, (u, 0, 1 - w)), (w, 0, 1)) * abs(J.det())
    vs = [([poly(F[1][cc], L1) for cc in range(2)], [poly(F[2][cc], L2) for cc in range(2)]) for F in V]
    gradv = lambda p: [sp.diff(p, X), sp.diff(p, Y)]
    ok = True
    for i in range(2):
        for j in range(i, 2):
            g1 = intT(sum(dp * dq for cc in range(2) for dp, dq in zip(gradv(vs[i][0][cc]), gradv(vs[j][0][cc]))), z, a, c)
            g2 = intT(sum(dp * dq for cc in range(2) for dp, dq in zip(gradv(vs[i][1][cc]), gradv(vs[j][1][cc]))), z, c, b)
            ge = sp.integrate(sum(vs[i][0][cc] * vs[j][0][cc] for cc in range(2)).subs(Y, 0), (X, 0, 1))
            gd = sp.integrate(sum(sp.diff(vs[i][0][cc], Y) * sp.diff(vs[j][0][cc], Y) for cc in range(2)).subs(Y, 0), (X, 0, 1))
            cf = [q.subs({sy: (x1 if sy.name == 'x1' else y1) for sy in q.free_symbols}) for q in CF["Q1"][(i, j)]]
            t2 = CF["gC"][i] * CF["gC"][j] * (x1**2 + 1) * sp.Rational(str(C.Ffun(Fr(str(x2)), Fr(str(y2)))))
            ok &= (sp.simplify(g1 - cf[0]) == 0) and (ge == cf[1]) and (gd == cf[2]) and (sp.simplify(g2 - t2) == 0)
    # pairing on the physical T1 (P_3 projection by direct Gram solve)
    mons = [X**p * Y**q for p in range(4) for q in range(4 - p)]
    Gm = sp.Matrix(len(mons), len(mons), lambda p, q: intT(mons[p] * mons[q], z, a, c))
    lz = -X - y1 * Y; la = X - x1 * Y
    P = sp.zeros(2, 5)
    for jj in range(5):
        H = sp.expand(lz**jj * la**(4 - jj))
        cc = Gm.LUsolve(sp.Matrix([intT(H * m, z, a, c) for m in mons]))
        eta = sp.expand((H - sum(cc[q] * mons[q] for q in range(len(mons)))).subs(Y, 0))
        for fi in range(2):
            P[fi, jj] = sp.integrate(eta * (-vs[fi][0][1].subs(Y, 0)), (X, 0, 1))
    Pref = sp.Matrix([[0, 0, sp.Rational(-1, 3780), sp.Rational(1, 2520), 0], [0, sp.Rational(-1, 2520), 0, sp.Rational(1, 2520), 0]])
    okP = (P - Pref).is_zero_matrix
    okG = all(p == q for p, q in zip(C.GWperp(Fr(str(x1)), Fr(str(y1))), C.GWperp_cf(Fr(str(x1)), Fr(str(y1)))))
    return ok, okP, okG
random.seed(11)
for _ in range(4):
    sh = [sp.Rational(random.randint(1, 20), 10), sp.Rational(random.randint(-3, 15), 10), sp.Rational(random.randint(1, 20), 10), sp.Rational(random.randint(-3, 15), 10)]
    print("shape", sh, " grams/closed forms/T2-part:", *run(*sh), flush=True)
