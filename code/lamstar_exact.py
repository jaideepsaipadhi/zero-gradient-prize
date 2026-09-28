"""Exact (rational-arithmetic) certificate that N_h is indefinite on divergence-free P4 fields below a penalty threshold.

Star: fan of m triangles at O=(0,0) with rational rim points; free boundary = the two segments of y=0 at O
(Nitsche boundary, outward normal (0,-1)); zero Dirichlet on the rim chords.
v is represented per triangle by P4 polynomials in global (x,y) with unknown rational coefficients.
Exact constraints: continuity across interior edges, zero on rim chords, div v = 0 on every triangle.
Exact quadratic forms: A = sum_T int (1/2) Dv:Dv,  B = int_{free} (d_n v).v,  C = int_{free} |v|^2.
Certificate: a rational v in the exact kernel with  2B - A - lam0*C > 0  (checked exactly).
By dilation invariance this certifies: N_h(v,v) < 0 whenever mu/rho_z <= lam0 for any mesh containing this star
(scaled) at a vertex of a straight piece of boundary.
"""
import sys, itertools
import numpy as np
import sympy as sy
from fractions import Fraction as Fr

x, y, t = sy.symbols("x y t")
MON = [(a, b) for a in range(5) for b in range(5 - a)]            # 15 monomials, degree <= 4


def tri_moment(T, a, b):
    """exact int_T x^a y^b over triangle with rational vertices"""
    P0, P1, P2 = [sy.Matrix(p) for p in T]
    u, w = sy.symbols("u w")
    X = P0 + u * (P1 - P0) + w * (P2 - P0)
    J = abs((P1 - P0).row_join(P2 - P0).det())
    integrand = sy.expand(X[0] ** a * X[1] ** b) * J
    return sy.integrate(sy.integrate(integrand, (w, 0, 1 - u)), (u, 0, 1))


def run(rim, lam0):
    O = (sy.Integer(0), sy.Integer(0))
    rim = [(sy.Rational(p[0]), sy.Rational(p[1])) for p in rim]
    tris = [(O, rim[i], rim[i + 1]) for i in range(len(rim) - 1)]
    m = len(tris)
    nvar = m * 2 * 15
    cs = sy.symbols(f"c0:{nvar}")

    def poly(ti, comp):
        base = (ti * 2 + comp) * 15
        return sum(cs[base + k] * x ** a * y ** b for k, (a, b) in enumerate(MON))

    V = [[poly(i, c) for c in range(2)] for i in range(m)]
    eqs = []

    def on_segment(expr, P, Q):
        sub = sy.expand(expr.subs({x: P[0] + t * (Q[0] - P[0]), y: P[1] + t * (Q[1] - P[1])}, simultaneous=True))
        return sy.Poly(sub, t).all_coeffs()

    for i in range(m - 1):                       # interior edge O - rim[i+1] shared by T_i and T_{i+1}
        for c in range(2):
            eqs += on_segment(V[i][c] - V[i + 1][c], O, rim[i + 1])
    for i in range(m):                           # rim chords: zero
        for c in range(2):
            eqs += on_segment(V[i][c], rim[i], rim[i + 1])
    for i in range(m):                           # pointwise div-free
        d = sy.expand(sy.diff(V[i][0], x) + sy.diff(V[i][1], y))
        eqs += sy.Poly(d, x, y).coeffs() if d != 0 else []
    Aeq = sy.Matrix([[sy.Poly(e, *cs).coeff_monomial(cv) for cv in cs] for e in eqs])
    ker = Aeq.nullspace()
    K = sy.Matrix.hstack(*ker)                   # nvar x k, rational
    print(f"star m={m}: {len(eqs)} exact constraints, divergence-free kernel dim = {K.shape[1]}", flush=True)

    # exact quadratic forms in the monomial coefficients
    def qform(expr_fn):
        pass

    # moments
    mom = []
    for T in tris:
        cache = {}
        for a in range(9):
            for b in range(9 - a):
                cache[(a, b)] = tri_moment(T, a, b)
        mom.append(cache)
    # A = sum_T int 2 u1x^2 + (u1y + u2x)^2 + 2 u2y^2
    Aform = sy.zeros(nvar, nvar)
    for i in range(m):
        u1, u2 = V[i]
        dens = sy.expand(2 * sy.diff(u1, x) ** 2 + (sy.diff(u1, y) + sy.diff(u2, x)) ** 2 + 2 * sy.diff(u2, y) ** 2)
        P = sy.Poly(dens, x, y)
        val = 0
        for (a, b), coef in zip(P.monoms(), P.coeffs()):
            val += coef * mom[i][(a, b)]
        H = sy.hessian(sy.expand(val), cs)
        Aform += H / 2
    # free boundary: y=0, segment from O to rim[0] (in T_0) and from rim[-1] to O (in T_{m-1}); n = (0,-1), d_n = -d_y
    Bval = 0; Cval = 0
    for (i, P, Q) in ((0, O, rim[0]), (m - 1, rim[-1], O)):
        u1, u2 = V[i]
        lo, hi = sorted([P[0], Q[0]])
        Bval += sy.integrate(sy.expand((-sy.diff(u1, y) * u1 - sy.diff(u2, y) * u2).subs(y, 0)), (x, lo, hi))
        Cval += sy.integrate(sy.expand((u1 ** 2 + u2 ** 2).subs(y, 0)), (x, lo, hi))
    Bform = sy.hessian(sy.expand(Bval), cs) / 2
    Cform = sy.hessian(sy.expand(Cval), cs) / 2
    Ak = K.T * Aform * K; Bk = K.T * Bform * K; Ck = K.T * Cform * K
    # floating-point maximizer of (2B - A)/C on the kernel, restricted to range(C)
    f = lambda M: np.array(M.evalf(30).tolist(), dtype=float)
    Af, Bf, Cf = f(Ak), f((Bk + Bk.T) / 2), f(Ck)
    Mf = 2 * Bf - Af
    ew, Q = np.linalg.eigh(Cf)
    keep = ew > 1e-12 * ew.max()
    R = Q[:, keep]; Nsp = Q[:, ~keep]
    Myy = R.T @ Mf @ R; Cyy = R.T @ Cf @ R
    if Nsp.shape[1]:
        Mkk = Nsp.T @ Mf @ Nsp; Mky = Nsp.T @ Mf @ R
        Seff = Myy - Mky.T @ np.linalg.solve(Mkk, Mky)
    else:
        Seff = Myy
    L = np.linalg.cholesky(Cyy); Li = np.linalg.inv(L)
    ev, EV = np.linalg.eigh(Li @ Seff @ Li.T)
    lam_float = ev[-1]
    yv = Li.T @ EV[:, -1]
    zf = R @ yv
    if Nsp.shape[1]:
        zf = zf + Nsp @ (-np.linalg.solve(Mkk, Mky @ yv))
    # rationalize and check EXACTLY
    zr = sy.Matrix([sy.Rational(Fr(val).limit_denominator(10 ** 6)) for val in zf / np.abs(zf).max()])
    num = (2 * zr.T * Bk * zr - zr.T * Ak * zr)[0]
    den = (zr.T * Ck * zr)[0]
    lam_exact = num / den
    ok = sy.simplify(num - lam0 * den) > 0
    print(f"  float lambda* = {lam_float:.6f};  exact Rayleigh quotient of rational witness = {float(lam_exact):.6f}")
    print(f"  EXACT check  2B - A - {lam0}*C > 0 :  {bool(ok)}   (value {float(num - lam0*den):.4e})")
    return lam_float, lam_exact, bool(ok)


if __name__ == "__main__":
    lam0 = sy.Rational(99, 10)
    run([(1, 0), (1, 1), (-1, 1), (-1, 0)], lam0)                 # m = 3
    run([(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0)], lam0)         # m = 4
