"""Exact (sympy, rational) local algebra for robust4.
Boundary triangle T = (0,0),(1,0),(a,r): e = [0,1] on y = 0, interior y > 0 (x tangential, y inward normal of Omega_h).
(1) z_T := Riesz rep in P_{k-1}(T) of r -> int_e r  has CONSTANT trace on e  (so the (D)-normal functional = const * flux).
(2) f_* in P_k[0,1]: orthogonal to the bubbles s(1-s)P_{k-2} and of zero mean (unique up to scale).
(3) K_T := {H in Hom_k : (I - pi_T)H |_e in span{1, f_*}}  -- the periodic-flat invisible set; dim, basis, H_*.
(4) Ker_T := span{n^k, l_a^k, l_b^k}; check K_T subset Ker_T.
(5) <f_*, L_k> = 0 for even k (L_k fields blind to K_T).
usage: python sym_kernel.py k a r
"""
import sys
import sympy as sp
x, y, s, u, w = sp.symbols('x y s u w')


def run(k, a, r):
    def intT(f):
        g = sp.expand(f.subs({x: u + w * a, y: w * r}, simultaneous=True) * r)
        return sp.integrate(sp.integrate(g, (u, 0, 1 - w)), (w, 0, 1))
    def on_e(f):
        return sp.expand(f.subs({y: 0}).subs(x, s))
    mons = [x**i * y**j for i in range(k) for j in range(k - i)]
    G = sp.Matrix(len(mons), len(mons), lambda i, j: intT(mons[i] * mons[j]))
    def proj(f):
        c = G.LUsolve(sp.Matrix([intT(f * m) for m in mons]))
        return sum(c[i] * mons[i] for i in range(len(mons)))
    c = G.LUsolve(sp.Matrix([sp.integrate(on_e(m), (s, 0, 1)) for m in mons]))
    zT = sum(c[i] * mons[i] for i in range(len(mons)))
    rho = sp.simplify(on_e(zT))
    # f_*
    cs = sp.symbols('c0:%d' % (k + 1)); p = sum(cs[i] * s**i for i in range(k + 1))
    eqs = [sp.integrate(p * s * (1 - s) * s**j, (s, 0, 1)) for j in range(k - 1)] + [sp.integrate(p, (s, 0, 1))]
    Mf = sp.Matrix([[sp.diff(q, c_) for c_ in cs] for q in eqs])
    nsf = Mf.nullspace(); assert len(nsf) == 1
    fstar = sp.expand(sum(nsf[0][i] * s**i for i in range(k + 1)))
    fstar = sp.expand(fstar / sp.Poly(fstar, s).LC())
    etas = [on_e(x**j * y**(k - j) - proj(x**j * y**(k - j))) for j in range(k + 1)]
    hs = sp.symbols('h0:%d' % (k + 1)); al = sp.symbols('al0:2')
    expr = sp.expand(sum(hs[j] * etas[j] for j in range(k + 1)) - al[0] - al[1] * fstar)
    eqs = sp.Poly(expr, s).all_coeffs()
    M = sp.Matrix([[sp.diff(q, v) for v in list(hs) + list(al)] for q in eqs])
    ns = M.nullspace()
    KT = [sp.expand(sum(v[j] * x**j * y**(k - j) for j in range(k + 1))) for v in ns]
    # Ker_T : normals of the three edges: e (y), (1,0)-(a,r): normal (r, a-1)... use line forms l = nx*x+ny*y
    la = r * x + (1 - a) * y          # normal to edge (1,0)-(a,r)
    lb = r * x - a * y                # normal to edge (0,0)-(a,r)
    ker = [y**k, sp.expand(la**k), sp.expand(lb**k)]
    def coeffs(P):
        P = sp.Poly(P, x, y); return [P.coeff_monomial(x**j * y**(k - j)) for j in range(k + 1)]
    Kmat = sp.Matrix([coeffs(P) for P in ker])
    inKer = [ (sp.Matrix.vstack(Kmat, sp.Matrix([coeffs(P)]))).rank() == Kmat.rank() for P in KT]
    Lk = sp.legendre(k, 2 * s - 1)
    return dict(rho=rho, fstar=fstar, dimK=len(ns), KT=KT, K_in_KerT=inKer,
                fstar_Lk=sp.integrate(fstar * Lk, (s, 0, 1)),
                fstar_antisym=sp.expand(fstar + fstar.subs(s, 1 - s)) == 0)


if __name__ == "__main__":
    k = int(sys.argv[1]); a = sp.Rational(sys.argv[2]); r = sp.Rational(sys.argv[3])
    for kk, vv in run(k, a, r).items():
        print(kk, ':', vv)
