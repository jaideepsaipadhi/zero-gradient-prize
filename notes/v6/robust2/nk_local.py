"""
Local analysis of the leading-order map  H (homogeneous degree-k) -> ((I - pi_T) H)|_e  tested on edge spaces.
Reference boundary triangle T = conv{(0,0),(1,0),(xi,zeta)}, e = [(0,0),(1,0)], x = tangential, y = normal (into T).
  P_e  = { p in P_k(e): p(0)=p(1)=0, int p = 0 }            (traces of the local div-free fields X_e)
  L_k  = shifted Legendre polynomial of degree k on e       (shape-independent pairing, Lemma L2 of REPORT.tex)
Exact rational arithmetic (sympy) for given shapes; float sweep over shapes for singular values.
usage: python nk_local.py k
"""
import sys, itertools
import sympy as sp
import numpy as np

t, x, y, xi, ze = sp.symbols('t x y xi zeta', real=True)


def tri_int(expr, xi_v, ze_v):
    # integral over T = conv{(0,0),(1,0),(xi,ze)} via affine map from reference (s,r): X = s*(1,0)+r*(xi,ze)
    s, r = sp.symbols('s r')
    X = s + r * xi_v; Y = r * ze_v
    f = sp.expand(expr.subs({x: X, y: Y}) * ze_v)
    return sp.integrate(sp.integrate(f, (s, 0, 1 - r)), (r, 0, 1))


def local_maps(k, xi_v, ze_v):
    Pk1 = [x**a * y**b for a in range(k) for b in range(k - a)]
    G = sp.Matrix(len(Pk1), len(Pk1), lambda i, j: tri_int(Pk1[i] * Pk1[j], xi_v, ze_v))
    Gi = G.inv()
    Hs = [x**(k - j) * y**j for j in range(k + 1)]
    res = []
    for H in Hs:
        rhs = sp.Matrix([tri_int(H * q, xi_v, ze_v) for q in Pk1])
        c = Gi * rhs
        piH = sum(c[i] * Pk1[i] for i in range(len(Pk1)))
        r = sp.expand((H - piH).subs({y: 0, x: t}))
        res.append(r)
    # P_e basis: d/dt [t^2 (1-t)^2 t^j], j = 0..k-3
    Pe = [sp.diff(t**2 * (1 - t)**2 * t**j, t) for j in range(k - 2)]
    Ge = sp.Matrix(len(Pe), len(Pe), lambda i, j: sp.integrate(Pe[i] * Pe[j], (t, 0, 1)))
    Lam = sp.Matrix(len(Pe), k + 1, lambda i, j: sp.integrate(Pe[i] * res[j], (t, 0, 1)))
    Lk = sp.legendre(k, 2 * t - 1)
    lk = [sp.integrate(Lk * r, (t, 0, 1)) for r in res]
    ints = [sp.integrate(r, (t, 0, 1)) for r in res]
    return res, Ge, Lam, lk, ints


def main():
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    print(f"k = {k}")
    shapes = [(sp.Rational(1, 2), sp.sqrt(3) / 2), (sp.Rational(1, 2), sp.Rational(1, 2)), (sp.Rational(1, 5), sp.Rational(3, 5)),
              (sp.Rational(9, 10), sp.Rational(2, 5)), (sp.Rational(-1, 5), sp.Rational(7, 10))]
    for (a, b) in shapes:
        res, Ge, Lam, lk, ints = local_maps(k, a, b)
        # normalised matrix  Ge^{-1/2} Lam  (float)
        Gf = np.array(Ge.evalf(30).tolist(), float); Lf = np.array(Lam.evalf(30).tolist(), float)
        w, V = np.linalg.eigh(Gf); Gm12 = V @ np.diag(w**-0.5) @ V.T
        M = Gm12 @ Lf
        sv = np.linalg.svd(M, compute_uv=False)
        _, _, Vt = np.linalg.svd(M)
        ker = Vt[len(sv):] if len(sv) < k + 1 else Vt[-1:]
        # visibility of the pure tangential monomial x^k: |M e_0|
        print(f" shape c=({a},{b}):  sv(P_e map) = {np.array2string(sv, precision=3)}")
        print(f"   |P_e-part of (I-pi)x^k| = {np.linalg.norm(M[:,0]):.3e};  kernel basis (coeffs of x^(k-j) y^j):")
        for kv in ker:
            print("     ", np.array2string(kv / np.abs(kv).max(), precision=3, suppress_small=True))
        print(f"   <(I-pi)H, L_k>_e  = {[sp.nsimplify(v) for v in lk]}   (exact; expected k!^2/(2k+1)! = {sp.factorial(k)**2/sp.factorial(2*k+1)} on x^k only)")
        print(f"   int_e (I-pi)H   = {[float(v) for v in ints]}")


if __name__ == "__main__":
    main()
