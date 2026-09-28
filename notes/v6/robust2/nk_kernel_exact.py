"""Exact rational check (fractions) of  ker Lam_T = span{(nu_w.x)^k}, w = a,b,c,  for rational shapes.
Uses (nu_w . x)^k  ~  lambda_w^k mod P_{k-1}; nu need not be normalised.  Prints the exact residual (should be 0)."""
import sys
from fractions import Fraction as Fr
from math import comb, factorial
import sympy as sp


def mono_int(a, b, xi, ze):
    # int_T x^a y^b, T = conv{(0,0),(1,0),(xi,ze)}; x = s + r xi, y = r ze, dA = ze ds dr
    tot = Fr(0)
    for i in range(a + 1):
        # (s + r xi)^a = sum C(a,i) s^(a-i) (r xi)^i
        p, q = a - i, i + b
        tot += comb(a, i) * xi**i * ze**b * Fr(factorial(p) * factorial(q), factorial(p + q + 2))
    return tot * ze


def check(k, xi, ze):
    Pk1 = [(a, b) for a in range(k) for b in range(k - a)]
    G = sp.Matrix(len(Pk1), len(Pk1), lambda i, j: mono_int(Pk1[i][0] + Pk1[j][0], Pk1[i][1] + Pk1[j][1], xi, ze))
    Gi = G.inv()
    V = [(Fr(0), Fr(0)), (Fr(1), Fr(0)), (xi, ze)]
    # P_e basis p_j = d/dt [t^2(1-t)^2 t^j] as coefficient lists in t
    t = sp.Symbol('t')
    Pe = [sp.Poly(sp.diff(t**2 * (1 - t)**2 * t**j, t), t) for j in range(k - 2)]
    out = []
    for w in range(3):
        P, Q = V[(w + 1) % 3], V[(w + 2) % 3]
        nu = (Q[1] - P[1], -(Q[0] - P[0]))
        c = [comb(k, j) * nu[0]**(k - j) * nu[1]**j for j in range(k + 1)]      # coefficient of x^(k-j) y^j
        rhs = sp.Matrix([sum(c[j] * mono_int(k - j + a, j + b, xi, ze) for j in range(k + 1)) for (a, b) in Pk1])
        cp = Gi * rhs
        # (I - pi) H on e (y = 0): H(t,0) = c[0] t^k ;  pi H(t,0) = sum_{b=0} cp * t^a
        f = c[0] * t**k - sum(cp[i] * t**Pk1[i][0] for i in range(len(Pk1)) if Pk1[i][1] == 0)
        out.append([sp.integrate(sp.expand(f * p.as_expr()), (t, 0, 1)) for p in Pe])
    return out


if __name__ == "__main__":
    k = int(sys.argv[1])
    for (xi, ze) in [(Fr(1, 2), Fr(1, 2)), (Fr(1, 5), Fr(3, 5)), (Fr(-2, 7), Fr(5, 9)), (Fr(13, 10), Fr(2, 3))]:
        r = check(k, xi, ze)
        print(f"k={k} shape ({xi},{ze}):", "all zero" if all(v == 0 for row in r for v in row) else r)
