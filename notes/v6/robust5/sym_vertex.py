"""Lemma K1 constants (exact): S_k(s) = J_k(1-s) + J_k(s) acts on P_k as  g -> alpha_k int g + alpha'_k (g(0) + g(1))."""
import sympy as sp
s = sp.symbols('s')
for k in range(2, 9):
    c = sp.symbols('c0:%d' % k)
    J = s**k + sum(c[i] * s**i for i in range(k))
    sol = sp.solve([sp.integrate(J * (1 - s) * s**j, (s, 0, 1)) for j in range(k)], c, dict=True)[0]
    J = sp.expand(J.subs(sol)); S = sp.expand(J.subs(s, 1 - s) + J)
    al, alp = sp.symbols('al alp')
    eqs = [sp.integrate(S * s**j, (s, 0, 1)) - (al * sp.Rational(1, j + 1) + alp * ((1 if j == 0 else 0) + 1)) for j in range(k + 1)]
    sl = sp.solve([eqs[0], eqs[2]], [al, alp])
    ok = all(sp.simplify(e.subs(sl)) == 0 for e in eqs)
    print("k=%d  alpha=%s  alpha'=%s  identity on all of P_k: %s" % (k, sl[al], sl[alp], ok))
