"""Constants of robust4 Lemma 2.3/2.4: J_k (monic, orthogonal to P_{k-1} w.r.t. weight 1-s on [0,1]),
d_k(s) := J_k(1-s) - J_k(s); verify <d_k, g> = beta_k (g(0) - g(1)) for all g in P_k[0,1], and d_k != 0."""
import sympy as sp
s = sp.symbols('s')
for k in range(2, 9):
    c = sp.symbols('c0:%d' % k)
    J = s**k + sum(c[i] * s**i for i in range(k))
    sol = sp.solve([sp.integrate(J * (1 - s) * s**j, (s, 0, 1)) for j in range(k)], c, dict=True)[0]
    J = sp.expand(J.subs(sol)); d = sp.expand(J.subs(s, 1 - s) - J)
    vals = [sp.integrate(d * s**j, (s, 0, 1)) for j in range(k + 1)]      # g = s^j: g(0)-g(1) = [j==0] - 1
    beta = vals[1] / (-1)                                                  # g = s: g(0)-g(1) = -1
    ok = all(sp.simplify(vals[j] - beta * ((1 if j == 0 else 0) - 1)) == 0 for j in range(k + 1))
    print(k, 'J_k(0) =', J.subs(s, 0), ' beta_k =', beta, ' identity holds:', ok, ' d_k nonzero:', d != 0)
