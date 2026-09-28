"""(f): exact check of the edge-representer constants of Lemma pd:lem:rep.

Reference triangle T = {(0,0),(1,0),(0,1)}, edge e = [0,1]x{0}.  r_k in P_{k-1}(T) is the L2(T)
representer of r -> int_e r.  Claims (proved in REPORT.tex, Prop. F):
   r_k(x,y) = 2 P_k'(1-2y)          (P_k Legendre; depends on y only)
   c_k := int_e r_k = ||r_k||_T^2 = k(k+1),   d_k := ||r_k||_{L2(e)}^2 = k^2 (k+1)^2 .
Exact rational arithmetic (sympy), full P_{k-1}(T) basis, k = 1..K.
usage: python3 chat_symbolic.py [K=12]
"""
import sys
import sympy as sp

x, y = sp.symbols('x y')
K = int(sys.argv[1]) if len(sys.argv) > 1 else 12


def tri_int(f):
    return sp.integrate(sp.integrate(f, (x, 0, 1 - y)), (y, 0, 1))


for k in range(1, K + 1):
    mons = [x**a * y**b for a in range(k) for b in range(k - a)]
    n = len(mons)
    # Gram matrix: int_T x^a y^b = a! b! / (a+b+2)!
    def mom(a, b):
        return sp.factorial(a) * sp.factorial(b) / sp.factorial(a + b + 2)
    ex = [(a, b) for a in range(k) for b in range(k - a)]
    G = sp.Matrix(n, n, lambda i, j: mom(ex[i][0] + ex[j][0], ex[i][1] + ex[j][1]))
    L = sp.Matrix(n, 1, lambda i, _: sp.Rational(1, ex[i][0] + 1) if ex[i][1] == 0 else 0)
    c = G.LUsolve(L)
    r = sp.expand(sum(c[i] * mons[i] for i in range(n)))
    ck = sp.integrate(r.subs(y, 0), (x, 0, 1))
    dk = sp.integrate(r.subs(y, 0)**2, (x, 0, 1))
    t = sp.symbols('t')
    claim = sp.expand(2 * sp.diff(sp.legendre(k, t), t).subs(t, 1 - 2 * y))
    ok = sp.simplify(r - claim) == 0
    print(f"k={k:2d}  dim={n:3d}  c_k={ck}  k(k+1)={k*(k+1)}  d_k={dk}  c_k^2={ck**2}  "
          f"r_k==2P_k'(1-2y): {ok}  r_k(0)={r.subs(y,0)}", flush=True)
    assert ck == k * (k + 1) and dk == ck**2 and ok
print("all assertions passed")
