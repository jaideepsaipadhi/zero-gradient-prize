"""Exact checks for (P^T) with k >= 5 (robust5 REPORT, Sec. 2.5):
 (a) the reference space U(T^) = {u in lam_z lam_a P_{k-2}^2 : (u, grad r)_{T^} = 0, r in P_{k-1}} has normal traces
     spanning exactly P_e = {p in P_k: p(0) = p(1) = 0, int p = 0}  (dim k-2);
 (b) the Riesz representer zeta^ in P_{k-1}(T^) of r -> int_e r has constant trace on e (so the normal parts of (M), (D)
     are automatic for fields vanishing at both ends of e with zero-mean trace)."""
import sympy as sp, flint
from fractions import Fraction as Fr
import fan
from fan import U_space, trmat
x, y = sp.symbols('x y')
for k in range(4, 9):
    F, B = U_space([(0, 0), (1, 0), (0, 1)], k)
    T = trmat(F, B)                        # Bernstein coefficients of traces (degree k)
    r = T.rank()
    # zero mean and zero endpoints for every trace
    ok = True
    for i in range(T.nrows()):
        row = [T[i, j] for j in range(k + 1)]
        ok &= (row[0] == 0 and row[k] == 0 and sum(row) == 0)
    mons = [x**i * y**j for i in range(k) for j in range(k - i)]
    intT = lambda f: sp.integrate(sp.integrate(sp.expand(f), (y, 0, 1 - x)), (x, 0, 1))
    Gm = sp.Matrix(len(mons), len(mons), lambda i, j: intT(mons[i] * mons[j]))
    cz = Gm.LUsolve(sp.Matrix([sp.integrate(m.subs(y, 0), (x, 0, 1)) for m in mons]))
    ztr = sp.expand(sum(cz[i] * mons[i] for i in range(len(mons))).subs(y, 0))
    print("k=%d: dim U(T^)=%d, trace rank=%d (dim P_e=%d), all traces in P_e: %s, zeta^ trace on e: %s"
          % (k, len(B), r, k - 2, ok, ztr), flush=True)
