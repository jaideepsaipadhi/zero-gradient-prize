"""Correct reference space U(T^) = {u in lam_z lam_a P_{k-2}^2 : (u, grad r)_T^ = 0 for all r in P_{k-1}},
T^=(0,0),(1,0),(0,1), (z,a,c) = vertices 0,1,2, e = {y=0}.  Monomial basis, exact sympy; no BB integration code.
Prints dim, trace rank, whether traces lie in P_e (zero at ends, zero mean), and a basis with prescribed traces."""
import sympy as sp, sys
x, y = sp.symbols('x y')
def intT(p): return sp.integrate(sp.integrate(sp.expand(p), (y, 0, 1 - x)), (x, 0, 1))
def run(k, verbose=False):
    lz, la = 1 - x - y, x
    mons = [x**p * y**q for p in range(k - 1) for q in range(k - 1 - p)]       # P_{k-2}
    amb = [(lz * la * m, 0) for m in mons] + [(0, lz * la * m) for m in mons]
    rs = [x**p * y**q for p in range(k) for q in range(k - p)][1:]            # P_{k-1} mod constants
    A = sp.Matrix([[intT(u0 * sp.diff(r, x) + u1 * sp.diff(r, y)) for (u0, u1) in amb] for r in rs])
    N = A.nullspace()
    fields = [(sp.expand(sum(c * a[0] for c, a in zip(v, amb))), sp.expand(sum(c * a[1] for c, a in zip(v, amb)))) for v in N]
    traces = [sp.expand(-u1.subs(y, 0)) for (u0, u1) in fields]
    s = x
    tmat = sp.Matrix([[sp.Poly(t, x).coeff_monomial(x**j) for j in range(k + 1)] for t in traces])
    means = [sp.integrate(t, (x, 0, 1)) for t in traces]
    print(f"k={k}: rank A = {A.rank()} (of {len(rs)} rows, {len(amb)} unknowns), dim U = {len(N)}, trace rank = {tmat.rank()},"
          f" dim P_e = {k-2}, trace means = {set(means)}")
    return fields, traces
if __name__ == "__main__":
    for k in range(4, 8 if len(sys.argv) < 2 else int(sys.argv[1])):
        run(k)
