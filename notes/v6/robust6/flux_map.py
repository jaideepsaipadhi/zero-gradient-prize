"""Reference computation for the vertex-field construction (task 2).  For a trace g in P_4(e^) (e^ = [0,1]x{0},
z at s=0, a at s=1) find p in lam_z lam_a P_2^2 with (p, grad r)_T^ = -<g, r>_e^ for all r in P_3 (solvable iff int g = 0),
and report the flux  Lambda(g) := int_e^ p.n  (well defined: p is unique mod U(T^), whose traces have zero mean)."""
import sympy as sp, sys
k = int(sys.argv[1]) if len(sys.argv) > 1 else 4
x, y, s = sp.symbols('x y s')
def intT(p): return sp.integrate(sp.integrate(sp.expand(p), (y, 0, 1 - x)), (x, 0, 1))
lz, la = 1 - x - y, x
mons = [x**p * y**q for p in range(k - 1) for q in range(k - 1 - p)]
amb = [(lz * la * m, sp.Integer(0)) for m in mons] + [(sp.Integer(0), lz * la * m) for m in mons]
rs = [x**p * y**q for p in range(k) for q in range(k - p)]
A = sp.Matrix([[intT(u0 * sp.diff(r, x) + u1 * sp.diff(r, y)) for (u0, u1) in amb] for r in rs])
flux_row = sp.Matrix([[sp.integrate((-u1).subs(y, 0), (x, 0, 1)) for (u0, u1) in amb]])
basis = [s**j for j in range(k + 1)]
cs = sp.symbols('g0:%d' % (k + 1))
g = sum(c * b for c, b in zip(cs, basis))
g = g.subs(cs[0], -sum(cs[j] / (j + 1) for j in range(1, k + 1)))   # int g = 0
cs = cs[1:]
rhs = sp.Matrix([-sp.integrate(g * r.subs(y, 0).subs(x, s), (s, 0, 1)) for r in rs])
# general solution
sol, params = A[1:, :].gauss_jordan_solve(rhs[1:, :])
flux = sp.simplify((flux_row * sol)[0])
print("solvability: row r=1 of rhs =", rhs[0])
print("k =", k, " flux Lambda(g) =", sp.expand(flux.subs({p: 0 for p in params})), "  free params appear:", any(flux.has(p) for p in params))
