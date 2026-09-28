"""Polynomial reproduction (patch) test: u = curl psi in P_k^3, p in P_{k-1}, on the box with (i) consistent
Nitsche + pressure traction on all faces, (ii) strong Dirichlet.  The discrete solution must equal the exact one
up to round-off.  python3 patch_test.py"""
import os, sys, numpy as np, sympy as sy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svn3d, exact3d as e

x, y, z = e.X
for k in (3, 4):
    psi = [y**2 * z * x + z**3, x**3 * y - z * x, y**2 * z**2] if k == 3 else [y**3 * z * x + z**4, x**4 - z * x * y**2, y**2 * z**3]
    u = e._curl(psi); p = x * y + z**2 - x if k == 3 else x * y * z + z**3 - x**2
    f = [-sum(sy.diff(u[i], v, 2) for v in e.X) + sy.diff(p, e.X[i]) for i in range(3)]
    Ul = sy.lambdify(e.X, u); Gl = sy.lambdify(e.X, [[sy.diff(u[i], v) for v in e.X] for i in range(3)]); Fl = sy.lambdify(e.X, f)
    vec = lambda fn: (lambda a, b, c: np.stack([np.broadcast_to(np.asarray(t, float), np.shape(a)) for t in fn(a, b, c)], -1))

    def gm(a, b, c):
        v = Gl(a, b, c)
        return np.stack([np.stack([np.broadcast_to(np.asarray(v[i][j], float), np.shape(a)) for j in range(3)], -1) for i in range(3)], -2)
    ex = dict(u=vec(Ul), gu=gm, f=vec(Fl))
    for N in (1, 2):
        S = svn3d.build("box", N, k, L=1.0)
        mu = 400.0 if k == 3 else 700.0
        sysm = svn3d.assemble(S, mu, ex["f"], [(S.fBox, ex["u"], "sym")], [(S.fBox, True)])
        uh, ph, hist, rr = svn3d.solve_ipm(S, sysm, None)
        E1 = svn3d.errors(S, uh, ex, mu)
        sysm = svn3d.assemble(S, mu, ex["f"], [], [])
        uh, ph, hist2, rr2 = svn3d.solve_ipm(S, sysm, ex["u"])
        E2 = svn3d.errors(S, uh, ex, mu)
        print(f"k={k} N={N}  Nitsche: H1 err {E1['H1']:.2e} L2 {E1['L2']:.2e} (ipm {len(hist)}, |div| {hist[-1]:.1e}) | "
              f"strong: H1 err {E2['H1']:.2e} L2 {E2['L2']:.2e} (|div| {hist2[-1]:.1e})", flush=True)
