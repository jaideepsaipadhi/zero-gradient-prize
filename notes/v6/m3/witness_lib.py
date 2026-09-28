"""Library: witness Hessians + kappa by quadrature. (was: Floating-point cross-check: the closed-form witness reproduces the star constants kappa_z of
code/lower_bound_tests.py (logs/lower_bound_stars.log), which were computed as the optimum over
the numerically computed space Y^1(omega_z) (dim 1).  Independent of code/m3_certify.py.
Run: python3 crosscheck_meshes.py      (< 1 min)"""
import sys, math, numpy as np, sympy as sp
sys.path.insert(0, '/home/claude/zero-gradient-prize/code')
import svn
a1, b1, a2, b2, a3, b3, x, y = sp.symbols('a1 b1 a2 b2 a3 b3 x y')
import os
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'witness_def.py')).read())
hs = [[sp.lambdify((a1, b1, a2, b2, a3, b3, x, y), sp.diff(f, *v), 'numpy')
       for v in ((x, x), (x, y), (y, y))] for f in PHIS]
g, w = np.polynomial.legendre.leggauss(6); g = (g + 1) / 2; w = w / 2
def tri_int(F, V0, V1, V2):
    # Duffy: x = V0 + u (V1-V0) + u v (V2-V1), jac = u |det|
    U, W = np.meshgrid(g, g, indexing='ij'); WU, WW = np.meshgrid(w, w, indexing='ij')
    Pp = V0[None, None] + U[..., None] * (V1 - V0) + (U * W)[..., None] * (V2 - V1)
    det = abs((V1 - V0)[0] * (V2 - V1)[1] - (V1 - V0)[1] * (V2 - V1)[0])
    return (F(Pp[..., 0], Pp[..., 1]) * U * WU * WW).sum() * det
def kappa(p):   # p: 4x2 array, z=0, fan ccw, p0 on +x axis.  The symbolic witness is written with p0=(1,0):
    p = np.asarray(p, float) / np.linalg.norm(p[0])   # kappa is scale invariant; M and ||grad w|| below are for |p0|=1
    c = (p[1][0], p[1][1], p[2][0], p[2][1], p[3][0], p[3][1])
    nrm2 = 0
    for j in range(3):
        F = lambda X, Y, j=j: (hs[j][0](*c, X, Y)**2 + 2 * hs[j][1](*c, X, Y)**2 + hs[j][2](*c, X, Y)**2)
        nrm2 += tri_int(F, np.zeros(2), p[j], p[j + 1])
    cr = lambda u, v: u[0] * v[1] - u[1] * v[0]
    A = cr(p[1], p[3]) * cr(p[2], p[3]); A3 = cr(p[0], p[1]) * cr(p[0], p[2])
    M = (np.linalg.norm(p[0])**5 * A + np.linalg.norm(p[3])**5 * A3) / 60
    return M, math.sqrt(nrm2)
