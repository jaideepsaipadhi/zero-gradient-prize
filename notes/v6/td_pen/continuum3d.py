"""3D half-space continuum constant for the witness v = curl(F e_y), exact (SymPy).
Half-space {z>0} = fluid, wall z=0, n=(0,0,-1) (points out of the fluid, into P_h), d_n v = -d_z v.
F = Y*phi(x/X)*phi(y/X)*g(z/Y), g(t)=t(1-t)^p, phi(s)=(1-s^2)^q  (F=0 for z>=Y, |x|>=X, |y|>=X).
v = curl(F e_y) = (-F_z, 0, F_x): divergence-free, v.n = 0 on z=0.
A = int 1/2 Dv:Dv (Dv = grad v + grad v^T), B = int_{z=0} (d_n v).v, C = int_{z=0} |v|^2.
Scaling: A,B ~ ell, C ~ ell^2, so Y*(2B-A)/C is scale-free (depends on X/Y, p, q only)."""
import sympy as sp, sys
x, y, z = sp.symbols("x y z", real=True)

def lam(M, p=6, q=6):
    Y = sp.Integer(1); X = sp.Rational(M) * Y
    g = lambda t: t * (1 - t) ** p
    phi = lambda s: (1 - s ** 2) ** q
    F = phi(x / X) * phi(y / X) * Y * g(z / Y)
    v = [-sp.diff(F, z), sp.Integer(0), sp.diff(F, x)]
    X_ = [x, y, z]
    G = [[sp.diff(v[i], X_[j]) for j in range(3)] for i in range(3)]
    dens = sum((G[i][j] + G[j][i]) ** 2 for i in range(3) for j in range(3)) / 2
    dens = sp.expand(dens)
    A = sp.integrate(dens, (z, 0, Y), (x, -X, X), (y, -X, X))
    Bd = sp.expand(sum(-sp.diff(v[i], z) * v[i] for i in range(3)).subs(z, 0))
    B = sp.integrate(Bd, (x, -X, X), (y, -X, X))
    C = sp.integrate(sp.expand(sum(vi ** 2 for vi in v).subs(z, 0)), (x, -X, X), (y, -X, X))
    return sp.nsimplify(Y * (2 * B - A) / C), A, B, C

if __name__ == "__main__":
    for M in ["2", "3", "4", "8"]:
        L, A, B, C = lam(M)
        print(f"X/Y={M}: Y(2B-A)/C = {L} = {float(L):.6f}   (A={float(A):.5g}, B={float(B):.5g}, C={float(C):.5g})")
    # 1D limit X/Y -> infinity: -2 g'(0) g''(0) - int g''^2  (g'(0)=1)
    t = sp.symbols("t"); g = t * (1 - t) ** 6
    lim = -2 * sp.diff(g, t).subs(t, 0) * sp.diff(g, t, 2).subs(t, 0) - sp.integrate(sp.diff(g, t, 2) ** 2, (t, 0, 1))
    print("X/Y -> inf limit:", lim, float(lim))
