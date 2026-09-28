"""Exact solutions (sympy -> numpy) for the 3D tests.  nu = 1:  -Lap u + grad p = f, div u = 0."""
import numpy as np
import sympy as sy

x, y, z = sy.symbols('x y z', real=True)
X = (x, y, z)


def _curl(psi):
    return [sy.diff(psi[2], y) - sy.diff(psi[1], z), sy.diff(psi[0], z) - sy.diff(psi[2], x),
            sy.diff(psi[1], x) - sy.diff(psi[0], y)]


def make(kind):
    if kind == "mms":               # smooth manufactured solution on the box (no hole)
        psi = [sy.sin(y + 0.3) * sy.cos(z) * sy.exp(x / 2), sy.cos(x - z) * sy.sin(y) / 2, sy.sin(x + y) * sy.cos(1.3 * z)]
        u = _curl(psi); p = sy.sin(x) * sy.cos(y + 0.2) * sy.exp(z / 3)
        f = [-sum(sy.diff(u[i], v, 2) for v in X) + sy.diff(p, X[i]) for i in range(3)]
    elif kind == "rot":             # (1 - r^-3) omega x x, p = 0 : exact Stokes, u = 0 and p const on |x| = 1
        om = (0.3, -0.5, 1.0)
        r3 = (x**2 + y**2 + z**2) ** sy.Rational(3, 2)
        c = [om[1] * z - om[2] * y, om[2] * x - om[0] * z, om[0] * y - om[1] * x]
        u = [(1 - 1 / r3) * ci for ci in c]; p = sy.Integer(0); f = [0, 0, 0]
    elif kind == "sph":             # Stokes flow past the fixed unit sphere, uniform stream U at infinity
        U = (1.0, 0.4, -0.3)
        r = sy.sqrt(x**2 + y**2 + z**2); Ux = U[0] * x + U[1] * y + U[2] * z
        u = [U[i] - sy.Rational(3, 4) * (U[i] / r + Ux * X[i] / r**3) - sy.Rational(1, 4) * (U[i] / r**3 - 3 * Ux * X[i] / r**5)
             for i in range(3)]
        p = -sy.Rational(3, 2) * Ux / r**3; f = [0, 0, 0]
    elif kind == "P":               # pure pressure: u = 0, f = grad p, p = the sphere-flow wall pressure field
        U = (1.0, 0.4, -0.3)
        r = sy.sqrt(x**2 + y**2 + z**2); Ux = U[0] * x + U[1] * y + U[2] * z
        u = [sy.Integer(0)] * 3; p = -sy.Rational(3, 2) * Ux / r**3
        f = [sy.diff(p, v) for v in X]
    else:
        raise ValueError(kind)
    gu = [[sy.diff(u[i], v) for v in X] for i in range(3)]
    Ul = sy.lambdify(X, u, 'numpy'); Gl = sy.lambdify(X, gu, 'numpy')
    Fl = sy.lambdify(X, f, 'numpy'); Pl = sy.lambdify(X, p, 'numpy')

    def vec(fn):
        def g(a, b, c):
            v = fn(a, b, c)
            return np.stack([np.broadcast_to(np.asarray(t, float), np.shape(a)) for t in v], -1)
        return g

    def gmat(a, b, c):
        v = Gl(a, b, c)
        return np.stack([np.stack([np.broadcast_to(np.asarray(v[i][j], float), np.shape(a)) for j in range(3)], -1)
                         for i in range(3)], -2)
    return dict(u=vec(Ul), gu=gmat, f=vec(Fl),
                p=lambda a, b, c: np.broadcast_to(np.asarray(Pl(a, b, c), float), np.shape(a)), sym=(u, p))


def check(kind):
    ex = make(kind); u, p = ex["sym"]
    div = sy.simplify(sum(sy.diff(u[i], X[i]) for i in range(3)))
    return div
