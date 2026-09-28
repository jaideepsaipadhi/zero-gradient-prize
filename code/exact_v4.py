"""Exact solutions for the v4 tests (sympy -> numpy), in the format of svn.make_exact:
dict(u, gu, f, p) with u(X,Y) -> (...,2), gu -> (...,2,2) [c,d] = d_d u_c, p -> (...).

  circle(kind, nu, conv, amp)   kind "B": psi = amp (r^2-1)^2 (x + y^2/2)/10, p = x^2 y - y + x/2
                                 kind "P": u = 0, same p.
                                 f = -nu Lap u + conv*(u.grad)u + grad p   (conv = 0: Stokes)
  levelset(curve, kind, nu, conv)  psi = phi^2 (x + y^2/2)/10 with phi a level-set function of the
                                 curve (phi = 0 on Gamma, grad phi != 0 there), so u = curl psi vanishes on
                                 Gamma together with psi and grad psi; p = x^2 y - y + x/2 (in Pi_h for k = 4):
                                   ellipse x^2/a^2 + y^2/b^2 = 1:   phi = sqrt(x^2/a^2 + y^2/b^2) - 1
                                   polar r = s(1 + eps cos(m t)):  phi = r - s - s*eps*Re((x+iy)^m)/r^m
  twodisk(c1, c2, a)             u = 0, p = c1 + (c2 - c1) S((x + a)/(2a)), S the C^3 smoothstep
                                 t^4 (35 - 84 t + 70 t^2 - 20 t^3) clamped to [0,1]; p = c1 for x <= -a,
                                 c2 for x >= a; f = grad p (piecewise polynomial).
"""
import numpy as np
import sympy as sy

x, y = sy.symbols('x y', real=True)


def _pack(u, p, nu, conv):
    f = [-nu * (sy.diff(u[i], x, 2) + sy.diff(u[i], y, 2))
         + conv * (u[0] * sy.diff(u[i], x) + u[1] * sy.diff(u[i], y))
         + sy.diff(p, [x, y][i]) for i in range(2)]
    gu = [[sy.diff(u[i], v) for v in (x, y)] for i in range(2)]
    U = sy.lambdify((x, y), u, 'numpy'); Fl = sy.lambdify((x, y), f, 'numpy')
    GU = sy.lambdify((x, y), gu, 'numpy'); P = sy.lambdify((x, y), p, 'numpy')

    def vec(fn):
        def g(X, Y):
            vals = fn(X, Y)
            return np.stack([np.broadcast_to(np.asarray(v, float), np.shape(X)) for v in vals], -1)
        return g

    def gmat(X, Y):
        vals = GU(X, Y)
        return np.stack([np.stack([np.broadcast_to(np.asarray(vals[i][j], float), np.shape(X))
                                   for j in range(2)], -1) for i in range(2)], -2)
    return dict(u=vec(U), f=vec(Fl), gu=gmat,
                p=lambda X, Y: np.broadcast_to(np.asarray(P(X, Y), float), np.shape(X)).copy())


P_B = x**2 * y - y + x / 2


def circle(kind, nu=1.0, conv=0, amp=1.0):
    r2 = x**2 + y**2
    if kind == "B":
        psi = amp * (r2 - 1)**2 * (x + y**2 / 2) / 10
        u = [sy.diff(psi, y), -sy.diff(psi, x)]
    elif kind == "P":
        u = [sy.Integer(0), sy.Integer(0)]
    else:
        raise ValueError(kind)
    return _pack(u, P_B, nu, conv)


def levelset(curve, kind, nu=1.0, conv=0):
    """kind "B" / "B<a>" (pressure amplitude a, velocity unchanged) or "P"."""
    pamp = 1
    if kind.startswith("B") and len(kind) > 1:
        pamp = sy.nsimplify(kind[1:]); kind = "B"
    cx, cy = curve.centre
    X, Y = x - cx, y - cy
    r = sy.sqrt(X**2 + Y**2)
    if curve.__class__.__name__ == "Ellipse":
        phi = sy.sqrt(X**2 / sy.nsimplify(curve.a)**2 + Y**2 / sy.nsimplify(curve.b)**2) - 1
    elif curve.__class__.__name__ == "PolarCurve":
        s, eps, m = sy.nsimplify(curve.s), sy.nsimplify(curve.eps), curve.m
        phi = r - s - s * eps * sy.re(sy.expand((X + sy.I * Y)**m)) / r**m
    elif curve.__class__.__name__ == "Circle":
        phi = r - sy.nsimplify(curve.r)
    else:
        raise ValueError
    if kind == "B":
        psi = phi**2 * (x + y**2 / 2) / 10
        u = [sy.diff(psi, y), -sy.diff(psi, x)]
    elif kind == "P":
        u = [sy.Integer(0), sy.Integer(0)]
    else:
        raise ValueError(kind)
    return _pack(u, pamp * P_B, nu, conv)


def twodisk(c1=1.0, c2=-1.0, a=0.3):
    def S(t):
        t = np.clip(t, 0, 1)
        return t**4 * (35 - 84 * t + 70 * t**2 - 20 * t**3)

    def dS(t):
        tc = np.clip(t, 0, 1)
        return np.where((t > 0) & (t < 1), 140 * tc**3 * (1 - tc)**3, 0.0)

    def p(X, Y):
        return c1 + (c2 - c1) * S((np.asarray(X) + a) / (2 * a)) + 0 * np.asarray(Y)

    def f(X, Y):
        fx = (c2 - c1) * dS((np.asarray(X) + a) / (2 * a)) / (2 * a)
        return np.stack([fx + 0 * np.asarray(Y), 0 * fx + 0 * np.asarray(Y)], -1)
    zero = lambda X, Y: np.zeros(np.shape(X) + (2,))
    return dict(u=zero, gu=lambda X, Y: np.zeros(np.shape(X) + (2, 2)), f=f, p=p)
