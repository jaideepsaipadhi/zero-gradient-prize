"""Check of the continuum identity used in Theorem A (half-plane, n = (0,-1)):
F = phi(x/X) * Y g(y/Y), g = t(1-t)^p, phi = (1-s^2)^q.   v = curl F = (F_y, -F_x).
 Y*(2B - A)/C = -2 g'(0) g''(0) - |g''|^2 - 2 (Y/X)^2 |g'|^2 |phi'|^2/|phi|^2 - (Y/X)^4 |g|^2 |phi''|^2/|phi|^2.
Computed symbolically-exactly (SymPy) from the definitions A = int 1/2 Dv:Dv, B = int_{y=0} (d_n v).v, C = int_{y=0}|v|^2."""
import sympy as sp
x, y = sp.symbols("x y", real=True)
p, q = 6, 6
X, Y = sp.Integer(2), sp.Integer(1)          # M = X/Y = 2
g = lambda t: t * (1 - t) ** p
phi = lambda s: (1 - s ** 2) ** q
F = phi(x / X) * Y * g(y / Y)
v1, v2 = sp.diff(F, y), -sp.diff(F, x)
D11, D22 = 2 * sp.diff(v1, x), 2 * sp.diff(v2, y)
D12 = sp.diff(v1, y) + sp.diff(v2, x)
dens = sp.expand((D11 ** 2 + 2 * D12 ** 2 + D22 ** 2) / 2)   # 1/2 Dv:Dv, Dv = grad v + grad v^T
A = sp.integrate(sp.integrate(dens, (y, 0, Y)), (x, -X, X))
Bd = sp.expand((-sp.diff(v1, y) * v1 - sp.diff(v2, y) * v2).subs(y, 0))
B = sp.integrate(Bd, (x, -X, X))
C = sp.integrate(sp.expand((v1 ** 2 + v2 ** 2).subs(y, 0)), (x, -X, X))
lam = sp.nsimplify((2 * B - A) / C)
t, s = sp.symbols("t s")
gg, pp = g(t), phi(s)
n = lambda f, a, b: sp.integrate(f ** 2, (a, *b))
pred = (-2 * sp.diff(gg, t).subs(t, 0) * sp.diff(gg, t, 2).subs(t, 0) - n(sp.diff(gg, t, 2), t, (0, 1))
        - 2 * (Y / X) ** 2 * n(sp.diff(gg, t), t, (0, 1)) * n(sp.diff(pp, s), s, (-1, 1)) / n(pp, s, (-1, 1))
        - (Y / X) ** 4 * n(gg, t, (0, 1)) * n(sp.diff(pp, s, 2), s, (-1, 1)) / n(pp, s, (-1, 1)))
print("Y*lambda(F) direct =", lam * Y, "=", float(lam * Y))
print("formula            =", sp.nsimplify(pred), "=", float(pred))
print("match:", sp.simplify(lam * Y - pred) == 0)
