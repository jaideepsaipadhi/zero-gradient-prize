"""Explicit constants for 2D Theorem A (penalty necessity via curl of an Argyris interpolant).

Part 1 (exact, rational): reference Argyris basis on That = conv{(0,0),(1,0),(0,1)} with DOFs
   vertex: f, f_x, f_y, f_xx, f_xy, f_yy   (18)
   edge  : grad f(m_e) . nu_e, nu_0=(1,1) on the hypotenuse [a1,a2], nu_1=(-1,0) on [a0,a2], nu_2=(0,-1) on [a0,a1].
 Rigorous sup bounds of all partial derivatives of order j=1,2 of every basis function over That via
 Bernstein coefficients (max|b_alpha| >= sup|p| on the simplex), exact Fractions.
Part 2: quintic-Hermite midpoint-derivative coefficients.
Part 3: profile sup norms (exact critical points) and the scaled derivative bounds m_j = Y^{j-1} M_j.
Part 4: continuum half-plane values Ahat, Bhat, Chat (Y=1, X=2), exact.
Writes consts.json.
"""
import sympy as sp, json, itertools, math
from fractions import Fraction as Fr

x, y, t, s = sp.symbols('x y t s', real=True)
mons = [x**i * y**(d - i) for d in range(6) for i in range(d + 1)]
V = [(0, 0), (1, 0), (0, 1)]
edges = [((1, 0), (0, 1), (1, 1)), ((0, 0), (0, 1), (-1, 0)), ((0, 0), (1, 0), (0, -1))]  # (p,q,nu), edge e_i opposite a_i


def dofs(f):
    out = []
    for (a, b) in V:
        for g in (f, sp.diff(f, x), sp.diff(f, y), sp.diff(f, x, 2), sp.diff(f, x, y), sp.diff(f, y, 2)):
            out.append(g.subs({x: a, y: b}))
    for (p, q, nu) in edges:
        m = (sp.Rational(p[0] + q[0], 2), sp.Rational(p[1] + q[1], 2))
        out.append((nu[0] * sp.diff(f, x) + nu[1] * sp.diff(f, y)).subs({x: m[0], y: m[1]}))
    return out


M = sp.Matrix([dofs(mm) for mm in mons]).T          # M[i,j] = dof_i(mon_j)
Minv = M.inv()
basis = [sp.expand(sum(Minv[j, i] * mons[j] for j in range(21))) for i in range(21)]
# sanity: dual
for i in range(21):
    d = dofs(basis[i])
    assert all(sp.simplify(d[k] - (1 if k == i else 0)) == 0 for k in range(21))
names = [f"v{v}_{n}" for v in range(3) for n in ("f", "x", "y", "xx", "xy", "yy")] + ["e0(hyp)", "e1", "e2"]


def bern_sup(p):
    """rigorous upper bound of sup_{That}|p| : max |Bernstein coefficient| (degree = total degree, >=1)."""
    p = sp.Poly(sp.expand(p), x, y)
    if p.is_zero:
        return Fr(0)
    n = max(p.total_degree(), 1)
    l0, l1, l2 = sp.symbols('l0 l1 l2')
    H = 0
    for (a, b), c in p.terms():
        H += c * l1**a * l2**b * (l0 + l1 + l2)**(n - a - b)
    H = sp.Poly(sp.expand(H), l0, l1, l2)
    best = Fr(0)
    for (i, j, k), c in H.terms():
        mult = sp.factorial(n) / (sp.factorial(i) * sp.factorial(j) * sp.factorial(k))
        b = abs(Fr(str(sp.nsimplify(c / mult))))
        best = max(best, b)
    return best


def sqrt_up(q):  # rigorous upper bound of sqrt of a Fraction, as float rounded up
    return math.sqrt(float(q)) * (1 + 1e-12)


norms = {}
for i, b in enumerate(basis):
    bx, by = sp.diff(b, x), sp.diff(b, y)
    s1 = [bern_sup(bx), bern_sup(by)]
    s2 = [bern_sup(sp.diff(b, x, 2)), bern_sup(sp.diff(b, x, y)), bern_sup(sp.diff(b, y, 2))]
    # operator norm of gradient <= sqrt(sup bx^2 + sup by^2); of Hessian <= Frobenius <= sqrt(sxx^2+2sxy^2+syy^2)
    n1 = sqrt_up(s1[0]**2 + s1[1]**2)
    n2 = sqrt_up(s2[0]**2 + 2 * s2[1]**2 + s2[2]**2)
    norms[names[i]] = (n1, n2)

# ---- Part 2: quintic Hermite on [0,1], derivative at 1/2 in terms of f(0),f'(0),f''(0),f(1),f'(1),f''(1)
cs = sp.symbols('c0:6')
P = sum(cs[i] * s**i for i in range(6))
f0, f1, g0, g1, k0, k1 = sp.symbols('f0 f1 g0 g1 k0 k1')
sol = sp.solve([P.subs(s, 0) - f0, P.subs(s, 1) - f1, sp.diff(P, s).subs(s, 0) - g0, sp.diff(P, s).subs(s, 1) - g1,
                sp.diff(P, s, 2).subs(s, 0) - k0, sp.diff(P, s, 2).subs(s, 1) - k1], cs)
Hd = sp.expand(sp.diff(P, s).subs(sol).subs(s, sp.Rational(1, 2)))
hc = {str(v): Hd.coeff(v) for v in (f0, f1, g0, g1, k0, k1)}
print("Hermite (Hf)'(1/2) =", Hd)
# bound |(Hf)'(1/2)| <= M6 h^6 * cH, with |f^{(m)}(end)| <= M6 h^6/(6-m)!  (endpoint at the Taylor centre: all zero,
# but we do not use that for edges through a0 except below)
fac = {0: Fr(1, 720), 1: Fr(1, 120), 2: Fr(1, 24)}
cH_full = sum(abs(Fr(str(hc[n]))) * fac[m] for n, m in (("f0", 0), ("f1", 0), ("g0", 1), ("g1", 1), ("k0", 2), ("k1", 2)))
cH_half = sum(abs(Fr(str(hc[n]))) * fac[m] for n, m in (("f1", 0), ("g1", 1), ("k1", 2)))  # edge from a0: data at a0 vanish

# ---- reference interpolation sums P_j (vertex DOFs of the Taylor remainder r, centre c=a0, so v0 DOFs vanish)
# |value|<=1/720, |first|<=1/120, |second|<=1/24 (times M6 h^6); edge DOF gamma_e bounded separately.
w = {"f": Fr(1, 720), "x": Fr(1, 120), "y": Fr(1, 120), "xx": Fr(1, 24), "xy": Fr(1, 24), "yy": Fr(1, 24)}
Pv = [0.0, 0.0]
for v in (1, 2):
    for n in ("f", "x", "y", "xx", "xy", "yy"):
        for j in (0, 1):
            Pv[j] += float(w[n]) * norms[f"v{v}_{n}"][j]
Ne = [[norms[e][j] for e in ("e0(hyp)", "e1", "e2")] for j in (0, 1)]

# ---- Part 3: profile
g = t * (1 - t)**6
phi = (1 - s**2)**6
tmin = sp.Rational(-1, 8)


def supabs(expr, var, a, b):
    expr = sp.expand(expr)
    cands = [a, b]
    d = sp.diff(expr, var)
    if d != 0:
        for r in sp.Poly(d, var).real_roots():
            if a <= r <= b:
                cands.append(r)
    return max(abs(float(sp.N(expr.subs(var, c), 30))) for c in cands) * (1 + 1e-12)


Phi = [supabs(sp.diff(phi, s, i), s, -1, 1) for i in range(7)]
G = [supabs(sp.diff(g, t, i), t, tmin, 1) for i in range(7)]


def cmax(i, j):  # max over unit u of |u1|^i |u2|^(j-i)
    if j == 0:
        return 1.0
    a = i / j
    return (a**(i / 2)) * ((1 - a)**((j - i) / 2)) if 0 < i < j else 1.0


def mj(j, XoverY=2.0):
    """Y^{j-1} * sup ||D^j F||_op, F = Y phi(x/X) g(y/Y); u-sup bounded termwise."""
    return sum(math.comb(j, i) * cmax(i, j) * XoverY**(-i) * Phi[i] * G[j - i] for i in range(j + 1))


m = {j: mj(j) for j in range(1, 7)}
print("Phi", Phi)
print("G", G)
print("m_j", m)

# ---- Part 4: continuum at Y=1, X=2
X, Y = 2, 1
F = phi.subs(s, x / X) * Y * g.subs(t, y / Y)
v1, v2 = sp.diff(F, y), -sp.diff(F, x)
dens = sp.expand(((2 * sp.diff(v1, x))**2 + 2 * (sp.diff(v1, y) + sp.diff(v2, x))**2 + (2 * sp.diff(v2, y))**2) / 2)
Ah = sp.integrate(sp.integrate(dens, (y, 0, 1)), (x, -2, 2))
Bh = sp.integrate(sp.expand((-sp.diff(v1, y) * v1 - sp.diff(v2, y) * v2).subs(y, 0)), (x, -2, 2))
Ch = sp.integrate(sp.expand((v1**2 + v2**2).subs(y, 0)), (x, -2, 2))
print("Ahat,Bhat,Chat =", Ah, Bh, Ch, float(Ah), float(Bh), float(Ch), "quot", float((2 * Bh - Ah) / Ch))

out = dict(norms=norms, Pv=Pv, Ne=Ne, cH_full=float(cH_full), cH_half=float(cH_half), hermite=str(Hd),
           Phi=Phi, G=G, m=m, Ahat=float(Ah), Bhat=float(Bh), Chat=float(Ch),
           Ahat_exact=str(Ah), Bhat_exact=str(Bh), Chat_exact=str(Ch), tmin=float(tmin))
json.dump(out, open(__file__.replace("argyris_consts.py", "consts.json"), "w"), indent=1)
for k_, v_ in norms.items():
    print(f"{k_:9s} |D1|<= {v_[0]:10.3f}  |D2|<= {v_[1]:10.3f}")
print("Pv", Pv, "Ne", Ne, "cH_full", float(cH_full), "cH_half", float(cH_half))
