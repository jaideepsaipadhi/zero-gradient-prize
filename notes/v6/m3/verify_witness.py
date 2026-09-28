"""Independent exact verification of the (M3) witness of notes/v4/m3_proof.tex (Definition m3:wit).

Written from the barycentric formulas only (no code shared with code/m3_certify.py).
Checks, as identities of rational functions in the free coordinates p1=(a1,b1), p2=(a2,b2), p3=(a3,b3)
(p0=(1,0), z=0 -- no loss by similarity):
  (1) phi, grad phi vanish on the fan boundary (rays [0,p0],[0,p3], outer edges p0p1,p1p2,p2p3);
  (2) phi is C^1 across the interior rays [0,p1], [0,p2];
  (3) int_{chord} d_nn phi = 0 on [0,p0], [0,p3];
  (4) int s^2 d_nn phi ds = |p0|^5 A/30 on [0,p0] and |p3|^5 A3/30 on [0,p3]  (A=D13 D23, A3=D01 D02);
  (5) A = r1 r2 r3^2 sin(a2+a3) sin a3 etc. (angular form, checked numerically at random points);
No condition beta_j != 0 or D03 != 0 is used: only D12 (> 0) appears in a denominator.
Run: python3 verify_witness.py        (about 1 minute)
"""
import time, random
import sympy as sp
from sympy import Rational as Q

t0 = time.time()
a1, b1, a2, b2, a3, b3, x, y, t = sp.symbols('a1 b1 a2 b2 a3 b3 x y t')
P = [sp.Matrix([1, 0]), sp.Matrix([a1, b1]), sp.Matrix([a2, b2]), sp.Matrix([a3, b3])]
cr = lambda u, v: u[0] * v[1] - u[1] * v[0]
D = {(i, j): cr(P[i], P[j]) for i in range(4) for j in range(4)}
X = sp.Matrix([x, y])
Z = sp.Matrix([0, 0])

def bary(A_, B_, C_):
    """barycentric coords (lam_A, lam_B, lam_C) of X in triangle A_,B_,C_ (rational in coords)"""
    area2 = cr(B_ - A_, C_ - A_)
    lA = cr(B_ - X, C_ - X) / area2
    lB = cr(C_ - X, A_ - X) / area2
    lC = cr(A_ - X, B_ - X) / area2
    return lA, lB, lC

A = D[1, 3] * D[2, 3]
A3 = D[0, 1] * D[0, 2]
# T1=(0,p0,p1), T2=(0,p1,p2), T3=(0,p2,p3); lambda0 at z, lambda1 at p_{j-1}, lambda2 at p_j
l = [bary(Z, P[j - 1], P[j]) for j in (1, 2, 3)]
L0, L1, L2 = l[0]
phi1 = A * D[0, 1]**2 * L0**2 * L2**2 * (L0 - 3 * L1)
L0, L1, L2 = l[2]
phi3 = A3 * D[2, 3]**2 * L0**2 * L1**2 * (L0 - 3 * L2)
L0, L1, L2 = l[1]
phi2 = L0**2 * (A * D[0, 1]**2 * L1**2 * (L0 + L2) + A3 * D[2, 3]**2 * L2**2 * (L0 + L1)
                 + 2 * A * D[0, 1] * D[0, 2] * L1 * L2 * (L0 + L1 + L2)
                 + A * D[0, 1] * (6 * D[1, 2] - 5 * D[0, 2] + 2 * D[0, 1]) * L1**2 * L2
                 + A3 * D[2, 3] * (6 * D[1, 2] - 5 * D[1, 3] + 2 * D[2, 3]) * L1 * L2**2)
phis = [phi1, phi2, phi3]

def grad(f):
    return [sp.diff(f, x), sp.diff(f, y)]

def on_seg(f, U, V):
    """f restricted to U + t (V-U), as a polynomial in t with rational-function coefficients"""
    return f.subs({x: U[0] + t * (V[0] - U[0]), y: U[1] + t * (V[1] - U[1])}, simultaneous=True)

def is_zero(e):
    e = sp.together(sp.expand(e))
    num = sp.numer(e)
    return sp.expand(num) == 0

res = {}
# (1) boundary
checks = [(phi1, Z, P[0], 'ray p0'), (phi3, Z, P[3], 'ray p3'),
          (phi1, P[0], P[1], 'edge p0p1'), (phi2, P[1], P[2], 'edge p1p2'), (phi3, P[2], P[3], 'edge p2p3')]
for f, U, V, name in checks:
    ok = is_zero(on_seg(f, U, V)) and all(is_zero(on_seg(g, U, V)) for g in grad(f))
    res['vanish+grad on ' + name] = ok
# (2) C^1 across interior rays
for j, (f, g) in ((1, (phi1, phi2)), (2, (phi2, phi3))):
    d = f - g
    ok = is_zero(on_seg(d, Z, P[j])) and all(is_zero(on_seg(gg, Z, P[j])) for gg in grad(d))
    res['C1 across ray p%d' % j] = ok
print('[%.1fs] boundary and C^1 checks:' % (time.time() - t0), res)

# (3),(4): d_nn phi on the chords, with unnormalised normal nu = (-Vy, Vx):  d_nn = nu^T H nu / |V|^2
def hess(f):
    return sp.Matrix([[sp.diff(f, x, x), sp.diff(f, x, y)], [sp.diff(f, x, y), sp.diff(f, y, y)]])

for f, V, Aval, name in ((phi1, P[0], A, 'p0'), (phi3, P[3], A3, 'p3')):
    nu = sp.Matrix([-V[1], V[0]])
    g = (nu.T * hess(f) * nu)[0, 0]          # = |V|^2 d_nn phi
    g = on_seg(g, Z, V)
    I0 = sp.integrate(sp.expand(sp.numer(sp.together(g))), (t, 0, 1)) / sp.denom(sp.together(g))
    I2 = sp.integrate(sp.expand(sp.numer(sp.together(g)) * (t - Q(1, 2))**2), (t, 0, 1)) / sp.denom(sp.together(g))
    V2 = V[0]**2 + V[1]**2
    # int_e d_nn phi ds = |V| * I0/|V|^2 ;  int_e s^2 d_nn phi ds = |V|^3 * I2/|V|^2 = |V| I2.
    # claim: |V| I2 = |V|^5 Aval/30  <=>  I2 = V2^2 Aval/30
    res['mean zero on ' + name] = is_zero(I0)
    res['moment on ' + name] = is_zero(I2 - V2**2 * Aval / 30)
print('[%.1fs] means/moments:' % (time.time() - t0), {k: v for k, v in res.items() if 'mean' in k or 'moment' in k})

# (5) angular form of A, A3 (numerical, random fans)
import math
worst = 0
for _ in range(200):
    th = [0]; r = [1] + [random.uniform(.3, 3) for _ in range(3)]
    for _j in range(3): th.append(th[-1] + random.uniform(.2, 1.4))
    pts = [(r[i] * math.cos(th[i]), r[i] * math.sin(th[i])) for i in range(4)]
    sub = {a1: pts[1][0], b1: pts[1][1], a2: pts[2][0], b2: pts[2][1], a3: pts[3][0], b3: pts[3][1]}
    al = [th[1] - th[0], th[2] - th[1], th[3] - th[2]]
    Aang = r[1] * r[2] * r[3]**2 * math.sin(al[1] + al[2]) * math.sin(al[2])
    A3ang = r[0]**2 * r[1] * r[2] * math.sin(al[0]) * math.sin(al[0] + al[1])
    worst = max(worst, abs(float(A.subs(sub)) - Aang), abs(float(A3.subs(sub)) - A3ang))
res['angular form of A, A3 (max abs err, 200 random fans)'] = worst
print('[%.1fs] ALL:' % (time.time() - t0))
for k, v in res.items(): print('   ', k, ':', v)
assert all(v is True for k, v in res.items() if not k.startswith('angular')) and worst < 1e-9
print('ALL CHECKS PASSED')
