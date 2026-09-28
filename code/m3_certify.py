"""
(M3) for k = 4 (hence k = 5, 6, and every k >= 4) on boundary stars: proof support and certificates.
Companion to notes/v4/m3_proof.tex (labels m3:*).  Log: logs/m3_certify.log.

    python3 m3_certify.py [symbolic|dimension|controls|meshes|certify|elementary|all]      (default: all)

Setting.  z = 0 is a vertex of Gamma_h, the triangles at z form a fan T_j = (0, p_{j-1}, p_j), j = 1..m,
counterclockwise, with [0,p_0] = e_z and [0,p_m] = e'_z the two chords.  We use a 3-triangle fan
F = (0; p0,p1,p2,p3): the whole star if m = 3 (then [0,p3] is the second chord), or the three triangles next
to e_z if m >= 4 (then [0,p3] is an interior edge).  D_ij := p_i x p_j (cross product),
    A := D13 D23,   A3 := D01 D02,
and in the barycentrics (l0 at z, l1 at the first, l2 at the second outer vertex) of each triangle
    phi_1 = A D01^2 l0^2 l2^2 (l0 - 3 l1)                                         on T1 = (0,p0,p1)
    phi_2 = l0^2 [ A D01^2 l1^2 (l0+l2) + A3 D23^2 l2^2 (l0+l1) + 2 A D01 D02 l1 l2 (l0+l1+l2)
                   + A D01 (6 D12 - 5 D02 + 2 D01) l1^2 l2 + A3 D23 (6 D12 - 5 D13 + 2 D23) l1 l2^2 ]  on T2
    phi_3 = A3 D23^2 l0^2 l1^2 (l0 - 3 l2)                                        on T3 = (0,p2,p3)
The witness is w = curl phi.  Sections:

  symbolic   phi is C^1, vanishes to first order on the boundary of F, and int d_nn phi = 0 on [0,p0] and
             [0,p3]: verified as identities of rational functions in the free coordinates p1, p2, p3
             (p0 = (1,0) by similarity).  Second moments = |p0|^5 A/30 and |p3|^5 A3/30.  The Plucker relation
             behind the construction.  Exact derivation of the Hessian-norm forms P1, P3, K (cotangent
             formula) and exact rational cross-checks of ||grad w||^2 against direct integration.
  dimension  exact rational linear algebra on the C^1-quintic space of random fans: dim Y(F) = 2m-3,
             dim Y^1(F) = 1 for m = 3 (so kappa_z = kappa_w exactly), the witness spans Y^1(F); the jump when
             the outer vertices are collinear (beta_j = 0).
  controls   negative controls: (a) an m = 3 star with alpha_1+alpha_2 > pi on which M changes sign, so (M3)
             fails somewhere on the segment (the hypothesis is needed); (b) C^2 (Argyris-type) stream functions
             give M = 0; (c) a tampered witness fails the C^1 check; (d) the interval certificate refuses a
             target above the true minimum.
  meshes     exact rational kappa_z (via the closed form, no linear solves) at every boundary vertex of the
             production meshes N = 8..1024, compared with logs/lower_bound_stars.log.
  certify    rigorous interval branch-and-bound (outward-rounded float64, Taylor-enclosed sin/cos) of a
             decomposed lower bound for kappa_w over the families F3(theta_min, delta) (m = 3) and
             Fs(theta_min, delta) (m >= 4), for the (theta_min, delta, kappa_0) listed in TARGETS.
  elementary a crude closed-form lower bound kappa_el(theta_min, delta) valid for every theta_min, delta.

One core, < 300 MB.  Everything except the [time] lines is deterministic.
"""
import os, sys, time, math, itertools
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
from fractions import Fraction as Fr
import numpy as np
import sympy as sy

HERE = os.path.dirname(os.path.abspath(__file__))
x, y, t = sy.symbols("x y t")
L = sy.symbols("l0:3")          # barycentrics
KS = sy.symbols("k0:3")         # cotangents of the angles at the vertices with barycentric l0, l1, l2

# certification targets: (theta_min [deg], delta [deg], family, kappa_0)
TARGETS = [(30, 10, "m3", Fr(2, 1000)), (30, 10, "sub", Fr(12, 10000)),
           (20, 15, "m3", Fr(5, 10000)), (20, 15, "sub", Fr(3, 10000))]


def T(msg, t0):
    print(f"[time] {msg}: {time.time() - t0:.1f}s", flush=True)


# =============================================================== closed-form witness
def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def bary(Q1, Q2):
    """barycentrics (l0 at 0, l1 at Q1, l2 at Q2) of the triangle (0, Q1, Q2) as functions of (x, y)"""
    d = cross(Q1, Q2)
    l1 = cross((x, y), Q2) / d
    l2 = cross(Q1, (x, y)) / d
    return 1 - l1 - l2, l1, l2


def witness_bary():
    """the three pieces as polynomials in (l0,l1,l2) with coefficient symbols for the D's"""
    D01, D02, D12, D13, D23 = sy.symbols("D01 D02 D12 D13 D23")
    A = D13 * D23; A3 = D01 * D02
    l0, l1, l2 = L
    f1 = A * D01**2 * l0**2 * l2**2 * (l0 - 3 * l1)
    f2 = l0**2 * (A * D01**2 * l1**2 * (l0 + l2) + A3 * D23**2 * l2**2 * (l0 + l1)
                  + 2 * A * D01 * D02 * l1 * l2 * (l0 + l1 + l2)
                  + A * D01 * (6 * D12 - 5 * D02 + 2 * D01) * l1**2 * l2
                  + A3 * D23 * (6 * D12 - 5 * D13 + 2 * D23) * l1 * l2**2)
    f3 = A3 * D23**2 * l0**2 * l1**2 * (l0 - 3 * l2)
    return (f1, f2, f3), (D01, D02, D12, D13, D23)


def witness_xy(P, tamper=False):
    """phi_1, phi_2, phi_3 as functions of (x, y) for the fan P = [p0,p1,p2,p3] (sympy numbers or symbols)"""
    (f1, f2, f3), Ds = witness_bary()
    D = lambda i, j: cross(P[i], P[j])
    sub = dict(zip(Ds, (D(0, 1), D(0, 2), D(1, 2), D(1, 3), D(2, 3))))
    if tamper:
        f2 = f2.subs(6, 7)          # replaces the coefficient 6 D12 by 7 D12
    out = []
    for f, (Q1, Q2) in zip((f1, f2, f3), ((P[0], P[1]), (P[1], P[2]), (P[2], P[3]))):
        out.append(f.subs(sub).subs(dict(zip(L, bary(Q1, Q2))), simultaneous=True))
    return out


def vanishes_on(f, Q, R=None):
    """f and grad f vanish identically on the ray through 0 and Q (R None) or on the line through Q and R"""
    s = {x: t * Q[0], y: t * Q[1]} if R is None else {x: Q[0] + t * (R[0] - Q[0]), y: Q[1] + t * (R[1] - Q[1])}
    for g in (f, sy.diff(f, x), sy.diff(f, y)):
        n, _ = sy.fraction(sy.together(g.subs(s, simultaneous=True)))
        if sy.expand(n) != 0:
            return False
    return True


def dnn_moments(f, Q):
    """(int_0^1 Q^perp.H.Q^perp dt,  int_0^1 (t-1/2)^2 Q^perp.H.Q^perp dt) along the ray [0,Q]; H = Hessian of f.
    With n = Q^perp/|Q| and s = |Q|(t-1/2):  int_e d_nn f = |Q|^{-1} * first,  int_e s^2 d_nn f = |Q| * second."""
    n = (-Q[1], Q[0])
    g = (sy.diff(f, x, 2) * n[0]**2 + 2 * sy.diff(f, x, y) * n[0] * n[1] + sy.diff(f, y, 2) * n[1]**2)
    g = sy.together(g.subs({x: t * Q[0], y: t * Q[1]}, simultaneous=True))
    return (sy.factor(sy.integrate(g, (t, 0, 1))), sy.factor(sy.integrate((t - sy.Rational(1, 2))**2 * g, (t, 0, 1))))


# =============================================================== Hessian norm on a triangle (cotangent formula)
def pform(f):
    """for f a polynomial in barycentrics: ||D^2 f||^2_{L^2(T)} = pform(f)(cot g0, cot g1, cot g2) / (2|T|)"""
    def g(a, b):                       # 2|T| grad l_a . grad l_b
        if a == b:
            return sum(KS[c] for c in range(3) if c != a)
        return -KS[3 - a - b]

    def bint(e):                       # int_T l^alpha / (2|T|)
        p = sy.Poly(sy.expand(e), *L); s = 0
        for ex, c in p.terms():
            s += c * sy.Rational(math.prod(math.factorial(i) for i in ex), math.factorial(sum(ex) + 2))
        return s
    d2 = [[sy.diff(f, L[a], L[b]) for b in range(3)] for a in range(3)]
    tot = 0
    for a, b, c, d in itertools.product(range(3), repeat=4):
        tot += g(a, c) * g(b, d) * bint(d2[a][b] * d2[c][d])
    return sy.expand(tot)


_FORMS = {}


def forms():
    """P1 (T1), P3 (T3) and the 5x5 matrix K (T2) of quadratic polynomials in the cotangents"""
    if _FORMS:
        return _FORMS["P1"], _FORMS["P3"], _FORMS["K"]
    l0, l1, l2 = L
    P1 = pform(l0**2 * l2**2 * (l0 - 3 * l1))
    P3 = pform(l0**2 * l1**2 * (l0 - 3 * l2))
    B2 = [l0**2 * l1**2 * (l0 + l2), l0**2 * l2**2 * (l0 + l1), l0**2 * l1 * l2 * (l0 + l1 + l2), l0**2 * l1**2 * l2,
          l0**2 * l1 * l2**2]
    c = sy.symbols("c0:5")
    PB = pform(sum(c[i] * B2[i] for i in range(5)))
    K = [[sy.expand(sy.diff(PB, c[i], c[j]) / 2) for j in range(5)] for i in range(5)]
    _FORMS.update(P1=P1, P3=P3, K=K)
    return P1, P3, K


def cot_at(V, B, C):
    """cot of the angle at V in the triangle (V,B,C) -- exact for rational input"""
    u = (B[0] - V[0], B[1] - V[1]); w = (C[0] - V[0], C[1] - V[1])
    return (u[0] * w[0] + u[1] * w[1]) / abs(cross(u, w))


def kappa_parts_exact(P):
    """exact (rational for rational P) pieces: A, A3, N1, N2, N3 with ||grad w||^2 = N1+N2+N3,
    and M(w) = +-(|p0|^5 A + |p3|^5 A3)/60 (two chords) or +-|p0|^5 A/60 (one chord)."""
    P1, P3, K = forms()
    D = lambda i, j: cross(P[i], P[j])
    O = (0, 0)
    A = D(1, 3) * D(2, 3); A3 = D(0, 1) * D(0, 2)
    cots = [(cot_at(O, P[j], P[j + 1]), cot_at(P[j], P[j + 1], O), cot_at(P[j + 1], O, P[j])) for j in range(3)]
    evf = lambda p, k: p.subs(dict(zip(KS, k)))
    N1 = A**2 * D(0, 1)**3 * evf(P1, cots[0])
    N3 = A3**2 * D(2, 3)**3 * evf(P3, cots[2])
    kv = [A * D(0, 1)**2, A3 * D(2, 3)**2, 2 * A * D(0, 1) * D(0, 2), A * D(0, 1) * (6 * D(1, 2) - 5 * D(0, 2) + 2 * D(0, 1)),
          A3 * D(2, 3) * (6 * D(1, 2) - 5 * D(1, 3) + 2 * D(2, 3))]
    Kv = [[evf(K[i][j], cots[1]) for j in range(5)] for i in range(5)]
    N2 = sum(kv[i] * Kv[i][j] * kv[j] for i in range(5) for j in range(5)) / D(1, 2)
    return A, A3, N1, N2, N3


# fast rational evaluation of the forms (coefficient tables)
_TAB = {}


def form_tables():
    if _TAB:
        return _TAB["t"]
    P1, P3, K = forms()

    def tab(p):
        return [(e, Fr(int(sy.numer(c)), int(sy.denom(c)))) for e, c in sy.Poly(p, *KS).terms()]
    _TAB["t"] = (tab(P1), tab(P3), [[tab(K[i][j]) for j in range(5)] for i in range(5)])
    return _TAB["t"]


def peval(tb, k):
    s = Fr(0)
    for e, c in tb:
        term = c
        for v, n in zip(k, e):
            if n:
                term *= v**n
        s += term
    return s


def kappa_parts_fr(P):
    """same as kappa_parts_exact, with Fractions only (fast)"""
    tP1, tP3, tK = form_tables()
    D = lambda i, j: cross(P[i], P[j])
    O = (Fr(0), Fr(0))
    A = D(1, 3) * D(2, 3); A3 = D(0, 1) * D(0, 2)
    cots = [(cot_at(O, P[j], P[j + 1]), cot_at(P[j], P[j + 1], O), cot_at(P[j + 1], O, P[j])) for j in range(3)]
    N1 = A**2 * D(0, 1)**3 * peval(tP1, cots[0])
    N3 = A3**2 * D(2, 3)**3 * peval(tP3, cots[2])
    kv = [A * D(0, 1)**2, A3 * D(2, 3)**2, 2 * A * D(0, 1) * D(0, 2), A * D(0, 1) * (6 * D(1, 2) - 5 * D(0, 2) + 2 * D(0, 1)),
          A3 * D(2, 3) * (6 * D(1, 2) - 5 * D(1, 3) + 2 * D(2, 3))]
    N2 = sum(kv[i] * peval(tK[i][j], cots[1]) * kv[j] for i in range(5) for j in range(5)) / D(1, 2)
    return A, A3, N1, N2, N3


def sqrt_bounds(q, bits=200):
    """rational lo <= sqrt(q) <= hi for a positive Fraction q"""
    s = 1 << bits
    n = q.numerator * s * s // q.denominator
    r = math.isqrt(n)
    return Fr(r, s), Fr(r + 1, s)


def kappa_lower_fr(P, two_chords, norm_by=0):
    """rigorous rational lower bound of kappa_w = |M(w)| / (|e|^2 ||grad w||), e = [0,p0] (norm_by=0) or [0,p3]"""
    A, A3, N1, N2, N3 = kappa_parts_fr(P)
    R0 = P[0][0]**2 + P[0][1]**2; R3 = P[3][0]**2 + P[3][1]**2
    r0 = sqrt_bounds(R0); r3 = sqrt_bounds(R3)
    terms = [(A, r0)] + ([(A3, r3)] if two_chords else [])
    lo = hi = Fr(0)
    for c, (rl, rh) in terms:
        a, b = c * rl**5, c * rh**5
        lo += min(a, b); hi += max(a, b)
    Mabs_lo = max(lo, -hi, Fr(0)) / 60         # |M| >= this
    E = R0 if norm_by == 0 else R3
    den = E**2 * (N1 + N2 + N3)
    k2 = Mabs_lo**2 / den                       # kappa^2 >= k2
    return k2, (A, A3)


# =============================================================== 1. symbolic verification
def cmd_symbolic():
    t0 = time.time()
    print("== 1. SYMBOLIC: the closed-form witness (free coordinates p1=(X1,Y1), p2=(X2,Y2), p3=(X3,Y3); p0=(1,0))")
    X1, Y1, X2, Y2, X3, Y3 = sy.symbols("X1 Y1 X2 Y2 X3 Y3")
    P = [(sy.Integer(1), sy.Integer(0)), (X1, Y1), (X2, Y2), (X3, Y3)]
    D = lambda i, j: cross(P[i], P[j])
    A = D(1, 3) * D(2, 3); A3 = D(0, 1) * D(0, 2)
    # Plucker relation: A N0^2 + c1 N1^2 + c2 N2^2 - A3 N3^2 = 0 with N_j(x) = p_j x (x,y)
    N = [cross(p, (x, y)) for p in P]
    c1 = -D(0, 2) * D(0, 3) * D(2, 3) / D(1, 2); c2 = D(0, 1) * D(0, 3) * D(1, 3) / D(1, 2)
    rel = sy.together(A * N[0]**2 + c1 * N[1]**2 + c2 * N[2]**2 - A3 * N[3]**2)
    print("   Plucker relation  A N0^2 - (D02 D03 D23/D12) N1^2 + (D01 D03 D13/D12) N2^2 - A3 N3^2 == 0 :",
          sy.expand(sy.fraction(rel)[0]) == 0)
    f1, f2, f3 = witness_xy(P)
    print("   C^1 across the ray [0,p1] (phi1-phi2 and its gradient vanish there):", vanishes_on(f1 - f2, P[1]))
    print("   C^1 across the ray [0,p2]:", vanishes_on(f2 - f3, P[2]))
    print("   phi = grad phi = 0 on [0,p0] and on [0,p3]:", vanishes_on(f1, P[0]), vanishes_on(f3, P[3]))
    print("   phi = grad phi = 0 on the outer edges p0p1, p1p2, p2p3:",
          vanishes_on(f1, P[0], P[1]), vanishes_on(f2, P[1], P[2]), vanishes_on(f3, P[2], P[3]))
    m0, s0 = dnn_moments(f1, P[0]); m3, s3 = dnn_moments(f3, P[3])
    R0 = P[0][0]**2 + P[0][1]**2; R3 = P[3][0]**2 + P[3][1]**2
    print("   int_{[0,p0]} d_nn phi = 0:", sy.simplify(m0) == 0, "   int_{[0,p3]} d_nn phi = 0:", sy.simplify(m3) == 0)
    # int_e s^2 d_nn phi ds = |Q| * second  (see dnn_moments) and must equal |Q|^5 A / 30:  second = |Q|^4 A/30
    print("   int_{[0,p0]} s^2 d_nn phi ds = |p0|^5 A/30 :", sy.simplify(s0 - R0**2 * A / 30) == 0)
    print("   int_{[0,p3]} s^2 d_nn phi ds = |p3|^5 A3/30:", sy.simplify(s3 - R3**2 * A3 / 30) == 0)
    T("symbolic identities", t0)
    # the forms
    t0 = time.time()
    P1, P3, K = forms()
    print("   Hessian-norm forms: ||D^2 f||^2_T = P_f(cot g0, cot g1, cot g2)/(2|T|), g_i the angle at the vertex of l_i")
    print("     P1 (T1, f = l0^2 l2^2 (l0-3 l1))  =", P1)
    print("     P3 (T3, f = l0^2 l1^2 (l0-3 l2))  =", P3)
    print("     K  (T2, basis l0^2 l1^2(l0+l2), l0^2 l2^2(l0+l1), l0^2 l1 l2(l0+l1+l2), l0^2 l1^2 l2, l0^2 l1 l2^2):")
    for i in range(5):
        print("       ", [str(K[i][j]) for j in range(5)])
    # exact cross-check of ||grad w||^2 = sum N_j and of M against direct integration, random rational fans
    u, w = sy.symbols("u w")
    rng = np.random.default_rng(20260928)
    print("   exact cross-checks on random rational fans (||grad w||^2 by direct integration of |D^2 phi|^2):")
    for trial in range(3):
        while True:
            ang = np.sort(rng.uniform(0.15, 3.1, 3)); rad = rng.uniform(0.5, 2.0, 3)
            if np.all(np.diff(np.concatenate([[0], ang])) > 0.3):
                break
        Pn = [(sy.Integer(1), sy.Integer(0))] + [(sy.Rational(round(r * math.cos(a), 3)), sy.Rational(round(r * math.sin(a), 3)))
                                                for a, r in zip(ang, rad)]
        fs = witness_xy(Pn)
        tot = 0
        for j, f in enumerate(fs):
            Q1, Q2 = Pn[j], Pn[j + 1]
            h = sy.diff(f, x, 2)**2 + 2 * sy.diff(f, x, y)**2 + sy.diff(f, y, 2)**2
            g = sy.expand(h.subs({x: u * Q1[0] + w * Q2[0], y: u * Q1[1] + w * Q2[1]}, simultaneous=True))
            tot += abs(cross(Q1, Q2)) * sy.integrate(sy.integrate(g, (w, 0, 1 - u)), (u, 0, 1))
        A_, A3_, N1, N2, N3 = kappa_parts_exact(Pn)
        print(f"     fan {trial}: direct = closed form: {sy.simplify(tot - (N1 + N2 + N3)) == 0}   "
              f"(||grad w||^2 = {float(tot):.6e}, A = {float(A_):.4f}, A3 = {float(A3_):.4f})")
    # reflection invariance of the closed form: mirror (x,y)->(x,-y) and reverse the order of the rays
    ok = True
    for trial in range(4):
        pts = [(Fr(1), Fr(0))] + [(Fr(int(a), 10), Fr(int(b), 10)) for a, b in ((3, 9), (-6, 8), (-11, 1))][:3]
        pts[1] = (pts[1][0] + Fr(trial, 7), pts[1][1]); pts[3] = (pts[3][0], pts[3][1] - Fr(trial, 11))
        if trial == 3:                                    # a star with beta_1 = 0 (p0, p1, p2 collinear)
            pts = [(Fr(1), Fr(0)), (Fr(1, 2), Fr(1)), (Fr(0), Fr(2)), (Fr(-1), Fr(1, 10))]
        mir = [(p[0], -p[1]) for p in reversed(pts)]
        A, A3, N1, N2, N3 = kappa_parts_fr(pts); Am, A3m, N1m, N2m, N3m = kappa_parts_fr(mir)
        ok &= (A == A3m and A3 == Am and N1 == N3m and N3 == N1m and N2 == N2m)
    print(f"   reflection invariance of the witness data (A<->A3, N1<->N3, N2 fixed), exact on 4 fans incl. beta_1 = 0: {ok}")
    T("forms and cross-checks", t0)


# =============================================================== 2. exact dimension counts
def fan_space(P, extra=None):
    """exact basis data of { phi C^1 piecewise quintic on the fan (0;P), phi = grad phi = 0 on its boundary }.
    Returns (unknown symbols, constraint matrix rows as sympy Matrix, list of per-triangle polynomials)."""
    m = len(P) - 1
    cs = []; polys = []
    mons = [(i, j) for i in range(6) for j in range(6 - i)]
    for k in range(m):
        c = sy.symbols(f"c{k}_0:{len(mons)}"); cs += list(c)
        polys.append(sum(ci * x**i * y**j for ci, (i, j) in zip(c, mons)))
    eqs = []

    def on_line(f, Q, R=None):
        s = {x: t * Q[0], y: t * Q[1]} if R is None else {x: Q[0] + t * (R[0] - Q[0]), y: Q[1] + t * (R[1] - Q[1])}
        for g in (f, sy.diff(f, x), sy.diff(f, y)):
            eqs.extend(sy.Poly(sy.expand(g.subs(s, simultaneous=True)), t).all_coeffs())
    for k in range(m):
        on_line(polys[k], P[k], P[k + 1])                  # outer edge
    on_line(polys[0], P[0]); on_line(polys[m - 1], P[m])   # the two boundary rays
    for k in range(m - 1):
        on_line(polys[k] - polys[k + 1], P[k + 1])         # C^1 across interior rays
    if extra:
        eqs.extend(extra(polys))
    Mx, _ = sy.linear_eq_to_matrix(eqs, cs)
    return cs, Mx, polys


def dnn_line_functionals(poly, Q):
    """linear functionals (in the coefficients) int d_nn phi and int (t-1/2)^2 d_nn phi along [0,Q] (unnormalised)"""
    n = (-Q[1], Q[0])
    g = sy.diff(poly, x, 2) * n[0]**2 + 2 * sy.diff(poly, x, y) * n[0] * n[1] + sy.diff(poly, y, 2) * n[1]**2
    g = sy.expand(g.subs({x: t * Q[0], y: t * Q[1]}, simultaneous=True))
    return sy.integrate(g, (t, 0, 1)), sy.integrate(sy.expand((t - sy.Rational(1, 2))**2 * g), (t, 0, 1))


def space_report(P, label, chords=(0, None)):
    """dim Y(F), dim Y^1(F) (means on the chord rays) and max |M| direction check against the witness"""
    m = len(P) - 1
    cs, Mx, polys = fan_space(P)
    ns = Mx.nullspace()
    dimY = len(ns)
    # mean constraints on chord rays: ray 0 always, ray m if two chords
    rays = [(0, P[0])] + ([(m - 1, P[m])] if chords[1] is not None else [])
    rows = [dnn_line_functionals(polys[k], Q)[0] for k, Q in rays]
    # evaluate on the nullspace basis
    Bm = sy.Matrix([[r.subs(dict(zip(cs, v))) for v in ns] for r in rows]) if ns else sy.zeros(len(rows), 0)
    sub = Bm.nullspace() if ns else []
    dimY1 = len(sub)
    return dimY, dimY1, (cs, ns, sub, polys)


def cmd_dimension():
    t0 = time.time()
    print("\n== 2. DIMENSION (exact rational linear algebra on C^1 piecewise quintics, phi = grad phi = 0 on the fan boundary)")
    R = sy.Rational
    fans = {
        "m=3 star, angle pi+0.1 (generic)": [(R(1), R(0)), (R(1, 5), R(9, 10)), (R(-7, 10), R(3, 4)), (R(-11, 10), R(-1, 9))],
        "m=3 star, angle < pi (convex boundary)": [(R(1), R(0)), (R(3, 10), R(4, 5)), (R(-1, 2), R(7, 10)), (R(-1), R(1, 10))],
        "m=3 star, straight boundary (angle = pi)": [(R(1), R(0)), (R(1, 2), R(1)), (R(-3, 5), R(4, 5)), (R(-6, 5), R(0))],
        "m=3 star, p0 p1 p2 collinear (beta_1 = 0)": [(R(1), R(0)), (R(1, 2), R(1)), (R(0), R(2)), (R(-1), R(1, 10))],
    }
    for name, P in fans.items():
        dimY, dimY1, (cs, ns, sub, polys) = space_report(P, name, chords=(0, 3))
        # witness membership: express the witness polynomials and test the constraints
        f = witness_xy(P)
        wvec = []
        mons = [(i, j) for i in range(6) for j in range(6 - i)]
        for k in range(3):
            pk = sy.Poly(sy.expand(sy.together(f[k])), x, y)
            wvec += [pk.coeff_monomial(x**i * y**j) for (i, j) in mons]
        _, Mx, _ = fan_space(P)
        inY = all(v == 0 for v in (Mx * sy.Matrix(wvec)))
        a0, _ = dnn_line_functionals(polys[0], P[0]); a3, _ = dnn_line_functionals(polys[2], P[3])
        means = [a0.subs(dict(zip(cs, wvec))), a3.subs(dict(zip(cs, wvec)))]
        Dd = lambda i, j: cross(P[i], P[j])
        print(f"   {name}: dim Y = {dimY} (2m-3 = 3), dim Y^1 (both chords) = {dimY1};  witness in Y: {inY}, "
              f"means {means};  A = {Dd(1,3)*Dd(2,3)}, A3 = {Dd(0,1)*Dd(0,2)}")
    # m = 4: single chord / two chords
    P4 = [(R(1), R(0)), (R(1, 2), R(4, 5)), (R(-1, 5), R(1)), (R(-4, 5), R(3, 5)), (R(-6, 5), R(-1, 10))]
    dimY, dimY1, _ = space_report(P4, "m=4", chords=(0, 4))
    print(f"   m=4 star: dim Y = {dimY} (2m-3 = 5), dim Y^1 (both chords) = {dimY1} (2m-5 = 3)")
    # one and two triangles: nothing
    P = [(R(1), R(0)), (R(1, 5), R(1))]
    cs, Mx, _ = fan_space(P)
    print(f"   fan of 1 triangle: dim Y = {len(Mx.nullspace())}")
    P = [(R(1), R(0)), (R(1, 5), R(1)), (R(-1), R(1, 5))]
    dimY, dimY1, _ = space_report(P, "m=2", chords=(0, 2))
    print(f"   fan of 2 triangles (m = 2, excluded by (M1')): dim Y = {dimY}, dim Y^1 = {dimY1}")
    T("dimension", t0)


# =============================================================== 3. negative controls
def cmd_controls():
    t0 = time.time()
    print("\n== 3. NEGATIVE CONTROLS")
    R = sy.Rational
    # (a) sign change of Xi = |p0|^5 A + |p3|^5 A3 along a segment of m=3 stars with alpha_1+alpha_2 > pi
    print("  (a) m = 3 stars with alpha_1 + alpha_2 > pi (outside the hypothesis): Xi = |p0|^5 A + |p3|^5 A3 changes sign.")
    print("      Since dim Y^1 = 1 there (beta_1 beta_2 != 0) and M = +-Xi/60 on the witness, no star-supported field works at the zero.")
    base = [(Fr(1), Fr(0)), (Fr(-1, 10), Fr(1)), (Fr(-1), Fr(-1, 20))]      # p0, p1, p2 with angle(p0,p2) > pi
    vals = []
    for lam in (Fr(3, 10), Fr(2), Fr(6)):
        p3 = (Fr(-2, 5) * lam, Fr(-9, 10) * lam)       # direction beyond p2, scaled
        P = [base[0], base[1], base[2], p3]
        A, A3, *_ = kappa_parts_fr(P)
        R0 = Fr(1); R3 = p3[0]**2 + p3[1]**2
        lo3, hi3 = sqrt_bounds(R3)
        xi_lo = A + min(A3 * lo3**5, A3 * hi3**5); xi_hi = A + max(A3 * lo3**5, A3 * hi3**5)
        D = lambda i, j: cross(P[i], P[j])
        vals.append((lam, xi_lo, xi_hi))
        print(f"      |p3| scale {float(lam):4.1f}: D01,D12,D23 = {float(D(0,1)):.3f},{float(D(1,2)):.3f},{float(D(2,3)):.3f} (>0), "
              f"D02 = {float(D(0,2)):.3f} (<0: alpha1+alpha2 > pi), Xi in [{float(xi_lo):.4e}, {float(xi_hi):.4e}]")
    sign_change = any(a[1] > 0 for a in vals) and any(a[2] < 0 for a in vals)
    print(f"      rigorous sign change of Xi on the segment: {sign_change}")
    # (b) Argyris-type (C^2 at z) stream functions: M = 0
    P = [(R(1), R(0)), (R(1, 5), R(9, 10)), (R(-7, 10), R(3, 4)), (R(-11, 10), R(-1, 9))]

    def c2_at_z(polys):
        e = []
        for k in range(len(polys) - 1):
            d = polys[k] - polys[k + 1]
            for dd in (sy.diff(d, x, 2), sy.diff(d, x, y), sy.diff(d, y, 2)):
                e.append(dd.subs({x: 0, y: 0}))
        return e
    cs, Mx, polys = fan_space(P, extra=c2_at_z)
    ns = Mx.nullspace()
    a0, s0 = dnn_line_functionals(polys[0], P[0]); a3, s3 = dnn_line_functionals(polys[2], P[3])
    vals = [(a0.subs(dict(zip(cs, v))), a3.subs(dict(zip(cs, v))), s0.subs(dict(zip(cs, v))), s3.subs(dict(zip(cs, v)))) for v in ns]
    # on the subspace with zero means, second moments:
    if ns:
        Bm = sy.Matrix([[v[0] for v in vals], [v[1] for v in vals]])
        sub = Bm.nullspace()
        mom = [sum(c * v[2] for c, v in zip(sv, vals)) for sv in sub] + [sum(c * v[3] for c, v in zip(sv, vals)) for sv in sub]
    else:
        sub, mom = [], []
    print(f"  (b) C^2-at-z (Argyris-type) subspace of the generic m=3 star: dim = {len(ns)}, with zero means dim = {len(sub)}, "
          f"so Y^1 contains no nonzero C^2-at-z field (second moments: {mom})")
    # (c) tampered witness
    X1, Y1, X2, Y2, X3, Y3 = sy.symbols("X1 Y1 X2 Y2 X3 Y3")
    PP = [(sy.Integer(1), sy.Integer(0)), (X1, Y1), (X2, Y2), (X3, Y3)]
    f1, f2, f3 = witness_xy(PP, tamper=True)
    print(f"  (c) tampered witness (6 D12 -> 7 D12 in phi_2): C^1 across [0,p1]: {vanishes_on(f1 - f2, PP[1])} (must be False)")
    T("controls", t0)
    return sign_change


# =============================================================== 4. production meshes
def cmd_meshes(Ns=(8, 12, 16, 20, 24, 32, 48, 64, 96, 128, 192, 256, 512, 1024)):
    t0 = time.time()
    sys.path.insert(0, HERE)
    import svn
    print("\n== 4. PRODUCTION MESHES: exact rational kappa_z (closed form; floats of svn.make_mesh read as exact rationals)")
    print("   kappa_z = |M(w)|/(|e_z|^2 ||grad w||), e_z = [c_i, c_{i+1}] (as in lower_bound_tests.py); since beta_1 beta_2 != 0")
    print("   is checked exactly, dim Y^1(star) = 1 and these are the exact star constants.  Also the min over both chords.")
    print("   N     stars  min kappa_z (rig. lower)   mean      max       min over both chords   all beta!=0  all A,A3>0")
    ref = {}
    try:
        for line in open(os.path.join(HERE, "..", "logs", "lower_bound_stars.log")):
            p = line.split()
            if p and p[0].isdigit():
                ref[int(p[0])] = float(p[3])
    except OSError:
        pass
    worst = 1.0
    for N in Ns:
        verts, tris, circ, outer = svn.make_mesh(N)
        V = [(Fr(float(a)), Fr(float(b))) for a, b in verts]
        ks = []; kboth = []; allbeta = True; pos = True
        for i in range(N):
            # circle vertex i is vertex id i; star = triangles 2i, 2i+1, 2i-2 ; e_z = [c_i, c_{i+1}]
            z = V[i]
            ip, im = (i + 1) % N, (i - 1) % N
            b_ip, b_i = N + ip, N + i                                  # layer-1 vertices
            order = [im, b_i, b_ip, ip]                                # counterclockwise around z; e_z = [0,p3]
            assert set(tuple(sorted(tr)) for tr in (tris[2 * i], tris[2 * i + 1], tris[(2 * i - 2) % (2 * N)])) == \
                {tuple(sorted((i, ip, b_ip))), tuple(sorted((i, b_ip, b_i))), tuple(sorted((im, i, b_i)))}
            P = [(V[k][0] - z[0], V[k][1] - z[1]) for k in order]
            D = lambda a, b: cross(P[a], P[b])
            assert D(0, 1) > 0 and D(1, 2) > 0 and D(2, 3) > 0 and D(0, 2) != 0 and D(1, 3) != 0
            for j in (1, 2):
                if cross((P[j - 1][0] - P[j][0], P[j - 1][1] - P[j][1]), (P[j + 1][0] - P[j][0], P[j + 1][1] - P[j][1])) == 0:
                    allbeta = False
            k2, (A, A3) = kappa_lower_fr(P, True, 1)                  # normalised by |e_z| = |p3|
            k2b, _ = kappa_lower_fr(P, True, 0)
            pos &= (A > 0 and A3 > 0)
            ks.append(k2); kboth.append(min(k2, k2b))
        kmin = min(ks)
        klo = math.sqrt(float(kmin)) * (1 - 1e-12)
        kb = math.sqrt(float(min(kboth))) * (1 - 1e-12)
        kf = [math.sqrt(float(q)) for q in ks]
        worst = min(worst, kb)
        rs = f"   (log: {ref[N]:.5f})" if N in ref else ""
        print(f"  {N:5d}  {N:5d}   {klo:.6f}                 {np.mean(kf):.5f}   {max(kf):.5f}   {kb:.6f}               "
              f"{allbeta}         {pos}{rs}", flush=True)
    print(f"   => (M3) holds on every boundary star of every production mesh, with kappa_0 >= {worst:.5f} for either choice of e_z.")
    T("meshes", t0)
    return worst


# =============================================================== 5. interval arithmetic and the certificate
INF = np.inf
_dn = lambda v: np.nextafter(v, -INF)
_up = lambda v: np.nextafter(v, INF)


class I:
    """vectorised interval arithmetic in float64: every correctly-rounded IEEE result is widened by one ulp"""
    __slots__ = ("lo", "hi")

    def __init__(s, lo, hi=None):
        lo = np.asarray(lo, float); s.lo = lo; s.hi = lo if hi is None else np.asarray(hi, float)

    @staticmethod
    def const(q):
        q = Fr(q); f = float(q)
        return I(f, f) if Fr(f) == q else I(_dn(f), _up(f))

    @staticmethod
    def _c(o):
        return o if isinstance(o, I) else I.const(o)

    def __add__(a, b):
        b = I._c(b); return I(_dn(a.lo + b.lo), _up(a.hi + b.hi))
    __radd__ = __add__

    def __neg__(a):
        return I(-a.hi, -a.lo)

    def __sub__(a, b):
        b = I._c(b); return I(_dn(a.lo - b.hi), _up(a.hi - b.lo))

    def __rsub__(a, b):
        return I._c(b) - a

    def __mul__(a, b):
        b = I._c(b)
        p = np.stack([a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi])
        return I(_dn(np.nanmin(p, 0)), _up(np.nanmax(p, 0)))
    __rmul__ = __mul__

    def inv(a):
        bad = (a.lo <= 0) & (a.hi >= 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            lo = _dn(1 / a.hi); hi = _up(1 / a.lo)
        return I(np.where(bad, -INF, lo), np.where(bad, INF, hi))

    def __truediv__(a, b):
        return a * I._c(b).inv()

    def sq(a):
        hi = np.maximum(_up(a.hi * a.hi), _up(a.lo * a.lo))
        lo = np.where((a.lo <= 0) & (a.hi >= 0), 0.0, np.minimum(_dn(a.lo * a.lo), _dn(a.hi * a.hi)))
        return I(np.maximum(lo, 0.0), hi)

    def __pow__(a, n):
        if n == 1:
            return a
        if n % 2 == 0:
            return (a**(n // 2)).sq()
        return a * (a**(n - 1))

    def __getitem__(a, k):
        return I(a.lo[k], a.hi[k])


PI = I(3.141592653589793, 3.1415926535897936)        # math.pi < pi < nextafter(math.pi)


def _taylor(xv, odd, J=20):
    """rigorous enclosure of sin (odd=1) or cos (odd=0) at the points xv, |xv| <= 7"""
    X = I(xv); X2 = X * X
    acc = I.const(Fr((-1)**J, math.factorial(2 * J + odd)))
    for j in range(J - 1, -1, -1):
        acc = acc * X2 + I.const(Fr((-1)**j, math.factorial(2 * j + odd)))
    if odd:
        acc = acc * X
    Rm = _up(np.abs(xv)**(2 * J + 2 + odd) / math.factorial(2 * J + 2 + odd)) * 1.0000001   # Lagrange remainder
    return I(_dn(acc.lo - Rm), _up(acc.hi + Rm))


HALFPI = PI * Fr(1, 2); THREEHALFPI = PI * Fr(3, 2)


def sin_iv(a):
    """sin over intervals contained in [0, 2 pi]"""
    sl = _taylor(a.lo, 1); sh = _taylor(a.hi, 1)
    lo = np.minimum(sl.lo, sh.lo); hi = np.maximum(sl.hi, sh.hi)
    hi = np.where((a.lo <= HALFPI.hi) & (a.hi >= HALFPI.lo), 1.0, hi)
    lo = np.where((a.lo <= THREEHALFPI.hi) & (a.hi >= THREEHALFPI.lo), -1.0, lo)
    return I(np.maximum(lo, -1.0), np.minimum(hi, 1.0))


def cot_dec(a):
    """cot over intervals inside (0, pi), where it is decreasing; [-inf, inf] otherwise"""
    bad = (a.lo <= 0) | (a.hi >= PI.lo)
    ch = _taylor(a.hi, 0) / _taylor(a.hi, 1); cl = _taylor(a.lo, 0) / _taylor(a.lo, 1)
    return I(np.where(bad, -INF, ch.lo), np.where(bad, INF, cl.hi))


def q_enclosure(a, b, fallback):
    """q = sin(b)/sin(a+b) on the box a x b.  dq/db = sin(a)/sin^2(a+b) > 0 and dq/da = sin(b) cos(g)/sin^2(g),
    g = pi-a-b, so q is monotone in both variables when g stays on one side of pi/2: then the exact range is
    attained at two corners (evaluated with point enclosures).  Otherwise the plain quotient is used."""
    g_lo = PI.lo - (a.hi + b.hi); g_hi = PI.hi - (a.lo + b.lo)          # conservative range of gamma
    acute = g_hi < HALFPI.lo                                             # cos g > 0: q increasing in a
    obtuse = g_lo > HALFPI.hi                                            # cos g < 0: q decreasing in a

    def qpt(av, bv):
        return sin_iv(I(bv)) / sin_iv(I(av) + I(bv))
    q_min = np.where(acute, qpt(a.lo, b.lo).lo, qpt(a.hi, b.lo).lo)
    q_max = np.where(acute, qpt(a.hi, b.hi).hi, qpt(a.lo, b.hi).hi)
    use = (acute | obtuse) & (a.lo > 0) & (b.lo > 0) & ((a.hi + b.hi) < PI.lo)
    return I(np.where(use, np.maximum(q_min, fallback.lo), fallback.lo), np.where(use, np.minimum(q_max, fallback.hi), fallback.hi))


def _ipoly(tb):
    def f(k):
        acc = None
        for e, c in tb:
            term = None
            for v, n in zip(k, e):
                if n:
                    term = v**n if term is None else term * (v**n)
            term = term * c
            acc = term if acc is None else acc + term
        return acc
    return f


def _terms(al, be, fam, SIN, COT, Q):
    """the pieces of (m3:dec) in interval arithmetic"""
    tP1, tP3, tK = form_tables()
    fP1, fP3 = _ipoly(tP1), _ipoly(tP3)
    s = [SIN(a) for a in al]
    ab = [al[j] + be[j] for j in range(3)]
    q = [Q(al[j], be[j], SIN(be[j]) / SIN(ab[j])) for j in range(3)]    # q_j = r_j/r_{j-1} = sin(beta)/sin(gamma)
    cot = [(COT(al[j]), COT(be[j]), -COT(ab[j])) for j in range(3)]    # cot(gamma) = -cot(alpha+beta)
    s12 = SIN(al[0] + al[1]); s23 = SIN(al[1] + al[2])
    rho = (s[0] * s12) / (s23 * s[2])
    K2 = [[_ipoly(tK[i][j])(cot[1]) for j in range(5)] for i in range(5)]

    def qf(v):
        acc = None
        for i in range(5):
            for j in range(i, 5):
                if v[i] is None or v[j] is None:
                    continue
                term = v[i] * K2[i][j] * v[j] * (1 if i == j else 2)
                acc = term if acc is None else acc + term
        return acc
    F1 = (q[0] * s[0])**3 * fP1(cot[0])
    ut = [s[0].sq(), None, 2 * s[0] * q[1] * s12, s[0] * (6 * q[0] * q[1] * s[1] - 5 * q[1] * s12 + 2 * s[0]), None]
    U = q[0].sq() * qf(ut) / (q[1] * s[1])
    wt = [None, q[1] * s[2].sq(), None, None, s[2] * (6 * s[1] / q[2] - 5 * s23 + 2 * q[1] * s[2])]
    W = rho.sq() * q[0].sq() * q[1] * qf(wt) / s[1]
    Ph3 = rho.sq() * (q[0] * q[1]).sq() * s[2]**3 * fP3(cot[2]) / q[2]
    r3 = q[0] * q[1] * q[2]
    b = 1 + r3**3 * rho if fam == "m3" else None
    return F1, U, W, Ph3, b, r3


def decomposed_bound(al, be, fam):
    """rigorous interval data of the decomposed bound (m3:dec).  Returns (boost, Tsum_hi, r3) with
         kappa_w >= boost.lo / sqrt(3600 * Tsum_hi)   (normalisation |e| = |p0|)."""
    F1, U, W, Ph3, b, r3 = _terms(al, be, fam, sin_iv, cot_dec, q_enclosure)
    su = _up(_up(np.sqrt(np.maximum(U.hi, 0.0))) + _up(np.sqrt(np.maximum(W.hi, 0.0))))
    Tsum = _up(_up(F1.hi + _up(su * su)) + Ph3.hi)
    boost = b if fam == "m3" else I(np.ones_like(Tsum))
    return boost, Tsum, r3


def decomposed_float(al, be, fam):
    """plain float value of the same bound (for searches and controls)"""
    b, Ts, r3 = decomposed_bound([I(np.array([v])) for v in al], [I(np.array([v])) for v in be], fam)
    return float(b.lo[0] / math.sqrt(3600 * Ts[0])), float(r3.lo[0])


def contract(B, fam, TM, DL, sweeps=2):
    """shrink each box to a box containing all its feasible points (outward rounding keeps it a superset):
         gamma_j >= td:  alpha_j + beta_j <= pi - td;   F3: alpha3 = S - alpha1 - alpha2 in [td, pi - 2 td];
         Fs: alpha1 + alpha2 + alpha3 <= pi + dd - td.   Empty boxes come back with lo > hi."""
    B = B.copy(); lo = B[:, :, 0]; hi = B[:, :, 1]
    cap = float((PI - TM).hi); th = float(TM.lo); top3 = float((PI - TM * 2).hi)
    for _ in range(sweeps):
        if fam == "m3":
            # alpha3 = S - a1 - a2 >= td ;  alpha3 <= pi - 2 td ;  gamma3: alpha3 + beta3 <= pi - td
            hi[:, 0] = np.minimum(hi[:, 0], _up(_up(hi[:, 2] - lo[:, 1]) - th))
            hi[:, 1] = np.minimum(hi[:, 1], _up(_up(hi[:, 2] - lo[:, 0]) - th))
            lo[:, 2] = np.maximum(lo[:, 2], _dn(_dn(lo[:, 0] + lo[:, 1]) + th))
            a3lo = _dn(_dn(lo[:, 2] - hi[:, 0]) - hi[:, 1])
            hi[:, 5] = np.minimum(hi[:, 5], _up(cap - a3lo))
            # alpha3 + beta3 <= pi - td  =>  S <= pi - td - beta3 + a1 + a2
            hi[:, 2] = np.minimum(hi[:, 2], _up(_up(_up(cap - lo[:, 5]) + hi[:, 0]) + hi[:, 1]))
            hi[:, 2] = np.minimum(hi[:, 2], _up(_up(hi[:, 0] + hi[:, 1]) + top3))
            for j in range(2):
                hi[:, 3 + j] = np.minimum(hi[:, 3 + j], _up(cap - lo[:, j]))
                hi[:, j] = np.minimum(hi[:, j], _up(cap - lo[:, 3 + j]))
        else:
            tot = float((PI + DL - TM).hi)
            for j in range(3):
                others = _dn(lo[:, (j + 1) % 3] + lo[:, (j + 2) % 3])
                hi[:, j] = np.minimum(hi[:, j], _up(tot - others))
                hi[:, 3 + j] = np.minimum(hi[:, 3 + j], _up(cap - lo[:, j]))
                hi[:, j] = np.minimum(hi[:, j], _up(cap - lo[:, 3 + j]))
    return B


def branch_and_bound(td, dd, fam, kap0, batch=20000, maxboxes=10**9, progress=True, min_width=1e-9):
    """Certify  kappa_w >= kap0  on the family (angles in degrees):
       F3(td,dd): 3 triangles, all 9 angles >= td, |alpha1+alpha2+alpha3 - 180| <= dd, and |p3| <= |p0|
                  (the case |p3| > |p0| is the mirror image; boxes with r3 > 1 are discarded as covered by symmetry)
       Fs(td,dd): 3 triangles, all 9 angles >= td, alpha1+alpha2+alpha3 <= 180 + dd - td.
       Coordinates: (alpha1, alpha2, alpha3 or S = alpha1+alpha2+alpha3 for F3, beta1, beta2, beta3) in radians."""
    TM = PI * Fr(td, 180); DL = PI * Fr(dd, 180)
    lo0 = float(TM.lo); hi0 = float((PI - TM * 2).hi)        # every angle lies in [td, 180-2td]
    boxes = np.array([[[lo0, hi0]] * 6]); depth = np.zeros(1, int)
    if fam == "m3":
        boxes[0, 2] = [float((PI - DL).lo), float((PI + DL).hi)]
    C = I.const(Fr(3600) * Fr(kap0)**2)
    nc = npr = nsym = tot = maxd = 0; t0 = time.time()
    while len(boxes):
        B = boxes[-batch:]; d = depth[-batch:]; boxes = boxes[:-batch]; depth = depth[:-batch]
        B = contract(B, fam, TM, DL)
        al = [I(B[:, j, 0], B[:, j, 1]) for j in range(3)]; be = [I(B[:, 3 + j, 0], B[:, 3 + j, 1]) for j in range(3)]
        prune = np.any(B[:, :, 0] > B[:, :, 1], axis=1)
        if fam == "m3":
            al[2] = al[2] - al[0] - al[1]                      # alpha3 = S - alpha1 - alpha2 (interval)
            prune |= al[2].hi < TM.lo
        else:
            S = al[0] + al[1] + al[2]; prune |= S.lo > (PI + DL - TM).hi
        for j in range(3):
            prune |= (al[j] + be[j]).lo > (PI - TM).hi         # gamma_j < td on the whole box
        ok = np.zeros(len(B), bool); sym = np.zeros(len(B), bool)
        idx = np.where(~prune)[0]
        if len(idx):
            with np.errstate(all="ignore"):
                boost, Ts, r3 = decomposed_bound([v[idx] for v in al], [v[idx] for v in be], fam)
                good = np.isfinite(Ts) & (boost.lo > 0) & (_dn(boost.lo * boost.lo) >= _up(C.hi * Ts))
                if fam == "m3":
                    sym[idx] = (r3.lo > 1.0) & ~good
                ok[idx] = good
        nc += ok.sum(); npr += prune.sum(); nsym += sym.sum(); tot += len(B)
        rest = ~(ok | prune | sym)
        if rest.any():
            R_ = B[rest]; dr = d[rest]; ar = np.arange(len(R_))
            wd = R_[:, :, 1] - R_[:, :, 0]
            if wd.max(1).min() < min_width:            # a tiny box that cannot be certified: refuted
                return False, dict(certified=int(nc), pruned=int(npr), symmetric=int(nsym), total=int(tot), depth=maxd,
                                   pending=int(len(boxes)), refuted_at_deg=np.degrees(R_[wd.max(1).argmin(), :, 0]).round(4).tolist())
            k = wd.argmax(1); mid = 0.5 * (R_[ar, k, 0] + R_[ar, k, 1])
            c1 = R_.copy(); c2 = R_.copy(); c1[ar, k, 1] = mid; c2[ar, k, 0] = mid
            boxes = np.concatenate([boxes, c1, c2]); depth = np.concatenate([depth, dr + 1, dr + 1])
            maxd = max(maxd, int(dr.max()) + 1)
        if progress and tot // 2000000 != (tot - len(B)) // 2000000:
            print(f"[time]    ... {tot} boxes, {len(boxes)} pending, {time.time()-t0:.0f}s", flush=True)
        if tot > maxboxes:
            return False, dict(certified=int(nc), pruned=int(npr), symmetric=int(nsym), total=int(tot), depth=maxd,
                               pending=int(len(boxes)))
    return True, dict(certified=int(nc), pruned=int(npr), symmetric=int(nsym), total=int(tot), depth=maxd, pending=0)


def float_search(td, dd, fam, n=20000, iters=60, seed=7):
    """float minimisation of the decomposed bound (for reporting and for the negative control): batch random
    sampling of the family, then a shrinking random local search, all with vectorised point evaluations"""
    rng = np.random.default_rng(seed); tm = math.radians(td); de = math.radians(dd)

    def vals(P):
        al, be = P[:, :3], P[:, 3:]
        g = np.min(np.minimum(np.minimum(al, be), math.pi - al - be), 1) - tm
        S = al.sum(1)
        g = np.minimum(g, de - np.abs(S - math.pi)) if fam == "m3" else np.minimum(g, math.pi + de - tm - S)
        with np.errstate(all="ignore"):
            b, Ts, r3 = decomposed_bound([I(al[:, j]) for j in range(3)], [I(be[:, j]) for j in range(3)], fam)
            v = b.lo / np.sqrt(3600 * Ts)
        v = np.where(g < 0, 1 - g, v)
        if fam == "m3":
            v = np.where(r3.lo > 1, 1 + r3.lo, v)
        return np.nan_to_num(v, nan=9.0)
    S = math.pi + rng.uniform(-de, de, n) if fam == "m3" else rng.uniform(3 * tm, math.pi + de - tm, n)
    al = rng.dirichlet([1, 1, 1], n) * (S - 3 * tm)[:, None] + tm
    be = tm + rng.uniform(0, 1, (n, 3)) * np.maximum(math.pi - al - 2 * tm, 0)
    P = np.concatenate([al, be], 1); v = vals(P)
    best = P[np.argsort(v)[:8]]
    for it in range(iters):
        r = 0.2 * 0.9**it
        cand = np.concatenate([best, (best[:, None, :] + rng.normal(0, r, (len(best), 250, 6))).reshape(-1, 6)])
        cv = vals(cand); best = cand[np.argsort(cv)[:8]]
    vb = vals(best[:1])[0]
    return float(vb), best[0]


def cmd_certify(targets=TARGETS):
    print("\n== 6. INTERVAL CERTIFICATE of the decomposed bound (m3:dec) on compact families")
    print("   F3(td,dd): m = 3 stars, all angles >= td, |Theta_z - 180| <= dd  [normalised by the longer chord]")
    print("   Fs(td,dd): the 3 triangles next to e_z of an m >= 4 star, all angles >= td, Theta_z <= 180 + dd")
    # consistency: the decomposed bound is below kappa_w (exact closed form) at random points of the families
    rng = np.random.default_rng(11); worst = np.inf
    for fam in ("m3", "sub"):
        tm = math.radians(20); de = math.radians(15); cnt = 0
        while cnt < 150:
            S = math.pi + rng.uniform(-de, de) if fam == "m3" else rng.uniform(3 * tm, math.pi + de - tm)
            al = rng.dirichlet([1, 1, 1]) * (S - 3 * tm) + tm
            if max(al) > math.pi - 2 * tm:
                continue
            be = tm + rng.uniform(0, 1, 3) * (math.pi - al - 2 * tm)
            th = np.cumsum([0, *al]); r = [1.0]
            for j in range(3):
                r.append(r[-1] * math.sin(be[j]) / math.sin(al[j] + be[j]))
            if fam == "m3" and r[3] > 1:
                continue
            P = [(Fr(r[j] * math.cos(th[j])), Fr(r[j] * math.sin(th[j]))) for j in range(4)]
            k2, _ = kappa_lower_fr(P, fam == "m3", 0)
            kd, _ = decomposed_float(al, be, fam)
            worst = min(worst, math.sqrt(float(k2)) / kd); cnt += 1
    print(f"   check: min over 300 random points of kappa_w / (decomposed bound) = {worst:.6f}  (must be >= 1 up to rounding)")
    res = []
    for td, dd, fam, k0 in targets:
        t0 = time.time()
        fmin, xm = float_search(td, dd, fam)
        ok, info = branch_and_bound(td, dd, fam, k0)
        name = "F3" if fam == "m3" else "Fs"
        print(f"   {name}({td},{dd}):  kappa_w >= {float(k0):.5f} CERTIFIED={ok}  {info};  "
              f"float min of the bound {fmin:.6f} at angles(deg) {np.degrees(xm).round(2).tolist()}", flush=True)
        T(f"certificate {name}({td},{dd})", t0)
        res.append((td, dd, fam, k0, ok))
    # negative control: a target 5% above the float minimum of the bound must be refused
    td, dd, fam = 30, 10, "sub"
    fmin, xm = float_search(td, dd, fam)
    k_bad = Fr(fmin * 1.05).limit_denominator(10**9)
    ok, info = branch_and_bound(td, dd, fam, k_bad, maxboxes=3 * 10**6, progress=False)
    print(f"   negative control: Fs(30,10) with kappa_0 = {float(k_bad):.6f} (5% above the float minimum {fmin:.6f} of the "
          f"bound): certified = {ok} (must be False), {info}")
    return res


# =============================================================== 6. elementary bound
def elementary(td, dd):
    """crude closed-form lower bound of kappa_w on F3(td,dd) (without the r3 <= 1 reduction) and Fs(td,dd):
       every factor of (m3:dec) bounded by sigma = sin td, c = cot td, sig* = min(sin 2td, sin(td-dd))."""
    import mpmath as mp
    mp.mp.dps = 30
    th = mp.radians(td); de = mp.radians(dd)
    sg = mp.sin(th); c = mp.cot(th); ss = min(mp.sin(2 * th), mp.sin(th - de))
    P1, P3, K = forms()
    n1 = lambda p: sum(abs(mp.mpf(int(sy.numer(cf))) / int(sy.denom(cf))) for _, cf in sy.Poly(p, *KS).terms())
    k1 = [[n1(K[i][j]) for j in range(5)] for i in range(5)]
    qmax = 1 / sg; rho = 1 / (sg * ss)
    F1 = qmax**3 * n1(P1) * c**2
    ub = [1, 0, 2 * qmax, 6 * qmax**2 + 5 * qmax + 2, 0]
    U = qmax**2 / (sg * sg) * sum(ub[i] * k1[i][j] * ub[j] for i in range(5) for j in range(5)) * c**2
    wb = [0, qmax, 0, 0, 6 * qmax + 5 + 2 * qmax]
    W = rho**2 * qmax**3 / sg * sum(wb[i] * k1[i][j] * wb[j] for i in range(5) for j in range(5)) * c**2
    Ph3 = rho**2 * qmax**4 * qmax * n1(P3) * c**2
    return 1 / mp.sqrt(3600 * (F1 + (mp.sqrt(U) + mp.sqrt(W))**2 + Ph3))


def cmd_elementary():
    print("\n== 5. ELEMENTARY closed-form bound kappa_el(td,dd) (valid on F3 and Fs for every td <= 60, dd < td)")
    for td, dd in ((10, 5), (20, 15), (20, 5), (30, 10), (45, 10)):
        print(f"   td = {td:2d}, dd = {dd:2d}:  kappa_el = {mp_fmt(elementary(td, dd))}")


def mp_fmt(v):
    return f"{float(v):.3e}"


# =============================================================== main
if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    t_all = time.time()
    print("m3_certify.py -- (M3) for k = 4 on boundary stars; see notes/v4/m3_proof.tex")
    if what in ("symbolic", "all"):
        cmd_symbolic()
    if what in ("dimension", "all"):
        cmd_dimension()
    if what in ("controls", "all"):
        cmd_controls()
    if what in ("meshes", "all"):
        cmd_meshes()
    if what in ("elementary", "all"):
        cmd_elementary()
    if what in ("certify", "all"):
        cmd_certify()
    T("total", t_all)
