"""
Rigorous interval certificate for (P^T), k = 4 (robust5 REPORT, Sec. 2-3).

  kappa_T^2 >= 1 / ( (3780^2 Q00 + 2520^2 Q11) * lambda_max(G_W) ),
  Qff = 13 [ g_ff(x1,y1) + gC_f^2 (x1^2+1) F(x2,y2) ] + (3 CI2 + CT2) e_ff(x1,y1) + 3 CT2 d_ff(x1,y1),
  CI2 = 20 s1,  CT2 = 2 s1 (D^2/pi^2 + R D/pi),  R^2 = max(x1^2+1, y1^2+1)/s1^2,  D^2 = max(1, R^2),  s1 = x1+y1,
  gC = (-1/20, -1/10),  F(x2,y2) = 10 p(x2,y2) / (7 s2^3),
  G_W = Frobenius Gram of m1 = l_z l_a^3, m2 = l_z^2 l_a^2 with l_a = x - x1 y, l_z = -x - y1 y.
(closed forms from sym_Q.py; cross-checked here against the generic exact pipeline pt_fields.py.)

Admissible shapes (M0_theta): every angle of T1 and of T2 >= theta.  In cotangent coordinates (x = cot at one base
vertex, y = cot at the other) this is  x, y <= cot(theta),  1 - x y <= cot(theta) (x + y),  x + y > 0; it implies
x, y >= -cot(2 theta) and x + y >= sin(theta).

Branch and bound (python-flint arb, outward-rounded ball arithmetic):
  step 1: F_max(theta) >= sup F over admissible (x2, y2)            (maximisation)
  step 2: kappa^2 >= kappa0^2 on all admissible (x1, y1) with F := F_max  (Q is increasing in F, and the bound is
          decreasing in Q, so this is valid).
A box is discarded only if it provably contains no admissible point; it is accepted only if the ball evaluation over
the WHOLE box proves the inequality; otherwise it is bisected along its longer side.

usage: python cert_PT.py THETA_DEG KAPPA0 [--check]
"""
import sys, time, random
from fractions import Fraction as Fr
import flint
from flint import arb

sys.setrecursionlimit(10000)


def gff(x1, y1):
    s1 = x1 + y1
    g00 = (4307*x1**4 - 953*x1**3*y1 + 1173*x1**2*y1**2 + 15670*x1**2 + 1427*x1*y1**3 + 10334*x1*y1 + 3817*y1**4
           + 12310*y1**2 + 8823) / (54880 * s1)
    g11 = (3047*x1**4 + 622*x1**3*y1 + 487*x1**2*y1**2 + 10399*x1**2 + 734*x1*y1**3 + 10236*x1*y1 + 2725*y1**4
           + 9643*y1**2 + 4903) / (13720 * s1)
    e00 = (981*x1**2 - 1230*x1*y1 + 925*y1**2 + 3136) / 123480
    e11 = (2881*x1**2 - 3646*x1*y1 + 2881*y1**2 + 9408) / 123480
    d00 = 3 * (2304*x1**4 - 2100*x1**3*y1 + 1062*x1**2*y1**2 + 8513*x1**2 + 2240*x1*y1**3 + 9130*x1*y1 + 3046*y1**4
               + 6889*y1**2) / 13720
    d11 = (9666*x1**4 - 5439*x1**3*y1 - 3428*x1**2*y1**2 + 36798*x1**2 + 4221*x1*y1**3 + 51644*x1*y1 + 14496*y1**4
           + 36798*y1**2) / 6860
    g01 = (3334*x1**4 + 153*x1**3*y1 + 193*x1**2*y1**2 + 11344*x1**2 + 1399*x1*y1**3 + 10824*x1*y1 + 2928*y1**4
           + 9286*y1**2 + 4903) / (27440 * s1)
    e01 = (2937*x1**2 - 3646*x1*y1 + 2825*y1**2 + 9408) / 246960
    d01 = (10422*x1**4 - 8064*x1**3*y1 - 1811*x1**2*y1**2 + 39234*x1**2 + 9786*x1*y1**3 + 51644*x1*y1 + 15063*y1**4
           + 34362*y1**2) / 13720
    return (g00, e00, d00), (g11, e11, d11), (g01, e01, d01)


def Ffun(x2, y2):
    s2 = x2 + y2
    p = 9*x2**4 + 27*x2**3*y2 + 39*x2**2*y2**2 + 4*x2**2 + 27*x2*y2**3 + 2*x2*y2 + 9*y2**4 + 4*y2**2 + 3
    return 10 * p / (7 * s2**3)


def GW(x1, y1):
    one = 1 + 0 * x1
    la = [-x1, one]        # coefficient lists [y-coef, x-coef]: l = a*y + b*x, poly as coefficients of x^j y^(1-j)
    lz = [-y1, -one]

    def mul(p, q):
        out = [0] * (len(p) + len(q) - 1)
        for i, a in enumerate(p):
            for j, b in enumerate(q):
                out[i + j] = out[i + j] + a * b
        return out
    m1 = mul(mul(mul(lz, la), la), la); m2 = mul(mul(mul(lz, lz), la), la)
    binom = [1, 4, 6, 4, 1]
    ip = lambda p, q: sum(p[j] * q[j] / binom[j] for j in range(5))
    return ip(m1, m1), ip(m1, m2), ip(m2, m2)


def GWperp_cf(x1, y1):
    """closed form of GWperp (sympy; cross-checked exactly in --check)"""
    s4 = (x1 + y1) ** 4
    Dl = (12*x1**4*y1**4 + 8*x1**4*y1**2 + 2*x1**4 - 16*x1**3*y1**3 - 4*x1**3*y1 + 8*x1**2*y1**4 + 18*x1**2*y1**2
          + 3*x1**2 - 4*x1*y1**3 - 6*x1*y1 + 2*y1**4 + 3*y1**2 + 2)
    a = s4 * (12*x1**4*y1**2 + 3*x1**4 - 12*x1**3*y1 + 12*x1**2*y1**2 + 8*x1**2 - 8*x1*y1 + 2*y1**2 + 3) / (4 * Dl)
    b = -s4 * (4*x1**4*y1**2 + x1**4 - 4*x1**3*y1**3 - 6*x1**3*y1 + 9*x1**2*y1**2 + 4*x1**2 - 2*x1*y1**3 - 6*x1*y1
               + 2*y1**2 + 2) / (2 * Dl)
    c = s4 * (4*x1**4*y1**2 + x1**4 - 8*x1**3*y1**3 - 8*x1**3*y1 + 4*x1**2*y1**4 + 18*x1**2*y1**2 + 6*x1**2
              - 8*x1*y1**3 - 12*x1*y1 + y1**4 + 6*y1**2 + 4) / (3 * Dl)
    return a, b, c


def GWperp(x1, y1):
    """Frobenius Gram of the projections of m1, m2 onto Ker_T^perp:  G_WW - B G_KK^{-1} B^T,
    B_ji = <m_j, l_i^4>_F = m_j(g_i),  (G_KK)_il = (g_i . g_l)^4,  g_z = (-1,-y1), g_a = (1,-x1), g_c = -g_z-g_a."""
    a, b, c = GW(x1, y1)
    gz = (-1 + 0 * x1, -y1); ga = (1 + 0 * x1, -x1); gc = (-gz[0] - ga[0], -gz[1] - ga[1])
    gs = (gz, ga, gc)
    dot = lambda u, w: u[0] * w[0] + u[1] * w[1]
    Bm = [[dot(gz, g) * dot(ga, g) ** 3 for g in gs], [dot(gz, g) ** 2 * dot(ga, g) ** 2 for g in gs]]
    K = [[dot(u, w) ** 4 for w in gs] for u in gs]
    # adjugate / det of the 3x3 K
    det = (K[0][0] * (K[1][1] * K[2][2] - K[1][2] * K[2][1]) - K[0][1] * (K[1][0] * K[2][2] - K[1][2] * K[2][0])
           + K[0][2] * (K[1][0] * K[2][1] - K[1][1] * K[2][0]))
    adj = [[K[(j + 1) % 3][(i + 1) % 3] * K[(j + 2) % 3][(i + 2) % 3] - K[(j + 1) % 3][(i + 2) % 3] * K[(j + 2) % 3][(i + 1) % 3]
            for j in range(3)] for i in range(3)]
    corr = [[sum(Bm[p][i] * adj[i][l] * Bm[q][l] for i in range(3) for l in range(3)) / det for q in range(2)] for p in range(2)]
    return a - corr[0][0], b - corr[0][1], c - corr[1][1]


# ---------------------------------------------------------------- Taylor-form enclosures of the closed forms
_POLYS = None


def polys():
    global _POLYS
    if _POLYS is None:
        import sympy as sp
        from taylor2 import P2
        X, Y = sp.symbols('x1 y1')
        (g00, e00, d00), (g11, e11, d11), (g01, e01, d01) = gff(X, Y)
        a, b, c = GWperp_cf(X, Y)
        s1 = X + Y
        ex = dict(G00=g00 * 54880 * s1, G11=g11 * 13720 * s1, G01=g01 * 27440 * s1, E00=e00, E11=e11, E01=e01,
                  D00=d00, D11=d11, D01=d01, S1=s1, XX=X**2 + 1, YY=Y**2 + 1)
        Dl = (12*X**4*Y**4 + 8*X**4*Y**2 + 2*X**4 - 16*X**3*Y**3 - 4*X**3*Y + 8*X**2*Y**4 + 18*X**2*Y**2
              + 3*X**2 - 4*X*Y**3 - 6*X*Y + 2*Y**4 + 3*Y**2 + 2)
        ex["DL"] = Dl
        ex["A"] = sp.cancel(a * Dl); ex["B"] = sp.cancel(b * Dl); ex["C"] = sp.cancel(c * Dl)
        _POLYS = {k: P2(sp.cancel(v), X, Y) for k, v in ex.items()}
    return _POLYS


def kappa2_box(x0, y0, rx, ry, Fmax):
    """ball enclosure of the kappa^2 lower-bound function over the box [x0 +- rx] x [y0 +- ry]  (all fmpq)"""
    P = polys(); E = {k: p.enclose(x0, y0, rx, ry) for k, p in P.items()}
    s1 = E["S1"]
    if s1.lower() <= 0 or E["DL"].lower() <= 0:
        return None
    g00, g11, g01 = E["G00"] / (54880 * s1), E["G11"] / (13720 * s1), E["G01"] / (27440 * s1)
    PI = arb.pi()
    CI2 = 20 * s1
    R2 = E["XX"].max(E["YY"]) / (s1 * s1)
    D2 = R2.max(arb(1))
    CT2 = 2 * s1 * (D2 / (PI * PI) + (R2 * D2).sqrt() / PI)
    x2p1 = E["XX"]
    Q00 = 13 * (g00 + x2p1 * Fmax / 400) + (3 * CI2 + CT2) * E["E00"] + 3 * CT2 * E["D00"]
    Q11 = 13 * (g11 + x2p1 * Fmax / 100) + (3 * CI2 + CT2) * E["E11"] + 3 * CT2 * E["D11"]
    Q01 = 13 * (g01 + x2p1 * Fmax / 200) + (3 * CI2 + CT2) * E["E01"] + 3 * CT2 * E["D01"]
    a = E["A"] / E["DL"]; b = E["B"] / E["DL"]; c = E["C"] / E["DL"]
    return 1 / (2520**2 * Q11 * a + 2 * 2520 * 3780 * Q01 * b + 3780**2 * Q00 * c)


_POLYS2 = {}


def polys2(Fmax):
    """S * DL = P0/s1 + CT2 * P1 with P0, P1 polynomials (Fmax a rational number)."""
    key = str(Fmax)
    if key not in _POLYS2:
        import sympy as sp
        from taylor2 import P2
        X, Y = sp.symbols('x1 y1')
        Fq = sp.Rational(str(Fmax))
        (g00, e00, d00), (g11, e11, d11), (g01, e01, d01) = gff(X, Y)
        s1 = X + Y; xx = X**2 + 1
        q0 = {(0, 0): 13 * (g00 + xx * Fq / 400) + 60 * s1 * e00, (1, 1): 13 * (g11 + xx * Fq / 100) + 60 * s1 * e11,
              (0, 1): 13 * (g01 + xx * Fq / 200) + 60 * s1 * e01}
        q1 = {(0, 0): e00 + 3 * d00, (1, 1): e11 + 3 * d11, (0, 1): e01 + 3 * d01}
        a, b, c = GWperp_cf(X, Y)
        Dl = (12*X**4*Y**4 + 8*X**4*Y**2 + 2*X**4 - 16*X**3*Y**3 - 4*X**3*Y + 8*X**2*Y**4 + 18*X**2*Y**2
              + 3*X**2 - 4*X*Y**3 - 6*X*Y + 2*Y**4 + 3*Y**2 + 2)
        A, B, C = sp.cancel(a * Dl), sp.cancel(b * Dl), sp.cancel(c * Dl)
        comb_ = lambda q: 2520**2 * q[(1, 1)] * A + 2 * 2520 * 3780 * q[(0, 1)] * B + 3780**2 * q[(0, 0)] * C
        P0 = sp.cancel(comb_(q0) * s1); P1 = sp.cancel(comb_(q1))
        _POLYS2[key] = dict(P0=P2(P0, X, Y), P1=P2(P1, X, Y), DL=P2(Dl, X, Y), S1=P2(s1, X, Y),
                            XX=P2(xx, X, Y), YY=P2(Y**2 + 1, X, Y))
    return _POLYS2[key]


def kappa2_box2(x0, y0, rx, ry, Fmax):
    P = polys2(Fmax); E = {k: p.enclose(x0, y0, rx, ry) for k, p in P.items()}
    s1 = E["S1"]
    if s1.lower() <= 0 or E["DL"].lower() <= 0:
        return None
    PI = arb.pi()
    R2 = E["XX"].max(E["YY"]) / (s1 * s1)
    D2 = R2.max(arb(1))
    CT2 = 2 * s1 * (D2 / (PI * PI) + (R2 * D2).sqrt() / PI)
    S = (E["P0"] / s1 + CT2 * E["P1"]) / E["DL"]
    return 1 / S


def kappa2(x1, y1, Fmax, PI, amax, asqrt):
    s1 = x1 + y1
    (g00, e00, d00), (g11, e11, d11), (g01, e01, d01) = gff(x1, y1)
    CI2 = 20 * s1
    R2 = amax(x1 * x1 + 1, y1 * y1 + 1) / (s1 * s1)
    D2 = amax(R2, 1 + 0 * R2)
    R = asqrt(R2); D = asqrt(D2)
    CT2 = 2 * s1 * (D2 / (PI * PI) + R * D / PI)
    x2p1 = x1 * x1 + 1
    Q00 = 13 * (g00 + x2p1 * Fmax / 400) + (3 * CI2 + CT2) * e00 + 3 * CT2 * d00
    Q11 = 13 * (g11 + x2p1 * Fmax / 100) + (3 * CI2 + CT2) * e11 + 3 * CT2 * d11
    Q01 = 13 * (g01 + x2p1 * Fmax / 200) + (3 * CI2 + CT2) * e01 + 3 * CT2 * d01
    a, b, c = GWperp_cf(x1, y1)
    # kappa^2 >= 1 / tr(P^-1 Q P^-T Gperp),  P^-1 = [[0,-2520],[-3780,0]]
    return 1 / (2520**2 * Q11 * a + 2 * 2520 * 3780 * Q01 * b + 3780**2 * Q00 * c)


# ---------------------------------------------------------------- cross-check of the closed forms (exact)
def check_closed_forms(n=6, seed=3):
    import pt_fields as PF
    random.seed(seed)
    for _ in range(n):
        p = [Fr(random.randint(-3, 20), 10) for _ in range(4)]
        if p[0] + p[1] <= 0 or p[2] + p[3] <= 0:
            continue
        G, V, info = PF.star_fields(*p, PF.NumExact)
        T1only = [{1: F[1], 2: [{a: 0 for a in PF.I4} for _ in range(2)]} for F in V]
        (g00, e00, d00), (g11, e11, d11), (g01, e01, d01) = gff(p[0], p[1])
        ok = [PF.gram_grad(T1only[0], T1only[1], G) == g01, PF.gram_e(V[0], V[1]) == e01, PF.gram_dn(V[0], V[1], G) == d01,
              PF.gram_grad(T1only[0], T1only[0], G) == g00, PF.gram_e(V[0], V[0]) == e00,
              PF.gram_dn(V[0], V[0], G) == d00, PF.gram_grad(T1only[1], T1only[1], G) == g11,
              PF.gram_e(V[1], V[1]) == e11, PF.gram_dn(V[1], V[1], G) == d11]
        full = [PF.gram_grad(V[f], V[f], G) - PF.gram_grad(T1only[f], T1only[f], G) for f in range(2)]
        gC = (Fr(-1, 20), Fr(-1, 10))
        ok += [full[f] == gC[f] ** 2 * (p[0] ** 2 + 1) * Ffun(p[2], p[3]) for f in range(2)]
        # G_W against pt_fields' construction from barycentric gradients
        g = G["g1"]
        lz = [g[0][1], g[0][0]]; la = [g[1][1], g[1][0]]
        m1 = PF.lin_mul(PF.lin_mul(PF.lin_mul(lz, la), la), la); m2 = PF.lin_mul(PF.lin_mul(PF.lin_mul(lz, lz), la), la)
        GWp = PF.frob_gram([m1, m2])
        a, b, c = GW(p[0], p[1])
        ok += [GWp[0][0] == a, GWp[0][1] == b, GWp[1][1] == c]
        ok += [u == w for u, w in zip(GWperp(p[0], p[1]), GWperp_cf(p[0], p[1]))]
        print("closed forms == pipeline at", [str(q) for q in p], ":", all(ok))
        assert all(ok)


# ---------------------------------------------------------------- branch and bound
class Box:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi):
        self.lo, self.hi = lo, hi

    def balls(self):
        return [arb(Fr(l + h, 2)) + arb(0, float(Fr(h - l, 2)) * (1 + 1e-15)) for l, h in zip(self.lo, self.hi)]

    def split(self):
        i = max(range(2), key=lambda j: self.hi[j] - self.lo[j])
        m = (self.lo[i] + self.hi[i]) / 2
        a_hi = list(self.hi); a_hi[i] = m; b_lo = list(self.lo); b_lo[i] = m
        return Box(list(self.lo), a_hi), Box(b_lo, list(self.hi))


def ball_of(lo, hi):
    """arb ball containing [lo, hi] (Fractions)"""
    mid = arb(lo.numerator) / lo.denominator / 2 + arb(hi.numerator) / hi.denominator / 2
    rad = (arb(hi.numerator) / hi.denominator - arb(lo.numerator) / lo.denominator) / 2
    return mid + arb(0, 1) * rad.upper()


def infeasible(x, y, C, SIN):
    """True only if the box provably contains no admissible (x, y)"""
    if (1 - x * y - C * (x + y)).lower() > 0:
        return True
    if (x + y).upper() < SIN.lower():
        return True
    if x.lower() > C.upper() or y.lower() > C.upper():      # an angle < theta
        return True
    return False


def amax(a, b):
    return a.max(b)


def asqrt(a):
    return a.sqrt()


def bb(theta_deg, func, accept, maxboxes=2_000_000, minwidth=Fr(1, 2 ** 22), label=""):
    th = arb.pi() * theta_deg / 180
    C = 1 / th.tan(); C2 = 1 / (2 * th).tan(); SIN = th.sin()
    import math
    lo = -Fr(math.ceil(float(C2.mid()) * 10**6) + 10, 10**6)      # <= -cot(2 theta)
    hi = Fr(math.ceil(float(C.mid()) * 10**6) + 10, 10**6)        # >= cot(theta)
    stack = [Box([lo, lo], [hi, hi])]
    stats = dict(accepted=0, discarded=0, split=0, worst=None)
    t0 = time.time()
    while stack:
        B = stack.pop()
        x = ball_of(B.lo[0], B.hi[0]); y = ball_of(B.lo[1], B.hi[1])
        if infeasible(x, y, C, SIN):
            stats["discarded"] += 1; continue
        if (x + y).lower() > 0:
            val = func(x, y)
            if accept(val):
                stats["accepted"] += 1
                if stats["worst"] is None or val.mid() < stats["worst"]:
                    stats["worst"] = val.mid()
                continue
        if B.hi[0] - B.lo[0] < minwidth and B.hi[1] - B.lo[1] < minwidth:
            raise RuntimeError("REFUTED/unresolved box %s %s value %s" % (B.lo, B.hi, func(x, y) if (x + y).lower() > 0 else None))
        stats["split"] += 1
        if stats["split"] > maxboxes:
            raise RuntimeError("too many boxes")
        stack.extend(B.split())
    stats["time_s"] = round(time.time() - t0, 1)
    return stats


def bb2(theta_deg, Fmax, kappa0, maxboxes=3_000_000, minwidth=Fr(1, 2 ** 24)):
    """branch and bound for  kappa^2 >= kappa0^2  (Taylor-form enclosures, exact rational boxes)"""
    import math
    from flint import fmpq
    th = arb.pi() * theta_deg / 180
    C = 1 / th.tan(); C2 = 1 / (2 * th).tan(); SIN = th.sin()
    lo = -Fr(math.ceil(float(C2.mid()) * 10**6) + 10, 10**6)
    hi = Fr(math.ceil(float(C.mid()) * 10**6) + 10, 10**6)
    k2 = (arb(kappa0) ** 2).upper()
    stack = [Box([lo, lo], [hi, hi])]
    st = dict(accepted=0, discarded=0, split=0, min_lower=None, min_width=None)
    t0 = time.time()
    while stack:
        B = stack.pop()
        x = ball_of(B.lo[0], B.hi[0]); y = ball_of(B.lo[1], B.hi[1])
        if infeasible(x, y, C, SIN):
            st["discarded"] += 1; continue
        q = lambda f: fmpq(f.numerator, f.denominator)
        v = kappa2_box2(q((B.lo[0] + B.hi[0]) / 2), q((B.lo[1] + B.hi[1]) / 2), q((B.hi[0] - B.lo[0]) / 2),
                        q((B.hi[1] - B.lo[1]) / 2), Fmax)
        if v is not None and v.is_finite() and v.lower() >= k2:
            st["accepted"] += 1
            lw = v.lower()
            if st["min_lower"] is None or lw < st["min_lower"]:
                st["min_lower"] = lw
            w = B.hi[0] - B.lo[0]
            if st["min_width"] is None or w < st["min_width"]:
                st["min_width"] = w
            continue
        if B.hi[0] - B.lo[0] < minwidth and B.hi[1] - B.lo[1] < minwidth:
            raise RuntimeError("REFUTED/unresolved box %s %s value %s" % (B.lo, B.hi, v))
        st["split"] += 1
        if st["split"] > maxboxes:
            raise RuntimeError("too many boxes")
        stack.extend(B.split())
    st["time_s"] = round(time.time() - t0, 1)
    st["min_lower"] = str(st["min_lower"]); st["min_width"] = str(st["min_width"])
    return st


if __name__ == "__main__":
    if "--check" in sys.argv:
        check_closed_forms()
    import math
    th = float(sys.argv[1]); kappa0 = sys.argv[2]
    C = 1 / math.tan(math.radians(th)); C2 = 1 / math.tan(math.radians(2 * th))
    n = 400; best = 0
    for i in range(n + 1):
        for j in range(n + 1):
            xx = -C2 + (C + C2) * i / n; yy = -C2 + (C + C2) * j / n
            if xx + yy > 0 and 1 - xx * yy <= C * (xx + yy) + 1e-12:
                best = max(best, Ffun(xx, yy))
    Fmax = math.ceil(best * 1.02)
    print("theta = %g deg: float max of F = %.6g; certify F <= %d" % (th, best, Fmax), flush=True)
    st = bb(th, Ffun, lambda v: v.upper() <= Fmax)
    print("  F certificate:", st, flush=True)
    kmin = 1e9
    for i in range(n + 1):
        for j in range(n + 1):
            xx = -C2 + (C + C2) * i / n; yy = -C2 + (C + C2) * j / n
            if xx + yy > 0 and 1 - xx * yy <= C * (xx + yy) + 1e-12 and xx <= C and yy <= C:
                kmin = min(kmin, kappa2(xx, yy, Fmax, math.pi, max, math.sqrt) ** 0.5)
    print("  float min of the kappa bound on a 401^2 grid: %.4g" % kmin, flush=True)
    st = bb2(th, Fmax, kappa0)
    print("  CERTIFIED: kappa_T >= %s on all admissible (T1, T2) with min angle >= %g deg:" % (kappa0, th), st, flush=True)
