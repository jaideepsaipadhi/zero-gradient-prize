"""
Explicit pressure-blind star fields for (P^T), k = 4 (robust5 REPORT, Sec. 2).

Two-triangle patch Q = T1 u T2 at the boundary vertex z:
  T1 = (z, a, c) = the boundary triangle T_{e_z}, e = [z, a];     T2 = (z, c, b) = its neighbour across [z, c].
Normalisation (similarity): z = (0,0), a = (1,0), interior of Omega_h above e (outward normal n = (0,-1)).
Shape parameters: cotangents  x1 = cot(angle of T1 at z), y1 = cot(angle of T1 at a),
                              x2 = cot(angle of T2 at z), y2 = cot(angle of T2 at c).
  c = (x1, 1)/(x1+y1),   b = (x2 c + R c)/(x2+y2),  R = rotation by +90 degrees.

Fields (Bernstein-Bezier coefficients, degree 4, per triangle, vector valued):
  v_f = Piola(u_f^) + gA curl psi_A + gC curl psi_C,     f = 1, 2,
  u_f^ : reference fields on T^ = (0,0),(1,0),(0,1) vanishing on the two edges other than e, pressure-blind on T^,
         normal traces B_1 - B_2 and B_1 - B_3 (quartic Bernstein on e); they span P_e.
  psi_A = 30 lam_z^2 lam_a^2 lam_c on T1 (BB coefficient d_{221} = 1), 0 on T2;
  psi_C : BB d^{T1}_{212} = 1, d^{T2}_{221} = beta_a   (beta = barycentrics of b w.r.t. T1)  -- C^1 across [z,c].
  gA, gC solve  int_e v.t = 0,  int_e d_n v . t = 0  (a triangular 2x2 system with nonzero diagonal).
The code is generic in the number type: exact (fractions.Fraction) or rigorous balls (flint.arb).
"""
from fractions import Fraction as Fr
from math import factorial as fact, comb
import itertools

K = 4


def midx(n):
    return [(a, b, n - a - b) for a in range(n, -1, -1) for b in range(n - a, -1, -1)]


I4, I3, I5 = midx(4), midx(3), midx(5)


def iz(x):
    """structurally zero (int/Fraction 0); balls are never skipped"""
    return isinstance(x, (int, Fr)) and x == 0


def mf(a):
    r = 1
    for x in a:
        r *= fact(x)
    return r


# ---------------------------------------------------------------- reference data (exact, computed once)
# reference fields u^_1, u^_2 on T^ (z,a,c) = (0,0),(1,0),(0,1): nonzero BB coefficients {(comp, alpha): value}
# (from fan.py: U(T^) = {u in lam_z lam_a P_2^2 : (u, grad r) = 0, r in P_3}, dim 3; basis printed by t3.py)
UREF = [
    {(0, (3, 1, 0)): Fr(5, 14), (0, (2, 2, 0)): Fr(11, 28), (0, (2, 1, 1)): Fr(-5, 28), (0, (1, 3, 0)): Fr(-15, 56),
     (0, (1, 2, 1)): Fr(-3, 56), (0, (1, 1, 2)): Fr(1, 28), (1, (3, 1, 0)): Fr(-1), (1, (2, 2, 0)): Fr(1)},
    {(0, (3, 1, 0)): Fr(5, 56), (0, (2, 2, 0)): Fr(25, 14), (0, (2, 1, 1)): Fr(-13, 56), (0, (1, 3, 0)): Fr(-51, 56),
     (0, (1, 2, 1)): Fr(-13, 56), (0, (1, 1, 2)): Fr(1, 14), (1, (3, 1, 0)): Fr(-1), (1, (1, 3, 0)): Fr(1)},
]


class NumExact:
    """Fraction arithmetic (no sqrt/pi: the constants C_I, C_T are not formed in exact mode)."""
    @staticmethod
    def q(p, r=1):
        return Fr(p, r)
    exact = True


class NumArb:
    exact = False

    def __init__(self):
        import flint
        self.flint = flint
        self.arb = flint.arb

    def q(self, p, r=1):
        if isinstance(p, Fr):
            return self.arb(p.numerator) / self.arb(p.denominator) / r
        return self.arb(p) / self.arb(r)


# ---------------------------------------------------------------- polynomial helpers in BB form
def int_T(area2, a, b, m, n):
    """int_T B^m_a B^n_b  (area2 = 2|T|), as a Fraction multiple of area2"""
    ab = tuple(x + y for x, y in zip(a, b))
    return Fr(fact(m) * fact(n) * mf(ab), mf(a) * mf(b) * fact(m + n + 2))


def int_e(a, b, m, n, live=(0, 1)):
    """int over the edge {third index = 0} (unit length) of B^m_a B^n_b"""
    i, j = live
    return Fr(fact(m) * fact(n) * fact(a[i] + b[i]) * fact(a[j] + b[j]),
              fact(a[i]) * fact(a[j]) * fact(b[i]) * fact(b[j]) * fact(m + n + 1))


def grad_bb(coef, g, n):
    """coef: dict alpha(|alpha|=n) -> value (scalar); g: 3 gradients (gx, gy) of barycentrics.
    returns (dx, dy): dicts beta(|beta|=n-1) -> value."""
    out = []
    for d in range(2):
        D = {}
        for be in midx(n - 1):
            s = 0
            for i in range(3):
                al = list(be); al[i] += 1; al = tuple(al)
                if al in coef:
                    s = s + coef[al] * g[i][d]
            D[be] = s * n
        out.append(D)
    return out


def geometry(x1, y1, x2, y2, N):
    one = N.q(1)
    s1 = x1 + y1
    c = (x1 / s1, one / s1)
    s2 = x2 + y2
    b = ((x2 * c[0] - c[1]) / s2, (x2 * c[1] + c[0]) / s2)
    z = (N.q(0), N.q(0)); a = (one, N.q(0))

    def bary_grads(P0, P1, P2):
        det = (P1[0] - P0[0]) * (P2[1] - P0[1]) - (P2[0] - P0[0]) * (P1[1] - P0[1])
        g = [((P1[1] - P2[1]) / det, (P2[0] - P1[0]) / det), ((P2[1] - P0[1]) / det, (P0[0] - P2[0]) / det),
             ((P0[1] - P1[1]) / det, (P1[0] - P0[0]) / det)]
        return g, det
    g1, det1 = bary_grads(z, a, c)
    g2, det2 = bary_grads(z, c, b)
    bc = b[1] / c[1]; ba = b[0] - bc * c[0]; bz = one - ba - bc
    return dict(c=c, b=b, g1=g1, g2=g2, det1=det1, det2=det2, beta=(bz, ba, bc))


def fields(x1, y1, x2, y2, N):
    G = geometry(x1, y1, x2, y2, N)
    c = G["c"]; det1 = G["det1"]
    J = ((N.q(1), c[0]), (N.q(0), c[1]))
    # Piola images of the reference fields (T1 only)
    U = []
    for ur in UREF:
        F = {1: [{al: N.q(0) for al in I4} for _ in range(2)], 2: [{al: N.q(0) for al in I4} for _ in range(2)]}
        for al in I4:
            u0 = ur.get((0, al), Fr(0)); u1 = ur.get((1, al), Fr(0))
            if u0 == 0 and u1 == 0:
                continue
            for d in range(2):
                F[1][d][al] = (J[d][0] * N.q(u0) + J[d][1] * N.q(u1)) / det1
        U.append(F)

    def curl_field(d1, d2):
        F = {}
        for t, coef, g in ((1, d1, G["g1"]), (2, d2, G["g2"])):
            dx, dy = grad_bb(coef, g, 5)
            F[t] = [{al: dy.get(al, 0) for al in I4}, {al: -dx.get(al, 0) for al in I4}]
        return F
    wA = curl_field({(2, 2, 1): N.q(1)}, {})
    wC = curl_field({(2, 1, 2): N.q(1)}, {(2, 2, 1): G["beta"][1]})
    return G, U, wA, wC


def Mt(F):
    """int_e v . t, t = (1,0): e = {lam_c = 0} of T1; int_e B^4 = 1/5"""
    t = sum(F[1][0][al] for al in I4 if al[2] == 0)
    return Fr(t, 5) if isinstance(t, (int, Fr)) else t / 5


def Dt(F, G, comp=0):
    """int_e (n . grad) v_comp,  n = (0,-1)"""
    g = G["g1"]; s = 0
    for be in I3:
        if be[2] != 0:
            continue
        for i in range(3):
            al = list(be); al[i] += 1; al = tuple(al)
            s = s + F[1][comp][al] * (-g[i][1])
    return s          # 4 * sum_i (n.g_i) v_{be+e_i} * (1/4)


def lin(F, coefs):
    """sum_j coef_j * F_j"""
    out = {1: [dict(), dict()], 2: [dict(), dict()]}
    for t in (1, 2):
        for d in range(2):
            for al in I4:
                s = 0
                for cf, Fj in coefs:
                    s = s + cf * Fj[t][d][al]
                out[t][d][al] = s
    return out


_N = NumExact


def star_fields(x1, y1, x2, y2, N):
    global _N
    _N = N
    G, U, wA, wC = fields(x1, y1, x2, y2, N)
    mA, mC = Mt(wA), Mt(wC)
    dA, dC = Dt(wA, G), Dt(wC, G)
    V = []
    for u in U:
        gA = -Mt(u) / mA                            # needs mC == 0 (checked exactly in verify)
        gC = -(Dt(u, G) + gA * dA) / dC
        V.append(lin(None, [(N.q(1), u), (gA, wA), (gC, wC)]))
    return G, V, dict(mA=mA, mC=mC, dA=dA, dC=dC)


# ---------------------------------------------------------------- quadratic forms
def gram_grad(F1, F2, G):
    """(grad v, grad w)_Q"""
    s = 0
    for t, g, det in ((1, G["g1"], G["det1"]), (2, G["g2"], G["det2"])):
        area2 = det                      # both triangles are counterclockwise
        for comp in range(2):
            d1 = grad_bb(F1[t][comp], g, 4); d2 = grad_bb(F2[t][comp], g, 4)
            for d in range(2):
                for a in I3:
                    if iz(d1[d][a]):
                        continue
                    for b in I3:
                        if iz(d2[d][b]):
                            continue
                        s = s + d1[d][a] * d2[d][b] * area2 * _N.q(int_T(1, a, b, 3, 3))
    return s


def gram_e(F1, F2):
    """(v, w)_e"""
    s = 0
    for comp in range(2):
        for a in I4:
            if a[2] or iz(F1[1][comp][a]):
                continue
            for b in I4:
                if b[2] or iz(F2[1][comp][b]):
                    continue
                s = s + F1[1][comp][a] * F2[1][comp][b] * _N.q(int_e(a, b, 4, 4))
    return s


def dn_on_e(F, G):
    """(n . grad) v on e as BB degree-3 coefficients per component (beta with beta_c = 0)"""
    g = G["g1"]; out = []
    for comp in range(2):
        dx, dy = grad_bb(F[1][comp], g, 4)
        out.append({be: -dy[be] for be in I3 if be[2] == 0})
    return out


def gram_dn(F1, F2, G):
    a1, a2 = dn_on_e(F1, G), dn_on_e(F2, G); s = 0
    for comp in range(2):
        for a, va in a1[comp].items():
            for b, vb in a2[comp].items():
                s = s + va * vb * _N.q(int_e(a, b, 3, 3))
    return s


# ---------------------------------------------------------------- the (P^T) lower bound
# reference pairing on the complement W = span{m_1 = l_z l_a^3, m_2 = l_z^2 l_a^2}  (ref_pairing.py; constant)
PW = ((Fr(0), Fr(-1, 3780)), (Fr(-1, 2520), Fr(0)))       # PW[f][j]
C3 = 10                                                   # c3const.py (exact: lambda_max = 20 on T^)


def lin_mul(p, q):
    """product of polynomials in (x, y) given as coefficient lists of x^j y^(deg-j), j = 0..deg"""
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] = out[i + j] + a * b
    return out


def frob_gram(polys):
    k = len(polys[0]) - 1
    return [[sum(p[j] * q[j] / comb(k, j) for j in range(k + 1)) for q in polys] for p in polys]


def kappa_bound(x1, y1, x2, y2, N, sqrt=None, pi=None, detail=False):
    """rigorous lower bound for kappa_T^2 (returns kappa^2 lower bound; N = NumArb for balls).
    sqrt, pi: functions/constants of the number type (needed for C_T)."""
    G, V, info = star_fields(x1, y1, x2, y2, N)
    c = G["c"]; h = c[1]
    # constants
    area1 = G["det1"] / 2
    CI2 = N.q(C3) / area1
    lc2 = c[0] * c[0] + c[1] * c[1]; lca2 = (c[0] - 1) * (c[0] - 1) + c[1] * c[1]
    R = sqrt(max_(lc2, lca2)); D = sqrt(max_(max_(lc2, lca2), N.q(1)))
    CT2 = 2 / h * (D * D / (pi * pi) + R * D / pi)
    # quadratic form
    Q = [[None, None], [None, None]]
    parts = {}
    for i in range(2):
        for j in range(i, 2):
            gg = gram_grad(V[i], V[j], G); ge = gram_e(V[i], V[j]); gd = gram_dn(V[i], V[j], G)
            Q[i][j] = Q[j][i] = 13 * gg + (3 * CI2 + CT2) * ge + 3 * CT2 * gd
            parts[(i, j)] = (gg, ge, gd)
    # pairing on W (constant) and Frobenius Gram of W
    g = G["g1"]
    lz = [g[0][1], g[0][0]]; la = [g[1][1], g[1][0]]        # coefficient lists [y-coef, x-coef] -> x^j y^(1-j) order
    m1 = lin_mul(lin_mul(lin_mul(lz, la), la), la); m2 = lin_mul(lin_mul(lin_mul(lz, lz), la), la)
    GW = frob_gram([m1, m2])
    trG = GW[0][0] + GW[1][1]; detG = GW[0][0] * GW[1][1] - GW[0][1] * GW[0][1]
    lmaxG = (trG + sqrt(max_(trG * trG - 4 * detG, N.q(0)))) / 2
    P = [[N.q(PW[f][j]) for j in range(2)] for f in range(2)]
    adj = [[Q[1][1], -Q[0][1]], [-Q[1][0], Q[0][0]]]
    # tr(P^T adj P)
    tr = 0
    for j in range(2):
        for f in range(2):
            for g_ in range(2):
                tr = tr + P[f][j] * adj[f][g_] * P[g_][j]
    detP2 = (P[0][0] * P[1][1] - P[0][1] * P[1][0]) ** 2
    k2 = detP2 / (tr * lmaxG)
    if detail:
        return k2, dict(Q=Q, GW=GW, CI2=CI2, CT2=CT2, parts=parts, info=info)
    return k2


def max_(a, b):
    """max for floats/balls (ball: the union upper envelope)"""
    try:
        import flint
        if isinstance(a, flint.arb) or isinstance(b, flint.arb):
            a = a if isinstance(a, flint.arb) else flint.arb(a)
            return a.max(b) if hasattr(a, "max") else _armax(a, b)
    except ImportError:
        pass
    return a if a >= b else b


def _armax(a, b):
    import flint
    lo = max(a.lower(), b.lower()) if hasattr(a, 'lower') else None
    return flint.arb.union(a, b)
