"""Exact (rational) verification of the local vertex-gradient claims in REPORT.md.

(A) First-order jet system at a wall vertex z of Gamma_h with m incident triangles:
    unknowns G_1..G_m (2x2, trace-free = div-free at z), constraints
    G_1 t_0 = 0, G_m t_m = 0 (v = 0 on the two Gamma_h edges), (G_{i+1}-G_i) s_i = 0 (C^0 across interior edges).
    Claim (Lemma 2): dim A_z = m + 1 - min(3, #distinct lines among the m+1 edges at z).
(B) Global polynomial realisability on the star: C^0 piecewise-P_k fields on the star, exactly div-free,
    zero on the two Gamma_h edges; rank of the map v -> (grad v|_T(z))_T.  Claim: equals dim A_z (k = 4,5,6).
(C) Same with v = 0 also on the outer boundary of the star (compactly supported in the star).
(D) gamma_k: min_{p in P_{k-1}, p(z)=0} ||1-p||_{L2(T)}^2 / |T| at a vertex z (affine invariant).
(E) Inscribed-polygon stars: distance of the constant tuple (tau (x) n, ..., tau (x) n) to A_z, as the chord
    angle 2t -> 0.  Claim (Lemma 3): = |c| exactly... O(t) for m >= 3; = sqrt(m)|c| for m <= 2 (lock).
All arithmetic over QQ (sympy DomainMatrix).
"""
from fractions import Fraction as Fr
import itertools, sys, time
from sympy import QQ, Rational
from sympy.polys.matrices import DomainMatrix

def F(e):
    return Fr(int(e.numerator), int(e.denominator))

def dm(rows, ncols):
    return DomainMatrix([[QQ(Fr(x).numerator, Fr(x).denominator) for x in r] for r in rows], (len(rows), ncols), QQ)

def nullspace(rows, n):
    if not rows:
        return [[1 if i == j else 0 for i in range(n)] for j in range(n)]
    M = dm(rows, n)
    N = M.nullspace()                 # rows of N span the kernel
    return [[F(N[i, j].element) for j in range(n)] for i in range(N.shape[0])] if N.shape[0] else []

def rank(rows, n):
    if not rows:
        return 0
    return dm(rows, n).rank()

def circ(t):                       # rational point on the unit circle, angle 2*atan(t)
    t = Fr(t); return (Fr(1) - t * t) / (1 + t * t), 2 * t / (1 + t * t)

# ---------------------------------------------------------------- (A) jet system
def jet_rows(rays):
    """rays: list of m+1 direction vectors (ccw), rays[0], rays[-1] on Gamma_h. unknowns: m blocks of 4
    (G[0][0],G[0][1],G[1][0],G[1][1]) with trace-free rows added."""
    m = len(rays) - 1; n = 4 * m; rows = []
    def G(i, a, b): return 4 * i + 2 * a + b
    for i in range(m):                                  # trace free
        r = [0] * n; r[G(i, 0, 0)] = 1; r[G(i, 1, 1)] = 1; rows.append(r)
    for (i, t) in ((0, rays[0]), (m - 1, rays[m])):     # G_i t = 0
        for a in range(2):
            r = [0] * n; r[G(i, a, 0)] = t[0]; r[G(i, a, 1)] = t[1]; rows.append(r)
    for i in range(m - 1):                              # (G_{i+1} - G_i) s_i = 0, s_i = rays[i+1]
        s = rays[i + 1]
        for a in range(2):
            r = [0] * n
            r[G(i + 1, a, 0)] += s[0]; r[G(i + 1, a, 1)] += s[1]
            r[G(i, a, 0)] -= s[0]; r[G(i, a, 1)] -= s[1]; rows.append(r)
    return rows, n

def nlines(rays):
    L = []
    for r in rays:
        if not any(r[0] * q[1] - r[1] * q[0] == 0 for q in L):
            L.append(r)
    return len(L)

# ---------------------------------------------------------------- (B),(C) polynomial star
def star_rank(z, P, k, compact=False):
    """z: vertex, P: m+1 outer points (ccw), triangles T_i = (z, P_{i-1}, P_i), i=1..m.
    Edges z-P_0 and z-P_m lie on Gamma_h."""
    m = len(P) - 1
    mons = [(a, b) for d in range(k + 1) for a in range(d + 1) for b in [d - a]]
    nm = len(mons); nb = 2 * nm; n = m * nb
    idx = {mn: j for j, mn in enumerate(mons)}
    def ev(p):                                           # monomials at p (global coords centred at z)
        x, y = p[0] - z[0], p[1] - z[1]
        return [x ** a * y ** b for (a, b) in mons]
    rows = []
    def pts(A, B):
        return [(A[0] + Fr(j, k) * (B[0] - A[0]), A[1] + Fr(j, k) * (B[1] - A[1])) for j in range(k + 1)]
    def zero_on(i, A, B):
        for p in pts(A, B):
            e = ev(p)
            for c in range(2):
                r = [0] * n
                for j in range(nm): r[i * nb + c * nm + j] = e[j]
                rows.append(r)
    zero_on(0, z, P[0]); zero_on(m - 1, z, P[m])
    for i in range(m - 1):                               # C^0 across z-P_{i+1}
        for p in pts(z, P[i + 1]):
            e = ev(p)
            for c in range(2):
                r = [0] * n
                for j in range(nm):
                    r[i * nb + c * nm + j] = e[j]; r[(i + 1) * nb + c * nm + j] = -e[j]
                rows.append(r)
    for i in range(m):                                   # div = 0 exactly
        for d in range(k):
            for a in range(d + 1):
                b = d - a; r = [0] * n
                r[i * nb + 0 * nm + idx[(a + 1, b)]] += a + 1       # d/dx of u1 x^{a+1} y^b
                r[i * nb + 1 * nm + idx[(a, b + 1)]] += b + 1       # d/dy of u2 x^a y^{b+1}
                rows.append(r)
    if compact:
        for i in range(m):
            zero_on(i, P[i], P[i + 1])
    K = nullspace(rows, n)
    # gradient at z in T_i: d u_c/dx = coeff of x (1,0), d/dy = coeff (0,1)
    gm = []
    for v in K:
        g = []
        for i in range(m):
            for c in range(2):
                g += [v[i * nb + c * nm + idx[(1, 0)]], v[i * nb + c * nm + idx[(0, 1)]]]
        gm.append(g)
    return len(K), rank(gm, 4 * m)

def generic_star(m, theta_t=Fr(1, 10)):
    """reflex fluid angle pi + theta at z = 0: rays from angle 0 ccw to pi + theta, rational directions."""
    z = (Fr(0), Fr(0))
    # rational directions at approx angles j*(pi+theta)/m via tan-half-angle, perturbed to avoid collinearity
    import math
    tot = math.pi + 2 * math.atan(float(theta_t))
    P = []
    for j in range(m + 1):
        ang = j * tot / m
        if j == m:
            P.append((Fr(-1), -theta_t))                 # angle pi + atan(theta_t)
            continue
        t = Fr(math.tan(ang / 2)).limit_denominator(97) if abs(math.cos(ang / 2)) > 1e-9 else None
        c = circ(t); P.append((c[0] * (1 + Fr(j, 7)), c[1] * (1 + Fr(j, 7))))
    return z, P

# ---------------------------------------------------------------- (D) gamma_k
def gamma_k(k):
    mons = [(a, b) for d in range(k) for a in range(d + 1) for b in [d - a]]
    from math import factorial
    def I(a, b): return Fr(factorial(a) * factorial(b), factorial(a + b + 2))
    M = DomainMatrix([[QQ(I(a1 + a2, b1 + b2).numerator, I(a1 + a2, b1 + b2).denominator)
                       for (a2, b2) in mons] for (a1, b1) in mons], (len(mons), len(mons)), QQ)
    e0 = DomainMatrix([[QQ(1 if mn == (0, 0) else 0)] for mn in mons], (len(mons), 1), QQ)
    x = M.lu_solve(e0)
    K00 = F(x[0, 0].element)                                 # e0^T M^{-1} e0
    return Fr(1) / K00 / Fr(1, 2)                         # min ||q||^2 / |T|, |T|=1/2

# ---------------------------------------------------------------- (E) inscribed stars
def inscribed_star(m, t):
    """z=(1,0) on the unit circle, polygon neighbours z1 = circ(t), z2 = circ(-t) (chord angle 2*atan t ~ 2t).
    Rays ccw: to z2, interior rays into the fluid (x > 1 side), to z1."""
    z = (Fr(1), Fr(0)); z1 = circ(t); z2 = circ(-t)
    dirs = {1: [], 2: [(Fr(1), Fr(0))], 3: [(Fr(1), Fr(-1, 2)), (Fr(1), Fr(1, 2))],
            4: [(Fr(1), Fr(-1)), (Fr(1), Fr(0)), (Fr(1), Fr(1))],
            5: [(Fr(1, 2), Fr(-1)), (Fr(1), Fr(-1, 3)), (Fr(1), Fr(1, 3)), (Fr(1, 2), Fr(1))]}[m]
    rays = [(z2[0] - z[0], z2[1] - z[1])] + dirs + [(z1[0] - z[0], z1[1] - z[1])]
    return rays

def dist_const(rays, c):
    rows, n = jet_rows(rays)
    K = nullspace(rows, n)
    m = n // 4
    target = [x for _ in range(m) for x in c]
    if not K:
        return sum(Fr(x) ** 2 for x in target), 0
    # least squares: min |B a - target|^2, B columns = K vectors
    nk = len(K)
    BtB = [[sum(K[i][r] * K[j][r] for r in range(n)) for j in range(nk)] for i in range(nk)]
    Btc = [[sum(K[i][r] * target[r] for r in range(n))] for i in range(nk)]
    a = dm(BtB, nk).lu_solve(dm(Btc, 1))
    a = [F(a[i, 0].element) for i in range(nk)]
    res = [sum(a[i] * K[i][r] for i in range(nk)) - target[r] for r in range(n)]
    return sum(Fr(x) ** 2 for x in res), nk

if __name__ == "__main__":
    parts = sys.argv[1] if len(sys.argv) > 1 else "ABDE"
    out = []
    def P(*s):
        line = " ".join(str(x) for x in s); print(line); out.append(line); sys.stdout.flush()
    if "A" in parts:
        P("(A) jet-system dimension at a corner vertex (reflex angle pi+theta), exact over QQ")
        for m in range(1, 7):
            z, Pp = generic_star(m)
            rays = [(p[0] - z[0], p[1] - z[1]) for p in Pp]
            rows, n = jet_rows(rays)
            d = n - rank(rows, n)
            P(f"  m={m}: lines={nlines(rays)}  dim A_z={d}  predicted={m + 1 - min(3, nlines(rays))}")
        # collinear case: m=2 with interior edge continuing e1 through z
        rays = [(Fr(1), Fr(0)), (Fr(-1), Fr(0)), (Fr(-1), Fr(-1, 10))]
        rows, n = jet_rows(rays); P(f"  m=2 collinear (e1 || s1, singular per (M1)): lines={nlines(rays)} dim A_z={n - rank(rows, n)}"
                                       f" predicted={3 - min(3, nlines(rays))}")
        rows_k = jet_rows(rays)[0]
        K = nullspace(rows_k, 8); P(f"    kernel basis (G_1 | G_2 entries): {K}")
    if "B" in parts:
        P("(B)/(C) C^0 P_k div-free fields on the star, zero on the two Gamma_h edges: dim of space, rank of vertex-gradient map")
        for k in (4, 5, 6):
            for m in range(1, 6 if k < 6 else 5):
                z, Pp = generic_star(m)
                t0 = time.time()
                dimB, rB = star_rank(z, Pp, k)
                dimC, rC = star_rank(z, Pp, k, compact=True)
                rays = [(p[0] - z[0], p[1] - z[1]) for p in Pp]
                P(f"  k={k} m={m}: free star: dim={dimB} grad-rank={rB} | compact in star: dim={dimC} grad-rank={rC} |"
                  f" dim A_z={m + 1 - min(3, nlines(rays))}   ({time.time() - t0:.1f}s)")
    if "D" in parts:
        P("(D) gamma_k^2 = min_{p in P_{k-1}, p(vertex)=0} ||1-p||^2_{L2(T)}/|T|  (affine invariant, exact)")
        for k in (4, 5, 6, 7, 8):
            g = gamma_k(k); P(f"  k={k}: gamma_k^2 = {g} = {float(g):.6g}")
    if "E" in parts:
        P("(E) inscribed star at z=(1,0), chord half-angle ~ t; c = tau(x)n = [[0,0],[1,0]];")
        P("    d^2 = min_{G in A_z} sum_T |G_T - c|_F^2 (exact), and d/t")
        c = [0, 0, 1, 0]
        for m in (2, 3, 4, 5):
            for t in (Fr(1, 10), Fr(1, 20), Fr(1, 40), Fr(1, 80), Fr(1, 160)):
                d2, nk = dist_const(inscribed_star(m, t), c)
                P(f"  m={m} t={t}: dim A_z={nk}  d^2={float(d2):.6e}  d={float(d2) ** .5:.6e}  d/t={float(d2) ** .5 / float(t):.6f}")
    with open(__file__.replace(".py", f"_{parts}.log"), "w") as f:
            f.write("\n".join(out) + "\n")
