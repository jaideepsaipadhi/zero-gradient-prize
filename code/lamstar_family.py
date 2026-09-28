"""
Broad penalty necessity (notes/v3/penalty_broad.tex).

PART 1 (bubble): the single-triangle witness.  Let T be a boundary triangle with edge e on Gamma_h, scaled so
that e = [(0,0),(1,0)], apex (a,b), b>0, fluid above.  With barycentrics l0 (apex), l1, l2 (base vertices),
    psi = l1^2 l2^2 (c0 l0 + c1 l1 + c2 l2),     v = curl psi  (P4, div-free, v = 0 on both interior edges),
so v (extended by 0) lies in Z_h for every k >= 4 and every mesh containing T.  The three forms
    A = int_T 1/2 Dv:Dv,   B = int_e (d_n v).v,   C = int_e |v|^2        (n = (0,-1))
are 3x3 matrices of rational functions of (a,b), derived here EXACTLY (SymPy).  lambda_T := max eig(2B-A, C).
Certificate: for all (a,b) in the region  Q = {apex angle >= 60 deg} cap {b >= s a, b >= s(1-a)} (s = 87/1000 <
tan 5 deg, i.e. both base angles >= 5 deg is inside Q) we prove   2B - A - 9 C > 0  at an explicit rational q per box,
by exact rational interval arithmetic on the polynomial  Phi_q(a,b) = 2520 b^3 q^T (2B - A - 9C) q  (b>0).

PART 2 (mesh): float lambda* for the actual boundary patches of the production meshes (paper Section 13;
chords NOT equal, h/h_Gamma ~ 3.3): one boundary triangle, one vertex star, and J full boundary layers (J <= 4),
vs the global gamma* (for N = 32, 64 the global value is taken from results/pod/thresh, not computed here).

PART 3 (levels): float values of lambda_T quoted in paper Section 8.1 (level-set minima and isosceles values
by apex angle; right triangles with the right angle at a base vertex; the tallest boundary triangles of the
production meshes).  (limitstar): exact check of the limit star S_{3/4} with lamstar_exact.run.

Usage: python lamstar_family.py [bubble|certify|mesh|limit|levels|limitstar|all]     (1 core, small memory)
"""
import os, sys, time, pickle
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import numpy as np
import sympy as sy
from fractions import Fraction as Fr

a, b = sy.symbols("a b", real=True)
HERE = os.path.dirname(os.path.abspath(__file__))


# ------------------------------------------------------------------ PART 1: exact forms
_FORMS = []


def bubble_forms():
    if _FORMS:
        return _FORMS[0]
    x, y, u, w = sy.symbols("x y u w", real=True)
    l0 = y / b
    l1 = 1 - x - (1 - a) * y / b
    l2 = x - a * y / b
    c = sy.symbols("c0:3")
    psi = l1 ** 2 * l2 ** 2 * (c[0] * l0 + c[1] * l1 + c[2] * l2)
    v = [sy.diff(psi, y), -sy.diff(psi, x)]

    def intT(f):
        g = sy.expand(f.subs({x: u + w * a, y: w * b}, simultaneous=True)) * b
        return sy.integrate(sy.integrate(g, (w, 0, 1 - u)), (u, 0, 1))
    ux, uy, vx, vy = sy.diff(v[0], x), sy.diff(v[0], y), sy.diff(v[1], x), sy.diff(v[1], y)
    Aval = intT(sy.expand(2 * ux ** 2 + (uy + vx) ** 2 + 2 * vy ** 2))          # = int 1/2 |Dv|^2
    Bval = sy.integrate(sy.expand((-sy.diff(v[0], y) * v[0] - sy.diff(v[1], y) * v[1]).subs(y, 0)), (x, 0, 1))
    Cval = sy.integrate(sy.expand((v[0] ** 2 + v[1] ** 2).subs(y, 0)), (x, 0, 1))
    H = lambda e: (sy.hessian(sy.expand(e), c) / 2).applyfunc(sy.factor)
    forms = tuple(H(e) for e in (Aval, Bval, Cval))
    _FORMS.append(forms)
    return forms


def lam_fn():
    MA, MB, MC = bubble_forms()
    f = sy.lambdify((a, b), [MA, MB, MC])
    import scipy.linalg as sl

    def lam(av, bv):
        A_, B_, C_ = [np.array(m, float) for m in f(av, bv)]
        M = 2 * B_ - A_
        ev, V = sl.eigh((M + M.T) / 2, C_)
        return ev[-1], V[:, -1]
    return lam


def apex_from_angles(al, be):
    s = np.sin(al + be)
    return np.sin(be) * np.cos(al) / s, np.sin(al) * np.sin(be) / s


def bubble():
    MA, MB, MC = bubble_forms()
    print("== single-triangle witness: exact forms (entries are rational in the apex (a,b))")
    for nm, M in (("A", MA), ("B", MB), ("C", MC)):
        for i in range(3):
            for j in range(i, 3):
                print(f"  {nm}[{i}{j}] = {M[i, j]}")
    print("  det C =", sy.factor(MC.det()))
    lam = lam_fn()
    deg = np.pi / 180
    print("\n== lambda_T(alpha, beta): rows = base angle alpha at z1, cols = base angle beta at z2 (degrees)")
    angs = [10, 20, 30, 40, 45, 50, 60, 70, 80, 90, 100, 110, 120]
    print("      " + "".join(f"{x:7d}" for x in angs))
    for al in angs:
        row = []
        for be in angs:
            if al + be >= 170:
                row.append("      -"); continue
            row.append(f"{lam(*apex_from_angles(al * deg, be * deg))[0]:7.2f}")
        print(f"{al:5d} " + "".join(row))
    # minimum over apex angle >= 60 deg, base angles >= 5 deg (float scan)
    best = (np.inf, None)
    for al in np.linspace(5, 115, 221):
        for be in np.linspace(5, 115, 221):
            if al + be <= 120 + 1e-9:
                l = lam(*apex_from_angles(al * deg, be * deg))[0]
                if l < best[0]:
                    best = (l, (al, be))
    print(f"\nfloat scan: min lambda_T over {{apex >= 60, base angles >= 5}} = {best[0]:.6f} at (alpha,beta) = {best[1]}")
    leq = lam(0.5, np.sqrt(3) / 2)[0]
    print(f"equilateral: lambda_T = {leq:.8f}")
    # the production-mesh boundary triangles are right triangles (beta = 90) with alpha = atan(radial/tangential spacing)
    print("right triangles beta=90 (production meshes (paper Section 13) have alpha in ~[43.7, 58.2] deg):")
    for al in (30, 35, 40, 43.7, 45, 50, 55, 58.2, 60):
        print(f"   alpha={al:5.1f}  lambda_T={lam(*apex_from_angles(al * deg, 90 * deg))[0]:8.4f}")
    # independent float check of the exact forms by Gauss quadrature of the same witness
    check_forms()


def check_forms(av=0.3, bv=0.7):
    MA, MB, MC = bubble_forms()
    ex = [np.array(M.subs({a: av, b: bv}), float) for M in (MA, MB, MC)]
    g, wg = np.polynomial.legendre.leggauss(12); g = (g + 1) / 2; wg = wg / 2
    Aq = np.zeros((3, 3)); Bq = np.zeros((3, 3)); Cq = np.zeros((3, 3))
    hstep = 1e-5

    def lam3(X, Y):
        l0 = Y / bv; l1 = 1 - X - (1 - av) * Y / bv; l2 = X - av * Y / bv
        return l0, l1, l2

    def vel(X, Y, i):
        def psi(X, Y):
            l = lam3(X, Y); return l[1] ** 2 * l[2] ** 2 * l[i]
        return np.array([(psi(X, Y + hstep) - psi(X, Y - hstep)) / (2 * hstep),
                         -(psi(X + hstep, Y) - psi(X - hstep, Y)) / (2 * hstep)])

    def grad(X, Y, i):
        return np.stack([(vel(X + hstep, Y, i) - vel(X - hstep, Y, i)) / (2 * hstep),
                         (vel(X, Y + hstep, i) - vel(X, Y - hstep, i)) / (2 * hstep)], 1)   # [comp, dir]
    for U_, wu in zip(g, wg):
        for W_, ww in zip(g, wg):
            uu, w2 = U_, W_ * (1 - U_); X = uu + w2 * av; Y = w2 * bv; wt = wu * ww * (1 - U_) * bv
            Gs = [grad(X, Y, i) for i in range(3)]; Ds = [Gm + Gm.T for Gm in Gs]
            for i in range(3):
                for j in range(3):
                    Aq[i, j] += 0.5 * np.sum(Ds[i] * Ds[j]) * wt
    for X, wx in zip(g, wg):
        vs = [vel(X, 0.0, i) for i in range(3)]
        dn = [-(vel(X, hstep, i) - vel(X, -hstep, i)) / (2 * hstep) for i in range(3)]
        for i in range(3):
            for j in range(3):
                Bq[i, j] += 0.5 * (dn[j] @ vs[i] + dn[i] @ vs[j]) * wx
                Cq[i, j] += vs[i] @ vs[j] * wx
    err = max(np.abs(Aq - ex[0]).max() / np.abs(ex[0]).max(), np.abs(Bq - (ex[1] + ex[1].T) / 2).max() / np.abs(ex[1]).max(),
              np.abs(Cq - ex[2]).max() / np.abs(ex[2]).max())
    print(f"independent quadrature/finite-difference check of the exact forms at (a,b)=({av},{bv}): rel. diff {err:.1e}")


# ------------------------------------------------------------------ exact interval certificate
def _ipow(lo, hi, k):
    if k == 0:
        return Fr(1), Fr(1)
    cands = [lo ** k, hi ** k]
    if lo < 0 < hi and k % 2 == 0:
        return Fr(0), max(cands)
    return min(cands), max(cands)


def poly_lower(P, box):
    """rigorous lower bound of the polynomial P (dict {(i,j): Fraction}) over box [a0,a1]x[b0,b1]"""
    (a0, a1), (b0, b1) = box
    lo = Fr(0)
    for (i, j), cf in P.items():
        al, ah = _ipow(a0, a1, i); bl, bh = _ipow(b0, b1, j)
        prods = [al * bl, al * bh, ah * bl, ah * bh]
        lo += cf * min(prods) if cf > 0 else cf * max(prods)
    return lo


SQ3_LO, SQ3_HI = Fr(173205, 100000), Fr(173206, 100000)      # rigorous enclosure of sqrt(3)
S_SLOPE = Fr(87, 1000)                                          # < tan(5 deg) = 0.087488...


def disk_params(apex_deg):
    """outward-rounded rational enclosures of the disk {apex angle >= apex_deg}: (R^2 upper, centre-height lo, hi)"""
    import mpmath as mp
    mp.mp.dps = 50
    ph = mp.radians(apex_deg)
    R2 = 1 / (4 * mp.sin(ph) ** 2); k = mp.cot(ph) / 2
    den = 10 ** 12
    return (Fr(int(mp.ceil(R2 * den)), den), Fr(int(mp.floor(k * den)), den), Fr(int(mp.ceil(k * den)), den))


def box_outside(box, disk=None):
    """conservative: True only if the box provably misses Q"""
    (a0, a1), (b0, b1) = box
    # b >= s*a and b >= s*(1-a):  fails on the whole box if b1 < s*max over box ... use min over box
    if b1 < S_SLOPE * a0 or b1 < S_SLOPE * (1 - a1):
        return True
    # apex angle >= 60  <=>  (a-1/2)^2 + (b - 1/(2 sqrt3))^2 <= 1/3 ; outside if min distance^2 > 1/3
    if disk is None:
        R2, cy_lo, cy_hi = Fr(1, 3), 1 / (2 * SQ3_HI), 1 / (2 * SQ3_LO)
    else:
        R2, cy_lo, cy_hi = disk
    dx = max(Fr(0), a0 - Fr(1, 2), Fr(1, 2) - a1)
    dy = max(Fr(0), b0 - cy_hi, cy_lo - b1)
    return dx * dx + dy * dy > R2


def certify(c=Fr(9), apex_deg=60, max_depth=16):
    MA, MB, MC = bubble_forms()
    lam = lam_fn()
    M = (2 * (MB + MB.T) / 2 - MA - sy.Rational(c.numerator, c.denominator) * MC) * 2520 * b ** 3
    M = M.applyfunc(sy.cancel)
    for e in M:
        assert not sy.denom(sy.together(e)).free_symbols, "clearing denominators failed"
    Mpoly = [[sy.Poly(M[i, j], a, b) for j in range(3)] for i in range(3)]
    t0 = time.time()
    disk = None if apex_deg == 60 else disk_params(apex_deg)
    R = float(disk[0]) ** 0.5 if disk else 3 ** -0.5
    k = float(disk[2]) if disk else 0.2887
    a_lo, a_hi = Fr(0.5 - R - 0.01).limit_denominator(100), Fr(0.5 + R + 0.01).limit_denominator(100)
    stack = [((a_lo, a_hi), (Fr(4, 100), Fr(k + R + 0.01).limit_denominator(100)), 0)]   # contains Q (b >= s/2 > 0.04)
    assert a_lo < 0.5 - R and a_hi > 0.5 + R and float(stack[0][1][1]) > k + R
    nbox = nskip = 0; maxd = 0
    while stack:
        A_, B_, d = stack.pop()
        box = (A_, B_)
        if box_outside(box, disk):
            nskip += 1; continue
        am, bm = float(sum(A_) / 2), float(sum(B_) / 2)
        _, vec = lam(am, bm)
        vec = vec / np.abs(vec).max()
        q = [Fr(x).limit_denominator(1000) for x in vec]
        P = {}
        for i in range(3):
            for j in range(3):
                for (ei, ej), cf in Mpoly[i][j].terms():
                    P[(ei, ej)] = P.get((ei, ej), Fr(0)) + Fr(int(sy.numer(cf)), int(sy.denom(cf))) * q[i] * q[j]
        if poly_lower(P, box) > 0:
            nbox += 1; maxd = max(maxd, d); continue
        if d >= max_depth:
            raise RuntimeError(f"certificate failed on box {[(float(x), float(y)) for x, y in box]}")
        am_, bm_ = sum(A_) / 2, sum(B_) / 2
        for aa in ((A_[0], am_), (am_, A_[1])):
            for bb in ((B_[0], bm_), (bm_, B_[1])):
                stack.append((aa, bb, d + 1))
    print(f"== EXACT certificate [apex >= {apex_deg} deg, base angles >= 5 deg]: 2B - A - {c}*C > 0 at an explicit rational "
          f"witness on every box meeting Q: "
          f"{nbox} certified boxes, {nskip} boxes provably outside Q, max depth {maxd}  ({time.time()-t0:.0f}s)")
    return nbox


# ------------------------------------------------------------------ PART 2: production-mesh patches
def mesh_patches(Ns=(8, 16, 32, 64)):
    sys.path.insert(0, HERE)
    import svn, penalty
    print("== production meshes (paper Section 13): local thresholds lambda* = (mu/h)* hG  of fields supported on boundary patches")
    print("   (lambda* of a patch is a LOWER bound for the global gamma* = mu* hG / h)")
    print("   tri/star = max over all boundary triangles / single-vertex stars; winW = first-layer window of W columns")
    print("   around theta = 0; J = first J full layers (J = m = N/4 is the whole mesh only for N <= 16)")
    for N in Ns:
        verts, tris, circ, outer = svn.make_mesh(N)
        m = N // 4
        res = {}
        # layer index of each triangle: tris are generated layer by layer, 2N per layer
        layer = np.repeat(np.arange(m), 2 * N)
        s4 = N // 4
        quad = np.tile(np.repeat(np.arange(N), 2), m)          # quad (column) index of each triangle
        patches = []
        # single boundary triangles and single-vertex stars, maximised over one symmetry octant of vertices
        for i in range(0, s4 // 2 + 1):
            vi = i * 1                                        # circle vertex i (vertex ids 0..N-1 are layer 0)
            st = np.array([vi in t for t in tris]) & (layer == 0)
            patches.append((f"star{i}", st))
            for tix in np.where(st)[0]:
                t1 = np.zeros(len(tris), bool); t1[tix] = True
                if sum(1 for vv in tris[tix] if vv in circ) == 2:
                    patches.append((f"tri{i}", t1))
        # first-layer windows of W columns centred at the side midpoint (vertex s4//2, direction theta = 0)
        c = s4 // 2
        for W in (2, 4, 8, 16, 32, 64):
            if W < N:
                cols = (np.arange(c - W // 2, c + W // 2) % N)
                patches.append((f"win{W}", (layer == 0) & np.isin(quad, cols)))
        for J in (1, 2, 3, 4):
            if J <= m:
                patches.append((f"J={J}", layer < J))
        for name, mask in patches:
            sub = tris[mask]
            used = np.unique(sub)
            remap = -np.ones(len(verts), int); remap[used] = np.arange(len(used))
            S = svn.Space(verts[used], remap[sub], {int(remap[i]) for i in circ if remap[i] >= 0},
                          {int(remap[i]) for i in outer if remap[i] >= 0})
            if name != "global":
                S.dnodes = patch_dirichlet(S)
            mu = penalty.threshold(S, lo=1e-3, hi=1e6, tol=1e-6)
            lam = mu / S.h * S.hGamma if mu > 1.001e-3 else 0.0     # threshold at the floor: no negative direction
            key = name.rstrip("0123456789") if name.startswith(("star", "tri")) else name
            res[key] = max(res.get(key, -1.0), lam)
        print(f"N={N:3d}  " + "  ".join(f"{k}: {v:8.4f}" for k, v in res.items()), flush=True)


def patch_dirichlet(S):
    """nodes on the boundary of the patch that are NOT on Gamma_h (zero Dirichlet there)"""
    from collections import Counter
    cnt = Counter()
    for t in S.tris:
        for (i, j) in ((0, 1), (1, 2), (2, 0)):
            cnt[tuple(sorted((t[i], t[j])))] += 1
    circ_nodes = {tuple(sorted((S.tris[t][u], S.tris[t][w]))) for (t, u, w) in S.bedges}
    bdry = {e for e, k in cnt.items() if k == 1 and e not in circ_nodes}
    import svn
    nodes = set()
    for t, tri in enumerate(S.tris):
        for k, (l1, l2, l0) in enumerate(svn.LBARY):
            lamv = {0: l0, 1: l1, 2: l2}
            nz = [vv for vv in range(3) if lamv[vv] > 0]
            if len(nz) <= 2:
                vs = tuple(sorted(tri[vv] for vv in nz))
                if len(vs) == 2 and vs in bdry:
                    nodes.add(S.ids[t, k])
                if len(vs) == 1 and any(vs[0] in e for e in bdry):
                    nodes.add(S.ids[t, k])
    return np.array(sorted(nodes))




# ------------------------------------------------------------------ PART 3: exact / float patches of arbitrary shape
from math import comb, factorial

MON4 = [(i, j) for i in range(5) for j in range(5 - i)]


def _seg_poly(P, Q):
    """for each monomial x^i y^j (deg<=4) the coefficient list (in t) of (P + t(Q-P))^(i,j), exact"""
    dx, dy = Q[0] - P[0], Q[1] - P[1]
    out = {}
    for (i, j) in MON4:
        cx = [Fr(comb(i, k)) * P[0] ** (i - k) * dx ** k for k in range(i + 1)]
        cy = [Fr(comb(j, k)) * P[1] ** (j - k) * dy ** k for k in range(j + 1)]
        c = [Fr(0)] * (i + j + 1)
        for p, u in enumerate(cx):
            for q, w in enumerate(cy):
                c[p + q] += u * w
        out[(i, j)] = c + [Fr(0)] * (5 - len(c))
    return out


def _tri_moments(T, maxdeg=8):
    """exact int_T x^p y^q for p+q<=maxdeg via barycentric integration"""
    (x0, y0), (x1, y1), (x2, y2) = T
    area2 = abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))
    mom = {}
    # x^p y^q = (sum l_i x_i)^p (sum l_i y_i)^q ; int l0^a l1^b l2^c = area2 a! b! c! / (a+b+c+2)!
    def expand(coef, n):
        res = {}
        for a0 in range(n + 1):
            for a1 in range(n + 1 - a0):
                a2 = n - a0 - a1
                m = Fr(factorial(n), factorial(a0) * factorial(a1) * factorial(a2))
                res[(a0, a1, a2)] = m * coef[0] ** a0 * coef[1] ** a1 * coef[2] ** a2
        return res
    for p in range(maxdeg + 1):
        ex = expand((x0, x1, x2), p)
        for q in range(maxdeg + 1 - p):
            ey = expand((y0, y1, y2), q)
            s = Fr(0)
            for (a0, a1, a2), u in ex.items():
                for (b0, b1, b2), w in ey.items():
                    c0, c1, c2 = a0 + b0, a1 + b1, a2 + b2
                    s += u * w * Fr(factorial(c0) * factorial(c1) * factorial(c2), factorial(c0 + c1 + c2 + 2))
            mom[(p, q)] = s * area2
    return mom


def patch_forms(tris, free_segs):
    """tris: list of triangles (rational vertex tuples). Edges on the patch boundary that are not in free_segs are
    Dirichlet (zero). free_segs: list of (P, Q) on the line y=0, fluid above (n = (0,-1)).
    Returns exact constraint rows, forms A, B, C (dicts of Fractions over the 30*len(tris) unknowns)."""
    nt = len(tris); nv = 30 * nt
    idx = lambda t, c, k: 30 * t + 15 * c + k
    rows = []
    edges = {}
    for t, T in enumerate(tris):
        for (i, j) in ((0, 1), (1, 2), (2, 0)):
            key = tuple(sorted((T[i], T[j])))
            edges.setdefault(key, []).append(t)
    free = {tuple(sorted(s)) for s in free_segs}
    for (P, Q), ts in edges.items():
        sp_ = _seg_poly(P, Q)
        if len(ts) == 2:
            for c in range(2):
                for d in range(5):
                    r = {}
                    for k, m in enumerate(MON4):
                        if sp_[m][d]:
                            r[idx(ts[0], c, k)] = sp_[m][d]; r[idx(ts[1], c, k)] = -sp_[m][d]
                    if r: rows.append(r)
        elif (P, Q) not in free:
            for c in range(2):
                for d in range(5):
                    r = {idx(ts[0], c, k): sp_[m][d] for k, m in enumerate(MON4) if sp_[m][d]}
                    if r: rows.append(r)
    for t in range(nt):                               # div v = 0 : coefficient of x^p y^q, p+q<=3
        for (p, q) in [(p, q) for p in range(4) for q in range(4 - p)]:
            r = {}
            for k, (i, j) in enumerate(MON4):
                if i >= 1 and (i - 1, j) == (p, q): r[idx(t, 0, k)] = r.get(idx(t, 0, k), 0) + Fr(i)
                if j >= 1 and (i, j - 1) == (p, q): r[idx(t, 1, k)] = r.get(idx(t, 1, k), 0) + Fr(j)
            if r: rows.append(r)
    # forms
    A = {}; B = {}; C = {}
    def add(D, i, j, v):
        if v: D[(i, j)] = D.get((i, j), Fr(0)) + v
    for t, T in enumerate(tris):
        mom = _tri_moments(T, 6)
        # derivative monomials: d/dx x^i y^j = i x^(i-1) y^j
        dx = [((i - 1, j), Fr(i)) if i else (None, 0) for (i, j) in MON4]
        dy = [((i, j - 1), Fr(j)) if j else (None, 0) for (i, j) in MON4]
        def I(m1, c1, m2, c2):
            if m1 is None or m2 is None: return Fr(0)
            return c1 * c2 * mom[(m1[0] + m2[0], m1[1] + m2[1])]
        # 1/2|Dv|^2 = 2 u_x^2 + (u_y + w_x)^2 + 2 w_y^2 ... times 1/2 * 2 -> the same density as lamstar_exact / 2 * 2
        for k in range(15):
            for l in range(15):
                ux = I(dx[k][0], dx[k][1], dx[l][0], dx[l][1]); uy = I(dy[k][0], dy[k][1], dy[l][0], dy[l][1])
                uxy = I(dy[k][0], dy[k][1], dx[l][0], dx[l][1])        # d_y(phi_k) d_x(phi_l)
                add(A, idx(t, 0, k), idx(t, 0, l), 2 * ux + uy)
                add(A, idx(t, 1, k), idx(t, 1, l), 2 * uy + ux)
                add(A, idx(t, 0, k), idx(t, 1, l), uxy)                # (u1_y)(u2_x) cross term
                add(A, idx(t, 1, l), idx(t, 0, k), uxy)
    for (P, Q) in free_segs:
        key = tuple(sorted((P, Q)))
        t = edges[key][0]
        lo, hi = sorted([P[0], Q[0]])
        # on y=0: v = sum c x^i (j=0 terms); d_n v = -d_y v = -(coef of j=1 terms) x^i
        for c in range(2):
            for k, (i, j) in enumerate(MON4):
                for l, (i2, j2) in enumerate(MON4):
                    if j == 0 and j2 == 0:
                        n_ = i + i2 + 1; add(C, idx(t, c, k), idx(t, c, l), Fr(hi ** n_ - lo ** n_, n_))
                    if j == 0 and j2 == 1:        # <d_n u, v>: v-term k, trial-derivative term l
                        n_ = i + i2 + 1; v = -Fr(hi ** n_ - lo ** n_, n_) / 2
                        add(B, idx(t, c, k), idx(t, c, l), v); add(B, idx(t, c, l), idx(t, c, k), v)
    return nv, rows, A, B, C


def patch_lambda(tris, free_segs, lam0=None, exact=False):
    from sympy.polys.matrices import DomainMatrix
    from sympy import QQ
    nv, rows, A, B, C = patch_forms(tris, free_segs)
    if not exact:
        R = np.zeros((len(rows), nv))
        for r, d in enumerate(rows):
            for k, v in d.items(): R[r, k] = float(v)
        U_, S_, Vt = np.linalg.svd(R)
        rk = int((S_ > 1e-11 * S_[0]).sum()); K = Vt[rk:].T
    else:
        dm = DomainMatrix([[QQ(r.get(k, 0).numerator, r.get(k, 0).denominator) if k in r else QQ(0) for k in range(nv)]
                           for r in rows], (len(rows), nv), QQ)
        ns = dm.nullspace().to_Matrix()             # rows = basis vectors
        Kx = [[Fr(int(sy.numer(e)), int(sy.denom(e))) for e in ns.row(i)] for i in range(ns.rows)]
        K = np.array([[float(e) for e in row] for row in Kx]).T
    def dense(D):
        M = np.zeros((nv, nv))
        for (i, j), v in D.items(): M[i, j] += float(v)
        return M
    Af, Bf, Cf = [K.T @ dense(D) @ K for D in (A, B, C)]
    Mf = 2 * Bf - Af
    ew, Q = np.linalg.eigh(Cf); keep = ew > 1e-12 * ew.max()
    Rr, Ns = Q[:, keep], Q[:, ~keep]
    Myy = Rr.T @ Mf @ Rr; Cyy = Rr.T @ Cf @ Rr
    if Ns.shape[1]:
        Mkk = Ns.T @ Mf @ Ns; Mky = Ns.T @ Mf @ Rr
        Seff = Myy - Mky.T @ np.linalg.solve(Mkk, Mky)
    else:
        Seff = Myy
    Lc = np.linalg.cholesky(Cyy); Li = np.linalg.inv(Lc)
    ev, EV = np.linalg.eigh(Li @ Seff @ Li.T)
    out = dict(dimZ=K.shape[1], lam_float=ev[-1])
    if exact:
        yv = Li.T @ EV[:, -1]; zf = Rr @ yv
        if Ns.shape[1]:
            zf = zf + Ns @ (-np.linalg.solve(Mkk, Mky @ yv))
        zr = [Fr(v).limit_denominator(10 ** 7) for v in zf / np.abs(zf).max()]
        vec = [sum(zr[b] * Kx[b][i] for b in range(len(zr))) for i in range(nv)]    # exact coefficient vector
        q = lambda D: sum(v * vec[i] * vec[j] for (i, j), v in D.items())
        # exact re-check of the constraints
        assert all(sum(v * vec[k] for k, v in r.items()) == 0 for r in rows)
        qa, qb, qc = q(A), q(B), q(C)
        out.update(rayleigh_exact=float((2 * qb - qa) / qc), certified=(2 * qb - qa - lam0 * qc) > 0, lam0=lam0)
    return out


def limit_geometry(W, layers=1, H=Fr(3, 4)):
    """production-mesh boundary layer at theta=0 in the limit N -> inf, scaled by the chord: columns [i,i+1] x [0,H],
    i = -W/2..W/2-1, quads split on the forward diagonal (as in svn.make_mesh)."""
    tris = []
    for jl in range(layers):
        for i in range(-W // 2, W // 2):
            p00, p10 = (Fr(i), H * jl), (Fr(i + 1), H * jl)
            p01, p11 = (Fr(i), H * (jl + 1)), (Fr(i + 1), H * (jl + 1))
            tris += [(p00, p10, p11), (p00, p11, p01)]
    free = [((Fr(i), Fr(0)), (Fr(i + 1), Fr(0))) for i in range(-W // 2, W // 2)]
    return tris, free


def limit_patches():
    print("== production-mesh limit geometry (theta = 0, N -> inf, chord-scaled; layer height 3/4 = (L-1)/2)")
    print("   (float: monomials in global coordinates, so W <= 32; at W = 64 the SVD rank decision fails)")
    t0 = time.time()
    star = [((Fr(-1), Fr(0)), (Fr(0), Fr(0)), (Fr(0), Fr(3, 4))),
            ((Fr(0), Fr(0)), (Fr(1), Fr(0)), (Fr(1), Fr(3, 4))),
            ((Fr(0), Fr(0)), (Fr(1), Fr(3, 4)), (Fr(0), Fr(3, 4)))]
    r = patch_lambda(star, [((Fr(-1), Fr(0)), (Fr(0), Fr(0))), ((Fr(0), Fr(0)), (Fr(1), Fr(0)))])
    print(f"  single-vertex star (cross-check of lamstar_exact.py): dimZ={r['dimZ']}  lambda*={r['lam_float']:.6f}")
    for layers in (1, 2, 3):
        for W in (2, 4, 8, 16, 32):
            if (layers == 2 and W > 32) or (layers == 3 and W > 16):
                continue
            tris, free = limit_geometry(W, layers)
            r = patch_lambda(tris, free)
            print(f"  layers={layers} window W={W:2d}: dimZ={r['dimZ']:4d}  lambda*={r['lam_float']:.6f}   ({time.time()-t0:.0f}s)",
                  flush=True)
    r = patch_lambda(star, [((Fr(-1), Fr(0)), (Fr(0), Fr(0))), ((Fr(0), Fr(0)), (Fr(1), Fr(0)))], lam0=Fr(17), exact=True)
    print(f"  EXACT star S_3/4: dimZ={r['dimZ']}, exact Rayleigh quotient of the rational witness={r['rayleigh_exact']:.6f};"
          f"  2B - A - 17 C > 0 exactly: {r['certified']}", flush=True)
    for W, lam0 in ((4, Fr(20)),):
        tris, free = limit_geometry(W, 1)
        r = patch_lambda(tris, free, lam0=lam0, exact=True)
        print(f"  EXACT window W={W}: dimZ={r['dimZ']}, float lambda*={r['lam_float']:.6f}, exact Rayleigh quotient of the "
              f"rational witness={r['rayleigh_exact']:.6f};  2B - A - {lam0} C > 0 exactly: {r['certified']}  "
              f"({time.time()-t0:.0f}s)", flush=True)


def levels():
    """Float values of lambda_T quoted in Section 8.1 (numerical observations, not certified)."""
    lam = lam_fn(); deg = np.pi / 180
    print("== lambda_T on the level sets {apex angle = phi, base angles >= 5 deg}: min over a scan, and isosceles value")
    for phi in (60, 50, 45, 40, 35, 30):
        best = (np.inf, None)
        for al in np.linspace(5, 175 - phi - 5, 2001):
            be = 180 - phi - al
            if be < 5:
                continue
            l = lam(*apex_from_angles(al * deg, be * deg))[0]
            if l < best[0]:
                best = (l, al)
        iso = lam(*apex_from_angles((90 - phi / 2) * deg, (90 - phi / 2) * deg))[0]
        print(f"  apex {phi:3d}:  min lambda_T = {best[0]:9.5f} at base angle {best[1]:7.3f}   isosceles: {iso:9.5f}")
    print("== right angle at a base vertex (beta = 90), apex angle = 90 - alpha")
    for ap in (30, 25, 20, 15, 10):
        al = 90 - ap
        print(f"  apex {ap:3d}:  lambda_T = {lam(*apex_from_angles(al * deg, 90 * deg))[0]:9.4f}")
    import svn
    print("== production meshes: boundary triangle with the smallest lambda_T (edge scaled to [(0,0),(1,0)])")
    for N in (16, 32, 64, 256):
        v, tr, circ, outer = svn.make_mesh(N)
        best = (np.inf, None)
        for T in tr:
            on = [i for i in T if i in circ]
            if len(on) != 2:
                continue
            P, Q = v[on[0]], v[on[1]]; X = v[[i for i in T if i not in on][0]]
            ex = Q - P; L = np.linalg.norm(ex); ex = ex / L; ey = np.array([-ex[1], ex[0]])
            y = X - P; a_, b_ = y @ ex / L, y @ ey / L
            if b_ < 0:
                ex = -ex; ey = np.array([-ex[1], ex[0]]); P = Q; y = X - P; a_, b_ = y @ ex / L, y @ ey / L
            l = lam(a_, b_)[0]
            if l < best[0]:
                best = (l, (a_, b_))
        print(f"  N={N:4d}: min lambda_T = {best[0]:8.4f} at apex ({best[1][0]:.3f}, {best[1][1]:.3f})")


def limitstar():
    import sympy
    import lamstar_exact
    print("limit star S_{3/4} at theta=0 (N->inf): rim (1,0),(1,3/4),(0,3/4),(-1,0); boundary y=0", flush=True)
    t0 = time.time()
    lamstar_exact.run([(1, 0), (1, sympy.Rational(3, 4)), (0, sympy.Rational(3, 4)), (-1, 0)], sympy.Integer(17))
    print("time", time.time() - t0)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode in ("bubble", "all"):
        bubble()
    if mode in ("certify", "all"):
        certify(Fr(9), 60)
        certify(Fr(93, 10), 60)
        certify(Fr(4), 45)
        certify(Fr(1), 35)
        print("negative controls (must FAIL: the bound exceeds the true minimum of lambda_T on Q):")
        for c_, ap in ((Fr(94, 10), 60), (Fr(45, 10), 45)):
            try:
                certify(c_, ap, max_depth=9); print("  UNEXPECTED success", c_, ap)
            except RuntimeError as e:
                print(f"  c={c_}, apex>={ap}: failed as expected: {e}", flush=True)
    if mode in ("mesh", "all"):
        mesh_patches()
    if mode in ("limit", "all"):
        limit_patches()
    if mode in ("levels", "all"):
        levels()
    if mode == "limitstar":
        limitstar()
