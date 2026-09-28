"""
Exact rational (fractions.Fraction) computation of the divergence-free space of a patch and of the
forms A, B, C on it; generic degree k.  Used for
  * NEGATIVE certificates:  A - 2 sym(B) + lam1*C  positive definite on the WHOLE divergence-free space
    (exact LDL^T, all pivots > 0)  ==>  (2B - A) < lam1*C for every nonzero field, i.e. lambda*_patch < lam1
    (and, since C >= 0, N_h(v,v) = A - 2B + (mu*l/h) C > 0 on the patch for every mu*l/h >= lam1).
    With lam1 <= 0 this says: no field supported in the patch is a witness, for any penalty mu >= 0.
  * POSITIVE certificates: an exact rational v in the space with 2B - A - lam0*C > 0.

Patch: rational vertices, triangles; Nitsche edges lie on the line y = 0 with outward normal n = (0,-1)
(fluid above); every other boundary edge of the patch carries v = 0.  Per triangle, each velocity
component is a polynomial of degree <= k in global (x, y) with unknown coefficients.  Constraints:
continuity across interior edges and vanishing on Dirichlet edges (as polynomial identities in the edge
parameter), and div v = 0 (as a polynomial identity).  All arithmetic exact.
"""
from fractions import Fraction as Fr
from math import comb, factorial
import itertools, sys, time
import numpy as np


def monos(k):
    return [(a, b) for a in range(k + 1) for b in range(k + 1 - a)]


# ------------------------------------------------------------------ exact polynomial helpers
def poly_on_edge(P, Q, k):
    """matrix E (k+1 x nm): coefficients in t (t^0..t^k) of x^a y^b along P + t (Q - P)."""
    mons = monos(k)
    dx, dy = Q[0] - P[0], Q[1] - P[1]
    rows = [[Fr(0)] * len(mons) for _ in range(k + 1)]
    for j, (a, b) in enumerate(mons):
        # (P0 + t dx)^a (P1 + t dy)^b
        for i in range(a + 1):
            ca = comb(a, i) * P[0] ** (a - i) * dx ** i
            if ca == 0:
                continue
            for l in range(b + 1):
                cb = comb(b, l) * P[1] ** (b - l) * dy ** l
                if cb:
                    rows[i + l][j] += ca * cb
    return rows


def tri_moments(T, maxdeg):
    """exact int_T x^a y^b for a+b <= maxdeg, via barycentric expansion."""
    (x0, y0), (x1, y1), (x2, y2) = T
    area2 = abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))   # = 2|T|
    # int_T l0^i l1^j l2^m = 2|T| i! j! m! / (i+j+m+2)!
    def bary_int(i, j, m):
        return Fr(area2 * factorial(i) * factorial(j) * factorial(m), factorial(i + j + m + 2))
    # power expansions: x^a = sum over multi-index of multinomial * x0^i x1^j x2^m l0^i l1^j l2^m
    def expand(c0, c1, c2, a):
        out = {}
        for i in range(a + 1):
            for j in range(a + 1 - i):
                m = a - i - j
                coef = factorial(a) // (factorial(i) * factorial(j) * factorial(m)) * c0 ** i * c1 ** j * c2 ** m
                if coef:
                    out[(i, j, m)] = out.get((i, j, m), 0) + coef
        return out
    xe = {a: expand(x0, x1, x2, a) for a in range(maxdeg + 1)}
    ye = {b: expand(y0, y1, y2, b) for b in range(maxdeg + 1)}
    mom = {}
    for a in range(maxdeg + 1):
        for b in range(maxdeg + 1 - a):
            s = Fr(0)
            for (i1, j1, m1), c1 in xe[a].items():
                for (i2, j2, m2), c2 in ye[b].items():
                    s += c1 * c2 * bary_int(i1 + i2, j1 + j2, m1 + m2)
            mom[(a, b)] = s
    return mom


def rref_nullspace(rows, n):
    """exact nullspace basis (list of vectors of length n) of the matrix given by rows."""
    M = [list(r) for r in rows if any(r)]
    piv_cols = []
    r = 0
    for c in range(n):
        p = None
        for i in range(r, len(M)):
            if M[i][c] != 0:
                p = i; break
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        inv = 1 / M[r][c]
        M[r] = [v * inv for v in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                Mi, Mr = M[i], M[r]
                M[i] = [Mi[j] - f * Mr[j] for j in range(n)]
        piv_cols.append(c)
        r += 1
        if r == len(M):
            break
    free = [c for c in range(n) if c not in set(piv_cols)]
    basis = []
    for fcol in free:
        v = [Fr(0)] * n
        v[fcol] = Fr(1)
        for i, pc in enumerate(piv_cols):
            v[pc] = -M[i][fcol]
        basis.append(v)
    return basis, len(piv_cols)


# ------------------------------------------------------------------ patch assembly
def assemble(P, tris, nit_edges, k):
    """P: list of rational points; tris: index triples; nit_edges: set of vertex pairs on y=0."""
    P = [(Fr(p[0]), Fr(p[1])) for p in P]
    mons = monos(k); nm = len(mons)
    midx = {mb: i for i, mb in enumerate(mons)}
    nT = len(tris)
    nvar = nT * 2 * nm
    var = lambda t, c, j: (t * 2 + c) * nm + j
    ecount = {}
    for ti, T in enumerate(tris):
        for e in ((T[0], T[1]), (T[1], T[2]), (T[0], T[2])):
            ecount.setdefault(tuple(sorted(e)), []).append(ti)
    nit = {tuple(sorted(e)) for e in nit_edges}
    rows = []
    for e, ts in ecount.items():
        E = poly_on_edge(P[e[0]], P[e[1]], k)
        if len(ts) == 2:
            t1, t2 = ts
            for c in range(2):
                for r in E:
                    row = [Fr(0)] * nvar
                    for j in range(nm):
                        if r[j]:
                            row[var(t1, c, j)] += r[j]; row[var(t2, c, j)] -= r[j]
                    rows.append(row)
        elif e not in nit:
            t1 = ts[0]
            for c in range(2):
                for r in E:
                    row = [Fr(0)] * nvar
                    for j in range(nm):
                        if r[j]:
                            row[var(t1, c, j)] += r[j]
                    rows.append(row)
        else:
            assert P[e[0]][1] == 0 and P[e[1]][1] == 0, "Nitsche edges must lie on y=0"
    # divergence: d/dx u + d/dy w = 0 coefficientwise
    for t in range(nT):
        for (a, b) in monos(k - 1):
            row = [Fr(0)] * nvar
            if (a + 1, b) in midx:
                row[var(t, 0, midx[(a + 1, b)])] += a + 1
            if (a, b + 1) in midx:
                row[var(t, 1, midx[(a, b + 1)])] += b + 1
            rows.append(row)
    basis, rank = rref_nullspace(rows, nvar)
    # per-triangle forms in monomial coefficients
    blocksA, blocksB, blocksC = [], [], []
    for t, T in enumerate(tris):
        pts = [P[i] for i in T]
        mom = tri_moments(pts, 2 * k)
        # derivative representation: d/dx x^a y^b = a x^{a-1} y^b
        def dx(j):
            a, b = mons[j]; return (a, (a - 1, b)) if a else (0, None)
        def dy(j):
            a, b = mons[j]; return (b, (a, b - 1)) if b else (0, None)
        def I(d1, d2):
            c1, m1 = d1; c2, m2 = d2
            if not c1 or not c2:
                return Fr(0)
            return c1 * c2 * mom[(m1[0] + m2[0], m1[1] + m2[1])]
        Gxx = [[I(dx(i), dx(j)) for j in range(nm)] for i in range(nm)]
        Gyy = [[I(dy(i), dy(j)) for j in range(nm)] for i in range(nm)]
        Gyx = [[I(dy(i), dx(j)) for j in range(nm)] for i in range(nm)]   # int d_y m_i d_x m_j
        n2 = 2 * nm
        A = [[Fr(0)] * n2 for _ in range(n2)]
        for i in range(nm):
            for j in range(nm):
                A[i][j] = 2 * Gxx[i][j] + Gyy[i][j]            # u-u
                A[nm + i][nm + j] = Gxx[i][j] + 2 * Gyy[i][j]  # w-w
                A[i][nm + j] = Gyx[i][j]                       # u_y w_x
                A[nm + j][i] = Gyx[i][j]
        B = [[Fr(0)] * n2 for _ in range(n2)]
        C = [[Fr(0)] * n2 for _ in range(n2)]
        for e in ((T[0], T[1]), (T[1], T[2]), (T[0], T[2])):
            if tuple(sorted(e)) in nit and len(ecount[tuple(sorted(e))]) == 1:
                xa, xb = sorted([P[e[0]][0], P[e[1]][0]])
                # on y=0 only monomials with b=0 survive; d_n = -d_y: -d_y x^a y^b at y=0 nonzero only for b=1
                def xint(p):
                    return (xb ** (p + 1) - xa ** (p + 1)) / (p + 1)
                for i, (a1, b1) in enumerate(mons):
                    if b1 != 0:
                        continue
                    for j, (a2, b2) in enumerate(mons):
                        if b2 == 0:
                            val = xint(a1 + a2)
                            C[i][j] += val; C[nm + i][nm + j] += val
                        if b2 == 1:
                            # B(u,v) = int (d_n u).v : bilinear [v-index i][u-index j], d_n x^a2 y = -x^a2
                            val = -xint(a1 + a2)
                            B[i][j] += val; B[nm + i][nm + j] += val
        blocksA.append(A); blocksB.append(B); blocksC.append(C)
    return dict(basis=basis, rank=rank, nvar=nvar, nm=nm, A=blocksA, B=blocksB, C=blocksC, nT=nT)


def reduce_forms(D):
    """forms restricted to the nullspace basis: returns d x d exact matrices A, Bsym, C."""
    basis = D["basis"]; d = len(basis); nm = D["nm"]
    A = [[Fr(0)] * d for _ in range(d)]; B = [[Fr(0)] * d for _ in range(d)]; C = [[Fr(0)] * d for _ in range(d)]
    for t in range(D["nT"]):
        sl = slice(t * 2 * nm, (t + 1) * 2 * nm)
        Nt = [b[sl] for b in basis]           # d x 2nm
        for M, Out, sym in ((D["A"][t], A, False), (D["B"][t], B, True), (D["C"][t], C, False)):
            if not any(any(r) for r in M):
                continue
            MN = [[sum(M[i][l] * Nt[q][l] for l in range(2 * nm) if M[i][l]) for i in range(2 * nm)] for q in range(d)]
            # MN[q][i] = (M N_q)_i
            for p in range(d):
                Np = Nt[p]
                nzp = [i for i in range(2 * nm) if Np[i]]
                if not nzp:
                    continue
                for q in range(d):
                    Out[p][q] += sum(Np[i] * MN[q][i] for i in nzp)
    Bs = [[(B[i][j] + B[j][i]) / 2 for j in range(d)] for i in range(d)]
    return A, Bs, C


def ldl_pd(M):
    """exact test of positive definiteness via LDL^T (Gaussian elimination without pivoting). Returns (bool, pivots)."""
    n = len(M); M = [row[:] for row in M]; piv = []
    for i in range(n):
        p = M[i][i]
        piv.append(p)
        if p <= 0:
            return False, piv
        for r in range(i + 1, n):
            if M[r][i] != 0:
                f = M[r][i] / p
                Mr, Mi = M[r], M[i]
                for c in range(i + 1, n):
                    Mr[c] -= f * Mi[c]
    return True, piv


def quad(M, v):
    return sum(v[i] * sum(M[i][j] * v[j] for j in range(len(v)) if M[i][j]) for i in range(len(v)) if v[i])


def float_top(A, Bs, C, lam0=None):
    """float lambda* by bisection on  t -> min eig of the pencil (A - 2B + tC, A)  (A is SPD on the space;
    the map is nondecreasing in t and lambda* is its zero).  This avoids splitting off ker C numerically
    (a float splitting of ker C misclassifies directions with tiny trace; cf. paper Section 9).
    Returns (lambda*, z) where z is the eigenvector of the most negative eigenvalue at t = lam0
    (or at t = lambda* - 1e-6 if lam0 is None) -- a witness candidate."""
    import scipy.linalg as sl
    Af = np.array(A, float); Bf = np.array(Bs, float); Cf = np.array(C, float)
    g = lambda t: sl.eigh(Af - 2 * Bf + t * Cf, Af, eigvals_only=True)[0]
    lo, hi = -1.0, 1.0
    while g(lo) > 0: lo *= 2
    while g(hi) < 0: hi *= 2
    for _ in range(200):
        mid = (lo + hi) / 2
        if g(mid) < 0: lo = mid
        else: hi = mid
        if hi - lo < 1e-12 * max(1, abs(hi)): break
    lam = (lo + hi) / 2
    t = float(lam0) if lam0 is not None else lam - 1e-6
    w, V = sl.eigh(Af - 2 * Bf + t * Cf, Af)
    return lam, V[:, 0]


def certify_negative(P, tris, nit, k, lam1=0, label=""):
    t0 = time.time()
    D = assemble(P, tris, nit, k)
    A, Bs, C = reduce_forms(D)
    d = len(A)
    lf, _ = float_top(A, Bs, C)
    M = [[A[i][j] - 2 * Bs[i][j] + Fr(lam1) * C[i][j] for j in range(d)] for i in range(d)]
    ok, piv = ldl_pd(M)
    print(f"[{label}] k={k}: constraints rank {D['rank']}, unknowns {D['nvar']}, div-free dim {d}; float lambda* = {lf:.6f}")
    print(f"   EXACT: A - 2B + ({lam1})C positive definite on the whole div-free space: {ok}"
          f"  (min pivot {float(min(piv)):.4e}, {time.time()-t0:.1f}s)", flush=True)
    return ok, lf, d


def certify_positive(P, tris, nit, k, lam0, label="", den=10 ** 6):
    t0 = time.time()
    D = assemble(P, tris, nit, k)
    A, Bs, C = reduce_forms(D)
    d = len(A)
    lf, z = float_top(A, Bs, C, lam0)
    z = z / np.abs(z).max()
    v = [Fr(float(x)).limit_denominator(den) for x in z]
    num = 2 * quad(Bs, v) - quad(A, v); cc = quad(C, v)
    lam0 = Fr(lam0)
    ok = num - lam0 * cc > 0 and cc > 0
    print(f"[{label}] k={k}: div-free dim {d}; float lambda* = {lf:.6f}; exact witness quotient = {float(num/cc):.6f};"
          f" EXACT 2B - A - {lam0}C > 0: {ok}  ({time.time()-t0:.1f}s)", flush=True)
    return ok, lf, float(num / cc), d


def fan_patch(rim):
    P = [(0, 0)] + list(rim); m = len(rim) - 1
    return P, [(0, i + 1, i + 2) for i in range(m)], {(0, 1), (0, m + 1)}


def strip_patch(W, J, H):
    idx = lambda i, j: j * (W + 1) + i
    P = [(i, j * H) for j in range(J + 1) for i in range(W + 1)]
    tris = []
    for j in range(J):
        for i in range(W):
            a, b, c, d = idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)
            tris += [(a, b, c), (a, c, d)]
    return P, tris, {(idx(i, 0), idx(i + 1, 0)) for i in range(W)}


def ct_patch(P, tris, nit):
    P = [(Fr(p[0]), Fr(p[1])) for p in P]; T2 = []
    for (a, b, c) in tris:
        g = ((P[a][0] + P[b][0] + P[c][0]) / 3, (P[a][1] + P[b][1] + P[c][1]) / 3)
        P.append(g); gi = len(P) - 1
        T2 += [(a, b, gi), (b, c, gi), (c, a, gi)]
    return P, T2, nit


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "regress"
    if mode == "regress":
        # reproduce the paper's S3 value 9.97854 (exact witness quotient) and equilateral-free check
        certify_positive(*fan_patch([(1, 0), (1, 1), (-1, 1), (-1, 0)]), 4, Fr(99, 10), "S3 regression")
