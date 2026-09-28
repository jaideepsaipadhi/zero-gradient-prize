"""Exact (rational-arithmetic) computations for the Clough--Tocher (barycentric) extension.

Part 1 (local surjectivity).  On the reference triangle T^ = conv{(0,0),(1,0),(0,1)} split at its barycentre,
  div : V_0(T^) = {continuous piecewise P_k on the split, vanishing on the boundary of T^}^2
        -> Q_0(T^) = {discontinuous piecewise P_{k-1} on the split, mean zero}
is checked to be ONTO for k = 2,...,KMAX by an exact rank computation.  (For general k this is the local
result behind Guzman--Neilan (SINUM 2018) in d = 2; for k = 2 it is Arnold--Qin.)  By the Piola transform the
statement transfers to every triangle, with constant depending only on shape regularity.

Part 2 (penalty-threshold certificates).  Star: fan of m macro triangles at O=(0,0) with rational rim points,
each macro triangle split at its barycentre.  Free (Nitsche) boundary = the two segments of y=0 at O, outward
normal (0,-1); zero Dirichlet data on the rim chords.  Velocity: continuous piecewise P_k on the refined star.
Exact constraints: continuity across all interior refined edges, zero on the rim chords, div v = 0 pointwise.
We report
  * dim V_star, dim Pi_star, rank of div : V_star -> Pi_star (exact; full rank <=> generic kernel dimension),
  * the exact divergence-free kernel dimension,
  * the floating-point maximiser lambda* of (2B - A)/C on the kernel,
  * a rational witness v in the exact kernel and its EXACT Rayleigh quotient (2B - A)/C,
  * the exact check 2B - A - lam0*C > 0 for a stated rational lam0.
Here A = sum_T int (1/2) Dv:Dv,  B = int_free (d_n v).v,  C = int_free |v|^2.
By dilation invariance, N_h(v,v) = A - 2B + (mu*l_z/h) C < 0 whenever mu*l_z/h < lam0.

Usage:  python3 lamstar_exact_ct.py  > ../logs/lamstar_exact_ct.log
"""
import sys, time
from fractions import Fraction as Fr
from math import factorial
import numpy as np
from sympy import QQ
from sympy.polys.matrices import DomainMatrix

KMAX_LOCAL = 6


# ---------------------------------------------------------------- exact polynomial helpers
def monos(deg):
    return [(a, b) for a in range(deg + 1) for b in range(deg + 1 - a)]


def tri_area2(T):
    (x0, y0), (x1, y1), (x2, y2) = T
    return abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))          # = 2|T|


def _expand_linear_power(coeffs, p):
    """(sum_i coeffs[i]*lam_i)^p as dict {alpha(3-tuple): Fraction}"""
    out = {(0, 0, 0): Fr(1)}
    for _ in range(p):
        new = {}
        for al, c in out.items():
            for i in range(3):
                be = list(al); be[i] += 1; be = tuple(be)
                new[be] = new.get(be, 0) + c * coeffs[i]
        out = new
    return out


def tri_moments(T, maxdeg):
    """exact {(a,b): int_T x^a y^b} for a+b <= maxdeg, via int lam^g = 2|T| g!/(|g|+2)!"""
    A2 = tri_area2(T)
    xs = [Fr(P[0]) for P in T]; ys = [Fr(P[1]) for P in T]
    px = [_expand_linear_power(xs, a) for a in range(maxdeg + 1)]
    py = [_expand_linear_power(ys, b) for b in range(maxdeg + 1)]
    mom = {}
    for a in range(maxdeg + 1):
        for b in range(maxdeg + 1 - a):
            s = Fr(0)
            for al, ca in px[a].items():
                for be, cb in py[b].items():
                    g = (al[0] + be[0], al[1] + be[1], al[2] + be[2])
                    s += ca * cb * Fr(factorial(g[0]) * factorial(g[1]) * factorial(g[2]), factorial(sum(g) + 2))
            mom[(a, b)] = A2 * s
    return mom


def ct_split(T):
    P0, P1, P2 = [(Fr(p[0]), Fr(p[1])) for p in T]
    g = ((P0[0] + P1[0] + P2[0]) / 3, (P0[1] + P1[1] + P2[1]) / 3)
    return [(P0, P1, g), (P1, P2, g), (P2, P0, g)]


def seg_points(P, Q, n):
    return [(P[0] + Fr(j, n) * (Q[0] - P[0]), P[1] + Fr(j, n) * (Q[1] - P[1])) for j in range(n + 1)]


def dm(rows, ncols):
    return DomainMatrix([[QQ(r.numerator, r.denominator) for r in row] for row in rows], (len(rows), ncols), QQ)


# ---------------------------------------------------------------- generic assembly on a refined mesh
class PwPoly:
    """unknown vector fields, one P_k polynomial per component per triangle, monomial basis in global coords"""

    def __init__(self, tris, k):
        self.tris, self.k = tris, k
        self.M = monos(k)
        self.nm = len(self.M)
        self.nvar = len(tris) * 2 * self.nm
        self.idx = {m: i for i, m in enumerate(self.M)}

    def base(self, t, c):
        return (t * 2 + c) * self.nm

    def eval_row(self, t, c, P):
        row = [Fr(0)] * self.nvar
        b = self.base(t, c)
        for i, (a, bb) in enumerate(self.M):
            row[b + i] = P[0] ** a * P[1] ** bb
        return row

    def edge_constraints(self, shared, zero_edges):
        rows = []
        for (t1, t2, P, Q) in shared:
            for X in seg_points(P, Q, self.k):
                for c in range(2):
                    r1 = self.eval_row(t1, c, X); r2 = self.eval_row(t2, c, X)
                    rows.append([u - v for u, v in zip(r1, r2)])
        for (t, P, Q) in zero_edges:
            for X in seg_points(P, Q, self.k):
                for c in range(2):
                    rows.append(self.eval_row(t, c, X))
        return rows

    def div_rows(self):
        """rows mapping coefficients -> coefficients of div on each triangle (monomials of degree <= k-1)"""
        rows = []
        for t in range(len(self.tris)):
            for (a, b) in monos(self.k - 1):
                row = [Fr(0)] * self.nvar
                row[self.base(t, 0) + self.idx[(a + 1, b)]] += a + 1
                row[self.base(t, 1) + self.idx[(a, b + 1)]] += b + 1
                rows.append(row)
        return rows


def shared_and_boundary_edges(tris):
    """shared refined edges (t1,t2,P,Q) and boundary edges (t,P,Q)"""
    emap = {}
    for t, T in enumerate(tris):
        for i in range(3):
            P, Q = T[i], T[(i + 1) % 3]
            key = frozenset([P, Q])
            emap.setdefault(key, []).append((t, P, Q))
    shared, bdry = [], []
    for key, lst in emap.items():
        if len(lst) == 2:
            shared.append((lst[0][0], lst[1][0], lst[0][1], lst[0][2]))
        else:
            bdry.append(lst[0])
    return shared, bdry


# ---------------------------------------------------------------- Part 1: local surjectivity
def local_surjectivity(k):
    T = [(0, 0), (1, 0), (0, 1)]
    tris = ct_split(T)
    S = PwPoly(tris, k)
    shared, bdry = shared_and_boundary_edges(tris)
    rows = S.edge_constraints(shared, bdry)
    Kc = dm(rows, S.nvar).nullspace().transpose()           # nvar x dimV0, basis of V_0
    dimV0 = Kc.shape[1]
    D = dm(S.div_rows(), S.nvar)
    rk = (D * Kc).rank()
    dimQ = 3 * len(monos(k - 1))
    nodes_int = 1 + 3 * (k - 1) + 3 * (k - 1) * (k - 2) // 2
    print(f"  k={k}: dim V_0(T^) = {dimV0} (node count {2*nodes_int}),  dim Q(T^) = {dimQ},"
          f"  rank div = {rk}  -> onto Q_0 (rank = dimQ-1): {rk == dimQ - 1};"
          f"  dim ker div = {dimV0 - rk}", flush=True)
    return rk == dimQ - 1


# ---------------------------------------------------------------- Part 2: stars
def exact_forms_value(S, tris, moms, free_segments, c):
    """exact A, B, C for a single coefficient vector c (list of Fractions)"""
    k = S.k
    A = Fr(0)
    for t in range(len(tris)):
        u = [{m: c[S.base(t, comp) + i] for i, m in enumerate(S.M)} for comp in range(2)]

        def d(p, dirn):
            out = {}
            for (a, b), v in p.items():
                if v == 0:
                    continue
                if dirn == 0 and a > 0:
                    out[(a - 1, b)] = out.get((a - 1, b), 0) + a * v
                if dirn == 1 and b > 0:
                    out[(a, b - 1)] = out.get((a, b - 1), 0) + b * v
            return out

        def mul(p, q):
            out = {}
            for (a, b), v in p.items():
                for (e, f), w in q.items():
                    out[(a + e, b + f)] = out.get((a + e, b + f), 0) + v * w
            return out

        def integ(p):
            return sum(v * moms[t][m] for m, v in p.items())

        u1x, u1y, u2x, u2y = d(u[0], 0), d(u[0], 1), d(u[1], 0), d(u[1], 1)
        s = {}
        for m, v in u1y.items():
            s[m] = s.get(m, 0) + v
        for m, v in u2x.items():
            s[m] = s.get(m, 0) + v
        A += 2 * integ(mul(u1x, u1x)) + integ(mul(s, s)) + 2 * integ(mul(u2y, u2y))
    B = Fr(0); C = Fr(0)
    for (t, lo, hi) in free_segments:       # on y = 0, x in [lo,hi]; n = (0,-1), d_n = -d_y
        for comp in range(2):
            cu = [c[S.base(t, comp) + i] for i in range(S.nm)]
            p0 = {}; p1 = {}       # u(x,0) and d_y u(x,0) as polys in x
            for i, (a, b) in enumerate(S.M):
                if b == 0:
                    p0[a] = p0.get(a, 0) + cu[i]
                if b == 1:
                    p1[a] = p1.get(a, 0) + cu[i]
            for a, v in p0.items():
                for e, w in p0.items():
                    C += v * w * (Fr(hi) ** (a + e + 1) - Fr(lo) ** (a + e + 1)) / (a + e + 1)
                for e, w in p1.items():
                    B += -v * w * (Fr(hi) ** (a + e + 1) - Fr(lo) ** (a + e + 1)) / (a + e + 1)
    return A, B, C


def exact_forms(S, tris, moms, free_segments):
    """exact Gram matrices (DomainMatrix over QQ) of A, B, C in the monomial-coefficient variables"""
    nv = S.nvar
    A = {}; B = {}; C = {}

    def add(Mx, i, j, v):
        if v:
            Mx[(i, j)] = Mx.get((i, j), 0) + v
    M = S.M
    dx = [(a, (a - 1, b)) if a > 0 else (0, None) for (a, b) in M]
    dy = [(b, (a, b - 1)) if b > 0 else (0, None) for (a, b) in M]
    for t in range(len(tris)):
        mo = moms[t]
        b1, b2 = S.base(t, 0), S.base(t, 1)

        def g(d1, i, d2, j):
            (ci, mi), (cj, mj) = d1[i], d2[j]
            if mi is None or mj is None:
                return 0
            return ci * cj * mo[(mi[0] + mj[0], mi[1] + mj[1])]
        for i in range(S.nm):
            for j in range(S.nm):
                gxx, gyy, gyx = g(dx, i, dx, j), g(dy, i, dy, j), g(dy, i, dx, j)
                # density 2 u1x^2 + (u1y + u2x)^2 + 2 u2y^2
                add(A, b1 + i, b1 + j, 2 * gxx + gyy)
                add(A, b2 + i, b2 + j, 2 * gyy + gxx)
                add(A, b1 + i, b2 + j, gyx)
                add(A, b2 + j, b1 + i, gyx)
    for (t, lo, hi) in free_segments:          # y = 0, n = (0,-1), d_n = -d_y
        for comp in range(2):
            b0 = S.base(t, comp)
            for i, (a, b) in enumerate(M):
                if b != 0:
                    continue
                for j, (e, f) in enumerate(M):
                    I = (Fr(hi) ** (a + e + 1) - Fr(lo) ** (a + e + 1)) / (a + e + 1)
                    if f == 0:
                        add(C, b0 + i, b0 + j, I)
                    if f == 1:
                        add(B, b0 + i, b0 + j, -I / 2)
                        add(B, b0 + j, b0 + i, -I / 2)

    def todm(Mx):
        rows = [[QQ(0)] * nv for _ in range(nv)]
        for (i, j), v in Mx.items():
            v = Fr(v)
            rows[i][j] = QQ(v.numerator, v.denominator)
        return DomainMatrix(rows, (nv, nv), QQ)
    return todm(A), todm(B), todm(C)


def tofloat(Mx):
    L = Mx.to_Matrix()
    return np.array([[float(L[i, j]) for j in range(L.shape[1])] for i in range(L.shape[0])])


def star(rim, k, lam0, label, split=True):
    t0 = time.time()
    O = (Fr(0), Fr(0))
    rim = [(Fr(p[0]), Fr(p[1])) for p in rim]
    macro = [(O, rim[i], rim[i + 1]) for i in range(len(rim) - 1)]
    m = len(macro)
    tris = []
    for T in macro:
        tris += ct_split(T) if split else [T]   # CT: (O, r_i, g), (r_i, r_{i+1}, g), (r_{i+1}, O, g)
    S = PwPoly(tris, k)
    shared, bdry = shared_and_boundary_edges(tris)
    rim_edges = [(t, P, Q) for (t, P, Q) in bdry if not (P[1] == 0 and Q[1] == 0)]
    free = [(t, P, Q) for (t, P, Q) in bdry if (P[1] == 0 and Q[1] == 0)]
    assert len(rim_edges) == m and len(free) == 2
    free_segments = [(t, min(P[0], Q[0]), max(P[0], Q[0])) for (t, P, Q) in free]
    # V_star: continuity + zero on rim
    rows = S.edge_constraints(shared, rim_edges)
    KV = dm(rows, S.nvar).nullspace().transpose()
    dimV = KV.shape[1]
    D = dm(S.div_rows(), S.nvar)
    DK = D * KV
    rkdiv = DK.rank()
    dimPi = D.shape[0]
    Kz = DK.nullspace().transpose()                       # coordinates in V_star basis
    K = KV * Kz                                            # nvar x kdim, exact kernel basis
    kdim = K.shape[1]
    print(f"[{label}{'' if split else ' UNSPLIT'}] k={k}, m={m} macro triangles ({len(tris)} elements): dim V_star = {dimV}, dim Pi_star = {dimPi},"
          f" rank div = {rkdiv} (onto: {rkdiv == dimPi}); div-free kernel dim = {kdim}"
          f" (= dimV - dimPi: {kdim == dimV - dimPi})", flush=True)
    moms = [tri_moments(T, 2 * k) for T in tris]
    Ad, Bd, Cd = exact_forms(S, tris, moms, free_segments)
    KT = K.transpose()
    Ak, Bk, Ck = KT * Ad * K, KT * Bd * K, KT * Cd * K       # exact kernel forms
    Mk = Bk * QQ(2) - Ak                                      # 2B - A (B is symmetrised already)
    # exact splitting of the kernel into ker C (fields with zero trace on the free boundary) and range C
    Nn = Ck.nullspace().transpose()                           # kdim x nn
    Rr = Ck.rref()[0]
    RrM = Rr.to_Matrix()
    rows_nz = [list(RrM.row(i)) for i in range(RrM.shape[0]) if any(v != 0 for v in RrM.row(i))]
    from sympy import Matrix as SM
    R = DomainMatrix.from_Matrix(SM(rows_nz)).convert_to(QQ).transpose()   # kdim x r, spans range(C)
    r = R.shape[1]; nn = Nn.shape[1]
    Myy = R.transpose() * Mk * R; Cyy = R.transpose() * Ck * R
    if nn:
        Mkk = Nn.transpose() * Mk * Nn; Mky = Nn.transpose() * Mk * R
        X = Mkk.lu_solve(Mky)                                 # exact  Mkk^{-1} Mky
        Seff = Myy - Mky.transpose() * X                      # exact Schur complement
    else:
        Seff = Myy
    Sf, Cf = tofloat(Seff), tofloat(Cyy)
    Lc = np.linalg.cholesky(Cf); Li = np.linalg.inv(Lc)
    ev, EV = np.linalg.eigh(Li @ Sf @ Li.T)
    lam_float = ev[-1]
    yv = Li.T @ EV[:, -1]; yv = yv / np.abs(yv).max()
    print(f"    exact split: rank C on kernel = {r}, dim(ker C) = {nn};"
          f" top eigenvalues of (2B-A)/C: {np.round(ev[::-1][:4], 6).tolist()}", flush=True)
    if lam0 is None:                                       # certify the float value truncated to 2 decimals
        lam0 = Fr(int(np.floor(lam_float * 100)), 100)
    for den in (10 ** 4, 10 ** 6, 10 ** 8, 10 ** 10):
        yr = [Fr(float(v)).limit_denominator(den) for v in yv]
        Y = DomainMatrix([[QQ(v.numerator, v.denominator)] for v in yr], (r, 1), QQ)
        z = R * Y
        if nn:
            z = z - Nn * (X * Y)                              # exact optimal ker-C component for this y
        cvec = K * z                                          # exact coefficient vector, exactly in the kernel
        cl = cvec.to_Matrix()
        c = [Fr(int(cl[i, 0].p), int(cl[i, 0].q)) for i in range(S.nvar)]
        A, B, C = exact_forms_value(S, tris, moms, free_segments, c)
        lam_exact = (2 * B - A) / C
        ok = (2 * B - A - lam0 * C) > 0
        if ok:
            break
    # independent exact re-check of the kernel membership of the witness
    res = D * cvec
    assert all(v == 0 for v in res.to_Matrix()), "witness not divergence-free"
    print(f"    float lambda* = {lam_float:.8f};  exact Rayleigh quotient of rational witness = {float(lam_exact):.8f}")
    print(f"    EXACT check  2B - A - ({lam0})*C > 0 :  {ok}   (value {float(2*B - A - lam0*C):.4e});"
          f"  witness denominators <= {den}; time {time.time()-t0:.1f}s", flush=True)
    return dict(k=k, m=m, dimV=dimV, dimPi=dimPi, rank=rkdiv, kdim=kdim, lam_float=lam_float,
                lam_exact=lam_exact, ok=ok)


STARS = {
    "S3": [(1, 0), (1, 1), (-1, 1), (-1, 0)],
    "S4": [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0)],
    "S4'": [(1, 0), (1, 1), (0, Fr(6, 5)), (-1, 1), (-1, 0)],
    "S2": [(1, 0), (0, 1), (-1, 0)],
}

if __name__ == "__main__":
    which = sys.argv[1:] or ["local", "validate", "2", "3", "4"]
    print("Clough--Tocher exact computations (code/lamstar_exact_ct.py)", flush=True)
    if "local" in which:
        print("Part 1: local surjectivity of div on the CT-split reference triangle", flush=True)
        for k in range(2, KMAX_LOCAL + 1):
            local_surjectivity(k)
    if "validate" in which:
        print("Validation: UNSPLIT fans with k = 4 must reproduce the paper's 9.97854 (S3) and 9.26747 (S4')", flush=True)
        star(STARS["S3"], 4, Fr(9), "S3", split=False)
        star(STARS["S4'"], 4, Fr(9), "S4'", split=False)
    lam0s = {2: None, 3: None, 4: None}
    for kk in (2, 3, 4):
        if str(kk) in which:
            print(f"Part 2: star certificates, k = {kk}", flush=True)
            for name in ("S2", "S3", "S4", "S4'"):
                star(STARS[name], kk, lam0s[kk], name)
