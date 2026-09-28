"""Patch library for the k=3 (any k) Alfeld (M3^3) question.  Independent of td_sharp/alfeld_patch.py.

A patch = list of macro-tetrahedra (4 points each, exact Fractions).  Each macro-tet is Alfeld-split
(barycentre B, sub-tet i = conv(B, face opposite vertex i)).  Fields: continuous piecewise P_k on the
split, zero on the patch boundary (union of macro faces not shared by two macro-tets), div-free.

Gamma-faces: macro faces lying on Gamma_h (they are on the patch boundary, so v=0 there).  For each,
with nu the inward unit normal of the macro-tet on F:
   cons_F(v) = int_F d_nu v dA                 (vector, tangent to F on div-free fields)
   mom_F(v)  = int_F Q_F d_nu v dA,  Q_F = 1/2 sum_{i<j} w^F_ij mu_i mu_j   (default w_ij = l_ij^2)
Total tangential moment: T (2x3 frame) . sum_F mom_F(v).
Everything is exact rational; the same matrices are evaluated mod p (rank certificates), or in floats
(singular values, energy-normalised constants).

Key algebra: v = sum_alpha c_alpha B_alpha (Bernstein, degree k, on a sub-tet with barycentrics lambda);
grad v = k sum_{|beta|=k-1} sum_j c_{beta+e_j} (x) grad lambda_j  B_beta.
On the face lambda_B = 0 where v = 0: grad v = G (x) grad lambda_B, G = k sum_{beta_B=0} c_{beta+e_B} B_beta,
so d_nu v = |grad lambda_B| G and  int_F f d_nu v dA = |grad lambda_B| |F| * (2 * int_ref f G),
with |grad lambda_B| |F| = |F|^2 / (3 vol(subtet))  (rational).
"""
from fractions import Fraction as Fr
from math import factorial
import itertools
import numpy as np

P_MOD = 2147483629  # prime < 2^31

def F3(p):
    return tuple(Fr(c) for c in p)

def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def add(a, b): return tuple(x + y for x, y in zip(a, b))
def scal(s, a): return tuple(s * x for x in a)
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])

def det3(a, b, c): return dot(a, cross(b, c))

def bary_grads(V):
    """gradients of barycentric coordinates of tet V[0..3] (exact)."""
    e1, e2, e3 = sub(V[1], V[0]), sub(V[2], V[0]), sub(V[3], V[0])
    D = det3(e1, e2, e3)
    assert D != 0, 'degenerate tet'
    g1 = scal(1 / D, cross(e2, e3)); g2 = scal(1 / D, cross(e3, e1)); g3 = scal(1 / D, cross(e1, e2))
    g0 = scal(-1, add(add(g1, g2), g3))
    return [g0, g1, g2, g3], D  # D = 6 vol (signed)

def mindex(k, n):
    return [a for a in itertools.product(range(k + 1), repeat=n) if sum(a) == k]

def multinom(a):
    r = factorial(sum(a))
    for x in a: r //= factorial(x)
    return r

def tri_mono_int(g):
    """int over ref triangle (area 1/2) of mu1^g1 mu2^g2 mu3^g3 = g!/(|g|+2)!."""
    return Fr(factorial(g[0]) * factorial(g[1]) * factorial(g[2]), factorial(sum(g) + 2))

class Patch:
    def __init__(self, macro, gfaces, k=3, weights=None):
        """macro: list of 4-tuples of points; gfaces: list of (macro index, (i,j,l) local vertex indices);
        weights: optional dict face_index -> (w12,w13,w23) (else l_ij^2)."""
        self.k = k
        self.macro = [[F3(p) for p in T] for T in macro]
        self.gfaces = gfaces
        # patch boundary faces
        cnt = {}
        for T in self.macro:
            for i in range(4):
                f = frozenset(T[j] for j in range(4) if j != i)
                cnt[f] = cnt.get(f, 0) + 1
        self.bfaces = [f for f, c in cnt.items() if c == 1]
        for (ti, loc) in gfaces:
            f = frozenset(self.macro[ti][j] for j in loc)
            assert cnt[f] == 1, 'Gamma face must be on patch boundary'
        # sub-tets
        self.subs = []
        for ti, T in enumerate(self.macro):
            B = scal(Fr(1, 4), add(add(T[0], T[1]), add(T[2], T[3])))
            for i in range(4):
                self.subs.append(([B] + [T[j] for j in range(4) if j != i], ti, i))
        # zero keys: lattice points of boundary faces
        zero = set()
        for f in self.bfaces:
            P = list(f)
            for b in mindex(k, 3):
                zero.add(tuple(sum(Fr(b[j], k) * P[j][c] for j in range(3)) for c in range(3)))
        self.key = {}
        self.sub_idx = []  # per sub-tet: dict alpha -> unknown node index or None
        for (V, ti, i) in self.subs:
            d = {}
            for a in mindex(k, 4):
                pt = tuple(sum(Fr(a[j], k) * V[j][c] for j in range(4)) for c in range(3))
                if pt in zero:
                    d[a] = None
                else:
                    if pt not in self.key: self.key[pt] = len(self.key)
                    d[a] = self.key[pt]
            self.sub_idx.append(d)
        self.nn = len(self.key); self.nunk = 3 * self.nn
        self.weights = weights or {}
        self._build()

    def _build(self):
        k = self.k
        rows = []  # list of dict col->Fr
        self.grads = []
        for (V, ti, i), d in zip(self.subs, self.sub_idx):
            G, D = bary_grads(V); self.grads.append((G, D))
            for b in mindex(k - 1, 4):
                r = {}
                for j in range(4):
                    a = tuple(b[m] + (1 if m == j else 0) for m in range(4))
                    n = d[a]
                    if n is None: continue
                    for c in range(3):
                        if G[j][c] != 0:
                            col = 3 * n + c
                            r[col] = r.get(col, 0) + k * G[j][c]
                rows.append(r)
        self.div_rows = rows
        # Gamma face functionals: each a list of 3 row-dicts (x,y,z components of the vector functional)
        self.cons = []; self.mom = []; self.fgeom = []
        for fi, (ti, loc) in enumerate(self.gfaces):
            T = self.macro[ti]
            opp = [j for j in range(4) if j not in loc][0]
            si = [s for s, (V, t, i) in enumerate(self.subs) if t == ti and i == opp][0]
            V = self.subs[si][0]  # V[0]=B, V[1..3] = macro verts except opp, in macro order
            P = V[1:]
            G, D = self.grads[si]
            area2 = dot(cross(sub(P[1], P[0]), sub(P[2], P[0])), cross(sub(P[1], P[0]), sub(P[2], P[0])))  # (2|F|)^2
            fac = (area2 / 4) / (abs(D) / 2)  # |F|^2/(3 vol) , vol = |D|/6
            l2 = {(0, 1): dot(sub(P[0], P[1]), sub(P[0], P[1])), (0, 2): dot(sub(P[0], P[2]), sub(P[0], P[2])),
                  (1, 2): dot(sub(P[1], P[2]), sub(P[1], P[2]))}
            w = self.weights.get(fi)
            if w is not None:
                # weights given w.r.t. the loc ordering (w12,w13,w23) -> map to P ordering
                # P is macro order of the 3 face vertices; loc may be permuted
                order = [j for j in range(4) if j != opp]  # macro indices in P order
                wl = {(0, 1): w[0], (0, 2): w[1], (1, 2): w[2]}
                pos = {loc[m]: m for m in range(3)}
                l2 = {}
                for (a_, b_) in [(0, 1), (0, 2), (1, 2)]:
                    A, Bq = sorted((pos[order[a_]], pos[order[b_]]))
                    l2[(a_, b_)] = Fr(wl[(A, Bq)])
            # Q = 1/2 sum l2 mu_a mu_b ; G_poly = k sum_{beta on F} c_{beta+e_B} B_beta(mu)
            crow = [dict() for _ in range(3)]; mrow = [dict() for _ in range(3)]
            d = self.sub_idx[si]
            for b in mindex(k - 1, 3):
                a = (1,) + b
                n = d[a]
                if n is None: continue
                cb = multinom(b)
                # int_ref B_b = cb * int mu^b ; int_ref Q B_b
                i0 = cb * tri_mono_int(b)
                i1 = 0
                for (A, Bq), lv in l2.items():
                    g = list(b); g[A] += 1; g[Bq] += 1
                    i1 += Fr(1, 2) * lv * cb * tri_mono_int(tuple(g))
                for c in range(3):
                    crow[c][3 * n + c] = fac * 2 * k * i0
                    mrow[c][3 * n + c] = fac * 2 * k * i1
            self.cons.append(crow); self.mom.append(mrow)
            self.fgeom.append(P)

    # ---------- matrices ----------
    def dense(self, rows, mode):
        n = self.nunk
        if mode == 'float':
            M = np.zeros((len(rows), n))
            for i, r in enumerate(rows):
                for c, v in r.items(): M[i, c] = float(v)
            return M
        M = np.zeros((len(rows), n), dtype=np.int64)
        for i, r in enumerate(rows):
            for c, v in r.items():
                M[i, c] = (v.numerator % P_MOD) * pow(v.denominator % P_MOD, P_MOD - 2, P_MOD) % P_MOD
        return M

    def mom_rows(self, frame, faces=None):
        """rows of frame . sum_F mom_F (frame: list of 3-vectors)."""
        faces = range(len(self.gfaces)) if faces is None else faces
        out = []
        for t in frame:
            t = F3(t); r = {}
            for fi in faces:
                for c in range(3):
                    if t[c] == 0: continue
                    for col, v in self.mom[fi][c].items():
                        r[col] = r.get(col, 0) + t[c] * v
            out.append(r)
        return out

    def cons_rows(self):
        out = []
        for crow in self.cons:
            out.extend(crow)
        return out

def rank_mod(M, p=P_MOD):
    M = M.copy() % p
    nr, nc = M.shape
    r = 0
    for c in range(nc):
        if r == nr: break
        piv = np.nonzero(M[r:, c])[0]
        if len(piv) == 0: continue
        pr = r + piv[0]
        if pr != r: M[[r, pr]] = M[[pr, r]]
        inv = pow(int(M[r, c]), p - 2, p)
        M[r] = (M[r] * inv) % p
        col = M[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr):
            M[nzr] = (M[nzr] - (col[nzr, None] * M[r][None, :]) % p) % p
        r += 1
    return r

def moment_rank_mod(P, frame, faces=None):
    """rank of the tangential total moment on Y^1 (mod p): rank([D;C;M]) - rank([D;C]).
    Returns (rank div, rank [D;C], rank M on Y^1).  If rank[D;C;M]-rank[D;C] = r mod p then
    rank over Q >= r provided rank_Q[D;C] = rank_p[D;C] (checked against the a-priori upper bound)."""
    D = P.dense(P.div_rows, 'mod'); C = P.dense(P.cons_rows(), 'mod'); M = P.dense(P.mom_rows(frame, faces), 'mod')
    rd = rank_mod(D); rdc = rank_mod(np.vstack([D, C])); rall = rank_mod(np.vstack([D, C, M]))
    return rd, rdc, rall - rdc

# ---------- float: energy-normalised constants ----------
def bern_gram(n, dim=3):
    """Gram matrix of degree-n Bernstein polys on a simplex of unit volume (dim=3)."""
    idx = mindex(n, dim + 1)
    G = np.zeros((len(idx), len(idx)))
    for i, a in enumerate(idx):
        for j, b in enumerate(idx):
            g = tuple(x + y for x, y in zip(a, b))
            num = multinom(a) * multinom(b) * factorial(dim)
            for x in g: num *= factorial(x)
            G[i, j] = num / factorial(2 * n + dim)
    return idx, G

def stiffness(P):
    """H^1 seminorm Gram matrix sum_K int |grad v|^2 on all unknowns (float)."""
    k = P.k
    idx, Gm = bern_gram(k - 1)
    pos = {a: i for i, a in enumerate(idx)}
    A = np.zeros((P.nunk, P.nunk))
    for (V, ti, i), d, (G, D) in zip(P.subs, P.sub_idx, P.grads):
        vol = abs(float(D)) / 6
        Gf = [np.array([float(c) for c in g]) for g in G]
        # grad v component c, derivative dir e: sum_beta sum_j k c_{beta+e_j} G_j[e] B_beta
        # build map: for each unknown (node,c) list of (beta index, e, coeff)
        for c in range(3):
            # matrix Bm[e][beta, node]
            nb = len(idx)
            cols = {}
            for b in idx:
                for j in range(4):
                    a = tuple(b[m] + (1 if m == j else 0) for m in range(4))
                    n = d[a]
                    if n is None: continue
                    cols.setdefault(n, []).append((pos[b], j))
            nodes = sorted(cols)
            if not nodes: continue
            for e in range(3):
                Bm = np.zeros((nb, len(nodes)))
                for jj, n in enumerate(nodes):
                    for (bi, j) in cols[n]:
                        Bm[bi, jj] += k * Gf[j][e]
                Kl = vol * Bm.T @ Gm @ Bm
                ii = np.array([3 * n + c for n in nodes])
                A[np.ix_(ii, ii)] += Kl
    return A

def null_float(M, tol=1e-9):
    if M.shape[0] == 0: return np.eye(M.shape[1])
    u, s, vt = np.linalg.svd(M, full_matrices=True)
    r = int((s > tol * s[0]).sum())
    return vt[r:].T

def kappa_float(P, frame, faces=None, return_all=False):
    """singular values of v -> frame.sum mom(v) on Y^1, with ||grad v||=1 normalisation (frame orthonormalised)."""
    fr = np.array([[float(c) for c in t] for t in frame]); fr, _ = np.linalg.qr(fr.T); fr = fr.T
    D = P.dense(P.div_rows, 'float'); C = P.dense(P.cons_rows(), 'float')
    Ny = null_float(D)
    Ny1 = Ny @ null_float(C @ Ny)
    A = stiffness(P)
    Gy = Ny1.T @ A @ Ny1
    L = np.linalg.cholesky(Gy + 0 * np.eye(len(Gy)))
    Mrows = []
    for t in fr:
        Mrows.append(P.dense(P.mom_rows([tuple(Fr(x).limit_denominator(10**12) for x in t)], faces), 'float')[0])
    M = np.array(Mrows) @ Ny1
    Ml = np.linalg.solve(L, M.T).T  # M L^{-T}
    s = np.linalg.svd(Ml, compute_uv=False)
    return (s, Ny1.shape[1]) if return_all else s

def edge_moment_rows(P, fi):
    """rows of the joint edge-moment map of Gamma face fi: for w in e_12,e_13,e_23, 3 comps each (9 rows;
    on div-free fields the normal comps vanish, so the rank is <= 6)."""
    rows = []
    for w in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:
        Q = Patch.__new__(Patch)
        rows.extend(_mom_face_rows(P, fi, w))
    return rows

def _mom_face_rows(P, fi, w):
    k = P.k
    ti, loc = P.gfaces[fi]
    T = P.macro[ti]
    opp = [j for j in range(4) if j not in loc][0]
    si = [s for s, (V, t, i) in enumerate(P.subs) if t == ti and i == opp][0]
    V = P.subs[si][0]; Pp = V[1:]
    G, D = P.grads[si]
    area2 = dot(cross(sub(Pp[1], Pp[0]), sub(Pp[2], Pp[0])), cross(sub(Pp[1], Pp[0]), sub(Pp[2], Pp[0])))
    fac = (area2 / 4) / (abs(D) / 2)
    order = [j for j in range(4) if j != opp]
    pos = {loc[m]: m for m in range(3)}
    wl = {(0, 1): w[0], (0, 2): w[1], (1, 2): w[2]}
    l2 = {}
    for (a_, b_) in [(0, 1), (0, 2), (1, 2)]:
        A, Bq = sorted((pos[order[a_]], pos[order[b_]]))
        l2[(a_, b_)] = Fr(wl[(A, Bq)])
    d = P.sub_idx[si]
    mrow = [dict() for _ in range(3)]
    for b in mindex(k - 1, 3):
        n = d[(1,) + b]
        if n is None: continue
        cb = multinom(b); i1 = 0
        for (A, Bq), lv in l2.items():
            g = list(b); g[A] += 1; g[Bq] += 1
            i1 += lv * cb * tri_mono_int(tuple(g))
        for c in range(3):
            mrow[c][3 * n + c] = fac * 2 * k * i1
    return mrow

def rank_on_Y1_mod(P, rows):
    D = P.dense(P.div_rows, 'mod'); C = P.dense(P.cons_rows(), 'mod'); M = P.dense(rows, 'mod')
    rdc = rank_mod(np.vstack([D, C])); rall = rank_mod(np.vstack([D, C, M]))
    return rdc, rall - rdc

def trace_rows(P, fi):
    """rows giving the Bernstein coefficients of G (d_nu v / |grad lambda_B|) on Gamma face fi:
    for each beta (|beta|=k-1 on the face, in macro vertex order) and comp c.  Returns (rows, betas)."""
    k = P.k
    ti, loc = P.gfaces[fi]
    opp = [j for j in range(4) if j not in loc][0]
    si = [s for s, (V, t, i) in enumerate(P.subs) if t == ti and i == opp][0]
    d = P.sub_idx[si]
    rows = []; betas = []
    for b in mindex(k - 1, 3):
        n = d[(1,) + b]
        for c in range(3):
            rows.append({3 * n + c: Fr(k)} if n is not None else {})
            betas.append((b, c))
    return rows, betas

def rank_on_Y_mod(P, rows):
    D = P.dense(P.div_rows, 'mod'); M = P.dense(rows, 'mod')
    rd = rank_mod(D); return rank_mod(np.vstack([D, M])) - rd
