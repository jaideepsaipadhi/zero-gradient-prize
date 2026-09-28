"""Single Alfeld-split tetrahedron witness for 3D penalty necessity (exact rational arithmetic).
T = conv(P0,P1,P2,P3), face F=(P1,P2,P3) in the plane z=0 (the Nitsche boundary), fluid in z>0, n=(0,0,-1).
Alfeld split: m = barycentre, S_i = conv(m, face opposite P_i); S_0 carries F.
Space W_k(T): continuous piecewise P_k (global monomials on each S_i), zero on the three faces of T other than F,
pointwise div-free.  Extended by zero it lies in Z_h for every mesh containing T (as long as T does not meet dR):
continuity across the other faces holds because v=0 there.
Forms: A = sum_S int 1/2 Dv:Dv,  B = int_F (d_n v).v, d_n=-d_z,  C = int_F |v|^2.
lambda_T := diam(F) * max_{v in W_k} (2B-A)/C.  If mu*diam(F)/h < lambda_T then N_h is indefinite on Z_h.
Certificate: a rational v in W_k with (2B-A-c*C)(v) > 0 exactly (c = certified lower bound)."""
from fractions import Fraction as Fr
from itertools import product
from math import factorial
import sys, numpy as np
from sympy.polys.matrices import DomainMatrix
from sympy import QQ

def monos(deg, nv=3):
    return [e for d in range(deg + 1) for e in product(range(d + 1), repeat=nv) if sum(e) == d]

def pmul(p, q):
    r = {}
    for a, ca in p.items():
        for b, cb in q.items():
            e = tuple(i + j for i, j in zip(a, b)); r[e] = r.get(e, 0) + ca * cb
    return {e: c for e, c in r.items() if c != 0}

def padd(p, q, s=1):
    r = dict(p)
    for e, c in q.items(): r[e] = r.get(e, 0) + s * c
    return {e: c for e, c in r.items() if c != 0}

def ppow(p, n, nv):
    r = {(0,) * nv: Fr(1)}
    for _ in range(n): r = pmul(r, p)
    return r

def subst(expo, lin):
    """x^expo with x_i = lin[i] (polys in params) -> poly in params"""
    nv = next(len(e) for p in lin for e in p) if any(lin) else 1
    assert all(len(e) == nv for p in lin for e in p)
    r = {(0,) * nv: Fr(1)}
    for i, n in enumerate(expo):
        if n: r = pmul(r, ppow(lin[i], n, nv))
    return r

def affine(o, vecs):
    """x_i = o_i + sum_j s_j vecs[j]_i as polys in s"""
    ns = len(vecs); lin = []
    for i in range(3):
        p = {(0,) * ns: Fr(o[i])}
        for j, v in enumerate(vecs):
            e = [0] * ns; e[j] = 1
            if v[i] != 0: p[tuple(e)] = Fr(v[i])
        lin.append({e: c for e, c in p.items() if c != 0})
    return lin

def ref_int(e):  # int over reference simplex of s^e
    n = len(e); return Fr(np.prod([factorial(a) for a in e]).item(), factorial(sum(e) + n))

def det3(a, b, c):
    return a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0])

def build(P, k):
    P = [tuple(Fr(c) for c in p) for p in P]
    m = tuple(sum(p[i] for p in P) / 4 for i in range(3))
    sub = [[m] + [P[j] for j in range(4) if j != i] for i in range(4)]
    M = monos(k); N = len(M); nun = 4 * 3 * N
    idx = lambda s, c, j: (s * 3 + c) * N + j
    rows = []
    def plane_rows(o, u, w, pieces):  # pieces: list of (sub, sign)
        lin = affine(o, [u, w]); ex = [subst(e, lin) for e in M]
        keys = sorted({kk for p in ex for kk in p})
        for c in range(3):
            for kk in keys:
                r = {}
                for (s, sg) in pieces:
                    for j in range(N):
                        if kk in ex[j]: r[idx(s, c, j)] = r.get(idx(s, c, j), 0) + sg * ex[j][kk]
                if r: rows.append(r)
    vsub = lambda a, b: tuple(x - y for x, y in zip(a, b))
    # continuity across internal faces conv(m,Pi,Pj), shared by S_k,S_l with {k,l} = complement
    for i in range(4):
        for j in range(i + 1, 4):
            kl = [t for t in range(4) if t not in (i, j)]
            plane_rows(m, vsub(P[i], m), vsub(P[j], m), [(kl[0], 1), (kl[1], -1)])
    # zero on outer faces of S_1,S_2,S_3 (faces of T through the apex P0)
    for s in (1, 2, 3):
        f = [P[j] for j in range(4) if j != s]
        plane_rows(f[0], vsub(f[1], f[0]), vsub(f[2], f[0]), [(s, 1)])
    # divergence
    Md = monos(k - 1); pos = {e: t for t, e in enumerate(Md)}
    for s in range(4):
        for t in range(len(Md)):
            r = {}
            for c in range(3):
                for j, e in enumerate(M):
                    if e[c] > 0:
                        ee = list(e); ee[c] -= 1
                        if pos[tuple(ee)] == t: r[idx(s, c, j)] = r.get(idx(s, c, j), 0) + e[c]
            if r: rows.append(r)
    Mat = DomainMatrix([[QQ(r.get(col, 0).numerator, r.get(col, 0).denominator) if col in r else QQ(0) for col in range(nun)] for r in rows], (len(rows), nun), QQ)
    K = Mat.nullspace().transpose()   # nun x d
    d = K.shape[1]
    Kl = [[Fr(int(K[i, j].element.numerator), int(K[i, j].element.denominator)) for j in range(d)] for i in range(nun)]
    # quadratic forms on full monomial basis
    def tet_table(V, deg):
        o = V[0]; vecs = [vsub(V[t], o) for t in (1, 2, 3)]
        vol6 = abs(det3(*vecs)); lin = affine(o, vecs); tab = {}
        for e in monos(deg):
            tab[e] = vol6 * sum(c * ref_int(ee) for ee, c in subst(e, lin).items())
        return tab
    def grad_basis(j):  # derivative of monomial j wrt each coord: list of (coef, exponent)
        e = M[j]; out = []
        for a in range(3):
            if e[a] == 0: out.append(None)
            else:
                ee = list(e); ee[a] -= 1; out.append((e[a], tuple(ee)))
        return out
    GB = [grad_basis(j) for j in range(N)]
    A = [[Fr(0)] * nun for _ in range(nun)]
    for s in range(4):
        tab = tet_table(sub[s], 2 * k - 2)
        # Dv:Dv/2 = sum_{a,b} (d_b v_a)^2 + sum_{a,b} d_b v_a d_a v_b
        for c1 in range(3):
            for j1 in range(N):
                I1 = idx(s, c1, j1)
                for c2 in range(3):
                    for j2 in range(N):
                        I2 = idx(s, c2, j2); val = Fr(0)
                        if c1 == c2:
                            for b in range(3):
                                g1, g2 = GB[j1][b], GB[j2][b]
                                if g1 and g2: val += g1[0] * g2[0] * tab[tuple(x + y for x, y in zip(g1[1], g2[1]))]
                        # cross term d_{c2} v_{c1} * d_{c1} v_{c2}
                        g1, g2 = GB[j1][c2], GB[j2][c1]
                        if g1 and g2: val += g1[0] * g2[0] * tab[tuple(x + y for x, y in zip(g1[1], g2[1]))]
                        if val: A[I1][I2] += val
    # face F in z=0: triangle P1,P2,P3
    o = P[1]; vecs = [vsub(P[2], o), vsub(P[3], o)]
    area2 = abs(vecs[0][0] * vecs[1][1] - vecs[0][1] * vecs[1][0]); lin = affine(o, vecs)
    ftab = {}
    for e in monos(2 * k):
        if e[2] == 0: ftab[e] = area2 * sum(c * ref_int(ee) for ee, c in subst(e, lin).items())
    B = [[Fr(0)] * nun for _ in range(nun)]; C = [[Fr(0)] * nun for _ in range(nun)]
    for c in range(3):
        for j1 in range(N):
            e1 = M[j1]
            if e1[2] > 0: continue  # vanishes on z=0
            for j2 in range(N):
                e2 = M[j2]
                ee = tuple(a + b for a, b in zip(e1, e2))
                if e2[2] == 0: C[idx(0, c, j1)][idx(0, c, j2)] += ftab[ee]
                if e2[2] == 1:  # -d_z(monomial j2) on z=0 times monomial j1
                    val = -ftab[(ee[0], ee[1], 0)]
                    B[idx(0, c, j1)][idx(0, c, j2)] += val / 2; B[idx(0, c, j2)][idx(0, c, j1)] += val / 2
    def red(Q):
        nz = [(i, j, Q[i][j]) for i in range(nun) for j in range(nun) if Q[i][j] != 0]
        R = [[Fr(0)] * d for _ in range(d)]
        # R = K^T Q K
        QK = {}
        for i, j, q in nz:
            row = QK.setdefault(i, [Fr(0)] * d)
            for b in range(d):
                if Kl[j][b]: row[b] += q * Kl[j][b]
        for i, row in QK.items():
            for a in range(d):
                if Kl[i][a]:
                    for b in range(d): R[a][b] += Kl[i][a] * row[b]
        return R
    diamF = max(float(np.linalg.norm(np.array(P[a], float) - np.array(P[b], float))) for a in (1, 2, 3) for b in (1, 2, 3))
    build.last = (Kl, M, sub, P, N)
    return d, len(rows), red(A), red(B), red(C), diamF

def lam_float(A, B, C):
    A, B, C = (np.array(X, float) for X in (A, B, C))
    Q = 2 * B - A
    lo, hi = -1e3, 1e4
    f = lambda t: np.linalg.eigvalsh(Q - t * C)[-1]
    if f(lo) <= 0: return None, None
    for _ in range(200):
        mid = (lo + hi) / 2
        (lo, hi) = (mid, hi) if f(mid) > 0 else (lo, mid)
    w, V = np.linalg.eigh(Q - (lo - 1e-9) * C)
    return lo, V[:, -1]

def certify(A, B, C, vec, c):
    q = [Fr(float(x)).limit_denominator(10 ** 6) for x in vec]; d = len(q)
    val = sum(q[a] * (2 * B[a][b] - A[a][b] - c * C[a][b]) * q[b] for a in range(d) for b in range(d))
    Cv = sum(q[a] * C[a][b] * q[b] for a in range(d) for b in range(d))
    Rq = sum(q[a] * (2 * B[a][b] - A[a][b]) * q[b] for a in range(d) for b in range(d)) / Cv
    return val > 0, Rq

TETS = {
 "near-regular": [(Fr(1, 2), Fr(7, 24), Fr(13, 16)), (0, 0, 0), (1, 0, 0), (Fr(1, 2), Fr(7, 8), 0)],
 "right-corner": [(0, 0, 1), (0, 0, 0), (1, 0, 0), (0, 1, 0)],
 "flat(h=1/2)":  [(Fr(1, 2), Fr(7, 24), Fr(1, 2)), (0, 0, 0), (1, 0, 0), (Fr(1, 2), Fr(7, 8), 0)],
 "flat(h=1/4)":  [(Fr(1, 2), Fr(7, 24), Fr(1, 4)), (0, 0, 0), (1, 0, 0), (Fr(1, 2), Fr(7, 8), 0)],
 "tall(h=2)":    [(Fr(1, 2), Fr(7, 24), 2), (0, 0, 0), (1, 0, 0), (Fr(1, 2), Fr(7, 8), 0)],
 "tall(h=3)":    [(Fr(1, 2), Fr(7, 24), 3), (0, 0, 0), (1, 0, 0), (Fr(1, 2), Fr(7, 8), 0)],
 "skew":         [(Fr(3, 2), Fr(1, 3), Fr(3, 4)), (0, 0, 0), (1, 0, 0), (Fr(1, 2), Fr(7, 8), 0)],
}

if __name__ == "__main__":
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    names = sys.argv[2:] or list(TETS)
    for nm in names:
        d, nr, A, B, C, dF = build(TETS[nm], k)
        L, vec = lam_float(A, B, C)
        if L is None:
            print(f"k={k} {nm}: dim={d}, pencil sup (2B-A)/C = -inf?"); continue
        lamT = L * dF
        c = Fr(float(L) * 0.999).limit_denominator(1000) if L > 0 else Fr(float(L) * 1.001).limit_denominator(1000)
        ok, Rq = certify(A, B, C, vec, c)
        print(f"k={k} {nm:14s} dim W={d:3d}  lam*(coords)={L:.6f}  diamF={dF:.4f}  lambda_T=diamF*lam*={lamT:.5f}  "
              f"exact cert 2B-A-({c})C>0 at rational witness: {ok} (witness quotient {float(Rq):.6f})", flush=True)
