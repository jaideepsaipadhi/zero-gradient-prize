"""
lambda* for the Gjerde-Scott Nitsche form on divergence-free P4 fields.

Reference star at a vertex x=0 of Gamma_h. Locally Gamma_h is the line y=0
(polygon angle -> pi as h->0), fluid in y>0, outward normal n=(0,-1).
The star is a fan of m triangles over the upper half-disk.

Space: continuous P4 vector fields on the star, zero on the outer chords,
FREE on the boundary line, with div v = 0 pointwise.

  A = 1/2 ||D(v)||^2,   B = <d_n v, v>_{Gamma},   C = ||v||^2_{Gamma}
Scale invariance (2D):  N_h(v,v) = A - 2B + (mu/rho) C  for a star of size h/rho.
  lambda* = sup (2B - A)/C.
If lambda* > 0, N_h is indefinite on Z_h whenever mu < rho*lambda*.
"""
import numpy as np, sys

DEG = 4
MONS = [(i, j) for i in range(DEG + 1) for j in range(DEG + 1 - i)]
g, w = np.polynomial.legendre.leggauss(10)


def lagr_nodes(T):
    P = []
    for i in range(DEG + 1):
        for j in range(DEG + 1 - i):
            k = DEG - i - j
            P.append((i * T[0] + j * T[1] + k * T[2]) / DEG)
    return np.array(P)


def mono(x, y):
    return np.array([x**a * y**b for (a, b) in MONS])


def dmono(x, y):
    dx = np.array([a * x**(a - 1) * y**b if a > 0 else 0.0 for (a, b) in MONS])
    dy = np.array([b * x**a * y**(b - 1) if b > 0 else 0.0 for (a, b) in MONS])
    return dx, dy


def triquad(T):
    pts, wts = [], []
    area = abs(np.cross(T[1] - T[0], T[2] - T[0]))
    for a, wa in zip(g, w):
        for b, wb in zip(g, w):
            u = (a + 1) / 2
            v = (b + 1) / 2 * (1 - u)
            pts.append(T[0] + u * (T[1] - T[0]) + v * (T[2] - T[0]))
            wts.append(wa * wb * (1 - u) / 4 * area)
    return np.array(pts), np.array(wts)


def linequad(P, Q):
    L = np.linalg.norm(Q - P)
    return np.array([P + (a + 1) / 2 * (Q - P) for a in g]), w * L / 2


def on_seg(p, A, B, tol=1e-9):
    if abs(np.cross(B - A, p - A)) > tol:
        return False
    t = np.dot(p - A, B - A) / np.dot(B - A, B - A)
    return -tol <= t <= 1 + tol


def run(angles):
    O = np.zeros(2)
    rim = [np.array([np.cos(t), np.sin(t)]) for t in angles]
    tris = [np.array([O, rim[i], rim[i + 1]]) for i in range(len(rim) - 1)]
    nodes, loc = {}, []
    for T in tris:
        N = lagr_nodes(T)
        V = np.array([mono(*p) for p in N])
        C = np.linalg.inv(V)          # column k = coeffs of basis fn k
        ids = []
        for p in N:
            key = (round(p[0], 9), round(p[1], 9))
            nodes.setdefault(key, len(nodes))
            ids.append(nodes[key])
        loc.append((T, N, C, ids))
    nn = len(nodes)
    nd = 2 * nn
    keys = [np.array(k) for k in nodes]
    fixed = {i for i, p in enumerate(keys)
             if any(on_seg(p, rim[j], rim[j + 1]) for j in range(len(rim) - 1))}
    free = [i for i in range(nn) if i not in fixed]

    A = np.zeros((nd, nd)); B = np.zeros((nd, nd)); Cm = np.zeros((nd, nd))
    Drows = []
    n = np.array([0.0, -1.0])
    for (T, N, C, ids) in loc:
        X, W = triquad(T)
        for x, wt in zip(X, W):
            dx, dy = dmono(*x)
            gx, gy = dx @ C, dy @ C
            # D(e_c phi) = grad(phi) e_c^T + e_c grad(phi)^T
            for a, ia in enumerate(ids):
                for c in range(2):
                    Ga = np.zeros((2, 2)); Ga[c] = [gx[a], gy[a]]
                    Da = Ga + Ga.T
                    for b, ib in enumerate(ids):
                        for d in range(2):
                            Gb = np.zeros((2, 2)); Gb[d] = [gx[b], gy[b]]
                            Db = Gb + Gb.T
                            A[2*ia+c, 2*ib+d] += 0.5 * np.sum(Da * Db) * wt
        # pointwise divergence-free: div is P3 on T; impose at 10 P3 nodes
        for i in range(4):
            for j in range(4 - i):
                k = 3 - i - j
                x = (i*T[0] + j*T[1] + k*T[2]) / 3
                dx, dy = dmono(*x)
                gx, gy = dx @ C, dy @ C
                row = np.zeros(nd)
                for a, ia in enumerate(ids):
                    row[2*ia] += gx[a]; row[2*ia+1] += gy[a]
                Drows.append(row)
        # boundary edges: those of T lying on y=0 (O-rim[0] and O-rim[-1])
        for P in (rim[0], rim[-1]):
            if any(np.allclose(P, T[k]) for k in (1, 2)):
                X2, W2 = linequad(O, P)
                for x, wt in zip(X2, W2):
                    m = mono(*x) @ C
                    dx, dy = dmono(*x)
                    dn = (dx*n[0] + dy*n[1]) @ C
                    for a, ia in enumerate(ids):
                        for b, ib in enumerate(ids):
                            for c in range(2):
                                B[2*ia+c, 2*ib+c] += dn[b] * m[a] * wt   # <d_n u, v>, u=b
                                Cm[2*ia+c, 2*ib+c] += m[a] * m[b] * wt
    fd = [2*i + c for i in free for c in range(2)]
    A = A[np.ix_(fd, fd)]; B = B[np.ix_(fd, fd)]; Cm = Cm[np.ix_(fd, fd)]
    D = np.array(Drows)[:, fd]
    U, S, Vt = np.linalg.svd(D)
    r = int((S > 1e-10 * S[0]).sum())
    Z = Vt[r:].T
    if Z.shape[1] == 0:
        return dict(dimZ=0)
    Bs = (B + B.T) / 2
    Ar, Br, Cr = Z.T @ A @ Z, Z.T @ Bs @ Z, Z.T @ Cm @ Z
    # restrict to complement of ker C (fields with zero trace): there 2B-A = -A <= 0
    # lambda* = max generalized eig of (2B - A, C) on range(C); if 2B-A>0 on ker C -> inf
    ev, Q = np.linalg.eigh(Cr)
    kerC = Q[:, ev < 1e-12 * ev.max()]
    ranC = Q[:, ev >= 1e-12 * ev.max()]
    M = 2*Br - Ar
    # Schur-complement: sup over v = y + k (y in ranC, k in kerC) of y^T.../y^T C y
    # maximize over k first: M_kk = -A_kk (neg def) -> k = -M_kk^{-1} M_ky
    Mkk = kerC.T @ M @ kerC; Mky = kerC.T @ M @ ranC; Myy = ranC.T @ M @ ranC
    Cyy = ranC.T @ Cr @ ranC
    if Mkk.size:
        Seff = Myy - Mky.T @ np.linalg.solve(Mkk, Mky)
    else:
        Seff = Myy
    L = np.linalg.cholesky(Cyy)
    Li = np.linalg.inv(L)
    lam = np.linalg.eigvalsh(Li @ Seff @ Li.T).max()
    # sanity: Nitsche trace-inverse constant  sup <d_n v, d_n v>/A  is finite
    return dict(dimZ=Z.shape[1], lam_star=lam)


if __name__ == "__main__":
    for m in (2, 3, 4, 5, 6):
        ang = np.linspace(0, np.pi, m + 1)
        print(f"m={m} uniform fan:", run(ang), flush=True)
    rng = np.random.default_rng(0)
    for trial in range(4):
        m = 4
        inner = np.sort(rng.uniform(0.3, np.pi - 0.3, m - 1))
        ang = np.concatenate([[0], inner, [np.pi]])
        print(f"m={m} random fan {np.round(ang,2)}:", run(ang), flush=True)
