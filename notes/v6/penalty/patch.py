"""
Float lambda* for a general triangulated patch (notes/v6/penalty).

Patch = vertices P (n x 2), triangles (index triples), a set of Nitsche edges (vertex pairs) with a
given outward normal each; every other boundary edge of the patch is Dirichlet (v = 0).  Space:
continuous P_k vector fields on the patch, zero on the Dirichlet edges, pointwise divergence-free.
  A = int 1/2 Dv:Dv,  B = int_{Nitsche} (d_n v).v,  C = int_{Nitsche} |v|^2
  lambda* = sup (2B - A)/C   (C > 0 part; on ker C, 2B - A = -A < 0).
Returns (lambda*, dim Z).  Uses the Schur complement onto range(C) as in code/lamstar.py.
"""
import numpy as np
import itertools

_cache = {}


def ref_data(k):
    if k in _cache:
        return _cache[k]
    mons = [(i, j) for i in range(k + 1) for j in range(k + 1 - i)]
    # lagrange nodes on reference triangle (0,0),(1,0),(0,1) indexed by barycentric multi-index
    nodes = []
    for i in range(k + 1):
        for j in range(k + 1 - i):
            nodes.append((k - i - j, i, j))      # (l0,l1,l2) * k ; point = (i/k, j/k)
    pts = np.array([[n[1] / k, n[2] / k] for n in nodes])
    V = np.array([[x ** a * y ** b for (a, b) in mons] for x, y in pts])
    Ci = np.linalg.inv(V)                          # coeffs of basis functions (columns)
    # quadrature on reference triangle (collapsed Gauss), exact to high degree
    g, w = np.polynomial.legendre.leggauss(k + 3)
    Q, W = [], []
    for u, wu in zip(g, w):
        for s, ws in zip(g, w):
            xx = (u + 1) / 2
            yy = (s + 1) / 2 * (1 - xx)
            Q.append((xx, yy)); W.append(wu * ws * (1 - xx) / 4)
    Q = np.array(Q); W = np.array(W)

    def evalmon(X):
        x, y = X[:, 0], X[:, 1]
        M = np.stack([x ** a * y ** b for (a, b) in mons], 1)
        Mx = np.stack([a * x ** max(a - 1, 0) * y ** b if a else 0 * x for (a, b) in mons], 1)
        My = np.stack([b * x ** a * y ** max(b - 1, 0) if b else 0 * x for (a, b) in mons], 1)
        return M @ Ci, Mx @ Ci, My @ Ci
    phiQ, phixQ, phiyQ = evalmon(Q)
    # divergence collocation points: P_{k-1} lattice (interior-shifted to be safe: use k-1 lattice scaled)
    dn = []
    for i in range(k):
        for j in range(k - i):
            dn.append(((i + 1 / 3) / (k + 0.0), (j + 1 / 3) / (k + 0.0)))
    dn = np.array(dn)
    _, dphix, dphiy = evalmon(dn)
    ge, we = np.polynomial.legendre.leggauss(k + 3)
    te = (ge + 1) / 2
    wedge = we / 2
    d = dict(k=k, nodes=nodes, Q=Q, W=W, phiQ=phiQ, phixQ=phixQ, phiyQ=phiyQ, dphix=dphix, dphiy=dphiy,
             evalmon=evalmon, te=te, we=wedge)
    _cache[k] = d
    return d


def lamstar(P, tris, nitsche, k=4, return_vec=False, tol=1e-10, c1vert=False):
    """nitsche: dict {(i,j): n (outward unit normal)}  (unordered vertex pair)."""
    R = ref_data(k)
    P = np.asarray(P, float)
    # global node numbering by combinatorial key
    key2id = {}
    tri_ids = []
    for T in tris:
        ids = []
        for (l0, l1, l2) in R["nodes"]:
            lam = {T[0]: l0, T[1]: l1, T[2]: l2}
            key = tuple(sorted((v, c) for v, c in lam.items() if c > 0))
            if key not in key2id:
                key2id[key] = len(key2id)
            ids.append(key2id[key])
        tri_ids.append(ids)
    nn = len(key2id)
    # edges count
    ecount = {}
    for T in tris:
        for a_, b_ in ((T[0], T[1]), (T[1], T[2]), (T[0], T[2])):
            e = tuple(sorted((a_, b_)))
            ecount[e] = ecount.get(e, 0) + 1
    bedges = [e for e, c in ecount.items() if c == 1]
    nit = {tuple(sorted(e)): np.asarray(n, float) for e, n in nitsche.items()}
    dir_edges = [e for e in bedges if e not in nit]
    fixed = set()
    for key, gid in key2id.items():
        vs = [v for v, c in key]
        if len(vs) <= 2:
            # node lies on edge (or is vertex); fixed if on a dirichlet edge
            for e in dir_edges:
                if set(vs) <= set(e):
                    fixed.add(gid)
    free = [i for i in range(nn) if i not in fixed]
    fmap = {g: i for i, g in enumerate(free)}
    nd = 2 * len(free)
    A = np.zeros((nd, nd)); B = np.zeros((nd, nd)); Cm = np.zeros((nd, nd))
    Drows = []
    for T, ids in zip(tris, tri_ids):
        X0, X1, X2 = P[T[0]], P[T[1]], P[T[2]]
        J = np.column_stack([X1 - X0, X2 - X0])
        detJ = abs(np.linalg.det(J))
        Jit = np.linalg.inv(J).T          # grad_x = Jit @ grad_ref
        gx = Jit[0, 0] * R["phixQ"] + Jit[0, 1] * R["phiyQ"]
        gy = Jit[1, 0] * R["phixQ"] + Jit[1, 1] * R["phiyQ"]
        W = R["W"] * detJ
        loc = [fmap.get(g, -1) for g in ids]
        nb = len(ids)
        # local A: v = sum_a (u_a e1 + w_a e2) phi_a
        # 1/2 Dv:Dv = 2 ux^2 + (uy+wx)^2 + 2 wy^2
        Gxx = (gx * W[:, None]).T @ gx
        Gyy = (gy * W[:, None]).T @ gy
        Gxy = (gx * W[:, None]).T @ gy
        Aloc = np.zeros((2 * nb, 2 * nb))
        Aloc[0::2, 0::2] = 2 * Gxx + Gyy
        Aloc[1::2, 1::2] = Gxx + 2 * Gyy
        Aloc[0::2, 1::2] = Gxy.T        # int uy_a * wx_b
        Aloc[1::2, 0::2] = Gxy
        dof = []
        for l in loc:
            dof += [2 * l, 2 * l + 1] if l >= 0 else [-1, -1]
        dof = np.array(dof)
        m = dof >= 0
        A[np.ix_(dof[m], dof[m])] += Aloc[np.ix_(m, m)]
        # divergence rows
        dx = Jit[0, 0] * R["dphix"] + Jit[0, 1] * R["dphiy"]
        dy = Jit[1, 0] * R["dphix"] + Jit[1, 1] * R["dphiy"]
        scale = np.sqrt(detJ)
        for r in range(dx.shape[0]):
            row = np.zeros(nd)
            for a_, l in enumerate(loc):
                if l >= 0:
                    row[2 * l] += dx[r, a_] * scale
                    row[2 * l + 1] += dy[r, a_] * scale
            Drows.append(row)
        # Nitsche edges of this triangle
        for (ia, ib) in ((0, 1), (1, 2), (0, 2)):
            e = tuple(sorted((T[ia], T[ib])))
            if e in nit:
                n = nit[e]
                Pa, Pb = P[T[ia]], P[T[ib]]
                L = np.linalg.norm(Pb - Pa)
                Xe = Pa[None, :] + R["te"][:, None] * (Pb - Pa)[None, :]
                # reference coords of Xe
                Xr = np.linalg.solve(J, (Xe - X0).T).T
                ph, phx, phy = R["evalmon"](Xr)
                gxe = Jit[0, 0] * phx + Jit[0, 1] * phy
                gye = Jit[1, 0] * phx + Jit[1, 1] * phy
                dnp = n[0] * gxe + n[1] * gye
                we = R["we"] * L
                Mb = (ph * we[:, None]).T @ dnp      # [a,b] = int phi_a d_n phi_b
                Mc = (ph * we[:, None]).T @ ph
                for c in range(2):
                    dd = np.array([2 * l + c if l >= 0 else -1 for l in loc])
                    mm = dd >= 0
                    B[np.ix_(dd[mm], dd[mm])] += Mb[np.ix_(mm, mm)]
                    Cm[np.ix_(dd[mm], dd[mm])] += Mc[np.ix_(mm, mm)]
    if c1vert:
        # Argyris subspace: grad v single-valued at every vertex (=> v = curl psi, psi C^1 P5 with C^2 at vertices)
        vt = {}
        for T, ids in zip(tris, tri_ids):
            for li in range(3):
                vt.setdefault(T[li], []).append((T, ids))
        for vtx, lst in vt.items():
            grads = []
            for (T, ids) in lst:
                X0, X1, X2 = P[T[0]], P[T[1]], P[T[2]]
                J = np.column_stack([X1 - X0, X2 - X0]); Jit = np.linalg.inv(J).T
                Xr = np.linalg.solve(J, (P[vtx] - X0))[None, :]
                ph, phx, phy = R["evalmon"](Xr)
                gxv = (Jit[0, 0] * phx + Jit[0, 1] * phy)[0]
                gyv = (Jit[1, 0] * phx + Jit[1, 1] * phy)[0]
                loc = [fmap.get(g, -1) for g in ids]
                rows4 = []
                for c in range(2):
                    for gg in (gxv, gyv):
                        row = np.zeros(nd)
                        for a_, l in enumerate(loc):
                            if l >= 0:
                                row[2 * l + c] += gg[a_]
                        rows4.append(row)
                grads.append(rows4)
            for q in range(1, len(grads)):
                for r in range(4):
                    Drows.append(grads[q][r] - grads[0][r])
    D = np.array(Drows)
    U, S, Vt = np.linalg.svd(D)
    r = int((S > tol * S[0]).sum())
    Z = Vt[r:].T
    if Z.shape[1] == 0:
        return (-np.inf, 0)
    Bs = (B + B.T) / 2
    Ar, Br, Cr = Z.T @ A @ Z, Z.T @ Bs @ Z, Z.T @ Cm @ Z
    ev, Qv = np.linalg.eigh(Cr)
    keep = ev >= 1e-11 * ev.max()
    kerC, ranC = Qv[:, ~keep], Qv[:, keep]
    M = 2 * Br - Ar
    Myy = ranC.T @ M @ ranC
    Cyy = ranC.T @ Cr @ ranC
    if kerC.shape[1]:
        Mkk = kerC.T @ M @ kerC; Mky = kerC.T @ M @ ranC
        Seff = Myy - Mky.T @ np.linalg.solve(Mkk, Mky)
    else:
        Seff = Myy
    Lc = np.linalg.cholesky(Cyy); Li = np.linalg.inv(Lc)
    evs, EV = np.linalg.eigh(Li @ Seff @ Li.T)
    lam = evs[-1]
    if return_vec:
        return lam, Z.shape[1], dict(A=Ar, B=Br, C=Cr, Z=Z, sv=S)
    return lam, Z.shape[1]


# ---------------------------------------------------------------- star builders
def fan(rim, k=4):
    """fan at O=(0,0); rim points listed from (x>0, y=0) counterclockwise to (x<0, y=0).
    Nitsche edges: O-rim[0] and O-rim[-1], normal (0,-1) (fluid above)."""
    P = [(0.0, 0.0)] + [tuple(p) for p in rim]
    m = len(rim) - 1
    tris = [(0, i + 1, i + 2) for i in range(m)]
    nit = {(0, 1): (0, -1), (0, m + 1): (0, -1)}
    return lamstar(P, tris, nit, k)


def fan_polar(angles, radii, k=4):
    rim = [(r * np.cos(t), r * np.sin(t)) for t, r in zip(angles, radii)]
    return fan(rim, k)


if __name__ == "__main__":
    # regression against the paper: S3, S4', uniform fans, S_{3/4}
    print("S3", fan([(1, 0), (1, 1), (-1, 1), (-1, 0)]))
    print("S4'", fan([(1, 0), (1, 1), (0, 1.2), (-1, 1), (-1, 0)]))
    print("S34", fan([(1, 0), (1, .75), (0, .75), (-1, 0)]))
    for m in (3, 4, 5):
        ang = np.linspace(0, np.pi, m + 1)
        print("uniform", m, fan_polar(ang, np.ones(m + 1)))
    # single triangle (a,b) = equilateral
    P = [(0, 0), (1, 0), (0.5, np.sqrt(3) / 2)]
    print("tri eq", lamstar(P, [(0, 1, 2)], {(0, 1): (0, -1)}))
