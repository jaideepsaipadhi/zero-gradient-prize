"""(IS3) constant on icosphere-Alfeld meshes, NUMERICAL.

Omega_h = shell between the icosphere polyhedron P_l (level l, vertices on the unit sphere) and R*P_l (R=2),
m layers of prisms (radii 2^(j/m)) split into 3 tets each (min-index rule), then Alfeld-refined.
beta_h = inf_{q in P_{k-1}^disc cap L^2_0} sup_{v in V_0} (div v, q)/(|grad v| |q|),  V_0 = P_k velocities vanishing on
all of dOmega_h (the space of (IS3)).  beta_h^{-2} = largest eigenvalue of M^{1/2} S^+ M^{1/2}, S = B A^{-1} B^T,
computed with ARPACK; each application is one sparse LU solve of the saddle-point system with a mean-zero multiplier.
usage: python3 is3_ico.py k level m [level m ...]
"""
import os, sys, time
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../code3d'))
import svn3d
import infsup3d  # noqa: F401  (documentation: same A, B as the dense code)

R_OUT = 2.0


def icosphere(level):
    t = (1 + 5 ** 0.5) / 2
    V = np.array([[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
                  [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]], float)
    F = [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6],
         [7, 1, 8], [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10],
         [8, 6, 7], [9, 8, 1]]
    V = list(V / np.linalg.norm(V, axis=1)[:, None]); F = [list(f) for f in F]
    for _ in range(level):
        cache = {}
        def mid(a, b):
            key = (min(a, b), max(a, b))
            if key not in cache:
                p = V[a] + V[b]; V.append(p / np.linalg.norm(p)); cache[key] = len(V) - 1
            return cache[key]
        NF = []
        for a, b, c in F:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            NF += [[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]
        F = NF
    return np.array(V), np.array(F)


def shell_mesh(level, m):
    S, tri = icosphere(level)
    ns = len(S)
    V = np.concatenate([S * R_OUT ** (j / m) for j in range(m + 1)])
    T = np.concatenate([svn3d._split_prisms(j * ns + np.arange(ns), (j + 1) * ns + np.arange(ns), tri) for j in range(m)])
    T = svn3d._orient(V, T)
    return V, T, tri


class IcoSpace(svn3d.Space):
    """svn3d.Space with the outer boundary R*P_l instead of the box."""
    def __init__(self, V, T, R):
        # reuse the parent constructor on a fake info, catching only the box/sphere classification
        self.R, self.V, self.T, self.info = R, V, T, dict(L=np.inf)
        k, nl, nt = R.k, R.nl, len(T)
        lam = np.array(R.LBARY)
        gv = np.broadcast_to(T[:, None, :], (nt, nl, 4)).astype(np.int64)
        lm = np.broadcast_to(lam[None], (nt, nl, 4)).astype(np.int64)
        gvm = np.where(lm > 0, gv, np.int64(1) << 40)
        o = np.argsort(gvm, axis=2)
        gk = np.take_along_axis(gvm, o, 2); lk = np.take_along_axis(np.where(lm > 0, lm, 0), o, 2)
        key = np.concatenate([gk, lk], 2).reshape(-1, 8)
        _, inv = np.unique(key, axis=0, return_inverse=True)
        self.ids = inv.reshape(nt, nl); self.nn = self.ids.max() + 1
        P = V[T]
        self.P0 = P[:, 0]
        self.J = np.stack([P[:, 1] - P[:, 0], P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]], 2)
        self.det = np.linalg.det(self.J); assert (self.det > 0).all()
        self.Jinv = np.linalg.inv(self.J)
        FL = [(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)]
        allf = np.concatenate([np.sort(T[:, f], 1) for f in FL])
        _, fi, cnt = np.unique(allf, axis=0, return_inverse=True, return_counts=True)
        bnd = cnt[fi.ravel()] == 1
        tt = np.tile(np.arange(nt), 4); ii = np.repeat(np.arange(4), nt)
        bt, bi = tt[bnd], ii[bnd]
        FLa = np.array(FL)
        FV = np.stack([V[T[bt, FLa[bi, j]]] for j in range(3)], 1)
        rr = np.linalg.norm(FV, axis=2)
        onsph = np.all(np.abs(rr - 1) < 1e-9, 1); onout = np.all(np.abs(rr - R_OUT) < 1e-9, 1)
        assert np.all(onsph ^ onout)
        self.fGamma = (bt[onsph], bi[onsph]); self.fBox = (bt[onout], bi[onout])
        gF = FV[onsph]
        self.hGamma = max(np.linalg.norm(gF[:, i] - gF[:, j], axis=1).max() for i, j in ((0, 1), (1, 2), (0, 2)))
        E = np.stack([np.linalg.norm(P[:, i] - P[:, j], axis=1) for i in range(4) for j in range(i + 1, 4)], 1)
        self.h = E.max()
        self.dofs = 3 * self.ids[:, :, None] + np.arange(3)[None, None, :]


def beta_h(S):
    R = S.R; nl = R.nl; nt = len(S.T); n3 = 3 * nl
    G = np.einsum('tim,tjn,ijab->tmnab', S.Jinv, S.Jinv, R.K) * S.det[:, None, None, None, None]
    tr = G[:, 0, 0] + G[:, 1, 1] + G[:, 2, 2]
    loc = np.zeros((nt, nl, 3, nl, 3))
    for c in range(3):
        loc[:, :, c, :, c] = tr
    D = S.dofs.reshape(nt, n3)
    A = sp.coo_matrix((loc.ravel(), (np.repeat(D, n3, 1).ravel(), np.tile(D, (1, n3)).ravel())), (3 * S.nn,) * 2).tocsr()
    B = svn3d.assemble(S, 0.0, lambda a, b, c: np.zeros(np.shape(a) + (3,)), [], [])["B"]
    lamb = np.array(R.LBARY); bn = set()
    bt = np.concatenate([S.fBox[0], S.fGamma[0]]); bi = np.concatenate([S.fBox[1], S.fGamma[1]])
    for i in range(4):
        ln = np.where(lamb[:, i] == 0)[0]
        bn.update(S.ids[bt[bi == i]][:, ln].ravel().tolist())
    bn = np.array(sorted(bn))
    dD = (3 * bn[:, None] + np.arange(3)).ravel(); free = np.setdiff1d(np.arange(3 * S.nn), dD)
    Af = A[free][:, free].tocsc(); Bf = B[:, free].tocsc()
    npr = B.shape[0]
    # block-diagonal pressure mass and its square root
    Lc = np.linalg.cholesky(R.Mp)                       # Mp_ref = L L^T ; M_t = det_t * Mp_ref
    Mh = sp.block_diag([Lc * np.sqrt(d) for d in S.det], format='csr')        # M = Mh Mh^T
    one = np.ones(npr)
    M1 = sp.block_diag([R.Mp * d for d in S.det], format='csr') @ one
    nv = Af.shape[0]
    K = sp.bmat([[Af, Bf.T, None], [Bf, None, sp.csc_matrix(M1[:, None])], [None, sp.csc_matrix(M1[None, :]), None]],
                format='csc')
    lu = spl.splu(K)
    w1 = Mh.T @ one; w1 /= np.linalg.norm(w1)          # direction of M^{1/2} 1 in the y variables (M^{1/2} := Mh^T)
    def op(y):
        y = y - w1 * (w1 @ y)
        b = Mh @ y                                      # b = M^{1/2} y  (Mh is the lower factor: M = Mh Mh^T)
        rhs = np.concatenate([np.zeros(nv), -b, [0.0]])
        x = lu.solve(rhs); p = x[nv:nv + npr]
        out = Mh.T @ p
        return out - w1 * (w1 @ out)
    O = spl.LinearOperator((npr, npr), matvec=op, dtype=float)
    ev = spl.eigsh(O, k=3, which='LA', tol=1e-8, return_eigenvectors=False)
    ev = np.sort(ev)[::-1]
    return 1 / np.sqrt(ev), nv, npr


if __name__ == "__main__":
    k = int(sys.argv[1])
    args = list(map(int, sys.argv[2:]))
    for level, m in zip(args[::2], args[1::2]):
        t0 = time.time()
        V, T, tri = shell_mesh(level, m)
        Va, Ta = svn3d.alfeld(V, T)
        S = IcoSpace(Va, Ta, svn3d.Ref(k))
        b, nv, npr = beta_h(S)
        print(f"icosphere level {level}, {m} layer(s), R={R_OUT}: macro tets={len(T)}, Gamma faces={len(tri)}, hG={S.hGamma:.4f}, "
              f"h={S.h:.4f}; k={k}: dofs V0={nv}, P={npr}; beta_h (3 smallest inf-sup values) = {np.array2string(b, precision=5)}; "
              f"1/beta_h = {1/b[0]:.3f}  ({time.time()-t0:.1f}s)", flush=True)
