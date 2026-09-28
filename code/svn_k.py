"""
Degree-general Scott-Vogelius-Nitsche solver: continuous P_k velocity, discontinuous P_{k-1}
pressure (all (k)(k+1)/2 moments per element), k >= 2.  Same forms as svn.py (which is the k = 4
special case and is reproduced to round-off):

  N_h(u,v) - (p, div v) + theta <p - pbar_Gamma(p), v.n>_{Gamma_h} = (f, v)      v = 0 on the square
  (q, div u) = 0                                                                   q in P_{k-1}-disc
  N_h(u,v) = 1/2 (Du, Dv) - <d_n u, v> - <u, d_n v> + (mu/h) <u, v>,  h = h_Omega (max diameter)

Meshes:  "std" = svn.make_mesh (Cavalcante-type structured mesh);
         "ct"  = Clough-Tocher / Alfeld barycentric refinement of the same macro mesh (each triangle
                 split into 3 at its barycentre).  Gamma_h edges, the outer square, and h_Omega (the
                 longest macro edge survives in a micro-triangle) are unchanged.

Quadrature: collapsed Gauss n = max(8, k+4) per direction on triangles, Gauss-Legendre
max(10, k+6) points on Gamma_h edges (identical to svn.py at k = 4).
"""
import os
os.environ.setdefault("PYPARDISO_MKL_RT", "/usr/local/lib/libmkl_rt.so.3")
import numpy as np
import scipy.sparse as sp
import svn                                   # reuse make_mesh and make_exact (unchanged)


# ---------------------------------------------------------------- reference element of degree k
def _mono(M, x, y):
    return np.stack([x**a * y**b for (a, b) in M], axis=-1)


def _dmono(M, x, y):
    dx = np.stack([a * x**max(a - 1, 0) * y**b if a > 0 else 0 * x for (a, b) in M], -1)
    dy = np.stack([b * x**a * y**max(b - 1, 0) if b > 0 else 0 * x for (a, b) in M], -1)
    return dx, dy


def tri_quad(n):
    g, w = np.polynomial.legendre.leggauss(n)
    g = (g + 1) / 2; w = w / 2
    U, V = np.meshgrid(g, g, indexing="ij"); WU, WV = np.meshgrid(w, w, indexing="ij")
    x = U.ravel(); y = (V * (1 - U)).ravel(); wt = (WU * WV * (1 - U)).ravel()
    return x, y, wt


class Ref:
    def __init__(self, k):
        assert k >= 2
        self.k = k
        self.MONV = [(a, b) for a in range(k + 1) for b in range(k + 1 - a)]
        self.MONP = [(a, b) for a in range(k) for b in range(k - a)]
        self.LBARY = [(i, j, k - i - j) for i in range(k + 1) for j in range(k + 1 - i)]
        self.nl = len(self.MONV); self.npl = len(self.MONP)
        self.REFNODES = np.array([[l1 / k, l2 / k] for (l1, l2, l0) in self.LBARY])
        self.COEF = np.linalg.inv(_mono(self.MONV, self.REFNODES[:, 0], self.REFNODES[:, 1]))
        self.QX, self.QY, self.QW = tri_quad(max(8, k + 4))
        self.PHI_Q = self.basis(self.QX, self.QY)
        self.QP_Q = _mono(self.MONP, self.QX, self.QY)
        gl, glw = np.polynomial.legendre.leggauss(max(10, k + 6))
        self.GL = (gl + 1) / 2; self.GLW = glw / 2

    def basis(self, x, y):
        return _mono(self.MONV, x, y) @ self.COEF

    def dbasis(self, x, y):
        dx, dy = _dmono(self.MONV, x, y)
        return dx @ self.COEF, dy @ self.COEF

    def pbasis(self, x, y):
        return _mono(self.MONP, x, y)


# ---------------------------------------------------------------- meshes
def barycentric_refine(verts, tris, circ, outer):
    """Alfeld / Clough-Tocher split: each triangle -> 3 triangles sharing its barycentre."""
    nv = len(verts)
    cen = verts[tris].mean(1)
    V = np.vstack([verts, cen])
    c = nv + np.arange(len(tris))
    T = np.concatenate([np.stack([tris[:, 0], tris[:, 1], c], 1),
                        np.stack([tris[:, 1], tris[:, 2], c], 1),
                        np.stack([tris[:, 2], tris[:, 0], c], 1)])
    P = V[T]
    ar = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    T[ar < 0] = T[ar < 0][:, [0, 2, 1]]
    return V, T, circ, outer


def make_mesh(N, mesh="std", L=2.5, layers=None):
    v, t, c, o = svn.make_mesh(N, L=L, layers=layers)
    if mesh == "ct":
        return barycentric_refine(v, t, c, o)
    assert mesh == "std"
    return v, t, c, o


class Space:
    def __init__(self, verts, tris, circ, outer, R):
        self.R = R
        self.verts, self.tris = verts, tris
        nt = len(tris); nl = R.nl
        edge_id = {}; ids = np.zeros((nt, nl), dtype=np.int64)
        nv = len(verts); nxt = nv
        for t in range(nt):
            g = tris[t]
            for kk, (l1, l2, l0) in enumerate(R.LBARY):
                lam = {0: l0, 1: l1, 2: l2}
                nz = [v for v in range(3) if lam[v] > 0]
                if len(nz) == 1:
                    ids[t, kk] = g[nz[0]]
                elif len(nz) == 2:
                    va, vb = nz; ga, gb = g[va], g[vb]
                    if ga > gb:
                        ga, gb, va, vb = gb, ga, vb, va
                    key = (ga, gb, lam[va])
                    if key not in edge_id:
                        edge_id[key] = nxt; nxt += 1
                    ids[t, kk] = edge_id[key]
                else:
                    ids[t, kk] = nxt; nxt += 1
        # vertices never referenced (none for these meshes) would be harmless zero rows
        self.nn = nxt; self.ids = ids
        P = verts[tris]
        self.P0 = P[:, 0]; self.J = np.stack([P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]], axis=2)
        self.detJ = np.linalg.det(self.J)
        self.Jinv = np.linalg.inv(self.J)
        xy = np.zeros((self.nn, 2))
        xy[ids.ravel()] = (self.P0[:, None, :] + np.einsum('tij,qj->tqi', self.J, R.REFNODES)).reshape(-1, 2)
        self.xy = xy
        self.hdiam = np.max(np.stack([np.linalg.norm(P[:, i] - P[:, (i + 1) % 3], axis=1)
                                      for i in range(3)], 1), 1)
        self.h = self.hdiam.max()
        be = []
        for t in range(nt):
            for (u, w) in [(0, 1), (1, 2), (2, 0)]:
                if tris[t, u] in circ and tris[t, w] in circ:
                    be.append((t, u, w))
        self.bedges = be
        L = np.abs(verts).max()
        on = (np.abs(np.abs(xy[:, 0]) - L) < 1e-10) | (np.abs(np.abs(xy[:, 1]) - L) < 1e-10)
        self.dnodes = np.where(on)[0]
        self.hGamma = max(np.linalg.norm(verts[tris[t, u]] - verts[tris[t, w]]) for (t, u, w) in be)

    def grads(self, x, y):
        dX, dY = self.R.dbasis(x, y)
        Ji = self.Jinv
        gx = np.einsum('t,qb->tqb', Ji[:, 0, 0], dX) + np.einsum('t,qb->tqb', Ji[:, 1, 0], dY)
        gy = np.einsum('t,qb->tqb', Ji[:, 0, 1], dX) + np.einsum('t,qb->tqb', Ji[:, 1, 1], dY)
        return gx, gy

    def phys(self, x, y):
        return self.P0[:, None, :] + np.einsum('tij,qj->tqi', self.J, np.stack([x, y], 1))

    def edge_data(self):
        R = self.R
        out = []
        Rr = np.array([[0, 0], [1, 0], [0, 1]], float)
        for (t, u, w) in self.bedges:
            ref = Rr[u][None, :] + R.GL[:, None] * (Rr[w] - Rr[u])[None, :]
            A = self.verts[self.tris[t, u]]; B = self.verts[self.tris[t, w]]
            L = np.linalg.norm(B - A)
            tvec = (B - A) / L
            n = np.array([tvec[1], -tvec[0]])
            C3 = self.verts[self.tris[t]].mean(0)
            if np.dot(n, C3 - A) > 0:
                n = -n
            X = A[None, :] + R.GL[:, None] * (B - A)[None, :]
            out.append((t, ref, X, R.GLW * L, n))
        return out


def dofs(ids):
    return np.stack([2 * ids, 2 * ids + 1], -1)


def build(N, k, mesh="std", L=2.5, layers=None):
    R = Ref(k)
    v, t, c, o = make_mesh(N, mesh, L=L, layers=layers)
    return Space(v, t, c, o, R)


def assemble(S, mu, f):
    R = S.R; nl, npl = R.nl, R.npl; n2 = 2 * nl
    nt = len(S.tris)
    gx, gy = S.grads(R.QX, R.QY)
    wq = R.QW[None, :] * np.abs(S.detJ)[:, None]
    Axx = np.einsum('tq,tqa,tqb->tab', wq, gx, gx); Ayy = np.einsum('tq,tqa,tqb->tab', wq, gy, gy)
    Axy = np.einsum('tq,tqa,tqb->tab', wq, gx, gy)
    loc = np.zeros((nt, nl, 2, nl, 2))
    loc[:, :, 0, :, 0] = 2 * Axx + Ayy
    loc[:, :, 1, :, 1] = Axx + 2 * Ayy
    loc[:, :, 0, :, 1] = np.transpose(Axy, (0, 2, 1))
    loc[:, :, 1, :, 0] = Axy
    D = dofs(S.ids).reshape(nt, n2)
    rows = np.repeat(D, n2, axis=1).ravel(); cols = np.tile(D, (1, n2)).ravel()
    ndof = 2 * S.nn
    A = sp.coo_matrix((loc.reshape(nt, n2, n2).ravel(), (rows, cols)), shape=(ndof, ndof)).tocsr()
    del loc
    Bl = np.zeros((nt, npl, nl, 2))
    Bl[:, :, :, 0] = np.einsum('tq,qj,tqa->tja', wq, R.QP_Q, gx)
    Bl[:, :, :, 1] = np.einsum('tq,qj,tqa->tja', wq, R.QP_Q, gy)
    prow = (np.arange(nt)[:, None] * npl + np.arange(npl)[None, :])
    rB = np.repeat(prow, n2, axis=1).ravel(); cB = np.tile(D, (1, npl)).ravel()
    B = sp.coo_matrix((Bl.reshape(nt, npl, n2).ravel(), (rB, cB)), shape=(npl * nt, ndof)).tocsr()
    del Bl, gx, gy
    Ml = np.einsum('tq,qi,qj->tij', wq, R.QP_Q, R.QP_Q)
    Minv = sp.block_diag([np.linalg.inv(m) for m in Ml], format='csr')
    X = S.phys(R.QX, R.QY)
    fv = f(X[..., 0], X[..., 1])
    Fl = np.einsum('tq,qa,tqc->tac', wq, R.PHI_Q, fv)
    F = np.zeros(ndof); np.add.at(F, D.ravel(), Fl.reshape(nt, n2).ravel())
    rN, cN, vN, vP = [], [], [], []
    rC, cC, vC = [], [], []
    gG = np.zeros(npl * nt); perim = 0.0
    for (t, ref, Xe, we, n) in S.edge_data():
        ph = R.basis(ref[:, 0], ref[:, 1])
        dX, dY = R.dbasis(ref[:, 0], ref[:, 1])
        Ji = S.Jinv[t]
        gxe = Ji[0, 0] * dX + Ji[1, 0] * dY; gye = Ji[0, 1] * dX + Ji[1, 1] * dY
        dn = n[0] * gxe + n[1] * gye
        Mphi = np.einsum('g,ga,gb->ab', we, ph, ph)
        Dn = np.einsum('g,ga,gb->ab', we, ph, dn)
        blk = -Dn - Dn.T + (mu / S.h) * Mphi
        d = dofs(S.ids[t])
        for c in range(2):
            dd = d[:, c]
            rN.append(np.repeat(dd, nl)); cN.append(np.tile(dd, nl))
            vN.append(blk.ravel()); vP.append(Mphi.ravel())
        q = R.pbasis(ref[:, 0], ref[:, 1])
        gG[prow[t]] += we @ q; perim += we.sum()
        Cl = np.einsum('g,gj,ga->ja', we, q, ph)
        for c in range(2):
            rC.append(np.repeat(prow[t], nl)); cC.append(np.tile(d[:, c], npl)); vC.append((Cl * n[c]).ravel())
    rN = np.concatenate(rN); cN = np.concatenate(cN)
    Nb = sp.coo_matrix((np.concatenate(vN), (rN, cN)), shape=(ndof, ndof)).tocsr()
    Pb = sp.coo_matrix((np.concatenate(vP), (rN, cN)), shape=(ndof, ndof)).tocsr() / S.h  # d/dmu of Nb
    C = sp.coo_matrix((np.concatenate(vC), (np.concatenate(rC), np.concatenate(cC))), shape=(npl * nt, ndof)).tocsr()
    return dict(A=A + Nb, Avol=A, Pen=Pb, B=B, C=C, Minv=Minv, F=F, Ml=Ml, gGamma=gG / perim)


def coupling(S, sysm, theta):
    B, C = sysm["B"], sysm["C"]
    if not theta:
        return B
    npl = S.R.npl
    e = np.zeros(B.shape[0]); e[0::npl] = 1.0          # MONP[0] = (0,0): constant mode
    phi = C.T @ e
    Cm = (C - sp.csr_matrix(sysm["gGamma"][:, None]) @ sp.csr_matrix(phi[None, :])).tocsr()
    return (B - Cm).tocsr()


def solve_lu(S, sysm, g, theta):
    """Direct solve of the saddle system with SuperLU partial pivoting (as svn.solve_lu)."""
    import scipy.sparse.linalg as spl
    A, B, F = sysm["A"], sysm["B"], sysm["F"]
    Bt = coupling(S, sysm, theta)
    ndof = A.shape[0]
    dD = dofs(S.dnodes).ravel(); free = np.setdiff1d(np.arange(ndof), dD)
    ug = np.zeros(ndof); ug[dD] = g(S.xy[S.dnodes, 0], S.xy[S.dnodes, 1]).ravel()
    K = sp.bmat([[A[free][:, free], -Bt[:, free].T], [B[:, free], None]], format="csc")
    rhs = np.concatenate([F[free] - A[free][:, dD] @ ug[dD], -B[:, dD] @ ug[dD]])
    lu = spl.splu(K, permc_spec="COLAMD", diag_pivot_thresh=1.0)
    x = lu.solve(rhs)
    for _ in range(3):
        x = x + lu.solve(rhs - K @ x)
    nnzLU = lu.L.nnz + lu.U.nnz
    del lu
    u = ug.copy(); u[free] = x[:len(free)]; p = x[len(free):]
    Bu = B @ u
    dres = np.sqrt(abs(Bu @ (sysm["Minv"] @ Bu)))
    relres = np.linalg.norm(K @ x - rhs) / np.linalg.norm(rhs)
    return u, p, [dres], relres, nnzLU


def solve_pardiso(S, sysm, g, theta, eps=1e-12, maxref=10):
    """Same saddle system as solve_lu, factorized with MKL PARDISO (unsymmetric, mtype 11, scaling +
    weighted matching, METIS ordering) plus iterative refinement.  Far less fill than SuperLU/COLAMD;
    used for the large cases.  The saddle matrix is ill-conditioned (a referee check found the
    divergence matrix B full row rank on the free dofs, smallest relative singular values 1.7e-2 ... 1.3e-4;
    an earlier comment here wrongly called it exactly singular).  Plain PARDISO + refinement diverged for
    theta = 0, so for robustness the pressure block is regularised by +eps*M (M = pressure mass matrix).  At eps = 1e-12 the velocity
    agrees with solve_lu to ~3e-10 (relative, H1) at k = 6, N = 32; relres is reported for the
    UNregularised system."""
    import pypardiso
    A, B, F = sysm["A"], sysm["B"], sysm["F"]
    Bt = coupling(S, sysm, theta)
    M = sp.block_diag(list(sysm["Ml"]), format="csr")
    ndof = A.shape[0]
    dD = dofs(S.dnodes).ravel(); free = np.setdiff1d(np.arange(ndof), dD)
    ug = np.zeros(ndof); ug[dD] = g(S.xy[S.dnodes, 0], S.xy[S.dnodes, 1]).ravel()
    Aff = A[free][:, free]; Btf = Bt[:, free]; Bf = B[:, free]
    K = sp.bmat([[Aff, -Btf.T], [Bf, eps * M]], format="csr"); K.sort_indices()
    rhs = np.concatenate([F[free] - A[free][:, dD] @ ug[dD], -B[:, dD] @ ug[dD]])
    slv = pypardiso.PyPardisoSolver(mtype=11)
    slv.set_iparm(11, 1); slv.set_iparm(13, 1)
    slv.factorize(K)
    x = slv.solve(K, rhs)
    nb = np.linalg.norm(rhs)
    for _ in range(maxref):
        res = rhs - K @ x
        if np.linalg.norm(res) < 1e-15 * nb:
            break
        x = x + slv.solve(K, res)
    slv.free_memory(everything=True)
    del K
    nf = len(free)
    u = ug.copy(); u[free] = x[:nf]; p = x[nf:]
    Bu = B @ u
    dres = np.sqrt(abs(Bu @ (sysm["Minv"] @ Bu)))
    r0 = np.concatenate([Aff @ x[:nf] - Btf.T @ p, Bf @ x[:nf]]) - rhs
    relres = np.linalg.norm(r0) / nb
    return u, p, [dres], relres, -1


def errors(S, u, ex, mu):
    R = S.R
    gx, gy = S.grads(R.QX, R.QY)
    wq = R.QW[None, :] * np.abs(S.detJ)[:, None]
    X = S.phys(R.QX, R.QY)
    D = dofs(S.ids)
    U = u[D]
    uhx = np.einsum('tqa,tac->tqc', gx, U); uhy = np.einsum('tqa,tac->tqc', gy, U)
    del gx, gy
    G = ex["gu"](X[..., 0], X[..., 1])
    e2 = ((G[..., 0] - uhx) ** 2 + (G[..., 1] - uhy) ** 2).sum(-1)
    H1 = np.sqrt((wq * e2).sum())
    uh = np.einsum('qa,tac->tqc', R.PHI_Q, U)
    L2 = np.sqrt((wq * ((ex["u"](X[..., 0], X[..., 1]) - uh) ** 2).sum(-1)).sum())
    bl2 = 0.0; bdn = 0.0
    for (t, ref, Xe, we, n) in S.edge_data():
        ph = R.basis(ref[:, 0], ref[:, 1]); dX, dY = R.dbasis(ref[:, 0], ref[:, 1])
        Ji = S.Jinv[t]
        gxe = Ji[0, 0] * dX + Ji[1, 0] * dY; gye = Ji[0, 1] * dX + Ji[1, 1] * dY
        Ut = u[dofs(S.ids[t])]
        uhe = ph @ Ut; dnh = (n[0] * gxe + n[1] * gye) @ Ut
        ue = ex["u"](Xe[:, 0], Xe[:, 1]); Ge = ex["gu"](Xe[:, 0], Xe[:, 1])
        dne = Ge[:, :, 0] * n[0] + Ge[:, :, 1] * n[1]
        bl2 += (we * ((ue - uhe) ** 2).sum(-1)).sum()
        bdn += (we * ((dne - dnh) ** 2).sum(-1)).sum()
    en = np.sqrt(H1**2 + (mu / S.h) * bl2 + S.h * bdn)
    return dict(H1=H1, L2=L2, bL2=np.sqrt(bl2), bdn=np.sqrt(bdn), energy=en)


make_exact = svn.make_exact
