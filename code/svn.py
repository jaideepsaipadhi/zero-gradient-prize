"""
P4 Scott-Vogelius-Nitsche solver on (square) minus (inscribed polygon of the unit circle).

Forms (theta = 0: Gjerde-Scott (27)-(28) as printed; theta = 1: pressure-consistent):
  N_h(u,v) - (p, div v) + theta <p, v.n>_{Gamma_h} = (f, v)      v = 0 on the square
  (q, div u) = 0                                                   q in P3-disc (all 10 moments)
  N_h(u,v) = 1/2 (Du, Dv) - <d_n u, v> - <u, d_n v> + (mu/h) <u, v>,  h = h_Omega (max diameter)
Strong Dirichlet on the square: u = Lagrange P4 interpolant of g.
Divergence enforced exactly by the iterated penalty method (Scott's IPM).
"""
import os
os.environ.setdefault("PYPARDISO_MKL_RT", "/usr/local/lib/libmkl_rt.so.3")
import numpy as np
import scipy.sparse as sp
try:
    import pypardiso
except Exception:
    pypardiso = None
import sympy as sy

# ---------------------------------------------------------------- reference data
DEG = 4
MON4 = [(a, b) for a in range(DEG + 1) for b in range(DEG + 1 - a)]          # 15
MON3 = [(a, b) for a in range(4) for b in range(4 - a)]                        # 10
LBARY = [(i, j, DEG - i - j) for i in range(DEG + 1) for j in range(DEG + 1 - i)]  # (l1,l2,l0)


def _mono(M, x, y):
    return np.stack([x**a * y**b for (a, b) in M], axis=-1)


def _dmono(M, x, y):
    dx = np.stack([a * x**max(a - 1, 0) * y**b if a > 0 else 0 * x for (a, b) in M], -1)
    dy = np.stack([b * x**a * y**max(b - 1, 0) if b > 0 else 0 * x for (a, b) in M], -1)
    return dx, dy


# local P4 nodes in reference coords: x = l1/4, y = l2/4
REFNODES = np.array([[l1 / DEG, l2 / DEG] for (l1, l2, l0) in LBARY])
COEF = np.linalg.inv(_mono(MON4, REFNODES[:, 0], REFNODES[:, 1]))   # column k = basis k


def basis(x, y):
    return _mono(MON4, x, y) @ COEF


def dbasis(x, y):
    dx, dy = _dmono(MON4, x, y)
    return dx @ COEF, dy @ COEF


def tri_quad(n=8):
    g, w = np.polynomial.legendre.leggauss(n)
    g = (g + 1) / 2; w = w / 2
    U, V = np.meshgrid(g, g, indexing="ij"); WU, WV = np.meshgrid(w, w, indexing="ij")
    x = U.ravel(); y = (V * (1 - U)).ravel(); wt = (WU * WV * (1 - U)).ravel()
    return x, y, wt


QX, QY, QW = tri_quad(8)
PHI_Q = basis(QX, QY)                     # (nq,15)
DPX_Q, DPY_Q = dbasis(QX, QY)             # reference derivatives
Q3_Q = _mono(MON3, QX, QY)                # pressure basis (monomials in ref coords)
GL, GLW = np.polynomial.legendre.leggauss(10)
GL = (GL + 1) / 2; GLW = GLW / 2


# ---------------------------------------------------------------- mesh
def make_mesh(N, L=2.5, layers=None):
    """Cavalcante-type structured mesh: N outer vertices on the square (N/4 per side),
    projected radially to the unit circle, m layers, quads split on the forward diagonal."""
    assert N % 4 == 0
    m = layers or N // 4
    s = N // 4
    pts = []
    corners = [(L, -L), (L, L), (-L, L), (-L, -L)]
    for c in range(4):
        A = np.array(corners[c]); B = np.array(corners[(c + 1) % 4])
        for k in range(s):
            pts.append(A + (B - A) * k / s)
    b = np.array(pts)
    a = b / np.linalg.norm(b, axis=1)[:, None]
    X = np.zeros((m + 1, N, 2))
    for j in range(m + 1):
        X[j] = (1 - j / m) * a + (j / m) * b
    vid = np.arange((m + 1) * N).reshape(m + 1, N)
    verts = X.reshape(-1, 2)
    tris = []
    for j in range(m):
        for i in range(N):
            i1 = (i + 1) % N
            tris.append([vid[j, i], vid[j, i1], vid[j + 1, i1]])
            tris.append([vid[j, i], vid[j + 1, i1], vid[j + 1, i]])
    tris = np.array(tris)
    # orient counterclockwise
    P = verts[tris]
    ar = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    tris[ar < 0] = tris[ar < 0][:, [0, 2, 1]]
    circ = set(vid[0].tolist()); outer = set(vid[m].tolist())
    return verts, tris, circ, outer


class Space:
    def __init__(self, verts, tris, circ, outer):
        self.verts, self.tris = verts, tris
        nt = len(tris)
        edge_id = {}; ids = np.zeros((nt, 15), dtype=np.int64)
        nv = len(verts); nxt = nv
        coords = [None] * 0
        for t in range(nt):
            g = tris[t]
            for k, (l1, l2, l0) in enumerate(LBARY):
                lam = {0: l0, 1: l1, 2: l2}
                nz = [v for v in range(3) if lam[v] > 0]
                if len(nz) == 1:
                    ids[t, k] = g[nz[0]]
                elif len(nz) == 2:
                    va, vb = nz; ga, gb = g[va], g[vb]
                    if ga > gb:
                        ga, gb, va, vb = gb, ga, vb, va
                    key = (ga, gb, lam[va])
                    if key not in edge_id:
                        edge_id[key] = nxt; nxt += 1
                    ids[t, k] = edge_id[key]
                else:
                    ids[t, k] = nxt; nxt += 1
        self.nn = nxt; self.ids = ids
        P = verts[tris]
        self.P0 = P[:, 0]; self.J = np.stack([P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]], axis=2)  # (nt,2,2) cols
        self.detJ = np.linalg.det(self.J)
        self.Jinv = np.linalg.inv(self.J)
        # node coordinates
        xy = np.zeros((self.nn, 2))
        for t in range(nt):
            xy[ids[t]] = self.P0[t] + REFNODES @ self.J[t].T
        self.xy = xy
        self.hdiam = np.max(np.stack([np.linalg.norm(P[:, i] - P[:, (i + 1) % 3], axis=1)
                                      for i in range(3)], 1), 1)
        self.h = self.hdiam.max()
        # boundary edges on Gamma_h: both vertices on circle
        be = []
        for t in range(nt):
            for (u, w) in [(0, 1), (1, 2), (2, 0)]:
                if tris[t, u] in circ and tris[t, w] in circ:
                    be.append((t, u, w))
        self.bedges = be
        # Dirichlet nodes: on the square (|x| or |y| == L)
        L = np.abs(verts).max()
        on = (np.abs(np.abs(xy[:, 0]) - L) < 1e-10) | (np.abs(np.abs(xy[:, 1]) - L) < 1e-10)
        self.dnodes = np.where(on)[0]
        self.hGamma = max(np.linalg.norm(verts[tris[t, u]] - verts[tris[t, w]]) for (t, u, w) in be)

    # physical gradients of basis at quadrature points: (nt,nq,15)
    def grads(self, x, y):
        dX, dY = dbasis(x, y)
        Ji = self.Jinv  # (nt,2,2): ref = Ji (X - P0); grad_x = Ji^T grad_ref
        gx = np.einsum('tk,qb->tqb', Ji[:, 0, :1].reshape(-1, 1), dX) + np.einsum('tk,qb->tqb', Ji[:, 1, :1].reshape(-1, 1), dY)
        gy = np.einsum('tk,qb->tqb', Ji[:, 0, 1:].reshape(-1, 1), dX) + np.einsum('tk,qb->tqb', Ji[:, 1, 1:].reshape(-1, 1), dY)
        return gx, gy

    def phys(self, x, y):
        return self.P0[:, None, :] + np.einsum('tij,qj->tqi', self.J, np.stack([x, y], 1))

    def edge_data(self):
        """for each Gamma_h edge: ref coords of Gauss points, physical points, weights, outward normal"""
        out = []
        for (t, u, w) in self.bedges:
            R = np.array([[0, 0], [1, 0], [0, 1]], float)
            ref = R[u][None, :] + GL[:, None] * (R[w] - R[u])[None, :]
            A = self.verts[self.tris[t, u]]; B = self.verts[self.tris[t, w]]
            L = np.linalg.norm(B - A)
            tvec = (B - A) / L
            n = np.array([tvec[1], -tvec[0]])
            C3 = self.verts[self.tris[t]].mean(0)
            if np.dot(n, C3 - A) > 0:      # outward of Omega_h = away from the element interior
                n = -n
            X = A[None, :] + GL[:, None] * (B - A)[None, :]
            out.append((t, ref, X, GLW * L, n))
        return out


def dofs(ids):
    return np.stack([2 * ids, 2 * ids + 1], -1)   # (...,15,2)


def assemble(S, mu, f, theta, want=("A", "B", "M", "C", "F")):
    nt = len(S.tris)
    gx, gy = S.grads(QX, QY)
    wq = QW[None, :] * np.abs(S.detJ)[:, None]                  # (nt,nq)
    # A: 1/2 D:D
    Axx = np.einsum('tq,tqa,tqb->tab', wq, gx, gx); Ayy = np.einsum('tq,tqa,tqb->tab', wq, gy, gy)
    Axy = np.einsum('tq,tqa,tqb->tab', wq, gx, gy)     # int phi_a,x phi_b,y
    loc = np.zeros((nt, 15, 2, 15, 2))
    loc[:, :, 0, :, 0] = 2 * Axx + Ayy
    loc[:, :, 1, :, 1] = Axx + 2 * Ayy
    loc[:, :, 0, :, 1] = np.transpose(Axy, (0, 2, 1))   # (test comp0 a, trial comp1 b): phi_a,y * phi_b,x ... see below
    loc[:, :, 1, :, 0] = Axy
    # 1/2 D(v e0):D(u e1) = v_y u_x  -> test a comp0, trial b comp1: int phi_a,y phi_b,x = Axy[b,a]
    D = dofs(S.ids).reshape(nt, 30)
    rows = np.repeat(D, 30, axis=1).ravel(); cols = np.tile(D, (1, 30)).ravel()
    ndof = 2 * S.nn
    A = sp.coo_matrix((loc.reshape(nt, 30, 30).ravel(), (rows, cols)), shape=(ndof, ndof)).tocsr()
    # B: (q_j, div v)
    Bl = np.zeros((nt, 10, 15, 2))
    Bl[:, :, :, 0] = np.einsum('tq,qj,tqa->tja', wq, Q3_Q, gx)
    Bl[:, :, :, 1] = np.einsum('tq,qj,tqa->tja', wq, Q3_Q, gy)
    prow = (np.arange(nt)[:, None] * 10 + np.arange(10)[None, :])
    rB = np.repeat(prow, 30, axis=1).ravel(); cB = np.tile(D, (1, 10)).ravel()
    B = sp.coo_matrix((Bl.reshape(nt, 10, 30).ravel(), (rB, cB)), shape=(10 * nt, ndof)).tocsr()
    Ml = np.einsum('tq,qi,qj->tij', wq, Q3_Q, Q3_Q)
    Minv = sp.block_diag([np.linalg.inv(m) for m in Ml], format='csr')
    # F
    X = S.phys(QX, QY)
    fv = f(X[..., 0], X[..., 1])                        # (nt,nq,2)
    Fl = np.einsum('tq,qa,tqc->tac', wq, PHI_Q, fv)
    F = np.zeros(ndof); np.add.at(F, D.ravel(), Fl.reshape(nt, 30).ravel())
    # Nitsche + C
    rN, cN, vN = [], [], []
    rC, cC, vC = [], [], []
    gG = np.zeros(10 * nt); perim = 0.0
    for (t, ref, Xe, we, n) in S.edge_data():
        ph = basis(ref[:, 0], ref[:, 1])                 # (ng,15)
        dX, dY = dbasis(ref[:, 0], ref[:, 1])
        Ji = S.Jinv[t]
        gxe = Ji[0, 0] * dX + Ji[1, 0] * dY; gye = Ji[0, 1] * dX + Ji[1, 1] * dY
        dn = n[0] * gxe + n[1] * gye
        Mphi = np.einsum('g,ga,gb->ab', we, ph, ph)
        Dn = np.einsum('g,ga,gb->ab', we, ph, dn)       # int phi_a * dn phi_b  (test a, trial b)
        blk = -Dn - Dn.T + (mu / S.h) * Mphi          # scalar block, same for both comps
        d = dofs(S.ids[t])                              # (15,2)
        for c in range(2):
            dd = d[:, c]
            rN.append(np.repeat(dd, 15)); cN.append(np.tile(dd, 15)); vN.append(blk.ravel())
        q3 = _mono(MON3, ref[:, 0], ref[:, 1])
        gG[prow[t]] += we @ q3; perim += we.sum()
        Cl = np.einsum('g,gj,ga->ja', we, q3, ph)        # int q_j phi_a
        for c in range(2):
            rC.append(np.repeat(prow[t], 15)); cC.append(np.tile(d[:, c], 10)); vC.append((Cl * n[c]).ravel())
    Nb = sp.coo_matrix((np.concatenate(vN), (np.concatenate(rN), np.concatenate(cN))), shape=(ndof, ndof)).tocsr()
    C = sp.coo_matrix((np.concatenate(vC), (np.concatenate(rC), np.concatenate(cC))), shape=(10 * nt, ndof)).tocsr()
    return dict(A=A + Nb, Avol=A, B=B, C=C, Minv=Minv, F=F, Ml=Ml, gGamma=gG / perim)


def coupling(sysm, theta):
    """returns Bt with momentum pressure term  -(p, div v) + theta <p - pbar_Gamma(p), v.n>  =  -p^T Bt v"""
    B, C = sysm["B"], sysm["C"]
    if not theta:
        return B
    npr0 = B.shape[0]
    e = np.zeros(npr0); e[0::10] = 1.0
    phi = C.T @ e
    Cm = (C - sp.csr_matrix(sysm["gGamma"][:, None]) @ sp.csr_matrix(phi[None, :])).tocsr()
    return (B - Cm).tocsr()


def solve(S, sysm, g, theta, r=1e4, tol=1e-10, maxit=200, verbose=False):
    A, B, Minv, F = sysm["A"], sysm["B"], sysm["Minv"], sysm["F"]
    Bt = coupling(sysm, theta).T.tocsr()
    K = (A + r * (Bt @ Minv @ B)).tocsr()
    ndof = A.shape[0]
    dD = dofs(S.dnodes).ravel()
    ug = np.zeros(ndof)
    gv = g(S.xy[S.dnodes, 0], S.xy[S.dnodes, 1])
    ug[dD] = gv.ravel()
    free = np.setdiff1d(np.arange(ndof), dD)
    Kff = K[free][:, free].tocsr(); Kfd = K[free][:, dD]
    slv = pypardiso.PyPardisoSolver()
    slv.factorize(Kff)
    w = np.zeros(B.shape[0]); u = ug.copy()
    hist = []
    for it in range(maxit):
        rhs = (F + Bt @ w)[free] - Kfd @ ug[dD]
        u = ug.copy(); u[free] = slv.solve(Kff, rhs)
        Bu = B @ u
        dres = np.sqrt(abs(Bu @ (Minv @ Bu)))
        hist.append(dres)
        w = w - r * (Minv @ Bu)
        if dres < tol:
            break
    slv.free_memory(everything=True)
    return u, w, hist


def solve_lu(S, sysm, g, theta):
    """Robust direct solve of the saddle system with SuperLU partial pivoting."""
    import scipy.sparse.linalg as spl
    A, B, F = sysm["A"], sysm["B"], sysm["F"]
    Bt = coupling(sysm, theta)
    ndof = A.shape[0]
    dD = dofs(S.dnodes).ravel(); free = np.setdiff1d(np.arange(ndof), dD)
    ug = np.zeros(ndof); ug[dD] = g(S.xy[S.dnodes, 0], S.xy[S.dnodes, 1]).ravel()
    K = sp.bmat([[A[free][:, free], -Bt[:, free].T], [B[:, free], None]], format="csc")
    rhs = np.concatenate([F[free] - A[free][:, dD] @ ug[dD], -B[:, dD] @ ug[dD]])
    lu = spl.splu(K, permc_spec="COLAMD", diag_pivot_thresh=1.0)
    x = lu.solve(rhs)
    for _ in range(3):                      # iterative refinement
        x = x + lu.solve(rhs - K @ x)
    u = ug.copy(); u[free] = x[:len(free)]; p = x[len(free):]
    Bu = B @ u
    dres = np.sqrt(abs(Bu @ (sysm["Minv"] @ Bu)))
    relres = np.linalg.norm(K @ x - rhs) / np.linalg.norm(rhs)
    return u, p, [dres], relres


def solve_al(S, sysm, g, theta, r=1e3, tol=1e-13, maxit=300):
    """GMRES on the augmented saddle system with the block-triangular augmented-Lagrangian
    preconditioner  P = [[K_r, -Bt^T], [0, -M/r]],  K_r = A + r Bt^T M^{-1} B  (factorized once)."""
    import scipy.sparse.linalg as spl
    A, B, Minv, F, Ml = sysm["A"], sysm["B"], sysm["Minv"], sysm["F"], sysm["Ml"]
    Bt = coupling(sysm, theta)
    ndof = A.shape[0]; npr = B.shape[0]
    dD = dofs(S.dnodes).ravel(); free = np.setdiff1d(np.arange(ndof), dD)
    ug = np.zeros(ndof); ug[dD] = g(S.xy[S.dnodes, 0], S.xy[S.dnodes, 1]).ravel()
    Af = A[free][:, free].tocsr(); Bf = B[:, free].tocsr(); Btf = Bt[:, free].tocsr()
    BtfT = Btf.T.tocsr()
    c = -(B[:, dD] @ ug[dD])
    Kr = (Af + r * (BtfT @ Minv @ Bf)).tocsr()
    M = sp.block_diag(list(Ml), format="csr")
    slv = pypardiso.PyPardisoSolver()
    slv.factorize(Kr)
    nf = len(free)
    rhs = np.concatenate([F[free] - A[free][:, dD] @ ug[dD] + r * (BtfT @ (Minv @ c)), c])

    def matvec(x):
        u, p = x[:nf], x[nf:]
        return np.concatenate([Kr @ u - BtfT @ p, Bf @ u])

    def prec(y):
        a, b = y[:nf], y[nf:]
        p = -r * (Minv @ b)
        u = slv.solve(Kr, a + BtfT @ p)
        return np.concatenate([u, p])

    Kop = spl.LinearOperator((nf + npr, nf + npr), matvec=matvec)
    Pop = spl.LinearOperator((nf + npr, nf + npr), matvec=prec)
    nit = [0]
    x, info = spl.gmres(Kop, rhs, M=Pop, rtol=tol, atol=0.0, restart=200, maxiter=maxit,
                        callback=lambda rk: nit.__setitem__(0, nit[0] + 1), callback_type="pr_norm")
    slv.free_memory(everything=True)
    u = ug.copy(); u[free] = x[:nf]; p = x[nf:]
    Bu = B @ u
    dres = np.sqrt(abs(Bu @ (Minv @ Bu)))
    res = np.linalg.norm(matvec(x) - rhs) / np.linalg.norm(rhs)
    return u, p, [dres] * max(1, nit[0]), dict(info=info, relres=res, its=nit[0])


def solve_direct(S, sysm, g, theta):
    """Direct saddle-point solve:  [A  -(B - theta C)^T ; B  0] [u;p] = [F;0],
    Dirichlet rows eliminated. For theta=1 the constant pressure mode is in the kernel
    (it drops out of the momentum row), so pressure dof 0 (element 0, constant mode) is pinned."""
    A, B, C, F = sysm["A"], sysm["B"], sysm["C"], sysm["F"]
    if theta:
        # <p - pbar_Gamma(p), v.n>:  C_mod = C - g (e^T C),  g_j = (1/|Gamma_h|) int_Gamma q_j,
        # e = constant-one pressure (coefficient 1 on each element's constant monomial)
        npr0 = B.shape[0]
        e = np.zeros(npr0); e[0::10] = 1.0
        phi = C.T @ e                                   # v -> <1, v.n>
        gvec = sysm["gGamma"]
        C = (C - sp.csr_matrix(np.outer(gvec, np.ones(1))) @ sp.csr_matrix(phi[None, :])).tocsr()
    Bt = (B - theta * C)
    ndof = A.shape[0]; npr = B.shape[0]
    dD = dofs(S.dnodes).ravel()
    ug = np.zeros(ndof)
    ug[dD] = g(S.xy[S.dnodes, 0], S.xy[S.dnodes, 1]).ravel()
    free = np.setdiff1d(np.arange(ndof), dD)
    pk = np.arange(npr)
    K = sp.bmat([[A[free][:, free], -Bt[:, free].T],
                 [B[:, free], None]], format='csr')
    rhs = np.concatenate([F[free] - A[free][:, dD] @ ug[dD], -B[:, dD] @ ug[dD]])
    slv = pypardiso.PyPardisoSolver(mtype=11)
    slv.set_iparm(11, 1)   # scaling (1-based iparm 11)
    slv.set_iparm(13, 1)   # weighted matching
    slv.factorize(K)
    x = slv.solve(K, rhs)
    nb = np.linalg.norm(rhs)
    for _ in range(30):
        res = rhs - K @ x
        if np.linalg.norm(res) < 1e-14 * nb:
            break
        x = x + slv.solve(K, res)
    slv.free_memory(everything=True)
    u = ug.copy(); u[free] = x[:len(free)]
    p = x[len(free):len(free) + npr].copy()
    Bu = B @ u
    dres = np.sqrt(abs(Bu @ (sysm["Minv"] @ Bu)))
    return u, p, [dres]


def errors(S, u, ex, mu):
    """H1-seminorm, boundary L2, boundary d_n, energy; ex = dict(u, gu) exact functions"""
    nt = len(S.tris)
    gx, gy = S.grads(QX, QY)
    wq = QW[None, :] * np.abs(S.detJ)[:, None]
    X = S.phys(QX, QY)
    D = dofs(S.ids)                                      # (nt,15,2)
    U = u[D]                                             # (nt,15,2)
    uhx = np.einsum('tqa,tac->tqc', gx, U); uhy = np.einsum('tqa,tac->tqc', gy, U)
    G = ex["gu"](X[..., 0], X[..., 1])                   # (nt,nq,2,2)  G[..., c, d] = d_d u_c
    e2 = ((G[..., 0] - uhx) ** 2 + (G[..., 1] - uhy) ** 2).sum(-1)
    H1 = np.sqrt((wq * e2).sum())
    uh = np.einsum('qa,tac->tqc', PHI_Q, U)
    L2 = np.sqrt((wq * ((ex["u"](X[..., 0], X[..., 1]) - uh) ** 2).sum(-1)).sum())
    bl2 = 0.0; bdn = 0.0; slipn = 0.0
    for (t, ref, Xe, we, n) in S.edge_data():
        ph = basis(ref[:, 0], ref[:, 1]); dX, dY = dbasis(ref[:, 0], ref[:, 1])
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


# ---------------------------------------------------------------- exact solutions
def make_exact(kind):
    x, y = sy.symbols('x y', real=True)
    r2 = x**2 + y**2
    if kind == "A":                       # Scott's shear flow, p = 0
        u = [(1 - 1 / r2) * (-y), (1 - 1 / r2) * x]; p = sy.Integer(0)
    elif kind.startswith("B"):            # non-constant pressure, u = curl psi, psi=(r^2-1)^2 phi
        amp = float(kind[1:]) if len(kind) > 1 else 1.0
        psi = (r2 - 1)**2 * (x + y**2 / 2) / 10
        u = [sy.diff(psi, y), -sy.diff(psi, x)]
        p = amp * (x**2 * y - y + x / 2)
    elif kind.startswith("P"):            # pure pressure: u = 0, f = grad p  (GS output = traction response exactly)
        amp = float(kind[1:]) if len(kind) > 1 else 1.0
        u = [sy.Integer(0), sy.Integer(0)]
        p = amp * (x**2 * y - y + x / 2)
    elif kind == "C":                     # same psi as B, NON-polynomial pressure (exercises p* - pi_h p*)
        psi = (r2 - 1)**2 * (x + y**2 / 2) / 10
        u = [sy.diff(psi, y), -sy.diff(psi, x)]
        p = sy.exp(x / 2) * sy.cos(y)
    else:
        raise ValueError
    f = [-(sy.diff(u[i], x, 2) + sy.diff(u[i], y, 2)) + sy.diff(p, [x, y][i]) for i in range(2)]
    f = [sy.simplify(fi) for fi in f]
    gu = [[sy.diff(u[i], v) for v in (x, y)] for i in range(2)]
    divu = sy.simplify(sy.diff(u[0], x) + sy.diff(u[1], y))
    assert divu == 0
    U = sy.lambdify((x, y), u, 'numpy'); Fl = sy.lambdify((x, y), f, 'numpy')
    GU = sy.lambdify((x, y), gu, 'numpy'); P = sy.lambdify((x, y), p, 'numpy')

    def vec(fn):
        def g(X, Y):
            vals = fn(X, Y)
            return np.stack([np.broadcast_to(np.asarray(v, float), np.shape(X)) for v in vals], -1)
        return g

    def gmat(X, Y):
        vals = GU(X, Y)
        return np.stack([np.stack([np.broadcast_to(np.asarray(vals[i][j], float), np.shape(X))
                                   for j in range(2)], -1) for i in range(2)], -2)
    return dict(u=vec(U), f=vec(Fl), gu=gmat, p=lambda X, Y: np.broadcast_to(np.asarray(P(X, Y), float), np.shape(X)))


def run(N, kind, mu, theta, r=1e4, verbose=False, method="ipm"):
    ex = make_exact(kind)
    verts, tris, circ, outer = make_mesh(N)
    S = Space(verts, tris, circ, outer)
    sysm = assemble(S, mu, ex["f"], theta)
    if method == "direct":
        u, w, hist = solve_direct(S, sysm, ex["u"], theta)
    elif method == "al":
        u, w, hist, info = solve_al(S, sysm, ex["u"], theta)
    else:
        u, w, hist = solve(S, sysm, ex["u"], theta, r=r)
    E = errors(S, u, ex, mu)
    E.update(N=N, kind=kind, mu=mu, theta=theta, h=S.h, hG=S.hGamma, ntri=len(tris),
             ndof=2 * S.nn, divres=hist[-1], its=len(hist))
    return E, S, u, w


if __name__ == "__main__":
    import sys, time
    t0 = time.time()
    for kind in ("A", "B"):
        for th in (0, 1):
            E, *_ = run(16, kind, 100.0, th)
            print({k: (float(f"{v:.4g}") if isinstance(v, float) else v) for k, v in E.items()}, time.time() - t0)
