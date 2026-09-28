"""
notes/v6/unifmu2: the periodic cell problem for the capture fraction phi(gamma)  (REPORT.tex, Section 4).

Rotationally symmetric version of the leak problem: Omega_h = (polygon of R_out) minus (regular N-gon inscribed in
the unit circle), N identical angular sectors, quads between polygonal layers split on the forward diagonal
(the same pattern as code/svn.py). Data p - pbar = e^{i m theta}. The discrete problem is invariant under
the rotation by alpha = 2 pi/N, so it decouples into Floquet modes; the data is pure mode m, and the
solution satisfies  u(R x) = e^{i m alpha} R u(x). We solve it on ONE sector (one wall edge) with this
quasi-periodic condition: a 1D-periodic strip cell problem with the exact polygonal curvature.

  X_gamma : N_h(X, v) = lam <g, v> on Z_h, lam = gamma/l (penalty per edge length l), g = -(p~ - pbar) n_e
  U       : exact annulus leak flow (Definition hc:def with Omega = {1 < r < R_out}), mode m (closed form)
  W       = X - U~,  A(gamma) = ||grad W||_{Omega_h} / (G l^{1/2})      (G = ||p - pbar||_{L2(Gamma)})
  phi_E(gamma) = ||grad W(gamma)|| / ||grad W(infinity)||   (energy capture fraction)
  phi_P(gamma) = <W, q>/<g - U~, q>, q = beta(s) tau_e       (pairing capture fraction of the unifmu REPORT)
  gamma_c : penalty below which N_h is not positive definite on Z_h (all Floquet phases), and on V_h.
Usage: python floquet_cell.py coerc N a            (threshold scan)
       python floquet_cell.py phi N a [m]          (capture fraction scan)
"""
import os, sys, json, time
for v_ in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v_, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, scipy.linalg as sla
import sympy as sy
import svn

ROUT = 2.5
NQ = 16
GQ, GQW = np.polynomial.legendre.leggauss(NQ); GQ = (GQ + 1) / 2; GQW = GQW / 2
R3 = np.array([[0, 0], [1, 0], [0, 1]], float)


def sector_mesh(N, a, rout=ROUT, J=None):
    """vertices v(j,i), i in {0,1} (angles 0, 2pi/N), radii r_0 = 1 < r_1 < ... < r_J = rout,
    r_{j+1} = r_j (1 + a*2*sin(pi/N)) (cells of aspect ~a), last layer adjusted. J given: truncated strip."""
    al = 2 * np.pi / N; ell = 2 * np.sin(np.pi / N)
    r = [1.0]
    while r[-1] < rout:
        r.append(r[-1] * (1 + a * ell))
        if J is not None and len(r) > J:
            break
    if J is None:
        r[-1] = rout
        if r[-1] - r[-2] < 0.5 * a * ell * r[-2]:
            r.pop(-2)
    r = np.array(r)
    verts = []
    for j, rj in enumerate(r):
        for i in range(2):
            verts.append([rj * np.cos(i * al), rj * np.sin(i * al)])
    verts = np.array(verts)
    vid = lambda j, i: 2 * j + i
    tris = []
    for j in range(len(r) - 1):
        tris.append([vid(j, 0), vid(j, 1), vid(j + 1, 1)])
        tris.append([vid(j, 0), vid(j + 1, 1), vid(j + 1, 0)])
    tris = np.array(tris)
    circ = {vid(0, 0), vid(0, 1)}; outer = {vid(len(r) - 1, 0), vid(len(r) - 1, 1)}
    return verts, tris, circ, outer, r, ell, al


class Cell:
    def __init__(self, N, a, J=None):
        verts, tris, circ, outer, r, ell, al = sector_mesh(N, a, J=J)
        S = svn.Space(verts, tris, circ, outer)
        S.h = 1.0                                   # so that assemble's penalty mu/h is lam itself
        self.S, self.N, self.a, self.ell, self.al, self.r = S, N, a, ell, al, r
        xy = S.xy
        ang = np.arctan2(xy[:, 1], xy[:, 0]); rad = np.hypot(xy[:, 0], xy[:, 1])
        tol = 1e-9
        left = np.where(np.abs(xy[:, 1]) < tol * np.maximum(1, rad))[0]
        right = np.where(np.abs(ang - al) < 1e-9)[0]
        Rm = np.array([[np.cos(al), -np.sin(al)], [np.sin(al), np.cos(al)]])
        # match right node x_R = R x_L
        pre = xy[right] @ Rm                          # R^{-1} x_R  (row vectors: x R^{-T} = x R)
        match = {}
        for k, node in enumerate(right):
            d = np.linalg.norm(xy[left] - pre[k], axis=1); j = np.argmin(d)
            assert d[j] < 1e-9
            match[node] = left[j]
        self.Rm, self.match = Rm, match
        # outer boundary nodes (Dirichlet 0)
        rmax = r[-1]
        outer_nodes = set()
        for t, g in enumerate(tris):
            for (u, w) in [(0, 1), (1, 2), (2, 0)]:
                if g[u] in outer and g[w] in outer:
                    for k, (l1, l2, l0) in enumerate(svn.LBARY):
                        if {0: l0, 1: l1, 2: l2}[3 - u - w] == 0:
                            outer_nodes.add(S.ids[t, k])
        self.outer_nodes = outer_nodes
        self.red_nodes = [n for n in range(S.nn) if n not in match and n not in outer_nodes]
        self.pos = {n: k for k, n in enumerate(self.red_nodes)}
        # wall edge data
        (t, ref, Xe, we, n) = S.edge_data()[0]
        u_, w_ = [(u, w) for (tt, u, w) in S.bedges][0]
        A = S.verts[S.tris[t, u_]]; B = S.verts[S.tris[t, w_]]
        if np.cross(np.append(A, 0), np.append(B, 0))[2] < 0:
            u_, w_ = w_, u_; A, B = B, A
        L = np.linalg.norm(B - A); tau = (B - A) / L
        ref = R3[u_][None, :] + GQ[:, None] * (R3[w_] - R3[u_])[None, :]
        ph = svn.basis(ref[:, 0], ref[:, 1]); dX, dY = svn.dbasis(ref[:, 0], ref[:, 1]); Ji = S.Jinv[t]
        gxe = Ji[0, 0] * dX + Ji[1, 0] * dY; gye = Ji[0, 1] * dX + Ji[1, 1] * dY
        self.wall = dict(t=t, ph=ph, dn=n[0] * gxe + n[1] * gye, n=n, tau=tau, L=L, w=GQW * L,
                         s=(GQ - 0.5) * L, X=A[None, :] + GQ[:, None] * (B - A)[None, :])
        sysm = svn.assemble(S, 1.0, lambda X, Y: np.zeros(np.shape(X) + (2,)), 0, want=("A", "B"))
        self.Avol = sysm["Avol"].tocsr(); self.Nb1 = (sysm["A"] - sysm["Avol"]).tocsr()   # penalty mu=1 included
        self.B = sysm["B"].tocsr(); self.Minv = sysm["Minv"]
        # split Nitsche part into consistency terms and penalty mass (penalty is linear in mu)
        sys2 = svn.assemble(S, 2.0, lambda X, Y: np.zeros(np.shape(X) + (2,)), 0, want=("A",))
        self.Mg = (sys2["A"] - sysm["A"]).tocsr()              # penalty mass  <u, v>_Gamma
        self.Cn = (self.Nb1 - self.Mg).tocsr()                 # -<dn u, v> - <u, dn v>

    def T(self, phase):
        """complex prolongation from reduced dofs to sector dofs for Floquet phase e^{i phase}"""
        S = self.S; rows, cols, vals = [], [], []
        om = np.exp(1j * phase)
        for n in range(S.nn):
            if n in self.outer_nodes:
                continue
            if n in self.match:
                k = self.pos[self.match[n]]
                for c in range(2):
                    for d in range(2):
                        rows.append(2 * n + c); cols.append(2 * k + d); vals.append(om * self.Rm[c, d])
            else:
                k = self.pos[n]
                for c in range(2):
                    rows.append(2 * n + c); cols.append(2 * k + c); vals.append(1.0 + 0j)
        return sp.coo_matrix((vals, (rows, cols)), shape=(2 * S.nn, 2 * len(self.red_nodes))).tocsr()

    def K(self, lam):
        return (self.Avol + self.Cn + lam * self.Mg).tocsr()


# ------------------------------------------------------------------ exact annulus leak flow, mode m
def leak_flow(m, rout=ROUT):
    x, y = sy.symbols('x y', real=True)
    r = sy.sqrt(x ** 2 + y ** 2)
    rs = sy.symbols('r', positive=True)
    if m == 1:
        basis = [rs, 1 / rs, rs ** 3, rs * sy.log(rs)]
    else:
        basis = [rs ** m, rs ** (-m), rs ** (m + 2), rs ** (2 - m)]
    c = sy.symbols('c0:4')
    F = sum(ci * b for ci, b in zip(c, basis))
    eqs = [F.subs(rs, 1) - 1, sy.diff(F, rs).subs(rs, 1), F.subs(rs, rout), sy.diff(F, rs).subs(rs, rout)]
    sol = sy.solve(eqs, c)
    F = F.subs(sol)
    Fr = F.subs(rs, r)
    eim = ((x + sy.I * y) / r) ** m
    Psi = Fr / (sy.I * m) * eim
    U = [sy.diff(Psi, y), -sy.diff(Psi, x)]
    G = [[sy.diff(U[cc], v) for v in (x, y)] for cc in range(2)]
    fu = sy.lambdify((x, y), U, "numpy"); fg = sy.lambdify((x, y), G, "numpy")

    def u(X, Y):
        a = fu(X, Y); return np.stack([np.broadcast_to(np.asarray(a[0], complex), np.shape(X)),
                                       np.broadcast_to(np.asarray(a[1], complex), np.shape(X))], -1)

    def gu(X, Y):
        a = fg(X, Y)
        return np.stack([np.stack([np.broadcast_to(np.asarray(a[cc][d], complex), np.shape(X)) for d in range(2)], -1)
                         for cc in range(2)], -2)
    return u, gu


def grad_err(cell, vec, gu):
    S = cell.S
    gx, gy = S.grads(svn.QX, svn.QY)
    wq = svn.QW[None, :] * np.abs(S.detJ)[:, None]
    X = S.phys(svn.QX, svn.QY)
    Vl = vec[svn.dofs(S.ids)]
    vx = np.einsum('tqa,tac->tqc', gx, Vl); vy = np.einsum('tqa,tac->tqc', gy, Vl)
    G = gu(X[..., 0], X[..., 1]) if gu is not None else 0
    e = np.stack([vx, vy], -1) - G
    return np.sqrt((wq * (np.abs(e) ** 2).sum((-1, -2))).sum())


# ------------------------------------------------------------------ coercivity threshold
def cmd_coerc(N, a, J=6):
    cell = Cell(N, a, J=J)
    out = []
    for phase in np.linspace(0, np.pi, 9):
        T = cell.T(phase).toarray()
        BT = cell.B @ T
        Z = sla.null_space(BT, rcond=1e-11)
        A0 = T.conj().T @ (cell.Avol + cell.Cn) @ T; Mg = T.conj().T @ cell.Mg @ T
        res = {}
        for name, Q in (("Z", Z), ("V", np.eye(T.shape[1]))):
            A0q = Q.conj().T @ A0 @ Q; Mq = Q.conj().T @ Mg @ Q
            A0q = (A0q + A0q.conj().T) / 2; Mq = (Mq + Mq.conj().T) / 2
            # threshold: smallest lam with A0 + lam M >= 0  <=>  lam >= max generalized eig of (-A0, M) on range(M)
            w, U = np.linalg.eigh(Mq); keep = w > 1e-12 * w.max()
            # restrict to complement of ker M: A0 must be PD on ker M (checked), then Schur
            Uk = U[:, keep]; U0 = U[:, ~keep]
            A00 = U0.conj().T @ A0q @ U0
            if A00.shape[0]:
                assert np.linalg.eigvalsh((A00 + A00.conj().T) / 2).min() > 0
                Ak0 = Uk.conj().T @ A0q @ U0; Akk = Uk.conj().T @ A0q @ Uk
                Sch = Akk - Ak0 @ np.linalg.solve(A00, Ak0.conj().T)
            else:
                Sch = Uk.conj().T @ A0q @ Uk
            Mk = np.diag(w[keep])
            ev = sla.eigh((Sch + Sch.conj().T) / 2, Mk, eigvals_only=True)
            res[name] = float(-ev.min() * cell.ell)                # gamma_c = lam_c * l
        out.append(dict(N=N, a=a, phase=float(phase), gamma_c_Z=res["Z"], gamma_c_V=res["V"]))
        print(json.dumps(out[-1]), flush=True)
    print(json.dumps(dict(N=N, a=a, gamma_c_Z_max=max(o["gamma_c_Z"] for o in out),
                          gamma_c_V_max=max(o["gamma_c_V"] for o in out))), flush=True)


# ------------------------------------------------------------------ capture fraction
def cmd_phi(N, a, m=2, gammas=None):
    t0 = time.time()
    cell = Cell(N, a); S = cell.S
    phase = m * cell.al
    T = cell.T(phase); TH = T.conj().T.tocsr()
    BT = (cell.B @ T).tocsr()
    u, gu = leak_flow(m)
    Wd = cell.wall
    # data g = -(p~ - pbar) n_e, p~ = e^{i m theta}
    th = np.arctan2(Wd["X"][:, 1], Wd["X"][:, 0])
    g = -np.exp(1j * m * th)[:, None] * Wd["n"][None, :]
    Ut = u(Wd["X"][:, 0], Wd["X"][:, 1])
    F = np.zeros(2 * S.nn, complex); d = svn.dofs(S.ids[Wd["t"]])
    for c in range(2):
        np.add.at(F, d[:, c], (Wd["w"][:, None] * Wd["ph"] * g[:, c:c + 1]).sum(0))
    beta = Wd["s"] * (Wd["L"] ** 2 / 4 - Wd["s"] ** 2) / Wd["L"] ** 2
    cq = np.exp(1j * m * np.arctan2(Wd["X"][NQ // 2, 1], Wd["X"][NQ // 2, 0]))
    q = cq * beta[:, None] * Wd["tau"][None, :]
    zq = (Wd["w"][:, None] * (g - Ut) * np.conj(q)).sum()
    Fr = TH @ F
    nB = BT.shape[0]
    gam_list = gammas or [0.5, 0.75, 1, 1.5, 2, 3, 5, 10, 20, 30, 50, 100, 300, 1e3, 3e3, 1e4, 1e10]
    Ntot = cell.N
    G = np.sqrt(2 * np.pi)                                   # ||e^{i m theta}||_{L2(Gamma)}
    res = []
    for gam in gam_list:
        lam = gam / cell.ell
        Kr = (TH @ cell.K(lam) @ T).tocsc()
        M = sp.bmat([[Kr, -BT.conj().T], [BT, None]], format="csc")
        rhs = np.concatenate([lam * Fr, np.zeros(nB)])
        lu = spl.splu(M, permc_spec="COLAMD", diag_pivot_thresh=1.0)
        x = lu.solve(rhs)
        for _ in range(2):
            x = x + lu.solve(rhs - M @ x)
        X = T @ x[:Kr.shape[0]]
        divres = np.abs(BT @ x[:Kr.shape[0]]).max()
        eW = grad_err(cell, X, gu)                            # sector
        gradW = np.sqrt(Ntot) * eW
        XW = Wd["ph"] @ X[d]
        W = XW - Ut
        phiP = ((Wd["w"][:, None] * W * np.conj(q)).sum() / zq)
        res.append(dict(N=N, a=a, m=m, gamma=gam, A=gradW / (G * np.sqrt(cell.ell)), phiP_re=phiP.real,
                        phiP_im=phiP.imag, divres=float(divres),
                        relres=float(np.linalg.norm(M @ x - rhs) / np.linalg.norm(rhs))))
    Ainf = res[-1]["A"]
    for rr in res:
        rr["phiE"] = rr["A"] / Ainf
        print(json.dumps({k: (float(f"{v:.6g}") if isinstance(v, (float, np.floating)) else v) for k, v in rr.items()}),
              flush=True)
    gradU = np.sqrt(Ntot) * grad_err(cell, np.zeros(2 * S.nn, complex), gu)
    print(json.dumps(dict(N=N, a=a, m=m, ell=cell.ell, layers=len(cell.r) - 1, ndof=int(T.shape[1]),
                          K_annulus=float(gradU / G), secs=round(time.time() - t0, 1))), flush=True)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "coerc":
        cmd_coerc(int(sys.argv[2]), float(sys.argv[3]))
    elif cmd == "phi":
        cmd_phi(int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else 2)
