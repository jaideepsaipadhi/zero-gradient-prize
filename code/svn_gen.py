"""
Scott-Vogelius + Nitsche on GENERAL meshes (mesh_general.py) with several curved components, and
steady Navier-Stokes.  Element code, quadrature, the volume/Nitsche assembly and the error norms
are reused from svn_k.py unchanged (svn_k.assemble, svn_k.errors); this file adds

  * GSpace: svn_k.Space for an arbitrary triangulation: Gamma_h edges with a component label,
    Dirichlet part of Sigma = listed box sides, optional natural (outflow) sides;
  * pressure couplings  -p^T Bt v  with Bt =
        "GS"   B                                  (printed method, theta = 0)
        "CNS"  B - (C - gbar (x) phi)             (theta = 1, GLOBAL boundary mean over all Gamma_h,i;
                                                   = svn_k.coupling(.., 1))
        "PER"  B - sum_i (C_i - g_i (x) phi_i)    (per-component means, Remark gd:rem:percomp)
        "NAIVE" B - C                             (theta = 1, no mean; singular when all of Sigma is
                                                   Dirichlet, the consistent method with an outflow)
    where C_i v.q = <q, v.n>_{Gamma_h,i}, phi_i(v) = int_{Gamma_h,i} v.n, g_i(q) = mean of q on Gamma_h,i;
  * steady Navier-Stokes  nu N_h(u,v) + c(u;u,v) - p^T Bt v = (f,v),  c = c1 (standard) or cs (skew),
    Newton's method (exact Jacobian), each linear step a regularised PARDISO saddle solve exactly
    as svn_k.solve_pardiso (pressure block +1e-12 M);
  * optional do-nothing outflow: the term -nu <(grad u)^T n, v>_out is added so that the natural
    boundary condition of nu*a_h is nu d_n u - p n = 0 (the Schafer-Turek condition), since
    a_h = 1/2 (Du, Dv) with D = grad + grad^T;
  * diagnostics: per-component fluxes, leak-profile misfit, pressure errors, drag/lift by the
    boundary traction on Gamma_h and by the volume (variational) formula.

Sign conventions: n on Gamma_h is the outward normal of Omega_h (points INTO the obstacle), so
Phi_i(u_h) = int u_h.n > 0 means fluid leaves Omega_h through Gamma_h,i.
"""
import os
os.environ.setdefault("PYPARDISO_MKL_RT", "/usr/local/lib/libmkl_rt.so.3")
import numpy as np
import scipy.sparse as sp
import svn_k
from svn_k import dofs


# ---------------------------------------------------------------- space on a general mesh
class GSpace(svn_k.Space):
    def __init__(self, M, R, dirichlet_sides=(1, 2, 3, 4), outflow_sides=()):
        self.R = R; self.M = M
        verts, tris = M.verts, M.tris
        self.verts, self.tris = verts, tris
        nt = len(tris); nl = R.nl
        edge_id = {}; ids = np.zeros((nt, nl), dtype=np.int64)
        nxt = len(verts)
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
        self.nn = nxt; self.ids = ids
        P = verts[tris]
        self.P0 = P[:, 0]; self.J = np.stack([P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]], axis=2)
        self.detJ = np.linalg.det(self.J); self.Jinv = np.linalg.inv(self.J)
        xy = np.zeros((self.nn, 2))
        xy[ids.ravel()] = (self.P0[:, None, :] + np.einsum('tij,qj->tqi', self.J, R.REFNODES)).reshape(-1, 2)
        self.xy = xy
        self.hdiam = np.max(np.stack([np.linalg.norm(P[:, i] - P[:, (i + 1) % 3], axis=1)
                                      for i in range(3)], 1), 1)
        self.h = self.hdiam.max()
        # edge -> (t, u, w)
        loc = {}
        for t in range(nt):
            for (u, w) in [(0, 1), (1, 2), (2, 0)]:
                loc[tuple(sorted((tris[t, u], tris[t, w])))] = (t, u, w)
        self.bedges, self.bcomp = [], []
        for (a, b, i) in M.gamma_edges:
            self.bedges.append(loc[tuple(sorted((a, b)))]); self.bcomp.append(i)
        self.bcomp = np.array(self.bcomp); self.ncomp = int(self.bcomp.max()) + 1
        self.hGamma = max(np.linalg.norm(verts[tris[t, u]] - verts[tris[t, w]]) for (t, u, w) in self.bedges)
        dn, self.out_edges = set(), []
        for (a, b, side) in M.sigma_edges:
            t, u, w = loc[tuple(sorted((a, b)))]
            if side in dirichlet_sides:
                opp = 3 - u - w                           # local vertex not on the edge
                for kk, (l1, l2, l0) in enumerate(R.LBARY):
                    lam = (l0, l1, l2)
                    if lam[opp] == 0:
                        dn.add(int(ids[t, kk]))
            elif side in outflow_sides:
                self.out_edges.append((t, u, w))
            else:
                raise ValueError(f"side {side} neither Dirichlet nor outflow")
        self.dnodes = np.array(sorted(dn))


def circle_space(N, k, L=2.5):
    """the paper's structured circle mesh (svn_k.build) with a component label."""
    S = svn_k.build(N, k, "std", L=L)
    S.bcomp = np.zeros(len(S.bedges), int); S.ncomp = 1; S.out_edges = []
    return S


# ---------------------------------------------------------------- per-component boundary data
def gamma_extras(S):
    """C_i (npres x ndof), phi_i (ndof), g_i (npres, mean functional on Gamma_h,i), |Gamma_h,i|."""
    R = S.R; nl, npl = R.nl, R.npl; nt = len(S.tris); ndof = 2 * S.nn
    prow = np.arange(nt)[:, None] * npl + np.arange(npl)[None, :]
    rows = [[] for _ in range(S.ncomp)]; cols = [[] for _ in range(S.ncomp)]; vals = [[] for _ in range(S.ncomp)]
    phi = np.zeros((S.ncomp, ndof)); g = np.zeros((S.ncomp, npl * nt)); lens = np.zeros(S.ncomp)
    for (t, ref, Xe, we, n), i in zip(S.edge_data(), S.bcomp):
        ph = R.basis(ref[:, 0], ref[:, 1]); q = R.pbasis(ref[:, 0], ref[:, 1])
        d = dofs(S.ids[t])
        Cl = np.einsum('g,gj,ga->ja', we, q, ph)
        for c in range(2):
            rows[i].append(np.repeat(prow[t], nl)); cols[i].append(np.tile(d[:, c], npl))
            vals[i].append((Cl * n[c]).ravel())
            np.add.at(phi[i], d[:, c], (we @ ph) * n[c])
        g[i, prow[t]] += we @ q; lens[i] += we.sum()
    C = [sp.coo_matrix((np.concatenate(vals[i]), (np.concatenate(rows[i]), np.concatenate(cols[i]))),
                       shape=(npl * nt, ndof)).tocsr() for i in range(S.ncomp)]
    return dict(C=C, phi=phi, g=g / lens[:, None], len=lens)


def coupling(sysm, ext, method):
    B = sysm["B"]
    if method == "GS":
        return B
    Ctot = sum(ext["C"])
    if method == "NAIVE":
        return (B - Ctot).tocsr()
    if method == "CNS":                                   # global mean (identical to svn_k.coupling)
        gbar = (ext["g"] * ext["len"][:, None]).sum(0) / ext["len"].sum()
        phit = ext["phi"].sum(0)
        return (B - (Ctot - sp.csr_matrix(gbar[:, None]) @ sp.csr_matrix(phit[None, :]))).tocsr()
    if method == "PER":
        Cm = Ctot
        for i in range(len(ext["C"])):
            Cm = Cm - sp.csr_matrix(ext["g"][i][:, None]) @ sp.csr_matrix(ext["phi"][i][None, :])
        return (B - Cm).tocsr()
    raise ValueError(method)


# ---------------------------------------------------------------- convection
def _elem(S):
    R = S.R
    gx, gy = S.grads(R.QX, R.QY)
    wq = R.QW[None, :] * np.abs(S.detJ)[:, None]
    return gx, gy, wq


def _uq(S, u, E):
    R = S.R
    gx, gy, wq = E
    U = u[dofs(S.ids)]                                    # (nt,nl,2)
    Uq = np.einsum('qa,tac->tqc', R.PHI_Q, U)
    Gq = np.stack([np.einsum('tqa,tac->tqc', gx, U), np.einsum('tqa,tac->tqc', gy, U)], -1)  # [c,d] = d_d u_c
    return Uq, Gq


def conv_parts(S, u, E=None):
    """element matrices of c1(w;.,.) linearisations at w = u:
      K1 (nt,nl,nl):      c1(u; d, v)   (same block for both components)
      K2 (nt,nl,2,nl,2):  c1(d; u, v)
      K3 (nt,nl,2,nl,2):  c1(d; v, u)   (test v, trial d)"""
    R = S.R
    E = E if E is not None else _elem(S)
    gx, gy, wq = E
    Uq, Gq = _uq(S, u, E)
    adv = Uq[..., 0:1] * gx + Uq[..., 1:2] * gy          # (nt,nq,nl): u.grad phi_b
    WP = wq[:, :, None] * R.PHI_Q[None]                   # (nt,nq,nl)
    K1 = np.einsum('tqa,tqb->tab', WP, adv, optimize=True)
    K2 = np.einsum('tqa,qb,tqcd->tacbd', WP, R.PHI_Q, Gq, optimize=True)
    gd = np.stack([gx, gy], -1)                           # (nt,nq,nl,d)
    WU = WP[:, :, :, None] * Uq[:, :, None, :]            # (nt,nq,nl_b,c)
    K3 = np.einsum('tqbc,tqad->tacbd', WU, gd, optimize=True)
    return K1, K2, K3


def conv_residual(S, u, form, E=None):
    """vector v -> c(u;u,v), without forming matrices."""
    R = S.R
    E = E if E is not None else _elem(S)
    gx, gy, wq = E
    Uq, Gq = _uq(S, u, E)
    ugu = np.einsum('tqd,tqcd->tqc', Uq, Gq)              # (u.grad)u
    loc = np.einsum('tq,qa,tqc->tac', wq, R.PHI_Q, ugu, optimize=True)
    if form == "cs":
        adv = Uq[..., 0:1] * gx + Uq[..., 1:2] * gy       # u.grad phi_a
        loc = 0.5 * loc - 0.5 * np.einsum('tq,tqa,tqc->tac', wq, adv, Uq, optimize=True)
    out = np.zeros(2 * S.nn)
    np.add.at(out, dofs(S.ids).ravel(), loc.ravel())
    return out


def _scatter(S, loc):
    nt = len(S.tris); n2 = 2 * S.R.nl; ndof = 2 * S.nn
    D = dofs(S.ids).reshape(nt, n2)
    rows = np.repeat(D, n2, axis=1).ravel(); cols = np.tile(D, (1, n2)).ravel()
    return sp.coo_matrix((loc.reshape(nt, n2, n2).ravel(), (rows, cols)), shape=(ndof, ndof)).tocsr()


def conv_system(S, u, form, E=None):
    """returns (N(u) residual vector, Jacobian matrix) for c = c1 or cs."""
    nl = S.R.nl
    K1, K2, K3 = conv_parts(S, u, E)
    K1f = np.zeros((len(S.tris), nl, 2, nl, 2))
    K1f[:, :, 0, :, 0] = K1; K1f[:, :, 1, :, 1] = K1
    M1 = _scatter(S, K1f); M2 = _scatter(S, K2)
    if form == "c1":
        return M1 @ u, (M1 + M2).tocsr()
    M3 = _scatter(S, K3)
    Nu = 0.5 * (M1 @ u - M1.T @ u)
    Jc = 0.5 * (M1 + M2) - 0.5 * (M1.T + M3)
    return Nu, Jc.tocsr()


def outflow_matrix(S):
    """-<(grad u)^T n, v>_out  (test comp c, trial comp d:  -int phi_a d_c phi_b n_d)."""
    R = S.R; nl = R.nl; ndof = 2 * S.nn
    rr, cc, vv = [], [], []
    Rr = np.array([[0, 0], [1, 0], [0, 1]], float)
    for (t, u, w) in S.out_edges:
        ref = Rr[u][None, :] + R.GL[:, None] * (Rr[w] - Rr[u])[None, :]
        A = S.verts[S.tris[t, u]]; B = S.verts[S.tris[t, w]]
        L = np.linalg.norm(B - A); tv = (B - A) / L; n = np.array([tv[1], -tv[0]])
        if np.dot(n, S.verts[S.tris[t]].mean(0) - A) > 0:
            n = -n
        we = R.GLW * L
        ph = R.basis(ref[:, 0], ref[:, 1]); dX, dY = R.dbasis(ref[:, 0], ref[:, 1])
        Ji = S.Jinv[t]
        g = [Ji[0, 0] * dX + Ji[1, 0] * dY, Ji[0, 1] * dX + Ji[1, 1] * dY]
        d = dofs(S.ids[t])
        for c in range(2):
            for dd in range(2):
                blk = -np.einsum('g,ga,gb->ab', we, ph, g[c]) * n[dd]
                rr.append(np.repeat(d[:, c], nl)); cc.append(np.tile(d[:, dd], nl)); vv.append(blk.ravel())
    if not rr:
        return sp.csr_matrix((ndof, ndof))
    return sp.coo_matrix((np.concatenate(vv), (np.concatenate(rr), np.concatenate(cc))), shape=(ndof, ndof)).tocsr()


# ---------------------------------------------------------------- solver
def _pardiso_saddle(Jff, Btf, Bf, M, eps):
    """factorise the REGULARISED saddle matrix [[J, -Bt^T], [B, eps M]] with PARDISO; returns
    (solver, K_reg, K0) where K0 is the unregularised matrix used for refinement residuals."""
    import pypardiso
    K = sp.bmat([[Jff, -Btf.T], [Bf, eps * M]], format="csr"); K.sort_indices()
    K0 = sp.bmat([[Jff, -Btf.T], [Bf, None]], format="csr")
    slv = pypardiso.PyPardisoSolver(mtype=11)
    slv.set_iparm(11, 1); slv.set_iparm(13, 1)
    slv.factorize(K)
    return slv, K, K0


def solve(S, sysm, ext, method, g, nu=1.0, conv=None, Tout=None, eps=1e-8, tol=1e-12, maxit=40,
          verbose=False, x0=None):
    """nu N_h(u,v) [+ c(u;u,v)] [+ Tout] - p^T Bt v = (f,v),  (q, div u) = 0,  u = g on Dirichlet dofs.
    Linear (conv=None): one factorisation + iterative refinement.  Nonlinear: Newton from the Stokes
    solution, exact Jacobian, stopping at ||R||/||R_0|| < tol or ||du|| < 1e-13 ||u||.
    Returns u, p, info (relres of the UNregularised system, divres, Newton history)."""
    A, B, F = sysm["A"], sysm["B"], sysm["F"]
    Bt = coupling(sysm, ext, method)
    Ms = sp.block_diag(list(sysm["Ml"]), format="csr")
    ndof = A.shape[0]
    dD = dofs(S.dnodes).ravel(); free = np.setdiff1d(np.arange(ndof), dD)
    u = np.zeros(ndof); u[dD] = g(S.xy[S.dnodes, 0], S.xy[S.dnodes, 1]).ravel()
    L = nu * A + (Tout if Tout is not None else 0)
    L = L.tocsr()
    Btf = Bt[:, free]; Bf = B[:, free]; nf = len(free)
    p = np.zeros(B.shape[0])
    E = _elem(S) if conv else None

    def residual(u, p):
        Ru = L @ u - Bt.T @ p - F
        if conv:
            Ru = Ru + conv_residual(S, u, conv, E)
        return np.concatenate([Ru[free], B @ u])

    hist = []
    nb = np.linalg.norm(residual(u, p))                   # scale: residual at (lift of g, 0)
    if x0 is not None:                                    # continuation: start Newton from x0
        u = x0[0].copy(); p = x0[1].copy()
        u[dD] = g(S.xy[S.dnodes, 0], S.xy[S.dnodes, 1]).ravel()
    else:
        # Stokes step (the whole solve when conv is None): factorise once, iterative refinement
        slv, K, K0 = _pardiso_saddle(L[free][:, free], Btf, Bf, Ms, eps)
        for it in range(25):
            r = -residual(u, p)
            if np.linalg.norm(r) < 1e-15 * nb:
                break
            dx = slv.solve(K, r)
            u[free] += dx[:nf]; p += dx[nf:]
            if conv:
                break                                     # nonlinear: Stokes guess only
        slv.free_memory(everything=True); del K, K0
    if conv:
        for it in range(maxit):
            R_ = residual(u, p); rn = np.linalg.norm(R_)
            hist.append(float(rn / nb))
            if verbose:
                print(f"  newton {it} |R|/|R0| = {rn/nb:.3e}", flush=True)
            if rn < tol * nb:
                break
            import time as _t; _t0 = _t.time()
            _, Jc = conv_system(S, u, conv, E)
            if verbose:
                print(f"    jacobian {_t.time()-_t0:.1f}s", flush=True)
            J = (L + Jc).tocsr()
            slv, K, K0 = _pardiso_saddle(J[free][:, free], Btf, Bf, Ms, eps)
            dx = slv.solve(K, -R_)
            for _ in range(4):                            # refine the linear solve (unregularised)
                rr = -R_ - K0 @ dx
                dx = dx + slv.solve(K, rr)
            linres = np.linalg.norm(-R_ - K0 @ dx) / rn
            if verbose:
                print(f"    linear residual {linres:.2e}", flush=True)
            slv.free_memory(everything=True); del K, K0
            lam = 1.0                                     # backtracking on ||R||
            for _ in range(8):
                un = u.copy(); un[free] += lam * dx[:nf]; pn = p + lam * dx[nf:]
                if np.linalg.norm(residual(un, pn)) < (1 - 1e-4 * lam) * rn:
                    break
                lam /= 2
            else:
                hist.append(float(rn / nb)); break         # no descent: give up (reported as not converged)
            u, p = un, pn
            if lam == 1.0 and np.linalg.norm(dx[:nf]) < 1e-14 * np.linalg.norm(u):
                hist.append(float(np.linalg.norm(residual(u, p)) / nb)); break
        converged = hist[-1] < 1e-9 if hist else True
    else:
        converged = True
    relres = float(np.linalg.norm(residual(u, p)) / nb)
    Bu = B @ u
    dres = float(np.sqrt(abs(Bu @ (sysm["Minv"] @ Bu))))
    return u, p, dict(relres=relres, divres=dres, newton=hist, newton_its=len(hist), converged=bool(converged))


# ---------------------------------------------------------------- diagnostics
def fluxes(ext, u):
    return ext["phi"] @ u


def ph_at(S, p, X):
    """values of the discontinuous pressure p at physical quadrature points X (nt,nq,2) -> (nt,nq)."""
    R = S.R
    P = p.reshape(len(S.tris), R.npl)
    return np.einsum('qj,tj->tq', R.QP_Q, P)


def pressure_errors(S, p, pex):
    R = S.R
    X = S.phys(R.QX, R.QY)
    wq = R.QW[None, :] * np.abs(S.detJ)[:, None]
    d = pex(X[..., 0], X[..., 1]) - ph_at(S, p, X)
    area = wq.sum(); c = (wq * d).sum() / area
    return dict(pL2_raw=float(np.sqrt((wq * d**2).sum())), pL2=float(np.sqrt((wq * (d - c)**2).sum())),
                p_offset=float(c))


def gamma_stats(S, pex, curves=None):
    """p-bar (global) and G = ||p - pbar||_{L2(Gamma)} on the EXACT curves (4000-pt trapezoid) if
    given, else on Gamma_h; plus per-component means."""
    if curves is not None:
        vals, ws, comp = [], [], []
        for i, c in enumerate(curves):
            t = np.linspace(0, 2 * np.pi, 4001)[:-1]
            X = c.pos(t); w = np.linalg.norm(c.dpos(t), axis=1) * 2 * np.pi / 4000
            vals.append(pex(X[:, 0], X[:, 1])); ws.append(w); comp.append(np.full(len(t), i))
    else:
        vals, ws, comp = [], [], []
        for (t, ref, Xe, we, n), i in zip(S.edge_data(), S.bcomp):
            vals.append(pex(Xe[:, 0], Xe[:, 1])); ws.append(we); comp.append(np.full(len(we), i))
    v = np.concatenate(vals); w = np.concatenate(ws); cm = np.concatenate(comp)
    pbar = (w * v).sum() / w.sum()
    G = np.sqrt((w * (v - pbar) ** 2).sum())
    per = [float((w[cm == i] * v[cm == i]).sum() / w[cm == i].sum()) for i in range(cm.max() + 1)]
    lens = [float(w[cm == i].sum()) for i in range(cm.max() + 1)]
    return dict(pbar=float(pbar), G=float(G), pbar_i=per, len_i=lens)


def leak_profile(S, u, pex, pbar, scale):
    """relative L2(Gamma_h) misfit of u_h.n against scale*(p - pbar), and tangential/normal ratio."""
    R = S.R; num = den = tn = nn = 0.0
    for (t, ref, Xe, we, n) in S.edge_data():
        ph = R.basis(ref[:, 0], ref[:, 1]); uh = ph @ u[dofs(S.ids[t])]
        un = uh @ n; ut = uh @ np.array([-n[1], n[0]])
        pr = scale * (pex(Xe[:, 0], Xe[:, 1]) - pbar)
        num += (we * (un - pr) ** 2).sum(); den += (we * pr ** 2).sum()
        tn += (we * ut ** 2).sum(); nn += (we * un ** 2).sum()
    return dict(leak_misfit=float(np.sqrt(num / den)) if den > 0 else float("nan"),
                tang_over_norm=float(np.sqrt(tn / nn)) if nn > 0 else float("nan"))


def forces(S, sysm, u, p, nu):
    """force exerted by the fluid on the obstacle(s) (sum over Gamma_h):
       boundary:  F = -int_{Gamma_h} (-p_h n + nu D(u_h) n) ds          (n out of Omega_h)
       volume  :  F_j = -omega(chi_j),  omega(v) = nu/2 (Du,Dv) + c1(u;u,v) - (p, div v) - (f, v),
                  chi_j = e_j at every node on Gamma_h, 0 elsewhere (Gjerde-Scott (7))."""
    R = S.R; nl = R.nl
    Fb = np.zeros(2)
    gnodes = set()
    Rr = np.array([[0, 0], [1, 0], [0, 1]], float)
    for (t, ref, Xe, we, n), (tt, a, b) in zip(S.edge_data(), S.bedges):
        ph = R.basis(ref[:, 0], ref[:, 1]); dX, dY = R.dbasis(ref[:, 0], ref[:, 1])
        Ji = S.Jinv[t]
        gxe = Ji[0, 0] * dX + Ji[1, 0] * dY; gye = Ji[0, 1] * dX + Ji[1, 1] * dY
        U = u[dofs(S.ids[t])]
        Gx = gxe @ U; Gy = gye @ U                           # (ng,2): d_x u_c, d_y u_c
        grad = np.stack([Gx, Gy], -1)                        # [g, c, d] = d_d u_c
        Dm = grad + np.transpose(grad, (0, 2, 1))
        pv = R.pbasis(ref[:, 0], ref[:, 1]) @ p[t * R.npl:(t + 1) * R.npl]
        trac = -pv[:, None] * n[None, :] + nu * np.einsum('gcd,d->gc', Dm, n)
        Fb -= (we[:, None] * trac).sum(0)
        opp = 3 - a - b
        for kk, (l1, l2, l0) in enumerate(R.LBARY):
            if (l0, l1, l2)[opp] == 0:
                gnodes.add(int(S.ids[t, kk]))
    gn = np.array(sorted(gnodes))
    Nu = conv_residual(S, u, "c1")
    Fv = np.zeros(2)
    for j in range(2):
        chi = np.zeros(2 * S.nn); chi[2 * gn + j] = 1.0
        om = nu * (chi @ (sysm["Avol"] @ u)) + chi @ Nu - p @ (sysm["B"] @ chi) - chi @ sysm["F"]
        Fv[j] = -om
    return dict(Fb=Fb, Fv=Fv)


def point_pressure(S, p, pts):
    """p_h at points (averaged over all triangles whose closure contains the point)."""
    out = []
    for P in pts:
        vals = []
        lam = np.einsum('tij,tj->ti', S.Jinv, np.asarray(P)[None, :] - S.P0)   # ref coords
        ok = (lam[:, 0] >= -1e-10) & (lam[:, 1] >= -1e-10) & (lam.sum(1) <= 1 + 1e-10)
        for t in np.where(ok)[0]:
            q = S.R.pbasis(lam[t:t + 1, 0], lam[t:t + 1, 1])[0]
            vals.append(q @ p[t * S.R.npl:(t + 1) * S.R.npl])
        out.append(float(np.mean(vals)) if vals else float("nan"))
    return out
