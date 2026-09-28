"""
notes/v6/unifmu2: gamma-uniform H^1 lower bound for the leak field (REPORT.tex, Theorem 1).

Builds, on the production meshes of code/svn.py, the three-triangle test fields
    v_e = curl psi_e,  psi_e = l0 l1 l2^2 (l1 - l0)  on T_e = [x0, x1, z]   (x0, x1 on Gamma_h)
                       psi_e = a0 m_x0^2 m_z^2 m_y0  on T'_0 = [x0, z, y0]   (neighbour across [x0, z])
                       psi_e = a1 m_x1^2 m_z^2 m_y1  on T'_1 = [x1, z, y1]   (neighbour across [x1, z])
(a0, a1 fixed by C^1 matching), checks that v_e in Y_h^1 (continuous P4, div-free, zero on Gamma_h and
outside the patch, d_n v_e . tau_e odd on e), and computes
    kappa_e = int_e s d_n v_e . tau_e ds / (|e| ||grad v_e||)                   (shape constant, (M3b))
    cert    = <U~, d_n v> / (||D v|| + (sum_e C_PT(T_e)^2 |e| ||d_n v||_e^2)^{1/2}),  v = sum_e c_e v_e
(certified lower bound for ||grad(X - U~)||, every mu; C_PT explicit, see REPORT) and, with GS solves,
the identity  a_h(W, v) - <W, d_n v> = <U~, d_n v> - a_h(U~, v)  (W = X - U~) and ||grad W||.
Usage:  python lb_patch.py N [N ...]        (MUS=... penalty list; NOSOLVE=1 skips the GS solves)
"""
import os, sys, json, time
for v_ in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v_, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
import numpy as np, scipy.sparse as sp
import svn, h1_limit

MUS = [float(m) for m in os.environ.get("MUS", "1e1,3e1,1e2,1e3,1e4,1e8").split(",")]
NORMU = h1_limit.G * 2.7565157
NQ = 12
GQ, GQW = np.polynomial.legendre.leggauss(NQ); GQ = (GQ + 1) / 2; GQW = GQW / 2
R3 = np.array([[0, 0], [1, 0], [0, 1]], float)


def pfun(X, Y):
    return X ** 2 * Y - Y + X / 2


def bary(P):
    """barycentric coords of triangle P (3x2): returns (grad lambda_k (3x2), function lam(x))"""
    M = np.vstack([P.T, np.ones(3)])            # [x;y;1] = M lam
    Mi = np.linalg.inv(M)                       # lam = Mi [x;y;1]
    G = Mi[:, :2]
    return G, (lambda X: (Mi @ np.vstack([X.T, np.ones(len(X))])).T)


def curl_poly(terms, P, X):
    """v = curl psi at points X, psi = sum_c coef * prod lam_k^p_k (terms: list of (coef, (p0,p1,p2)))"""
    G, lamf = bary(P)
    L = lamf(X)
    gx = np.zeros(len(X)); gy = np.zeros(len(X))
    for coef, p in terms:
        for k in range(3):
            if p[k] == 0:
                continue
            q = list(p); q[k] -= 1
            f = coef * p[k] * np.prod([L[:, j] ** q[j] for j in range(3)], axis=0)
            gx += f * G[k, 0]; gy += f * G[k, 1]
    return np.stack([gy, -gx], -1)             # curl psi = (psi_y, -psi_x)


def patch_fields(S):
    """for every Gamma_h edge: list of (triangle, terms, local vertex order) defining psi_e"""
    tris = S.tris; V = S.verts
    tri_of = {}
    for t, g in enumerate(tris):
        for a in range(3):
            for b in range(a + 1, 3):
                tri_of.setdefault(tuple(sorted((g[a], g[b]))), []).append(t)
    out = []
    for (t, u, w) in S.bedges:
        g = tris[t]; x0, x1 = g[u], g[w]; z = g[3 - u - w]
        pieces = [(t, [x0, x1, z], [(1.0, (1, 2, 2)), (-1.0, (2, 1, 2))])]   # l0 l1^2 l2^2 - l0^2 l1 l2^2
        Pe = V[[x0, x1, z]]; Ge, _ = bary(Pe)
        for (xa, sgn, kE) in ((x0, -1.0, 1), (x1, +1.0, 0)):
            # on edge [xa, z] the T_e piece is sgn * l_a^2 l_E l_z^2 with l_E = lambda_kE (vanishing there)
            nb = [s for s in tri_of[tuple(sorted((xa, z)))] if s != t]
            assert len(nb) == 1
            s = nb[0]; y = [q for q in tris[s] if q not in (xa, z)][0]
            assert y not in S.circ_set, "neighbour triangle has a wall edge: (M1') violated"
            Gs, _ = bary(V[[xa, z, y]])
            # match gradients on [xa,z]: sgn * grad l_kE = alpha * grad m_y
            alpha = sgn * (Ge[kE] @ Gs[2]) / (Gs[2] @ Gs[2])
            pieces.append((s, [xa, z, y], [(alpha, (2, 2, 1))]))
        out.append(dict(t=t, x0=x0, x1=x1, z=z, pieces=pieces))
    return out


def field_vector(S, patch):
    """dof vector of v_e (Lagrange interpolation, exact since v_e is P4) + continuity defect"""
    vec = np.zeros(2 * S.nn); seen = {}
    defect = 0.0
    for (s, verts, terms) in patch["pieces"]:
        P = S.verts[verts]
        Xn = S.P0[s] + svn.REFNODES @ S.J[s].T          # 15 Lagrange nodes of triangle s
        vals = curl_poly(terms, P, Xn)
        for k, node in enumerate(S.ids[s]):
            if node in seen:
                defect = max(defect, np.abs(seen[node] - vals[k]).max())
            seen[node] = vals[k]
            vec[2 * node:2 * node + 2] = vals[k]
    # nodes of the patch boundary that are shared with triangles outside the patch must carry 0
    patch_tris = {p[0] for p in patch["pieces"]}
    for t in range(len(S.tris)):
        if t in patch_tris:
            continue
        for node in S.ids[t]:
            if node in seen:
                defect = max(defect, np.abs(seen[node]).max())
    return vec, defect


def edge_tab(S):
    E = []
    for (t, u, w) in S.bedges:
        A = S.verts[S.tris[t, u]]; B = S.verts[S.tris[t, w]]
        if np.cross(A, B) < 0:
            u, w = w, u; A, B = B, A
        L = np.linalg.norm(B - A); tau = (B - A) / L
        n = -np.array([tau[1], -tau[0]])                       # outward of Omega_h (into the polygon)
        ref = R3[u][None, :] + GQ[:, None] * (R3[w] - R3[u])[None, :]
        ph = svn.basis(ref[:, 0], ref[:, 1]); dX, dY = svn.dbasis(ref[:, 0], ref[:, 1]); Ji = S.Jinv[t]
        gxe = Ji[0, 0] * dX + Ji[1, 0] * dY; gye = Ji[0, 1] * dX + Ji[1, 1] * dY
        X = A[None, :] + GQ[:, None] * (B - A)[None, :]
        # explicit trace-Poincare constant: H ||w - wbar_T||_e^2 <= 2 ||w-wbar||_T^2 + 2 d ||w-wbar||_T ||grad w||_T
        P = S.verts[S.tris[t]]; area = abs(S.detJ[t]) / 2; H = 2 * area / L
        d = max(np.linalg.norm(P[i] - P[j]) for i in range(3) for j in range(3))
        CPT = np.sqrt((2 * (d / np.pi) ** 2 + 2 * d * d / np.pi) / (H * L))
        E.append(dict(t=t, L=L, tau=tau, n=n, ph=ph, dn=n[0] * gxe + n[1] * gye, w=GQW * L,
                      s=(GQ - 0.5) * L, X=X, CPT=CPT))
    return E


def dn_on_edges(S, E, vec):
    return [Ee["dn"] @ vec[svn.dofs(S.ids[Ee["t"]])] for Ee in E]


def grad_gram(S):
    nt = len(S.tris)
    gx, gy = S.grads(svn.QX, svn.QY)
    wq = svn.QW[None, :] * np.abs(S.detJ)[:, None]
    K = np.einsum('tq,tqa,tqb->tab', wq, gx, gx) + np.einsum('tq,tqa,tqb->tab', wq, gy, gy)
    loc = np.zeros((nt, 15, 2, 15, 2)); loc[:, :, 0, :, 0] = K; loc[:, :, 1, :, 1] = K
    D = svn.dofs(S.ids).reshape(nt, 30)
    rows = np.repeat(D, 30, axis=1).ravel(); cols = np.tile(D, (1, 30)).ravel()
    return sp.coo_matrix((loc.reshape(nt, 30, 30).ravel(), (rows, cols)), shape=(2 * S.nn, 2 * S.nn)).tocsr()


def aU(S, Uf, vec):
    """a_h(U~, v) = 1/2 (D U~, D v) by quadrature"""
    gx, gy = S.grads(svn.QX, svn.QY)
    wq = svn.QW[None, :] * np.abs(S.detJ)[:, None]
    X = S.phys(svn.QX, svn.QY)
    G = Uf.gu(X[..., 0], X[..., 1])                    # G[...,c,d] = d_d U_c
    Vl = vec[svn.dofs(S.ids)]                          # (nt,15,2)
    vx = np.einsum('tqa,tac->tqc', gx, Vl); vy = np.einsum('tqa,tac->tqc', gy, Vl)   # d_x v_c, d_y v_c
    Gv = np.stack([vx, vy], -1)                        # [..., c, d]
    DU = G + np.swapaxes(G, -1, -2); Dv = Gv + np.swapaxes(Gv, -1, -2)
    return 0.5 * (wq * (DU * Dv).sum((-1, -2))).sum()


def Dnorm(S, vec):
    gx, gy = S.grads(svn.QX, svn.QY)
    wq = svn.QW[None, :] * np.abs(S.detJ)[:, None]
    Vl = vec[svn.dofs(S.ids)]
    vx = np.einsum('tqa,tac->tqc', gx, Vl); vy = np.einsum('tqa,tac->tqc', gy, Vl)
    Gv = np.stack([vx, vy], -1); Dv = Gv + np.swapaxes(Gv, -1, -2)
    return np.sqrt((wq * (Dv ** 2).sum((-1, -2))).sum())


def run(N):
    t0 = time.time()
    S = svn.Space(*svn.make_mesh(N))
    S.circ_set = set(svn.make_mesh(N)[2])
    E = edge_tab(S); tE = {Ee["t"]: i for i, Ee in enumerate(E)}
    H = grad_gram(S)
    sysA = svn.assemble(S, 1.0, lambda X, Y: np.zeros(np.shape(X) + (2,)), 0, want=("A", "B"))
    B = sysA["B"]; Minv = sysA["Minv"]; Avol = sysA["Avol"]
    Uf = h1_limit.LimitField()
    gam_nodes = set()
    for Ee in E:
        pass
    import strong_bc
    gnodes = strong_bc.gamma_nodes(S)
    # Ubar.tau_e = c_e s + E_e : c_e by least squares in s (as in unifmu.py)
    Ut = [Uf.u(Ee["X"][:, 0], Ee["X"][:, 1]) for Ee in E]
    ce = np.array([(Ee["w"] * Ee["s"] * (U @ Ee["tau"])).sum() / (Ee["L"] ** 3 / 12) for Ee, U in zip(E, Ut)])
    Eres = max(np.abs(U @ Ee["tau"] - c * Ee["s"]).max() for Ee, U, c in zip(E, Ut, ce))
    hG = S.hGamma
    patches = patch_fields(S)
    kap = []; defects = []; divres = 0.0; gres = 0.0; oddres = 0.0; nres = 0.0
    vtot = np.zeros(2 * S.nn)
    for P, c in zip(patches, ce):
        vec, dfc = field_vector(S, P)
        defects.append(dfc)
        nrm = np.sqrt(vec @ (H @ vec)); vec /= nrm
        Bv = B @ vec; divres = max(divres, np.sqrt(abs(Bv @ (Minv @ Bv))))
        gres = max(gres, np.abs(vec[svn.dofs(np.array(sorted(gnodes))).ravel()]).max())
        i = tE[P["t"]]; Ee = E[i]
        dn = Ee["dn"] @ vec[svn.dofs(S.ids[Ee["t"]])]
        dt = dn @ Ee["tau"]; dnn = dn @ Ee["n"]
        oddres = max(oddres, abs((Ee["w"] * dt).sum()) / np.sqrt((Ee["w"] * dt ** 2).sum() * Ee["L"]))
        nres = max(nres, np.abs(dnn).max() / np.abs(dt).max())
        m = (Ee["w"] * Ee["s"] * dt).sum()
        kap.append(abs(m) / Ee["L"])
        vtot += c * np.sign(m) * vec
    kap = np.array(kap)
    # certificate
    dn = dn_on_edges(S, E, vtot)
    ell = sum((Ee["w"][:, None] * U * d).sum() for Ee, U, d in zip(E, Ut, dn))
    ell_lead = sum(c * (Ee["w"] * Ee["s"] * (d @ Ee["tau"])).sum() for Ee, c, d in zip(E, ce, dn))
    gradv = np.sqrt(vtot @ (H @ vtot)); Dv = Dnorm(S, vtot)
    Bnd = np.sqrt(sum(Ee["CPT"] ** 2 * Ee["L"] * (Ee["w"][:, None] * d ** 2).sum() for Ee, d in zip(E, dn)))
    aUv = aU(S, Uf, vtot)
    cert = (ell - abs(aUv)) / (Dv + Bnd)
    G2 = sum(c ** 2 * Ee["L"] for Ee, c in zip(E, ce))
    rec = dict(N=N, h=S.h, hG=hG, ntri_patch=3, max_C1_defect=max(defects), max_div=divres, max_on_Gamma=gres,
               max_rel_edge_mean=oddres, max_rel_normal=nres,
               kappa_min=kap.min(), kappa_mean=kap.mean(), kappa_max=kap.max(),
               Eres_over_hG2=Eres / hG ** 2, G2_riemann=G2, G2=h1_limit.G ** 2,
               ell=ell, ell_lead=ell_lead, aUv=aUv, gradv=gradv, Dv=Dv, trace_term=Bnd,
               CPT_max=max(Ee["CPT"] for Ee in E),
               cert=cert, cert_over_hG12=cert / np.sqrt(hG), cert_rel=cert / NORMU / np.sqrt(hG))
    print(json.dumps({k: float(f"{v:.6g}") if isinstance(v, (float, np.floating)) else v for k, v in rec.items()}),
          flush=True)
    if os.environ.get("NOSOLVE"):
        return
    exP = svn.make_exact("P")
    exU = dict(u=Uf.u, gu=Uf.gu)
    for mu in MUS:
        t1 = time.time()
        sysm = svn.assemble(S, mu, exP["f"], 0)
        uh, ph, dres, relres = svn.solve_lu(S, sysm, exP["u"], 0)
        X = -(mu / S.h) * uh
        Er = svn.errors(S, X, exU, mu)
        gW = Er["H1"]
        # identity check: a_h(X,v) - <X, d_n v> should vanish (v in Y_h); then
        # a_h(W,v) - <W,d_n v> = -a_h(U~,v) + <U~, d_n v>
        aXv = vtot @ (Avol @ X)                     # Avol assembles 1/2 (D.,D.)
        dnX = sum((Ee["w"][:, None] * (Ee["ph"] @ X[svn.dofs(S.ids[Ee["t"]])]) * d).sum() for Ee, d in zip(E, dn))
        lhsW = (aXv - aUv) - (dnX - ell)
        r = dict(N=N, mu=mu, gamma=mu * hG / S.h, gradW=gW, gradW_over_hG12=gW / np.sqrt(hG), d=gW / NORMU,
                 Yrow_residual=(aXv - dnX) / ell, identity_lhs=lhsW, identity_rhs=ell - aUv,
                 Lambda_eff=abs(lhsW) / (gW * gradv), cert_over_gradW=cert / gW, secs=time.time() - t1)
        print(json.dumps({k: float(f"{v:.6g}") if isinstance(v, (float, np.floating)) else v for k, v in r.items()}),
              flush=True)
    print(f"# N={N} {time.time() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    for N in (int(a) for a in sys.argv[1:]):
        run(N)
