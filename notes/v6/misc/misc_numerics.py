"""notes/v6/misc numerics (k = 4, production meshes of code/svn.py, CNS* = theta 1, GS = theta 0).

For each N:
 (a)  Stokes-vs-unconstrained comparison (Prop. A1 of REPORT.tex): u^L in V_h^R solves the Nitsche
      problem WITHOUT divergence constraint,  N_h(u^L, v) = (f~, v) - B*(pi_h p*, v)  for all v,
      same data on dR.  Reports ||e_p - mean||, ||div u^L||, ||grad(u_h - u^L)||, energy norm of u_h - u^L.
 (a') smooth volume moments (e_p - mean, g), g in {x, y, x^2 - y^2}  (Lemma A3: O(hG^2)).
 (b)  boundary pressure moments P_phi = <e_p - pbar_G(e_p), phi(x*)>_{Gamma_h} and normal-slip moments
      S_phi = <e_u . n, phi(x*)>; checks (mu/h) S_phi + P_phi = O(h^{3/2}) (Prop. B1).
 (e)  J^chi_h for psi = e1, e2, x^perp (CNS* and GS) vs J^N_h and the predicted torque defect
      (h/mu) int_Gamma sigma(u,p) n . d_n psi  (Prop. E1).
usage: python3 misc_numerics.py KIND MU N1,N2,...     (output: one JSON line per (N, method))
"""
import os, sys, json, time
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../code"))
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import svn

KIND = sys.argv[1]; MU = float(sys.argv[2]); NS = [int(s) for s in sys.argv[3].split(",")]


def make_exact_Q(amp):
    """test Q: u = curl((r^2-1)^2/4) = ((r^2-1) y, -(r^2-1) x) in P_3, p = amp (x^2 y - y + x/2) in P_3.
    Interior discretisation is exact (u in P_4, p in P_3, data on dR interpolated exactly): every error
    comes from the polygonal boundary.  Wall shear W = -2 != 0."""
    import sympy as sy
    x, y = sy.symbols('x y', real=True)
    r2 = x**2 + y**2
    if SKEW:   # test S: psi = (r^2-1)^2 (1+x)/4, u = curl psi in P_4 (no rotational symmetry)
        psi = (r2 - 1)**2 * (1 + x) / 4
        u = [sy.expand(sy.diff(psi, y)), sy.expand(-sy.diff(psi, x))]
    else:
        u = [(r2 - 1) * y, -(r2 - 1) * x]
    p = amp * (x**2 * y - y + x / 2)
    f = [sy.expand(-(sy.diff(u[i], x, 2) + sy.diff(u[i], y, 2)) + sy.diff(p, [x, y][i])) for i in range(2)]
    gu = [[sy.diff(u[i], v) for v in (x, y)] for i in range(2)]
    assert sy.simplify(sy.diff(u[0], x) + sy.diff(u[1], y)) == 0
    U = sy.lambdify((x, y), u, 'numpy'); Fl = sy.lambdify((x, y), f, 'numpy')
    GU = sy.lambdify((x, y), gu, 'numpy'); P = sy.lambdify((x, y), p, 'numpy')

    def vec(fn):
        def g(X, Y):
            return np.stack([np.broadcast_to(np.asarray(v, float), np.shape(X)) for v in fn(X, Y)], -1)
        return g

    def gmat(X, Y):
        vals = GU(X, Y)
        return np.stack([np.stack([np.broadcast_to(np.asarray(vals[i][j], float), np.shape(X))
                                   for j in range(2)], -1) for i in range(2)], -2)
    return dict(u=vec(U), f=vec(Fl), gu=gmat,
                p=lambda X, Y: np.broadcast_to(np.asarray(P(X, Y), float), np.shape(X)))


SKEW = KIND.startswith("S")
ex = make_exact_Q(float(KIND[1:] or 1)) if KIND[0] in "QS" else svn.make_exact(KIND)
PSIS = {"e1": (lambda X: np.stack([np.ones_like(X[:, 0]), 0 * X[:, 0]], -1), np.zeros((2, 2))),
        "e2": (lambda X: np.stack([0 * X[:, 0], np.ones_like(X[:, 0])], -1), np.zeros((2, 2))),
        "rot": (lambda X: np.stack([-X[:, 1], X[:, 0]], -1), np.array([[0., -1.], [1., 0.]]))}
PHIS = {"cos1": lambda th: np.cos(th), "sin1": lambda th: np.sin(th), "cos2": lambda th: np.cos(2 * th),
        "sin3": lambda th: np.sin(3 * th), "cos4": lambda th: np.cos(4 * th)}
GS_ = {"x": lambda X, Y: X, "y": lambda X, Y: Y, "x2y2": lambda X, Y: X**2 - Y**2}


def exact_J(psi, dpsi_mat=None, n_quad=4000):
    """J_psi = int_Gamma sigma(u,p) n . psi, n = -x (into D), sigma = grad u + grad u^T - p I.
    With dpsi_mat given, integrates sigma n . (dpsi_mat n) instead."""
    th = (np.arange(n_quad) + 0.5) * 2 * np.pi / n_quad
    X = np.stack([np.cos(th), np.sin(th)], -1); n = -X
    G = ex["gu"](X[:, 0], X[:, 1]); P = ex["p"](X[:, 0], X[:, 1])
    sig = G + np.transpose(G, (0, 2, 1)) - P[:, None, None] * np.eye(2)[None]
    t = np.einsum('qij,qj->qi', sig, n)
    w = psi(X) if dpsi_mat is None else n @ dpsi_mat.T
    return float((t * w).sum(-1).sum() * 2 * np.pi / n_quad)


def run(N):
    t0 = time.time()
    verts, tris, circ, outer = svn.make_mesh(N)
    S = svn.Space(verts, tris, circ, outer)
    nt = len(tris); h = S.h
    sysm = svn.assemble(S, MU, ex["f"], None)
    A, Avol, B, F = sysm["A"], sysm["Avol"], sysm["B"], sysm["F"]
    ndof = A.shape[0]
    dD = svn.dofs(S.dnodes).ravel(); free = np.setdiff1d(np.arange(ndof), dD)
    ug = np.zeros(ndof); ug[dD] = ex["u"](S.xy[S.dnodes, 0], S.xy[S.dnodes, 1]).ravel()
    X = S.phys(svn.QX, svn.QY); wq = svn.QW[None, :] * np.abs(S.detJ)[:, None]
    gx, gy = S.grads(svn.QX, svn.QY); Dm = svn.dofs(S.ids)
    ED = S.edge_data()
    # p* = p~ - pbar_Gamma(p~) and its L2 projection
    num = per = 0.0
    for (t, ref, Xe, we, n) in ED:
        num += (we * ex["p"](Xe[:, 0], Xe[:, 1])).sum(); per += we.sum()
    pbarex = num / per
    pst = ex["p"](X[..., 0], X[..., 1]) - pbarex
    rhsP = np.einsum('tq,qj,tq->tj', wq, svn.Q3_Q, pst)
    pip = np.linalg.solve(sysm["Ml"], rhsP[..., None])[..., 0].ravel()

    def pvals(p):
        return np.einsum('qj,tj->tq', svn.Q3_Q, p.reshape(nt, 10))

    def pedge(p, t, ref):
        return svn._mono(svn.MON3, ref[:, 0], ref[:, 1]) @ p[10 * t:10 * t + 10]

    def edge_fields(u, t, ref, n):
        ph = svn.basis(ref[:, 0], ref[:, 1]); dX, dY = svn.dbasis(ref[:, 0], ref[:, 1]); Ji = S.Jinv[t]
        gxe = Ji[0, 0] * dX + Ji[1, 0] * dY; gye = Ji[0, 1] * dX + Ji[1, 1] * dY
        Ut = u[svn.dofs(S.ids[t])]
        return ph @ Ut, (n[0] * gxe + n[1] * gye) @ Ut, gxe @ Ut, gye @ Ut

    recs = []
    for theta in (1, 0):
        Bt = svn.coupling(sysm, theta)
        K = sp.bmat([[A[free][:, free], -Bt[:, free].T], [B[:, free], None]], format="csc")
        lu = spl.splu(K, permc_spec="COLAMD", diag_pivot_thresh=1.0)

        def solve(rhs):
            x = lu.solve(rhs)
            for _ in range(3):
                x = x + lu.solve(rhs - K @ x)
            return x
        rhs = np.concatenate([F[free] - A[free][:, dD] @ ug[dD], -B[:, dD] @ ug[dD]])
        x = solve(rhs)
        u = ug.copy(); u[free] = x[:len(free)]; p = x[len(free):]
        rec = dict(N=N, kind=KIND, mu=MU, method=("GS", "CNS*")[theta], h=h, hG=S.hGamma)
        # ---------------- pressure error and moments
        ep = pst - pvals(p)
        area = wq.sum(); epm = (wq * ep).sum() / area
        rec["ep_L2opt"] = float(np.sqrt((wq * (ep - epm) ** 2).sum()))
        for gname, g in GS_.items():
            gv = g(X[..., 0], X[..., 1])
            rec["vol_" + gname] = float((wq * (ep - epm) * gv).sum())
        # boundary: p̄_Γ(e_p), P_phi, S_phi
        epG = 0.0
        eds = []
        for (t, ref, Xe, we, n) in ED:
            epe = ex["p"](Xe[:, 0], Xe[:, 1]) - pbarex - pedge(p, t, ref)
            uhe, dnh, _, _ = edge_fields(u, t, ref, n)
            eue = ex["u"](Xe[:, 0], Xe[:, 1]) - uhe
            eds.append((Xe, we, n, epe, eue, uhe, dnh, t, ref))
            epG += (we * epe).sum()
        epG /= per
        for pname, phi in PHIS.items():
            Pm = Sm = 0.0
            for (Xe, we, n, epe, eue, uhe, dnh, t, ref) in eds:
                ph_ = phi(np.arctan2(Xe[:, 1], Xe[:, 0]))
                Pm += (we * (epe - epG) * ph_).sum(); Sm += (we * (eue @ n) * ph_).sum()
            rec["P_" + pname] = float(Pm); rec["S_" + pname] = float(Sm)
            rec["muS+P_" + pname] = float(MU / h * Sm + Pm)
        # ---------------- (a) unconstrained comparison (CNS* coupling for p*)
        if theta == 1:
            rhsL = (F + Bt.T @ pip)[free] - A[free][:, dD] @ ug[dD]
            uL = ug.copy(); uL[free] = spl.spsolve(A[free][:, free].tocsc(), rhsL)
            d = u - uL
            U = uL[Dm]
            divL = np.einsum('tqa,ta->tq', gx, U[..., 0]) + np.einsum('tqa,ta->tq', gy, U[..., 1])
            rec["divL"] = float(np.sqrt((wq * divL ** 2).sum()))
            Dd = d[Dm]
            dx_ = np.einsum('tqa,tac->tqc', gx, Dd); dy_ = np.einsum('tqa,tac->tqc', gy, Dd)
            rec["dH1"] = float(np.sqrt((wq * (dx_ ** 2 + dy_ ** 2).sum(-1)).sum()))
            b2 = bn2 = 0.0
            for (t, ref, Xe, we, n) in ED:
                de, dnd, _, _ = edge_fields(d, t, ref, n)
                b2 += (we * (de ** 2).sum(-1)).sum(); bn2 += (we * (dnd ** 2).sum(-1)).sum()
            rec["dEnergy"] = float(np.sqrt(rec["dH1"] ** 2 + MU / h * b2 + h * bn2))
        # ---------------- (e) J^N and J^chi for psi in {e1, e2, rot}
        # omega_h(w) = a_h(u_h, w) - (p_h, div w) - (f~, w)
        def omega(w):
            return float(w @ (Avol @ u) - p @ (B @ w) - F @ w)
        for name, (psi, Jm) in PSIS.items():
            Jex = exact_J(psi)
            # J^N via flux identity: int t_h . psi + <u_h, d_n psi>
            JN = 0.0
            lam = sum((we * pedge(p, t, ref)).sum() for (Xe, we, n, epe, eue, uhe, dnh, t, ref) in eds) / per
            for (Xe, we, n, epe, eue, uhe, dnh, t, ref) in eds:
                th_ = dnh - theta * (pedge(p, t, ref) - lam)[:, None] * n[None, :] - MU / h * uhe
                JN += (we * (th_ * psi(Xe)).sum(-1)).sum() + (we * (uhe @ (Jm @ n))).sum()
            # chi_h: N_h(chi, w) + pressure = l_psi(w) = -<psi, d_n w> + mu/h <psi, w>, zero data on dR
            l = np.zeros(ndof)
            for (t, ref, Xe, we, n) in ED:
                ph = svn.basis(ref[:, 0], ref[:, 1]); dX, dY = svn.dbasis(ref[:, 0], ref[:, 1]); Ji = S.Jinv[t]
                dn = n[0] * (Ji[0, 0] * dX + Ji[1, 0] * dY) + n[1] * (Ji[0, 1] * dX + Ji[1, 1] * dY)
                ps = psi(Xe); dd = svn.dofs(S.ids[t])
                for c in range(2):
                    np.add.at(l, dd[:, c], np.einsum('g,g,ga->a', we, ps[:, c], -dn + MU / h * ph))
            xc = solve(np.concatenate([l[free], np.zeros(B.shape[0])]))
            chi = np.zeros(ndof); chi[free] = xc[:len(free)]
            Jchi = omega(chi)
            Dpred = (h / MU) * exact_J(None, Jm) if name == "rot" else 0.0
            if name == "rot":
                # zeta^1: N_h(z, w) + pressure = -<d_n psi, w> (the torque-specific part of v_psi - chi_h)
                l1 = np.zeros(ndof)
                for (t, ref, Xe, we, n) in ED:
                    ph = svn.basis(ref[:, 0], ref[:, 1]); dd = svn.dofs(S.ids[t]); Jn = Jm @ n
                    for c in range(2):
                        np.add.at(l1, dd[:, c], -Jn[c] * np.einsum('g,ga->a', we, ph))
                xz = solve(np.concatenate([l1[free], np.zeros(B.shape[0])]))
                z1 = np.zeros(ndof); z1[free] = xz[:len(free)]
                rec["-omega(zeta1)_rot"] = -omega(z1)
            rec["J_" + name] = Jex; rec["JN_" + name] = float(JN); rec["Jchi_" + name] = Jchi
            rec["Jchi-JN_" + name] = Jchi - float(JN); rec["pred_" + name] = Dpred
        rec["secs"] = time.time() - t0
        recs.append(rec)
        del lu
        print(json.dumps({k: (float(f"{v:.6g}") if isinstance(v, float) else v) for k, v in rec.items()}), flush=True)
    return recs


for N in NS:
    run(N)
