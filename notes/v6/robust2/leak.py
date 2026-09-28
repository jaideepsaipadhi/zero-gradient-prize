"""
Residual-leak law and lower-bound certificates for the CNS* pressure part W_*(phi)  (REPORT.tex, Thms R1-R4).

For the hydrostatic data (grad phi, g = 0), nu = 1, CNS* (theta = 1):
  eta   = phi - pi_T phi on Gamma_h (T = boundary triangle of the edge)
  Y     = Pi_T(eta n): L2(Gamma_h)-projection of eta*n onto the trace space T_h of V_h on Gamma_h
          (continuous piecewise-P_k vector fields) with zero net normal flux
  leak  = |W - (h/mu) Y|_Gamma / ((h/mu) |Y|_Gamma)          (Thm R2 predicts O(1/gamma))
  H1n   = |grad W| / (hG^{1/2} (h/mu) |eta|_Gamma)            (Thms R3/R4: bounded above and below)
  osc   = |(I-P0) W|_Gamma / |W|_Gamma  (P0 = edgewise mean)
  Lk    = (even k) L_k-mode certificate: v with v.n = g(m_e) L_k on every Gamma_h edge, g = grad^k phi[tau^k];
          D_Lk = <eta, v.n>/|||v|||  (Thm R1: |||W||| >= D_Lk / C)
usage: python leak.py N phi mu1,mu2,... [k]      phi in {sin, cubic, r1_4, r2_3}
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../code"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../robust"))
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import svn_k as K
from numpy.polynomial import legendre as Lg


def _r(X, Y):
    return np.sqrt(X**2 + Y**2)


PHIS = {
    "sin":   (lambda X, Y: np.sin(2 * X) * np.cos(Y),
              lambda X, Y: np.stack([2 * np.cos(2 * X) * np.cos(Y), -np.sin(2 * X) * np.sin(Y)], -1)),
    "cubic": (lambda X, Y: X**3 - 3 * X * Y**2,
              lambda X, Y: np.stack([3 * X**2 - 3 * Y**2, -6 * X * Y], -1)),
    "r1_4":  (lambda X, Y: (_r(X, Y) - 1)**4,
              lambda X, Y: np.stack([4 * (_r(X, Y) - 1)**3 * X / _r(X, Y), 4 * (_r(X, Y) - 1)**3 * Y / _r(X, Y)], -1)),
    "r2_3":  (lambda X, Y: (X**2 + Y**2 - 1)**3,
              lambda X, Y: np.stack([6 * X * (X**2 + Y**2 - 1)**2, 6 * Y * (X**2 + Y**2 - 1)**2], -1)),
}


def dtau_k(name, X, Y, tau, k):
    """grad^k phi (x)[tau,...,tau] (k-th directional derivative), by exact formulas or finite differences."""
    if name == "sin":
        out = 0
        for a in (np.array([2.0, 1.0]), np.array([2.0, -1.0])):
            s = a @ tau; arg = a[0] * X + a[1] * Y
            out = out + 0.5 * s**k * [np.sin, np.cos, lambda z: -np.sin(z), lambda z: -np.cos(z)][k % 4](arg)
        return out
    phi = PHIS[name][0]
    # high-order central difference of t -> phi(x + t tau)
    hh = 2e-2
    from math import comb
    return sum((-1)**j * comb(k, j) * phi(X + (k / 2 - j) * hh * tau[0], Y + (k / 2 - j) * hh * tau[1]) for j in range(k + 1)) / hh**k


def solve(S, A, Bt, B, F):
    ndof = A.shape[0]
    dD = K.dofs(S.dnodes).ravel(); free = np.setdiff1d(np.arange(ndof), dD)
    Km = sp.bmat([[A[free][:, free], -Bt[:, free].T], [B[:, free], None]], format="csc")
    rhs = np.concatenate([F[free], np.zeros(B.shape[0])])
    lu = spl.splu(Km, permc_spec="COLAMD", diag_pivot_thresh=1.0)
    x = lu.solve(rhs)
    for _ in range(3):
        x = x + lu.solve(rhs - Km @ x)
    rel = np.linalg.norm(Km @ x - rhs) / max(np.linalg.norm(rhs), 1e-300)
    u = np.zeros(ndof); u[free] = x[:len(free)]
    return u, x[len(free):], rel


ZERO = dict(u=lambda X, Y: np.zeros(np.shape(X) + (2,)), gu=lambda X, Y: np.zeros(np.shape(X) + (2, 2)))


def edge_info(S, phi):
    """per Gamma_h edge: GL points, weights, normal, basis values, eta."""
    R = S.R; out = []
    for (t, ref, Xe, we, n) in S.edge_data():
        ph = R.basis(ref[:, 0], ref[:, 1])
        Xq = S.phys(R.QX, R.QY)[t]; wq = R.QW * abs(S.detJ[t])
        M = np.einsum('q,qi,qj->ij', wq, R.QP_Q, R.QP_Q)
        cpi = np.linalg.solve(M, R.QP_Q.T @ (wq * phi(Xq[:, 0], Xq[:, 1])))
        eta = phi(Xe[:, 0], Xe[:, 1]) - R.pbasis(ref[:, 0], ref[:, 1]) @ cpi
        out.append(dict(t=t, X=Xe, w=we, n=n, ph=ph, eta=eta, d=K.dofs(S.ids[t])))
    return out


def trace_projection(S, sysm, E):
    """Y = Pi_T(eta n) with zero net flux; returns global dof vector."""
    Mg = (sysm["Pen"] * S.h).tocsr()
    ndof = Mg.shape[0]
    b = np.zeros(ndof); c = np.zeros(ndof)
    for e in E:
        for comp in range(2):
            np.add.at(b, e["d"][:, comp], e["ph"].T @ (e["w"] * e["eta"] * e["n"][comp]))
            np.add.at(c, e["d"][:, comp], e["ph"].T @ (e["w"] * e["n"][comp]))
    idx = np.where(np.abs(Mg).sum(1).A.ravel() > 0)[0]
    Mr = Mg[idx][:, idx]
    Kk = sp.bmat([[Mr, sp.csr_matrix(c[idx][:, None])], [sp.csr_matrix(c[idx][None, :]), None]], format="csc")
    sol = spl.spsolve(Kk, np.concatenate([b[idx], [0.0]]))
    y = np.zeros(ndof); y[idx] = sol[:-1]
    return y


def gnorms(E, u):
    """|u|_Gamma, |(I-P0)u|_Gamma, |u.n - mean|, at GL points."""
    tot = osc = 0.0
    for e in E:
        U = e["ph"] @ u[e["d"]]                       # (g, 2)
        tot += e["w"] @ (U**2).sum(1)
        m = (e["w"] @ U) / e["w"].sum()
        osc += e["w"] @ ((U - m)**2).sum(1)
    return np.sqrt(tot), np.sqrt(osc)


def Lk_field(S, E, name, k):
    """(even k) v in V_R with v.n = G_e L_k on each Gamma_h edge (exactly), G_e = grad^k phi(m_e)[tau_e^k]; see REPORT Thm R1."""
    assert k % 2 == 0
    R = S.R; ndof = 2 * S.nn; v = np.zeros(ndof)
    third = {}
    edges = []
    for (e, (t, u, w)) in zip(E, S.bedges):
        A = S.verts[S.tris[t, u]]; B = S.verts[S.tris[t, w]]
        tau = (B - A) / np.linalg.norm(B - A)
        m = (A + B) / 2
        G = dtau_k(name, np.array(m[0]), np.array(m[1]), tau, k)
        edges.append((t, u, w, A, B, tau, e["n"], float(G)))
    # vertex values: v(z).n_e = G_e for both edges at z
    vert = {}
    for (t, u, w, A, B, tau, n, G) in edges:
        for vid in (S.tris[t, u], S.tris[t, w]):
            vert.setdefault(vid, []).append((n, G))
    vval = {}
    for vid, lst in vert.items():
        assert len(lst) == 2
        Mx = np.array([lst[0][0], lst[1][0]]); rhs = np.array([lst[0][1], lst[1][1]])
        vval[vid] = np.linalg.solve(Mx, rhs)
    for (t, u, w, A, B, tau, n, G) in edges:
        th = 3 - u - w
        L = np.linalg.norm(B - A)
        ta = vval[S.tris[t, u]] @ tau; tb = vval[S.tris[t, w]] @ tau
        for kk, (l1, l2, l0) in enumerate(R.LBARY):
            lam = {0: l0, 1: l1, 2: l2}
            if lam[th] != 0:
                continue
            gid = S.ids[t, kk]
            X = S.xy[gid]
            s = np.dot(X - A, tau) / L
            if lam[u] == k:
                val = vval[S.tris[t, u]]
            elif lam[w] == k:
                val = vval[S.tris[t, w]]
            else:
                Lks = Lg.legval(2 * s - 1, [0] * k + [1])
                val = G * Lks * n + ((1 - s) * ta + s * tb) * tau
            v[2 * gid] = val[0]; v[2 * gid + 1] = val[1]
    return v


def run(N, name, mus, k=4):
    phi, gphi = PHIS[name]
    t0 = time.time()
    S = K.build(N, k)
    E = edge_info(S, phi)
    eta_n = np.sqrt(sum(e["w"] @ e["eta"]**2 for e in E))
    sysm0 = K.assemble(S, mus[0], gphi)
    Y = trace_projection(S, sysm0, E)
    Yn, Yosc = gnorms(E, Y)
    res = []
    vLk = Lk_field(S, E, name, k) if (k % 2 == 0 and name != "cubic") else None
    for mu in mus:
        sysm = K.assemble(S, mu, gphi) if mu != mus[0] else sysm0
        Bt = K.coupling(S, sysm, 1)
        u, p, rel = solve(S, sysm["A"], Bt, sysm["B"], sysm["F"])
        Er = K.errors(S, u, ZERO, mu)
        eps = S.h / mu
        Wn, Wosc = gnorms(E, u)
        dn, _ = gnorms(E, u - eps * Y)
        gam = mu * S.hGamma / S.h
        rec = dict(N=N, k=k, phi=name, mu=mu, gamma=gam, h=S.h, hG=S.hGamma, relres=rel,
                   H1=Er["H1"], energy=Er["energy"], Wgam=Wn, Wosc=Wosc, eta=eta_n, Y=Yn, Yosc=Yosc,
                   leak=dn / (eps * Yn), H1n=Er["H1"] / (S.hGamma**-0.5 * eps * eta_n), En=Er["energy"] / (eps**0.5 * eta_n),
                   H1_over_osc=Er["H1"] / (S.hGamma**-0.5 * Wosc),
                   E_n=Er["energy"] / (S.hGamma**(k + 0.5) * gam**-0.5), H1_hk=Er["H1"] / (S.hGamma**(k + 0.5) / gam))
        if vLk is not None:
            Ev = K.errors(S, vLk, ZERO, mu)
            pair = sum(e["w"] @ (e["eta"] * ((e["ph"] @ vLk[e["d"]]) @ e["n"])) for e in E)
            flux = sum(e["w"] @ ((e["ph"] @ vLk[e["d"]]) @ e["n"]) for e in E)
            rec.update(D_Lk=pair / Ev["energy"], D_Lk_n=pair / Ev["energy"] / (S.hGamma**(k + 0.5) * (1 + gam)**-0.5),
                       Lk_flux=flux, W_over_DLk=Er["energy"] / (pair / Ev["energy"]))
        res.append(rec)
        print(json.dumps(rec), flush=True)
    print(f"# N={N} {name} done in {time.time()-t0:.1f}s", flush=True)
    return res


def cert_only(N, name, mu=100.0, k=4):
    """D_Lk without solving (cheap for large N)."""
    phi, gphi = PHIS[name]
    S = K.build(N, k); E = edge_info(S, phi)
    v = Lk_field(S, E, name, k)
    Ev = K.errors(S, v, ZERO, mu)
    pair = sum(e["w"] @ (e["eta"] * ((e["ph"] @ v[e["d"]]) @ e["n"])) for e in E)
    gam = mu * S.hGamma / S.h
    # Riemann-sum prediction: kappa_k sum_e |e|^{k+1} g_e^2
    from math import factorial
    kap = factorial(k) / factorial(2 * k + 1)   # <phi,L_k>_e = |e|^(k+1) k!/(2k+1)! avg_w(D_tau^k phi)
    pred = 0.0; gsq = 0.0
    for (t, u, w) in S.bedges:
        A = S.verts[S.tris[t, u]]; B = S.verts[S.tris[t, w]]; L = np.linalg.norm(B - A); tau = (B - A) / L; m = (A + B) / 2
        g = float(dtau_k(name, np.array(m[0]), np.array(m[1]), tau, k)); pred += kap * L**(k + 1) * g * g; gsq += L * g * g
    rec = dict(N=N, phi=name, mu=mu, gamma=gam, hG=S.hGamma, pair=pair, pair_pred=pred, vH1=Ev["H1"], vbL2=Ev["bL2"],
               venergy=Ev["energy"], D_Lk=pair / Ev["energy"], D_Lk_n=pair / Ev["energy"] / (S.hGamma**(k + 0.5) * (1 + gam)**-0.5),
               vnorm_n=Ev["energy"] / ((1 + gam)**0.5 * S.hGamma**-0.5 * np.sqrt(gsq)))
    print(json.dumps(rec), flush=True)


def certY(N, name, mu=100.0, k=4):
    """D_Y = <eta, v.n>/|||v||| for v = nodal (boundary-layer) lift of Y = Pi_T0(eta n)  (Thm R1(i)); no solve."""
    phi, gphi = PHIS[name]
    S = K.build(N, k); E = edge_info(S, phi)
    sysm = K.assemble(S, mu, gphi)
    Y = trace_projection(S, sysm, E)          # dof vector: nonzero only on Gamma_h nodes = nodal lift
    Ev = K.errors(S, Y, ZERO, mu)
    pair = sum(e["w"] @ (e["eta"] * ((e["ph"] @ Y[e["d"]]) @ e["n"])) for e in E)
    Yn, Yosc = gnorms(E, Y)
    eta_n = np.sqrt(sum(e["w"] @ e["eta"]**2 for e in E))
    gam = mu * S.hGamma / S.h
    rec = dict(N=N, phi=name, mu=mu, gamma=gam, hG=S.hGamma, eta=eta_n, Y=Yn, Yosc=Yosc, pair_over_Y2=pair / Yn**2,
               D_Y=pair / Ev["energy"], lift_const=Ev["energy"] / ((1 + gam)**0.5 * S.hGamma**-0.5 * Yn),
               D_Y_over_eps_eta=pair / Ev["energy"] / ((S.h / mu)**0.5 * eta_n))
    print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "certY":
        for N in sys.argv[2].split(","):
            certY(int(N), sys.argv[3], float(sys.argv[4]) if len(sys.argv) > 4 else 100.0)
        sys.exit()
    if sys.argv[1] == "cert":
        for N in sys.argv[2].split(","):
            cert_only(int(N), sys.argv[3])
        sys.exit()
    N = int(sys.argv[1]); name = sys.argv[2]; mus = [float(s) for s in sys.argv[3].split(",")]
    k = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    run(N, name, mus, k)
