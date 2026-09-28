"""
Pressure-robustness experiments for GS, CNS* and CNS*-R (boundary reconstruction, REPORT.md Sec. R3).

  nu-explicit system:  nu * [a_h - <dn u,v> - <u,dn v>] + P(nu) <u,v>  - B*(p,v) = (f, Tv)
     penalty "scaled":   P = nu*mu/h         (the paper's nu N_h, Sec. 11)
     penalty "unscaled": P = muhat/h, fixed  (= nu N_h with mu = muhat/nu)
  Tv = v (GS, CNS*)  or  Tv = Rv = v - sum_e b_e(v) (CNS*-R).

Hydrostatic test: u = 0, f = grad(phi) on Omega_h, g = 0.  Exact velocity is 0.

usage: python hydro.py hydro  N1,N2,..  [k]      -> hydrostatic table (phi = sin(2x)cos(y) and x^3-3xy^2)
       python hydro.py hydrolite N1,..  [k]      -> sin only; scaled penalty at nu=1 only, unscaled at all nu
       python hydro.py rates  N1,N2,..  [k]      -> CNS* vs CNS*-R on the paper's tests A, B1, C (nu = 1)
"""
import sys, os, time, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../code"))
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import svn_k as K


PHIS = {
    "sin":   (lambda X, Y: np.sin(2 * X) * np.cos(Y),
              lambda X, Y: np.stack([2 * np.cos(2 * X) * np.cos(Y), -np.sin(2 * X) * np.sin(Y)], -1)),
    "cubic": (lambda X, Y: X**3 - 3 * X * Y**2,
              lambda X, Y: np.stack([3 * X**2 - 3 * Y**2, -6 * X * Y], -1)),
}


def recon_correction(S, f):
    """Vector c with c_i = (f, b(phi_i)) summed over boundary triangles, where for v in V_h and the
    boundary triangle T with Gamma_h edge e, b = b_e(v) in P_k(T)^2 is the L2(T)-minimal field with
      b.n = v.n on e,  b.n = 0 on the two other edges   (moments against P_k on each edge),
      (div b, r)_T = <v.n, r>_e  for all r in P_{k-1}(T).
    The modified load is F - c, i.e. (f, v - b(v)) = (f, Rv)."""
    R = S.R; k = R.k; nl = R.nl
    ndof = 2 * S.nn
    c = np.zeros(ndof)
    # monomial exponents for P_k and P_{k-1}
    Mk = [(a, b) for a in range(k + 1) for b in range(k + 1 - a)]
    Mk1 = [(a, b) for a in range(k) for b in range(k - a)]
    gl, glw = np.polynomial.legendre.leggauss(k + 6)
    gl = (gl + 1) / 2; glw = glw / 2
    for (t, ref, Xe, we, n) in S.edge_data():
        V = S.verts[S.tris[t]]
        x0 = V.mean(0); hT = S.hdiam[t]
        loc = lambda X: (X - x0) / hT
        mono = lambda X, M: np.stack([loc(X)[..., 0]**a * loc(X)[..., 1]**b for (a, b) in M], -1)

        def dmono(X, M):
            Z = loc(X)
            dx = np.stack([a * Z[..., 0]**max(a - 1, 0) * Z[..., 1]**b if a else 0 * Z[..., 0] for (a, b) in M], -1) / hT
            dy = np.stack([b * Z[..., 0]**a * Z[..., 1]**max(b - 1, 0) if b else 0 * Z[..., 0] for (a, b) in M], -1) / hT
            return dx, dy
        # b = sum_j (c1_j m_j, c2_j m_j): coefficient vector of length 2*nk
        nk = len(Mk)
        Xq = S.phys(R.QX, R.QY)[t]; wq = R.QW * abs(S.detJ[t])
        mq = mono(Xq, Mk)
        Gm = np.einsum('q,qi,qj->ij', wq, mq, mq)
        Gram = np.block([[Gm, 0 * Gm], [0 * Gm, Gm]])
        rows = []; rhs_rows = []   # constraint rows acting on b-coeffs; rhs as linear map of local v dofs (2*nl)
        # identify local edges
        lv = [(0, 1), (1, 2), (2, 0)]
        A_e = S.verts[S.tris[t]]
        # boundary edge endpoints in edge_data order: Xe runs from vertex u to w
        for (i0, i1) in lv:
            P0, P1 = A_e[i0], A_e[i1]
            L = np.linalg.norm(P1 - P0); tv = (P1 - P0) / L
            nn = np.array([tv[1], -tv[0]])
            if np.dot(nn, V.mean(0) - P0) > 0:
                nn = -nn                                  # outward from T
            Xg = P0[None, :] + gl[:, None] * (P1 - P0)[None, :]
            wg = glw * L
            q = np.polynomial.legendre.legvander(2 * gl - 1, k)   # P_k on edge
            mg = mono(Xg, Mk)
            Cn = np.einsum('g,gq,gj->qj', wg, q, mg)
            rows.append(np.hstack([Cn * nn[0], Cn * nn[1]]))
            is_bd = np.allclose(sorted([tuple(P0), tuple(P1)]), sorted([tuple(Xe[0] - (Xe[1]-Xe[0])*0), tuple(Xe[-1])]), atol=0) if False else None
            # decide if this local edge is the Gamma_h edge: both endpoints on the unit circle
            onG = abs(np.linalg.norm(P0) - 1) < 1e-12 and abs(np.linalg.norm(P1) - 1) < 1e-12
            if onG:
                # map local velocity dofs -> <v.n_T, q>_e ; v.n_T with n_T outward from T (= -n of Omega_h? no:
                # n (edge_data) points away from T too, into the polygon), consistent below.
                refg = np.linalg.solve(S.J[t], (Xg - S.P0[t]).T).T
                ph = R.basis(refg[:, 0], refg[:, 1])            # (g, nl)
                Dv = np.einsum('g,gq,ga->qa', wg, q, ph)
                rhs_rows.append(np.hstack([Dv * nn[0], Dv * nn[1]]))
                ge = (Xg, wg, refg, nn)
            else:
                rhs_rows.append(np.zeros((k + 1, 2 * nl)))
        # divergence moments
        dx, dy = dmono(Xq, Mk)
        r = mono(Xq, Mk1)
        Cd = np.hstack([np.einsum('q,qr,qj->rj', wq, r, dx), np.einsum('q,qr,qj->rj', wq, r, dy)])
        rows.append(Cd)
        Xg, wg, refg, nn = ge
        ph = R.basis(refg[:, 0], refg[:, 1]); rg = mono(Xg, Mk1)
        Dd = np.einsum('g,gr,ga->ra', wg, rg, ph)
        rhs_rows.append(np.hstack([Dd * nn[0], Dd * nn[1]]))
        C = np.vstack(rows); Dm = np.vstack(rhs_rows)
        # min-L2 solution: b = Gram^{-1} C^T (C Gram^{-1} C^T)^+ D v
        Gi = np.linalg.inv(Gram)
        Sm = Gi @ C.T @ np.linalg.pinv(C @ Gi @ C.T, rcond=1e-12) @ Dm        # (2nk, 2nl)
        res = np.linalg.norm(C @ Sm - Dm) / max(1.0, np.linalg.norm(Dm))
        assert res < 1e-8, res
        fq = f(Xq[:, 0], Xq[:, 1])
        gb = np.hstack([wq @ (fq[:, 0:1] * mq), wq @ (fq[:, 1:2] * mq)])   # (f, basis of b)
        cl = gb @ Sm                                                          # (2nl,) in [x-dofs, y-dofs] order
        d = K.dofs(S.ids[t])                                                  # (nl, 2)
        np.add.at(c, d[:, 0], cl[:nl]); np.add.at(c, d[:, 1], cl[nl:])
    return c


def solve(S, A, Bt, B, F, g=None):
    ndof = A.shape[0]
    dD = K.dofs(S.dnodes).ravel(); free = np.setdiff1d(np.arange(ndof), dD)
    ug = np.zeros(ndof)
    if g is not None:
        ug[dD] = g(S.xy[S.dnodes, 0], S.xy[S.dnodes, 1]).ravel()
    Kmat = sp.bmat([[A[free][:, free], -Bt[:, free].T], [B[:, free], None]], format="csc")
    rhs = np.concatenate([F[free] - A[free][:, dD] @ ug[dD], -B[:, dD] @ ug[dD]])
    solver = os.environ.get("HYDRO_SOLVER", "superlu")
    if solver in ("pardiso", "pardiso_dry"):
        # as svn_k.solve_pardiso: pressure block regularised by 1e-12 * (pressure mass); this sets a
        # floor ~1e-11 (relative) on quantities that are exactly zero (CNS*-R and cubic CNS* on the hydrostatic test)
        nf = len(free)
        Kreg = (Kmat + sp.block_diag((sp.csr_matrix((nf, nf)), 1e-12 * PRESSURE_MASS[0]), format="csc")).tocsr()
        Kreg.sort_indices()
        if solver == "pardiso":
            import pypardiso
            slv = pypardiso.PyPardisoSolver(mtype=11); slv.set_iparm(11, 1); slv.set_iparm(13, 1)
            slv.factorize(Kreg)
            sol = lambda r: slv.solve(Kreg, r)
        else:                                   # dry run of the regularised system with SuperLU (testing only)
            lu = spl.splu(Kreg.tocsc(), permc_spec="COLAMD", diag_pivot_thresh=1.0)
            sol = lu.solve
        x = sol(rhs)
        for _ in range(10):
            x = x + sol(rhs - Kmat @ x)
        if solver == "pardiso":
            slv.free_memory(everything=True)
    else:
        lu = spl.splu(Kmat, permc_spec="COLAMD", diag_pivot_thresh=1.0)
        x = lu.solve(rhs)
        for _ in range(3):
            x = x + lu.solve(rhs - Kmat @ x)
    rel = np.linalg.norm(Kmat @ x - rhs) / max(np.linalg.norm(rhs), 1e-300)
    u = ug.copy(); u[free] = x[:len(free)]
    return u, x[len(free):], rel


PRESSURE_MASS = []      # filled by hydro()/rates() for the PARDISO path


ZERO = dict(u=lambda X, Y: np.zeros(np.shape(X) + (2,)), gu=lambda X, Y: np.zeros(np.shape(X) + (2, 2)))


def hydro(Ns, k=4, mu=100.0, muhat=100.0, nus=(1.0, 1e-2, 1e-4, 1e-6), phis=("sin", "cubic"),
          pens=("scaled", "unscaled"), scaled_nus=None):
    out = []
    for N in Ns:
        t0 = time.time()
        S = K.build(N, k)
        for phiname in phis:
            phi, gphi = PHIS[phiname]
            sysm = K.assemble(S, mu, gphi)
            PRESSURE_MASS[:] = [sp.block_diag(list(sysm["Ml"]), format="csr")]
            Avol, Pen, A = sysm["Avol"], sysm["Pen"], sysm["A"]
            Ncons = A - mu * Pen                       # a_h - <dn u,v> - <u,dn v>
            F0 = sysm["F"]
            cR = recon_correction(S, gphi)
            B = sysm["B"]
            for meth, theta, F in [("GS", 0, F0), ("CNS*", 1, F0), ("CNS*-R", 1, F0 - cR)]:
                Bt = K.coupling(S, sysm, theta)
                for pen in pens:
                    for nu in (scaled_nus or nus) if pen == "scaled" else nus:
                        if pen == "scaled":
                            Anu = nu * A; mueff = mu
                        else:
                            Anu = nu * Ncons + muhat * Pen; mueff = muhat / nu
                        u, p, rel = solve(S, Anu.tocsr(), Bt, B, F)
                        E = K.errors(S, u, ZERO, mueff)
                        Bu = B @ u
                        div = np.sqrt(abs(Bu @ (sysm["Minv"] @ Bu)))
                        rec = dict(N=N, k=k, phi=phiname, method=meth, pen=pen, nu=nu, mu=mueff,
                                   h=S.h, hG=S.hGamma, H1=E["H1"], bL2=E["bL2"], energy=E["energy"],
                                   relres=rel, div=div)
                        out.append(rec)
                        print(json.dumps(rec), flush=True)
        print(f"# N={N} done in {time.time()-t0:.1f}s", flush=True)
    return out


def rates(Ns, k=4, mu=100.0, kinds=("A", "B1", "C")):
    for N in Ns:
        S = K.build(N, k)
        for kind in kinds:
            ex = K.make_exact(kind)
            sysm = K.assemble(S, mu, ex["f"])
            PRESSURE_MASS[:] = [sp.block_diag(list(sysm["Ml"]), format="csr")]
            cR = recon_correction(S, ex["f"])
            for meth, F in [("CNS*", sysm["F"]), ("CNS*-R", sysm["F"] - cR)]:
                Bt = K.coupling(S, sysm, 1)
                u, p, rel = solve(S, sysm["A"], Bt, sysm["B"], F, g=ex["u"])
                E = K.errors(S, u, ex, mu)
                rec = dict(N=N, k=k, kind=kind, method=meth, mu=mu, h=S.h, hG=S.hGamma,
                           H1=E["H1"], energy=E["energy"], bL2=E["bL2"], relres=rel)
                print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    mode = sys.argv[1]; Ns = [int(s) for s in sys.argv[2].split(",")]
    k = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    if mode == "hydro":
        hydro(Ns, k)
    elif mode == "hydrolite":      # scaled penalty only at nu = 1 (nu-dependence is exact 1/nu scaling), unscaled at all nu
        hydro(Ns, k, phis=("sin",), scaled_nus=(1.0,))
    else:
        rates(Ns, k)
