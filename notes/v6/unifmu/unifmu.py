"""
notes/v6/unifmu: the leak limit (mu/h) delta_R -> U uniformly in mu (REPORT.tex).

Test P (u = 0, f = grad p), so u~ - u_h = delta_R exactly and X := (mu/h) delta_R = -(mu/h) u_h solves
    N_h(X, v) = lam <g, v>   for all v in Z_h,   lam = mu/h,   g = -(p~ - pbar) n_e  on Gamma_h.
Key identity (REPORT, Lemma 1): with T_h = traces of Z_h = {continuous pw P_k on Gamma_h, zero flux}
and t* = Pi_T g (L2(Gamma_h)-projection), N_h(X, v) = lam <t*, v> on Z_h.

For each N this script computes
  X(mu)    GS solution, mu in MUS                                  d_h = |grad(X - U~)| / |grad U|
  X^D      discrete Stokes, strong Dirichlet data t* on Gamma_h     (REPORT: the comparison field)
  X^U      discrete Stokes, strong Dirichlet data Pi_T U~           (no sawtooth: smooth reference)
  trace quantities: |t*-g|, |X - t*|, |zeta|, and the sawtooth capture fraction
        phi = -<X - U~, q_h> / <zeta, q_h>   (-> 1 in the Dirichlet limit; REPORT Thm 3)
  lower-bound certificate  <zeta, q_h> / ||q_h||_{H^-1/2}  <=  |Pi_T zeta|_{H^1/2}   (Fourier, arclength)
  |X - U~|_{H^1/2(Gamma_h)}  (<= C_T |grad(X - U~)|)
Usage: python unifmu.py N [N ...]   (ONLYGEO=1: trace/certificate quantities only, no solves;
       MUS=1e2,1e8: penalty list)     (1 core; N=96 ~ 6-8 min total, each solve < 2 min, < 3 GB)
"""
import os, sys, time, json
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
import svn, h1_limit, strong_bc

MUS = [float(m) for m in os.environ.get("MUS", "1e2,1e3,1e4,1e8").split(",")]
NORMU = h1_limit.G * 2.7565157            # ||grad U||_{L2(Omega)}
NQ = 40                                   # Gauss points per edge for trace quantities
GQ, GQW = np.polynomial.legendre.leggauss(NQ); GQ = (GQ + 1) / 2; GQW = GQW / 2
R3 = np.array([[0, 0], [1, 0], [0, 1]], float)


def pfun(X, Y):                           # test P pressure; its mean over Gamma is 0
    return X ** 2 * Y - Y + X / 2


def edges_ccw(S):
    """Gamma_h edges ordered counterclockwise, each oriented A -> B counterclockwise."""
    out = []
    for (t, u, w) in S.bedges:
        A = S.verts[S.tris[t, u]]; B = S.verts[S.tris[t, w]]
        if np.cross(A, B) < 0:
            u, w = w, u; A, B = B, A
        L = np.linalg.norm(B - A); tau = (B - A) / L
        n = np.array([tau[1], -tau[0]])              # right normal of a ccw polygon = outward of P_h
        n = -n                                       # outward of Omega_h = into P_h
        mid = (A + B) / 2
        out.append((np.arctan2(mid[1], mid[0]), t, u, w, A, B, L, tau, n))
    out.sort(key=lambda r: r[0])
    return out


class Trace:
    def __init__(self, S):
        self.S = S
        self.E = edges_ccw(S)
        # sanity: n points into the polygon (towards the origin)
        for (_, t, u, w, A, B, L, tau, n) in self.E:
            assert np.dot(n, (A + B) / 2) < 0
        self.ph = []; self.X = []; self.w = []; self.sig = []; self.s = []
        acc = 0.0
        for (_, t, u, w, A, B, L, tau, n) in self.E:
            ref = R3[u][None, :] + GQ[:, None] * (R3[w] - R3[u])[None, :]
            self.ph.append(svn.basis(ref[:, 0], ref[:, 1]))
            self.X.append(A[None, :] + GQ[:, None] * (B - A)[None, :])
            self.w.append(GQW * L); self.sig.append(acc + GQ * L); self.s.append((GQ - 0.5) * L)
            acc += L
        self.Lam = acc
        self.gnodes = strong_bc.gamma_nodes(S)

    def ev(self, vec):
        """values of a global P4 vector field (dof vector) at the trace quadrature points: list of (NQ,2)"""
        return [ph @ vec[svn.dofs(self.S.ids[e[1]])] for ph, e in zip(self.ph, self.E)]

    def ip(self, F1, F2):
        return sum((w[:, None] * a * b).sum() for w, a, b in zip(self.w, F1, F2))

    def nrm(self, F):
        return np.sqrt(self.ip(F, F))

    def project(self, F):
        """L2(Gamma_h) projection of the edgewise field F onto T_h (continuous P4, zero flux): dof vector"""
        S = self.S; nd = 2 * S.nn
        rows, cols, vals = [], [], []; b = np.zeros(nd); c = np.zeros(nd)
        for ph, e, w, f in zip(self.ph, self.E, self.w, F):
            n = e[8]; d = svn.dofs(S.ids[e[1]])
            Mloc = np.einsum('g,ga,gb->ab', w, ph, ph)
            for k in range(2):
                rows.append(np.repeat(d[:, k], 15)); cols.append(np.tile(d[:, k], 15)); vals.append(Mloc.ravel())
                np.add.at(b, d[:, k], (w[:, None] * ph * f[:, k:k + 1]).sum(0))
                np.add.at(c, d[:, k], (w[:, None] * ph).sum(0) * n[k])
        M = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nd, nd)).tocsr()
        gd = svn.dofs(self.gnodes).ravel()
        Mg = M[gd][:, gd]; cg = c[gd]; bg = b[gd]
        K = sp.bmat([[Mg, sp.csr_matrix(cg[:, None])], [sp.csr_matrix(cg[None, :]), None]], format="csc")
        x = spl.spsolve(K, np.concatenate([bg, [0.0]]))
        out = np.zeros(nd); out[gd] = x[:-1]
        return out

    def fourier(self, F, nmax):
        """complex Fourier coefficients (per component) on the arclength circle of length Lam"""
        sig = np.concatenate(self.sig); w = np.concatenate(self.w); f = np.concatenate(F)
        n = np.arange(-nmax, nmax + 1)
        E = np.exp(-2j * np.pi * np.outer(n, sig) / self.Lam)
        return n, (E * w[None, :]) @ f / self.Lam          # (2nmax+1, 2)

    def hs(self, F, s, nmax):
        """homogeneous H^s norm (mean removed): Lam * sum_{n!=0} |2 pi n/Lam|^{2s} |F_n|^2"""
        n, c = self.fourier(F, nmax)
        k = np.abs(2 * np.pi * n / self.Lam); m = n != 0
        return np.sqrt(self.Lam * np.sum(k[m, None] ** (2 * s) * np.abs(c[m]) ** 2))


def dirichlet(S, tvec):
    """discrete Stokes (a_h = 1/2(D.,D.), SV pair), strong data tvec on Gamma_h nodes, 0 on dR, f = 0"""
    zero = lambda X, Y: np.zeros(np.shape(X) + (2,))
    sysm = svn.assemble(S, 1.0, zero, 0, want=("A", "B"))
    A, B = sysm["Avol"], sysm["B"]
    ndof = A.shape[0]; npr = B.shape[0]
    gn = strong_bc.gamma_nodes(S)
    dG = svn.dofs(gn).ravel(); dS = svn.dofs(S.dnodes).ravel(); dD = np.union1d(dS, dG)
    ug = np.zeros(ndof); ug[dG] = tvec[dG]
    free = np.setdiff1d(np.arange(ndof), dD)
    keep = np.arange(1, npr); Bk = B[keep]
    K = sp.bmat([[A[free][:, free], -Bk[:, free].T], [Bk[:, free], None]], format="csc")
    rhs = np.concatenate([-A[free][:, dD] @ ug[dD], -Bk[:, dD] @ ug[dD]])
    lu = spl.splu(K, permc_spec="COLAMD", diag_pivot_thresh=1.0)
    x = lu.solve(rhs)
    for _ in range(3):
        x = x + lu.solve(rhs - K @ x)
    del lu
    u = ug.copy(); u[free] = x[:len(free)]
    Bu = B @ u
    return u, float(np.sqrt(abs(Bu @ (sysm["Minv"] @ Bu)))), float(np.sum(B[0::10] @ u))


def run(N):
    t0 = time.time()
    exP = svn.make_exact("P")
    Uf = h1_limit.LimitField()
    exU = dict(u=Uf.u, gu=Uf.gu)
    S = svn.Space(*svn.make_mesh(N))
    T = Trace(S)
    lam_over_mu = 1.0 / S.h
    # edgewise data
    g = [-(pfun(X[:, 0], X[:, 1]))[:, None] * e[8][None, :] for X, e in zip(T.X, T.E)]
    Ut = [Uf.u(X[:, 0], X[:, 1]) for X in T.X]
    zeta = [gg - uu for gg, uu in zip(g, Ut)]                   # zeta = g - U~  (= -paper's zeta)
    # sawtooth test field q_h = c_e beta(s) tau_e in T_h, c_e = LS slope of zeta.tau_e in s
    q = []
    for z, s, e in zip(zeta, T.s, T.E):
        L = e[6]; tau = e[7]
        ce = (GQW * L * s * (z @ tau)).sum() / (L ** 3 / 12)
        beta = s * (L ** 2 / 4 - s ** 2) / L ** 2
        q.append(ce * beta[:, None] * tau[None, :])
    nmax = 16 * N
    tstar = T.project(g); tU = T.project(Ut); tz = T.project(zeta)
    Tt = T.ev(tstar); TtU = T.ev(tU); Ttz = T.ev(tz)
    rec = dict(N=N, h=S.h, hG=S.hGamma, Lam=T.Lam,
               zeta=T.nrm(zeta), tstar_minus_g=T.nrm([a - b for a, b in zip(Tt, g)]),
               PiTzeta=T.nrm(Ttz), PiTU_minus_U=T.nrm([a - b for a, b in zip(TtU, Ut)]),
               zeta_q=T.ip(zeta, q), q_L2=T.nrm(q), q_Hm12=T.hs(q, -0.5, nmax),
               PiTzeta_H12=T.hs(Ttz, 0.5, nmax), zeta_q_on_PiTzeta=T.ip(Ttz, q))
    rec["cert_H12"] = rec["zeta_q"] / rec["q_Hm12"]           # <= |Pi_T zeta|_{H^1/2}
    print(json.dumps({k: float(f"{v:.6g}") if isinstance(v, float) else v for k, v in rec.items()}), flush=True)
    if os.environ.get("ONLYGEO"):
        return
    # discrete Dirichlet comparison fields
    for name, tv in (("XD", tstar), ("XU", tU)):
        t1 = time.time()
        u, dres, flux = dirichlet(S, tv)
        E = svn.errors(S, u, exU, 1.0)
        r = dict(N=N, field=name, d=E["H1"] / NORMU, divres=dres, flux=flux, secs=time.time() - t1)
        if name == "XD":
            XD = u
        else:
            XU = u
        print(json.dumps({k: float(f"{v:.6g}") if isinstance(v, float) else v for k, v in r.items()}), flush=True)
    Esaw = svn.errors(S, XD - XU, dict(u=lambda X, Y: np.zeros(np.shape(X) + (2,)),
                                          gu=lambda X, Y: np.zeros(np.shape(X) + (2, 2))), 1.0)
    print(json.dumps(dict(N=N, field="XD-XU", d=float(f"{Esaw['H1'] / NORMU:.6g}"))), flush=True)
    for mu in MUS:
        t1 = time.time()
        sysm = svn.assemble(S, mu, exP["f"], 0)
        uh, ph, dres, relres = svn.solve_lu(S, sysm, exP["u"], 0)
        X = -(mu / S.h) * uh
        E = svn.errors(S, X, exU, mu)
        ED = svn.errors(S, X - XD, dict(u=lambda X_, Y_: np.zeros(np.shape(X_) + (2,)),
                                        gu=lambda X_, Y_: np.zeros(np.shape(X_) + (2, 2))), mu)
        TX = T.ev(X)
        W = [a - b for a, b in zip(TX, Ut)]
        r = dict(N=N, field="GS", mu=mu, gamma=mu * S.hGamma / S.h, gammahG=mu * S.hGamma ** 2 / S.h,
                 d=E["H1"] / NORMU, d_XD=ED["H1"] / NORMU,
                 X_minus_tstar=T.nrm([a - b for a, b in zip(TX, Tt)]), X_minus_g=T.nrm([a - b for a, b in zip(TX, g)]),
                 phi=-T.ip(W, q) / T.ip(zeta, q) * (-1),     # zeta here = g - U~, so capture = <W,q>/<g-U~,q>
                 W_H12=T.hs(W, 0.5, nmax), W_L2=T.nrm(W), relres=relres, divres=float(dres[0]),
                 secs=time.time() - t1)
        print(json.dumps({k: float(f"{v:.6g}") if isinstance(v, float) else v for k, v in r.items()}), flush=True)
    print(f"# N={N} total {time.time() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    for N in (int(a) for a in sys.argv[1:]):
        run(N)
