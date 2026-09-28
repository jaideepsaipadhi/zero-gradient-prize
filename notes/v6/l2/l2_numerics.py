"""
L2 velocity errors for GS and CNS* (notes/v6/l2/REPORT.md), k = 4 unless given, standard meshes (svn.make_mesh).

Modes
  slip  N [k] [mu]  : CNS* (theta=1) on tests A and B1 with THREE choices of Nitsche wall data on Gamma_h
                      (one LU factorisation, several right-hand sides):
                        zero  : the method as defined (u_h = 0 weakly on Gamma_h)
                        exact : g = u~|Gamma_h          (removes the whole trace inconsistency of u~)
                        mean  : g = (|e|^2/12) d_n u(x*) on each chord e (only the chord-MEAN slip removed)
                      Prints ||u~ - u_h||_{L2(Omega_h)} and its ratio to hG^2 for each; the difference zero-mean
                      is the discrete response V_h^FE to the mean-slip data (prediction: e_u = V_h + o(hG^2)).
                      Also prints the GS L2 error for A (G = 0).
  leak  N [mu]      : GS (theta=0) on test P: ||u~-u_h||_{L2}, its ratio to (h/mu)||U||_{L2(Omega)}, and the
                      direct distance ||(mu/h)(u~-u_h) - U~||_{L2(Omega_h)} / ||U||_{L2(Omega)},
                      U = leak flow of Definition hc:def (series solution from code/h1_limit.py).
  unorm             : ||U||_{L2(Omega)} for test P by polar Gauss quadrature (and ||grad U|| as a check).

Weak data g on Gamma_h enter N_h as -<u-g, d_n v> + (mu/h)<u-g, v>, i.e. RHS += -<g, d_n v> + (mu/h)<g, v>.
Single-threaded; N = 96 takes ~1-2 min and < 2 GB.  Usage: python3 l2_numerics.py slip 64
"""
import os, sys, time, json
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../../code"))
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import svn_k


def nitsche_data_rhs(S, mu, gfun):
    """vector of -<g, d_n v> + (mu/h)<g, v> over Gamma_h; gfun(X, n) -> (nq, 2) data at edge points X."""
    R = S.R; ndof = 2 * S.nn
    b = np.zeros(ndof)
    for (t, ref, Xe, we, n) in S.edge_data():
        ph = R.basis(ref[:, 0], ref[:, 1]); dX, dY = R.dbasis(ref[:, 0], ref[:, 1])
        Ji = S.Jinv[t]
        gxe = Ji[0, 0] * dX + Ji[1, 0] * dY; gye = Ji[0, 1] * dX + Ji[1, 1] * dY
        dn = n[0] * gxe + n[1] * gye
        g = gfun(Xe, n)
        d = svn_k.dofs(S.ids[t])
        for c in range(2):
            loc = -(we * g[:, c]) @ dn + (mu / S.h) * (we * g[:, c]) @ ph
            np.add.at(b, d[:, c], loc)
    return b


class Saddle:
    def __init__(self, S, sysm, theta):
        A, B = sysm["A"], sysm["B"]
        Bt = svn_k.coupling(S, sysm, theta)
        ndof = A.shape[0]
        self.dD = svn_k.dofs(S.dnodes).ravel(); self.free = np.setdiff1d(np.arange(ndof), self.dD)
        f = self.free
        self.K = sp.bmat([[A[f][:, f], -Bt[:, f].T], [B[:, f], None]], format="csc")
        self.lu = spl.splu(self.K, permc_spec="COLAMD", diag_pivot_thresh=1.0)
        self.A, self.B, self.ndof = A, B, ndof

    def solve(self, F, S, g):
        ug = np.zeros(self.ndof); ug[self.dD] = g(S.xy[S.dnodes, 0], S.xy[S.dnodes, 1]).ravel()
        f = self.free
        rhs = np.concatenate([F[f] - self.A[f][:, self.dD] @ ug[self.dD], -self.B[:, self.dD] @ ug[self.dD]])
        x = self.lu.solve(rhs)
        for _ in range(3):
            x = x + self.lu.solve(rhs - self.K @ x)
        u = ug.copy(); u[f] = x[:len(f)]
        return u, np.linalg.norm(self.K @ x - rhs) / max(np.linalg.norm(rhs), 1e-300)


def dnu_at_proj(ex, X):
    """d_n u(x*) with x* = X/|X| and n(x*) = -x* (outward normal of Omega, pointing into the disk)."""
    r = np.hypot(X[:, 0], X[:, 1]); xs = X / r[:, None]
    Gm = ex["gu"](xs[:, 0], xs[:, 1])            # [q, c, d] = d_d u_c
    return np.einsum('qcd,qd->qc', Gm, -xs)


def slip(N, k=4, mu=100.0):
    out = []
    for kind in ("A", "B1"):
        ex = svn_k.make_exact(kind)
        t0 = time.time()
        S = svn_k.build(N, k)
        sysm = svn_k.assemble(S, mu, ex["f"])
        rows = {}
        for theta in (1, 0):
            sad = Saddle(S, sysm, theta)
            lens = {}
            data = {"zero": None,
                    "exact": lambda X, n: ex["u"](X[:, 0], X[:, 1]),
                    "mean": None}
            for name in (("zero", "exact", "mean") if theta == 1 else ("zero",)):
                F = sysm["F"].copy()
                if name == "exact":
                    F = F + nitsche_data_rhs(S, mu, data["exact"])
                elif name == "mean":
                    # g = (|e|^2/12) d_n u(x*) per chord; |e| = sum of the edge quadrature weights
                    b = np.zeros(2 * S.nn); R = S.R
                    for (t, ref, Xe, we, n) in S.edge_data():
                        ph = R.basis(ref[:, 0], ref[:, 1]); dX, dY = R.dbasis(ref[:, 0], ref[:, 1])
                        Ji = S.Jinv[t]
                        dn = n[0] * (Ji[0, 0] * dX + Ji[1, 0] * dY) + n[1] * (Ji[0, 1] * dX + Ji[1, 1] * dY)
                        g = (we.sum() ** 2 / 12.0) * dnu_at_proj(ex, Xe)
                        d = svn_k.dofs(S.ids[t])
                        for c in range(2):
                            np.add.at(b, d[:, c], -(we * g[:, c]) @ dn + (mu / S.h) * (we * g[:, c]) @ ph)
                    F = F + b
                u, rr = sad.solve(F, S, ex["u"])
                E = svn_k.errors(S, u, ex, mu)
                rows[(theta, name)] = u
                rec = dict(N=N, k=k, mu=mu, kind=kind, method="CNS*" if theta else "GS", data=name,
                           L2=E["L2"], H1=E["H1"], bL2=E["bL2"], h=S.h, hG=S.hGamma, relres=rr)
                out.append(rec)
                print(f"{kind:3s} N={N:4d} {rec['method']:5s} data={name:5s}  L2={E['L2']:.4e}  L2/hG^2={E['L2']/S.hGamma**2:.5f}"
                      f"  H1={E['H1']:.4e}  slip={E['bL2']:.3e}  res={rr:.1e}", flush=True)
            del sad
        # the mean-slip flow V_h^FE := u_h(mean) - u_h(zero)  (discrete response to data g_m, zero forcing)
        V = rows[(1, "mean")] - rows[(1, "zero")]
        zero = dict(u=lambda X, Y: np.zeros(np.shape(X) + (2,)), gu=lambda X, Y: np.zeros(np.shape(X) + (2, 2)))
        EV = svn_k.errors(S, V, zero, mu)
        print(f"{kind:3s} N={N:4d} ||V_h^FE||_L2={EV['L2']:.4e}  /hG^2={EV['L2']/S.hGamma**2:.5f}   "
              f"({time.time()-t0:.0f}s)", flush=True)
        out.append(dict(N=N, k=k, mu=mu, kind=kind, method="V_h^FE", L2=EV["L2"], H1=EV["H1"], h=S.h, hG=S.hGamma))
    return out


def unorm(M=63, q=96):
    import h1_limit as HL
    Uf = HL.LimitField(M)
    g, w = np.polynomial.legendre.leggauss(q)
    L = HL.L; l2 = 0.0; h1 = 0.0
    for kk in range(8):
        a0 = -np.pi / 4 + kk * np.pi / 4; a1 = a0 + np.pi / 4
        th = a0 + (g + 1) / 2 * (a1 - a0); wt = w * (a1 - a0) / 2
        Rt = L / np.maximum(np.abs(np.cos(th)), np.abs(np.sin(th)))
        for tj, wj, Rj in zip(th, wt, Rt):
            rr = 1 + (g + 1) / 2 * (Rj - 1); ww = w * (Rj - 1) / 2
            X, Y = rr * np.cos(tj), rr * np.sin(tj)
            U = Uf.u(X, Y); GU = Uf.gu(X, Y)
            l2 += wj * np.sum(ww * rr * (U ** 2).sum(-1))
            h1 += wj * np.sum(ww * rr * (GU ** 2).sum((-1, -2)))
    print(f"||U||_L2(Omega) = {np.sqrt(l2):.8f}   ||grad U|| = {np.sqrt(h1):.8f}   G = {HL.G:.8f}   "
          f"||U||/G = {np.sqrt(l2)/HL.G:.8f}", flush=True)
    return np.sqrt(l2), Uf


def leak(N, mu=100.0, k=4):
    Unrm, Uf = unorm()
    ex = svn_k.make_exact("P1")
    S = svn_k.build(N, k)
    sysm = svn_k.assemble(S, mu, ex["f"])
    t0 = time.time()
    u, p, dres, relres, _ = svn_k.solve_lu(S, sysm, ex["u"], 0)
    E = svn_k.errors(S, u, ex, mu)
    s = S.h / mu
    exU = dict(u=lambda X, Y: -s * Uf.u(X, Y), gu=lambda X, Y: -s * Uf.gu(X, Y))   # e_u = -u_h ~ (h/mu) U
    D = svn_k.errors(S, u, exU, mu)
    rec = dict(N=N, mu=mu, k=k, L2=E["L2"], H1=E["H1"], ratio=E["L2"] / (s * Unrm),
               distL2=D["L2"] / (s * Unrm), distH1rel=D["H1"] / (s * 4.57024450), h=S.h, hG=S.hGamma, relres=relres)
    print(f"P  N={N:4d} mu={mu:g} GS  L2={E['L2']:.4e}  L2/((h/mu)||U||)={rec['ratio']:.5f}  "
          f"||(mu/h)e-U||_L2/||U||={rec['distL2']:.4f}  (H1 dist {rec['distH1rel']:.4f})  ({time.time()-t0:.0f}s)", flush=True)
    return [rec]


def moments(N, k=4, mu=100.0):
    """Chord means of e_u.n_e and e_u.tau_e (and of u~.tau_e) for CNS*, data zero and mean; Fourier
    coefficients (m <= 6) of the chord means against the midpoint angle, normalised by hG^2."""
    for kind in ("A", "B1"):
        ex = svn_k.make_exact(kind)
        S = svn_k.build(N, k)
        sysm = svn_k.assemble(S, mu, ex["f"])
        sad = Saddle(S, sysm, 1)
        R = S.R
        ed = S.edge_data()
        b = np.zeros(2 * S.nn)
        for (t, ref, Xe, we, n) in ed:
            ph = R.basis(ref[:, 0], ref[:, 1]); dX, dY = R.dbasis(ref[:, 0], ref[:, 1])
            Ji = S.Jinv[t]
            dn = n[0] * (Ji[0, 0] * dX + Ji[1, 0] * dY) + n[1] * (Ji[0, 1] * dX + Ji[1, 1] * dY)
            g = (we.sum() ** 2 / 12.0) * dnu_at_proj(ex, Xe)
            d = svn_k.dofs(S.ids[t])
            for c in range(2):
                np.add.at(b, d[:, c], -(we * g[:, c]) @ dn + (mu / S.h) * (we * g[:, c]) @ ph)
        for name, F in (("zero", sysm["F"]), ("mean", sysm["F"] + b)):
            u, _ = sad.solve(F.copy(), S, ex["u"])
            th, L, mn, mt, ut = [], [], [], [], []
            for (t, ref, Xe, we, n) in ed:
                ph = R.basis(ref[:, 0], ref[:, 1])
                uh = ph @ u[svn_k.dofs(S.ids[t])]
                eu = ex["u"](Xe[:, 0], Xe[:, 1]) - uh
                tau = np.array([-n[1], n[0]])
                X0 = Xe.mean(0); th.append(np.arctan2(X0[1], X0[0])); L.append(we.sum())
                if np.cross(X0, tau) < 0:      # orient tau counterclockwise (e_theta)
                    tau = -tau
                mn.append((we @ (eu @ n)) / we.sum()); mt.append((we @ (eu @ tau)) / we.sum())
                ut.append((we @ (ex["u"](Xe[:, 0], Xe[:, 1]) @ tau)) / we.sum())
            th, L, mn, mt, ut = map(np.array, (th, L, mn, mt, ut))
            hG2 = S.hGamma ** 2
            def four(v):
                return max(abs(np.sum(L * v * f(m * th))) / np.pi for m in range(0, 7) for f in (np.cos, np.sin))
            print(f"{kind:3s} N={N:4d} mu={mu:g} data={name:5s}  max_m|F_m(e.n mean)|/hG^2={four(mn)/hG2:.5f}  "
                  f"max_m|F_m(e.tau mean)|/hG^2={four(mt)/hG2:.5f}  max_m|F_m(u~.tau mean)|/hG^2={four(ut)/hG2:.5f}"
                  f"  L2(e.n means)/hG^2={np.sqrt(np.sum(L*mn**2))/hG2:.5f}  L2(e.tau means - u~.tau means)/hG^2="
                  f"{np.sqrt(np.sum(L*(mt-ut)**2))/hG2:.5f}", flush=True)

def pair(N, k=4, mu=100.0):
    """Test A, CNS* and GS: the weak identity of Theorem L3,  (e_u, phi~)_{Omega_h} = -Lambda_h + O(hG^3),
    with the pressure-free dual z = curl zeta, zeta = chi(r) (r-1)^2 / 2, phi~ = -Delta z (so W_z = 1), and
    -Lambda_h = -sum_e |e|^3/12 W_u W_z = + sum_e |e|^3/6  (W_u = -2 for test A)."""
    import sympy as sy
    x, y = sy.symbols('x y', real=True)
    r = sy.sqrt(x**2 + y**2)
    t = (r - 1.5) / 0.7                                  # chi = 1 for r <= 1.5, 0 for r >= 2.2 (C^infinity)
    bump = lambda s: sy.exp(-1 / s)
    chi = sy.Piecewise((1, r <= 1.5), (0, r >= 2.2),
                       (bump(1 - t) / (bump(1 - t) + bump(t)), True))
    zeta = chi * (r - 1)**2 / 2
    z = [sy.diff(zeta, y), -sy.diff(zeta, x)]
    phi = [-(sy.diff(zi, x, 2) + sy.diff(zi, y, 2)) for zi in z]
    PH = sy.lambdify((x, y), phi, 'numpy')
    ex = svn_k.make_exact("A")
    S = svn_k.build(N, k)
    sysm = svn_k.assemble(S, mu, ex["f"])
    R = S.R
    X = S.phys(R.QX, R.QY); wq = R.QW[None, :] * np.abs(S.detJ)[:, None]
    ph = PH(X[..., 0], X[..., 1])
    ph = np.stack([np.broadcast_to(np.nan_to_num(np.asarray(v, float)), X[..., 0].shape) for v in ph], -1)
    lam = sum(we.sum() ** 3 / 6.0 for (_, _, _, we, _) in S.edge_data())
    for theta in (1, 0):
        sad = Saddle(S, sysm, theta)
        u, _ = sad.solve(sysm["F"].copy(), S, ex["u"])
        uh = np.einsum('qa,tac->tqc', R.PHI_Q, u[svn_k.dofs(S.ids)])
        eu = ex["u"](X[..., 0], X[..., 1]) - uh
        val = (wq * (eu * ph).sum(-1)).sum()
        print(f"A N={N:4d} {'CNS*' if theta else 'GS  '} (e_u,phi)={val:.6e}  -Lambda_h={lam:.6e}  "
              f"ratio={val/lam:.5f}  (diff)/hG^3={(val-lam)/S.hGamma**3:.4f}", flush=True)
        del sad


if __name__ == "__main__":
    mode = sys.argv[1]
    res = []
    if mode == "slip":
        res = slip(int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 4,
                   float(sys.argv[4]) if len(sys.argv) > 4 else 100.0)
    elif mode == "leak":
        res = leak(int(sys.argv[2]), float(sys.argv[3]) if len(sys.argv) > 3 else 100.0)
    elif mode == "unorm":
        unorm()
    elif mode == "pair":
        pair(int(sys.argv[2]))
    elif mode == "moments":
        moments(int(sys.argv[2]), 4, float(sys.argv[3]) if len(sys.argv) > 3 else 100.0)
    if res:
        with open(os.path.join(HERE, "l2_results.jsonl"), "a") as fh:
            for r in res:
                fh.write(json.dumps(r) + "\n")
