"""
notes/v7/nsconst: the Navier--Stokes constants of Theorem n2:L for the paper's NS test (REPORT.tex, Section 2).

Test (code/job_v4.py ns, results/v4/report_v4.txt part B): Omega = (-2.5,2.5)^2 minus the closed unit disk,
u = curl psi, psi = (r^2-1)^2 (x + y^2/2)/10 (test B), nu in {1, 0.1}; f is manufactured, so u is the exact solution for
every nu.  Discretisation: isoparametric Taylor--Hood (iso_th.py) on the exact domain, production-mesh topology.

  python consts.py check N k          validation (manufactured Stokes solution; beta_O = nu for u = 0)
  python consts.py CF N k             Friedrichs constants: H^1_0(Omega), and zero only on the square (Neumann on Gamma)
  python consts.py beta N k nu        beta_O = inf-sup of the Oseen operator L on V (H^1_0 seminorm on both sides)
  python consts.py CR N k M           C_R lower-bound sequence: sup over g in P_M(box)^2 (fine solves) of Q(g)/||g||^2
  python consts.py Uinf               U_inf = max(||u||_inf, ||grad u||_inf) on the closed box (test B)
"""
import os, sys, json, time
for v_ in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v_, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, scipy.linalg as sla
import iso_th
import sympy as sy

L = 2.5


def exact_B():
    x, y = sy.symbols("x y", real=True)
    psi = (x ** 2 + y ** 2 - 1) ** 2 * (x + y ** 2 / 2) / 10
    u = [sy.diff(psi, y), -sy.diff(psi, x)]
    gu = [[sy.diff(u[c], v) for v in (x, y)] for c in range(2)]
    fu = sy.lambdify((x, y), u, "numpy"); fg = sy.lambdify((x, y), gu, "numpy")

    def U(X, Y):
        a = fu(X, Y); return np.stack([np.broadcast_to(np.asarray(v, float), np.shape(X)) for v in a], -1)

    def G(X, Y):
        a = fg(X, Y)
        return np.stack([np.stack([np.broadcast_to(np.asarray(a[c][d], float), np.shape(X)) for d in range(2)], -1)
                         for c in range(2)], -2)
    return U, G


class Saddle:
    """[[K, -B^T, 0], [B, 0, m], [0, m^T, 0]] on free velocity dofs; K may be nonsymmetric."""
    def __init__(self, I, K, B, m):
        self.I = I; f = I.vfree
        self.Kf = K[f][:, f].tocsc(); self.Bf = B[:, f].tocsc()
        nu_, np_ = self.Kf.shape[0], B.shape[0]
        mm = sp.csc_matrix(m[:, None])
        self.M = sp.bmat([[self.Kf, -self.Bf.T, None], [self.Bf, None, mm], [None, mm.T, None]], format="csc")
        self.lu = spl.splu(self.M, permc_spec="COLAMD")
        self.nu_, self.np_ = nu_, np_

    def solve(self, F, trans=False):
        """F: velocity load(s) on free dofs, (nf,) or (nf, r). returns velocity (nf,...) and pressure (np,...)"""
        F2 = F if F.ndim == 2 else F[:, None]
        rhs = np.zeros((self.M.shape[0], F2.shape[1])); rhs[:self.nu_] = F2
        x = self.lu.solve(rhs, trans="T" if trans else "N")
        r = rhs - (self.M.T @ x if trans else self.M @ x)
        x = x + self.lu.solve(r, trans="T" if trans else "N")
        u, p = x[:self.nu_], x[self.nu_:self.nu_ + self.np_]
        return (u[:, 0], p[:, 0]) if F.ndim == 1 else (u, p)


def build(N, k):
    I = iso_th.Iso(N, k)
    A = I.vec(I.scalar_stiff()); Mv = I.vec(I.scalar_mass()); B = I.div(); m = I.pmean()
    return I, A, Mv, B, m


# ------------------------------------------------------------------ validation
def cmd_check(N, k):
    I, A, Mv, B, m = build(N, k)
    x, y = sy.symbols("x y", real=True)
    # manufactured: phi = curl((r^2-1)^2 (x^2-L^2)^2 (y^2-L^2)^2 * sin(x+2y)/100), pi = x y^2 - mean
    s = (x ** 2 + y ** 2 - 1) ** 2 * (x ** 2 - L ** 2) ** 2 * (y ** 2 - L ** 2) ** 2 * sy.sin(x + 2 * y) / 100
    ph = [sy.diff(s, y), -sy.diff(s, x)]; pi = x * y ** 2 + sy.cos(y)
    g = [-sy.diff(ph[c], x, 2) - sy.diff(ph[c], y, 2) + sy.diff(pi, [x, y][c]) for c in range(2)]
    fg = sy.lambdify((x, y), g, "numpy"); fp = sy.lambdify((x, y), ph, "numpy")
    fgp = sy.lambdify((x, y), [[sy.diff(ph[c], v) for v in (x, y)] for c in range(2)], "numpy")
    G = lambda X, Y: np.stack([np.asarray(v, float) * np.ones_like(X) for v in fg(X, Y)], -1)
    F = I.load(G)[I.vfree]
    Sd = Saddle(I, A, B, m)
    u, p = Sd.solve(F)
    U = np.zeros(2 * I.nu); U[I.vfree] = u
    Uq = np.stack([np.einsum('tqa,ta->tq', I.gx, U[2 * I.S.ids + c]) for c in range(2)], -1)
    Vq = np.stack([np.einsum('tqa,ta->tq', I.gy, U[2 * I.S.ids + c]) for c in range(2)], -1)
    X = I.Xq; gp = fgp(X[..., 0], X[..., 1])
    ex = np.stack([np.asarray(gp[c][0], float) * np.ones_like(X[..., 0]) for c in range(2)], -1)
    ey = np.stack([np.asarray(gp[c][1], float) * np.ones_like(X[..., 0]) for c in range(2)], -1)
    err = np.sqrt((I.wq[..., None] * ((Uq - ex) ** 2 + (Vq - ey) ** 2)).sum())
    nrm = np.sqrt((I.wq[..., None] * (ex ** 2 + ey ** 2)).sum())
    # beta for u = 0 must be nu = 1
    b = beta_value(I, A, B, m, None, 1.0)
    print(json.dumps(dict(N=N, k=k, nfree=len(I.vfree), np=I.np_, H1_relerr=err / nrm, beta_u0=b)), flush=True)


# ------------------------------------------------------------------ Friedrichs
def cmd_CF(N, k):
    I = iso_th.Iso(N, k)
    K = I.scalar_stiff(); M = I.scalar_mass()
    out = dict(N=N, k=k)
    for name, fixed in (("CF0", I.dnodes), ("CF_R", I.S.dnodes)):
        f = np.setdiff1d(np.arange(I.nu), fixed)
        w = spl.eigsh(K[f][:, f].tocsc(), k=1, M=M[f][:, f].tocsc(), sigma=0, which="LM", return_eigenvectors=False)
        out[name] = float(1 / np.sqrt(w.min())); out["lam1_" + name] = float(w.min())
    print(json.dumps(out), flush=True)
    return out


# ------------------------------------------------------------------ Oseen inf-sup
def beta_value(I, A, B, m, conv, nu, nev=1):
    K = nu * A + (conv if conv is not None else 0 * A)
    SL = Saddle(I, K, B, m)
    Af = A[I.vfree][:, I.vfree].tocsr()

    def op(y):
        x, _ = SL.solve(Af @ y)
        z, _ = SL.solve(Af @ x, trans=True)
        return z
    n = len(I.vfree)
    Op = spl.LinearOperator((n, n), matvec=op, dtype=float)
    w = spl.eigs(Op, k=max(nev, 1), which="LM", return_eigenvectors=False, tol=1e-10, maxiter=5000)
    tau = np.max(np.abs(w))
    return float(1 / np.sqrt(tau))


def cmd_beta(N, k, nu):
    t0 = time.time()
    I, A, Mv, B, m = build(N, k)
    U, G = exact_B()
    C = I.oseen_conv(U, G)
    b = beta_value(I, A, B, m, C, nu)
    print(json.dumps(dict(N=N, k=k, nu=nu, beta_O=b, b_O=b / nu, ndof=len(I.vfree), secs=round(time.time() - t0, 1))),
          flush=True)


# ------------------------------------------------------------------ Stokes regularity constant (lower-bound sequence)
def leg_basis(M):
    from numpy.polynomial import legendre as Lg
    idx = [(i, j) for i in range(M + 1) for j in range(M + 1 - i)]

    def ev(X, Y):
        xs, ys = X / L, Y / L
        Px = np.stack([Lg.legval(xs, np.eye(M + 1)[i]) for i in range(M + 1)], -1)
        Py = np.stack([Lg.legval(ys, np.eye(M + 1)[j]) for j in range(M + 1)], -1)
        return np.stack([Px[..., i] * Py[..., j] for (i, j) in idx], -1)      # (..., nb)
    return idx, ev


def cmd_CR(N, k, Ms):
    t0 = time.time()
    I, A, Mv, B, m = build(N, k)
    Sd = Saddle(I, A, B, m)
    f = I.vfree
    ids = I.S.ids; pids = I.Sp.ids
    circ = I.circle_edges()
    area = I.wq.sum()
    for M in Ms:
        idx, ev = leg_basis(M)
        nb = len(idx)
        Pq = ev(I.Xq[..., 0], I.Xq[..., 1])                         # (nt,nq,nb)
        # vector basis: component c, scalar poly j  -> 2 nb functions
        loc = np.einsum('tq,qa,tqj->taj', I.wq, I.R.PHI_Q, Pq)       # (nt,nl,nb)
        F = np.zeros((2 * I.nu, 2 * nb))
        for c in range(2):
            np.add.at(F, (2 * ids + c).ravel(), np.concatenate([np.zeros((ids.size, c * nb)),
                      loc.reshape(-1, nb), np.zeros((ids.size, (1 - c) * nb))], 1))
        u, p = Sd.solve(F[f])
        Uf = np.zeros((2 * I.nu, 2 * nb)); Uf[f] = u
        # quadrature values
        Ux = [np.einsum('tqa,taj->tqj', I.gx, Uf[2 * ids + c]) for c in range(2)]
        Uy = [np.einsum('tqa,taj->tqj', I.gy, Uf[2 * ids + c]) for c in range(2)]
        Uv = [np.einsum('qa,taj->tqj', I.R.PHI_Q, Uf[2 * ids + c]) for c in range(2)]
        Pv = np.einsum('qa,taj->tqj', I.PP, p[pids]); Px = np.einsum('tqa,taj->tqj', I.pgx, p[pids])
        Py = np.einsum('tqa,taj->tqj', I.pgy, p[pids])
        pbar = np.einsum('tq,tqj->j', I.wq, Pv) / area
        Pv = Pv - pbar
        gq = [np.concatenate([Pq * (c == 0), Pq * (c == 1)], -1) for c in range(2)]   # g components (nt,nq,2nb)
        w = I.wq
        ip = lambda a, b: np.einsum('tq,tqj,tqk->jk', w, a, b, optimize=True)
        Gram = ip(gq[0], gq[0]) + ip(gq[1], gq[1])
        parts = dict(
            L2=ip(Uv[0], Uv[0]) + ip(Uv[1], Uv[1]),
            H1=sum(ip(a, a) for a in Ux + Uy),
            Lap=ip(Px - gq[0], Px - gq[0]) + ip(Py - gq[1], Py - gq[1]),     # ||Delta phi||^2 = ||grad pi - g||^2
            pL2=ip(Pv, Pv),
            pH1=ip(Px, Px) + ip(Py, Py))
        dn = np.zeros((2 * nb, 2 * nb))
        for (t, X, wg, n, gxe, gye) in circ:
            for c in range(2):
                vals = (n[:, 0:1] * gxe + n[:, 1:2] * gye) @ Uf[2 * ids[t] + c]    # (ng, 2nb)
                dn += np.einsum('g,gj,gk->jk', wg, vals, vals)
        parts["dnGamma"] = dn                                                   # Grisvard: ||D^2 phi||^2 = Lap + dnGamma
        Q = sum(parts.values())
        # generalized eigen, on a well-conditioned subspace of the Gram matrix
        gw, gv = np.linalg.eigh(Gram); keep = gw > 1e-11 * gw.max()
        T = gv[:, keep] / np.sqrt(gw[keep])
        Qs = T.T @ Q @ T; Qs = (Qs + Qs.T) / 2
        ev_, evec = np.linalg.eigh(Qs); top = ev_[-1]; x = T @ evec[:, -1]
        split = {kk: float(x @ v @ x) for kk, v in parts.items()}
        # sum form sup (||phi||_H2 + ||pi||_H1)/||g||: (a+b)^2 <= (1+s)a^2 + (1+1/s)b^2, minimise the sup over s
        Qphi = T.T @ (parts["L2"] + parts["H1"] + parts["Lap"] + parts["dnGamma"]) @ T
        Qpi = T.T @ (parts["pL2"] + parts["pH1"]) @ T
        lam_s = lambda s: np.linalg.eigvalsh((1 + s) * Qphi + (1 + 1 / s) * Qpi)[-1]
        sg = np.geomspace(1e-3, 1e3, 121); ls = [lam_s(s) for s in sg]
        print(json.dumps(dict(N=N, k=k, M=M, nb=2 * nb, rank=int(keep.sum()), CRquad=float(np.sqrt(top)),
                              CR_sum_bound=float(np.sqrt(min(ls))), s_opt=float(sg[int(np.argmin(ls))]),
                              CR_sqrt2_bound=float(np.sqrt(2 * top)), split_at_max=split,
                              secs=round(time.time() - t0, 1))), flush=True)


def cmd_Uinf():
    U, G = exact_B()
    s = np.linspace(-L, L, 2001); X, Y = np.meshgrid(s, s)
    msk = X ** 2 + Y ** 2 >= 1
    u = U(X, Y); g = G(X, Y)
    un = np.linalg.norm(u, axis=-1)[msk].max()
    gf = np.sqrt((g ** 2).sum((-1, -2)))[msk].max()
    g2 = np.linalg.norm(g, ord=2, axis=(-2, -1))[msk].max()
    # collar inside the disk, r in (1 - r0, 1], r0 = 0.1
    msk2 = (X ** 2 + Y ** 2 < 1) & (X ** 2 + Y ** 2 > 0.81)
    print(json.dumps(dict(u_inf=float(un), gradu_inf_frob=float(gf), gradu_inf_op=float(g2),
                          collar_u=float(np.linalg.norm(u, axis=-1)[msk2].max()),
                          collar_gradu=float(np.sqrt((g ** 2).sum((-1, -2)))[msk2].max()))), flush=True)


if __name__ == "__main__":
    c = sys.argv[1]
    if c == "check":
        cmd_check(int(sys.argv[2]), int(sys.argv[3]))
    elif c == "CF":
        cmd_CF(int(sys.argv[2]), int(sys.argv[3]))
    elif c == "beta":
        cmd_beta(int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]))
    elif c == "CR":
        cmd_CR(int(sys.argv[2]), int(sys.argv[3]), [int(v) for v in sys.argv[4].split(",")])
    elif c == "Uinf":
        cmd_Uinf()
