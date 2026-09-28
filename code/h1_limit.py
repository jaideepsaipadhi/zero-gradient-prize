"""
The H^1 leak constant of the printed method GS (notes/v3/h1_constant.tex).

Claim (hc:thm): (mu/h)(u~ - u_h) -> U in H^1, where U is the Stokes solution
    -Delta U + grad P = 0,  div U = 0  in Omega = (-L,L)^2 minus unit disk,
    U = (p - pbar) e_r  on Gamma (= -(p-pbar) n),   U = 0 on dR.
Hence ||grad(u~ - u_h)|| mu/(h G) -> K := ||grad U||_{L2(Omega)} / G.

Test P: p = x^2 y - y + x/2, so on Gamma  p - pbar = -(3/4) sin t + (1/4) sin 3t + (1/2) cos t,
stream function data Psi(1,t) = (3/4) cos t - (1/12) cos 3t + (1/2) sin t, dPsi/dr(1,t) = 0,
Psi = 0, dPsi/dnu = 0 on dR  (the constant on dR is 0 by the x -> -x antisymmetry of each class).

Three independent computations:
  annulus  : exact (40-digit) Stokes energies of concentric annuli 1<r<R; by domain monotonicity of the
             Dirichlet (minimum) principle,  K(R = L sqrt 2) <= K <= K(R = L)   (rigorous bracket).
  mps      : method of particular solutions (Michell series, r=1 conditions exact, least squares on dR),
             energy by Gauss quadrature on 8 polar sectors.  Converges to ~1e-9.
  fem      : the paper's GS solver (code/svn.py) on test P for mu up to 1e8 and N up to 128; prints the
             normalised H1 error, and the DIRECT distance || (mu/h)(u~-u_h) - U ||_{H1(Omega_h)} / ||grad U||.
Usage: python h1_limit.py [annulus|mps|fem|all]      (1 core, < 2 GB)
"""
import os, sys, time
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import numpy as np

L = 2.5
G2 = 7 * np.pi / 8                 # ||p - pbar||^2_{L2(Gamma)}
G = np.sqrt(G2)
A1C, A1S, A3C = 0.75, 0.5, -1.0 / 12   # Psi(1,t) = A1C cos t + A1S sin t + A3C cos 3t


# ------------------------------------------------------------------ annulus bracket
def annulus():
    import mpmath as mp
    mp.mp.dps = 40

    def E(n, R):
        if n == 1:
            B = [lambda r: r, lambda r: 1 / r, lambda r: r ** 3, lambda r: r * mp.log(r)]
        else:
            B = [lambda r, n=n: r ** n, lambda r, n=n: r ** -n, lambda r, n=n: r ** (n + 2), lambda r, n=n: r ** (2 - n)]
        rows = [[g(1) for g in B], [mp.diff(g, 1) for g in B], [g(R) for g in B], [mp.diff(g, R) for g in B]]
        k = mp.lu_solve(mp.matrix(rows), mp.matrix([1, 0, 0, 0]))
        f = lambda r: sum(k[i] * B[i](r) for i in range(4))

        def dens(r):     # |Hess Psi|^2 r, Psi = f(r) cos(n t), theta-integrated / pi
            f0 = f(r); f1 = mp.diff(f, r); f2 = mp.diff(f, r, 2)
            return (f2 ** 2 + 2 * n ** 2 * ((f1 * r - f0) / r ** 2) ** 2 + (f1 / r - n ** 2 * f0 / r ** 2) ** 2) * r
        return mp.pi * mp.quad(dens, [1, R])

    print("== annulus 1<r<R: Stokes energies E_n (Psi = f(r) cos nt, f(1)=1) and K(R)")
    out = {}
    for R in (mp.mpf(L), mp.mpf(L) * mp.sqrt(2), 10, 100, 10 ** 4):
        E1, E3 = E(1, R), E(3, R)
        tot = (A1C ** 2 + A1S ** 2) * E1 + A3C ** 2 * E3
        K = mp.sqrt(tot / (7 * mp.pi / 8))
        out[float(R)] = float(K)
        print(f"R={mp.nstr(R, 8):>10}  E1={mp.nstr(E1, 14)}  E3={mp.nstr(E3, 14)}  K={mp.nstr(K, 14)}", flush=True)
    print(f"rigorous bracket for the square: {out[L * np.sqrt(2)]:.10f} <= K <= {out[L]:.10f}")
    return out


# ------------------------------------------------------------------ Michell-series (MPS) solution
RC = L * np.sqrt(2)


def _terms(n):
    if n == 1:
        s = 1.0 / RC
        return [(lambda r: r * s, lambda r: s + 0 * r, lambda r: 0 * r),
                (lambda r: 1 / r, lambda r: -1 / r ** 2, lambda r: 2 / r ** 3),
                (lambda r: (r * s) ** 3, lambda r: 3 * s ** 3 * r ** 2, lambda r: 6 * s ** 3 * r),
                (lambda r: r * np.log(r), lambda r: np.log(r) + 1, lambda r: 1 / r)]
    F = []
    for p, sc in ((n, RC ** -n), (-n, 1.0), (n + 2, RC ** -(n + 2)), (2 - n, 1.0)):
        F.append((lambda r, p=p, sc=sc: sc * r ** p, lambda r, p=p, sc=sc: sc * p * r ** (p - 1),
                  lambda r, p=p, sc=sc: sc * p * (p - 1) * r ** (p - 2)))
    return F


class Series:
    """Psi = sum_{n odd <= M} f_n(r) cos(n t) for data Psi(1,t)=cos(n0 t), Psi_r(1,t)=0, clamped on dR."""

    def __init__(self, n0, M=63, npts=4000):
        self.modes = list(range(1, M + 1, 2))
        info = []
        for n in self.modes:
            F = _terms(n)
            Aeq = np.array([[f[0](1.0) for f in F], [f[1](1.0) for f in F]])
            info.append((n, F, np.linalg.inv(Aeq[:, [1, 3]]), Aeq[:, [0, 2]]))
        self.info = info
        self.d = [1.0 if n == n0 else 0.0 for n in self.modes]
        s = np.cos(np.pi * (np.arange(npts) + 0.5) / npts) * L
        pts, nrm = [], []
        for (X, Y, nx, ny) in ((L + 0 * s, s, 1, 0), (-L + 0 * s, s, -1, 0), (s, L + 0 * s, 0, 1), (s, -L + 0 * s, 0, -1)):
            pts.append(np.stack([X, Y], 1)); nrm.append(np.tile([nx, ny], (npts, 1)))
        pts = np.concatenate(pts); nrm = np.concatenate(nrm)
        r = np.hypot(pts[:, 0], pts[:, 1]); t = np.arctan2(pts[:, 1], pts[:, 0])
        nr = nrm[:, 0] * np.cos(t) + nrm[:, 1] * np.sin(t); nt = -nrm[:, 0] * np.sin(t) + nrm[:, 1] * np.cos(t)
        nu = 2 * len(self.modes)

        def residual(x):
            self._set(x)
            P, Pr, Pt, *_ = self.polar(r, t)
            return np.concatenate([P, (Pr * nr + Pt / r * nt) * L])
        r0 = residual(np.zeros(nu))
        J = np.stack([residual(np.eye(nu)[i]) - r0 for i in range(nu)], 1)
        x, *_ , sv = np.linalg.lstsq(J, -r0, rcond=None)
        self.resid = np.abs(residual(x)).max()
        self.cond = sv[0] / sv[-1]

    def _set(self, x):
        cf = []
        for k, ((n, F, Adi, Af), d) in enumerate(zip(self.info, self.d)):
            fr = x[2 * k:2 * k + 2]
            dep = Adi @ (np.array([d, 0.0]) - Af @ fr)
            cf.append((n, F, np.array([fr[0], dep[0], fr[1], dep[1]])))
        self.cf = cf

    def polar(self, r, t):
        """Psi, Psi_r, Psi_t, and orthonormal-frame Hessian (H_rr, H_rt, H_tt)"""
        P = Pr = Pt = Prr = Prt = Ptt = 0
        for n, F, c in self.cf:
            f = sum(ci * Fi[0](r) for ci, Fi in zip(c, F)); f1 = sum(ci * Fi[1](r) for ci, Fi in zip(c, F))
            f2 = sum(ci * Fi[2](r) for ci, Fi in zip(c, F))
            C, S = np.cos(n * t), np.sin(n * t)
            P = P + f * C; Pr = Pr + f1 * C; Pt = Pt - n * f * S
            Prr = Prr + f2 * C; Prt = Prt - n * f1 * S; Ptt = Ptt - n * n * f * C
        return P, Pr, Pt, Prr, Prt / r - Pt / r ** 2, Pr / r + Ptt / r ** 2


def energy(s1, s2, q=64):
    """int_Omega Hess(Psi1):Hess(Psi2), Omega = square minus disk, 8 polar sectors x Gauss(q x q)"""
    g, w = np.polynomial.legendre.leggauss(q)
    tot = 0.0
    for k in range(8):
        a0 = -np.pi / 4 + k * np.pi / 4; a1 = a0 + np.pi / 4
        th = a0 + (g + 1) / 2 * (a1 - a0); wt = w * (a1 - a0) / 2
        Rt = L / np.maximum(np.abs(np.cos(th)), np.abs(np.sin(th)))
        for tj, wj, Rj in zip(th, wt, Rt):
            rr = 1 + (g + 1) / 2 * (Rj - 1); ww = w * (Rj - 1) / 2
            H1 = s1.polar(rr, tj + 0 * rr)[3:]; H2 = s2.polar(rr, tj + 0 * rr)[3:]
            tot += wj * np.sum(ww * rr * (H1[0] * H2[0] + 2 * H1[1] * H2[1] + H1[2] * H2[2]))
    return tot


def mps(Ms=(15, 31, 47, 63, 79)):
    print("== Michell-series solution in the square minus disk")
    K = None
    for M in Ms:
        s1, s3 = Series(1, M), Series(3, M)
        for q in (48, 96):
            E11, E13, E33 = energy(s1, s1, q), energy(s1, s3, q), energy(s3, s3, q)
            tot = (A1C ** 2 + A1S ** 2) * E11 + 2 * A1C * A3C * E13 + A3C ** 2 * E33
            K = np.sqrt(tot / G2)
            print(f"M={M:3d} quad={q:3d} bdry-resid={max(s1.resid, s3.resid):.1e} cond={s1.cond:.1e} "
                  f"E11={E11:.10f} E13={E13:.10f} E33={E33:.10f}  ||grad U||={np.sqrt(tot):.10f}  K={K:.10f}", flush=True)
    return K


class LimitField:
    """U = curl Psi, Psi = A1C Psi_1 + A1S Psi_1(r, t - pi/2) + A3C Psi_3; Cartesian values and gradients."""

    def __init__(self, M=63):
        self.s1, self.s3 = Series(1, M), Series(3, M)

    def _cart(self, X, Y):
        r = np.hypot(X, Y); t = np.arctan2(Y, X)
        Pr = Pt = Hrr = Hrt = Htt = 0
        for s, a, sh in ((self.s1, A1C, 0.0), (self.s1, A1S, -np.pi / 2), (self.s3, A3C, 0.0)):
            _, pr, pt, hrr, hrt, htt = s.polar(r, t + sh)
            Pr = Pr + a * pr; Pt = Pt + a * pt; Hrr = Hrr + a * hrr; Hrt = Hrt + a * hrt; Htt = Htt + a * htt
        c, sn = np.cos(t), np.sin(t)
        Px = c * Pr - sn * Pt / r; Py = sn * Pr + c * Pt / r
        Pxx = c * c * Hrr - 2 * c * sn * Hrt + sn * sn * Htt
        Pyy = sn * sn * Hrr + 2 * c * sn * Hrt + c * c * Htt
        Pxy = c * sn * (Hrr - Htt) + (c * c - sn * sn) * Hrt
        return Px, Py, Pxx, Pyy, Pxy

    def u(self, X, Y):
        Px, Py, *_ = self._cart(X, Y)
        return np.stack([Py, -Px], -1)

    def gu(self, X, Y):                       # G[..., c, d] = d_d U_c
        _, _, Pxx, Pyy, Pxy = self._cart(X, Y)
        return np.stack([np.stack([Pxy, Pyy], -1), np.stack([-Pxx, -Pxy], -1)], -2)


def fem(Ns=(16, 24, 32, 48, 64, 96), mus=(1e2, 1e3, 1e4, 1e6, 1e8)):
    import svn
    ex_P = svn.make_exact("P")
    Uf = LimitField()
    # sanity: trace of U on Gamma equals (p - pbar) e_r
    t = np.linspace(0, 2 * np.pi, 7)[:-1]
    Ut = Uf.u(np.cos(t), np.sin(t)); pd = np.cos(t) ** 2 * np.sin(t) - np.sin(t) + np.cos(t) / 2
    print("== FEM (code/svn.py, method GS, test P); trace check max|U - (p-pbar)e_r| on Gamma =",
          f"{np.abs(Ut - pd[:, None] * np.stack([np.cos(t), np.sin(t)], 1)).max():.1e}")
    print(f"{'N':>4} {'mu':>7} {'h/mu':>9} {'gamma*hG':>9} {'K_h=|grad e|mu/(hG)':>20} {'slip':>8} "
          f"{'dist=|grad((mu/h)e-U)|/|grad U|':>32}", flush=True)
    for N in Ns:
        for mu in mus:
            t0 = time.time()
            S = svn.Space(*svn.make_mesh(N))
            sysm = svn.assemble(S, mu, ex_P["f"], 0)
            u, p, dres, relres = svn.solve_lu(S, sysm, ex_P["u"], 0)   # SuperLU + refinement (as in the paper)
            E = svn.errors(S, u, ex_P, mu)
            E.update(h=S.h, hG=S.hGamma)
            h = E["h"]
            ex = dict(u=lambda X, Y, s=h / mu: -s * Uf.u(X, Y), gu=lambda X, Y, s=h / mu: -s * Uf.gu(X, Y))
            D = svn.errors(S, u, ex, mu)                     # || grad(u_h + (h/mu) U) ||   (e = -u_h)
            Kh = E["H1"] * mu / (h * G)
            dist = D["H1"] * mu / h / (G * 2.7565157)
            gam_hG = mu * E["hG"] ** 2 / h
            print(f"{N:4d} {mu:7.0e} {h/mu:9.2e} {gam_hG:9.2e} {Kh:20.6f} {E['bL2']*mu/(h*G):8.5f} {dist:32.5f}"
                  f"   (res {relres:.1e}, div {dres[0]:.1e}, {time.time()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode in ("annulus", "all"):
        annulus()
    if mode in ("mps", "all"):
        mps()
    if mode in ("fem", "all"):
        fem()
