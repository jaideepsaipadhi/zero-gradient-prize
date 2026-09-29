"""notes/v7/mu100: is CNS* coercive / near-critical / indefinite at the production penalty mu = 100 (k = 4)?

Production meshes: svn.make_mesh(N) (L = 2.5, N/4 layers), P4 / P3-disc, Dirichlet on the square.
Free velocity dofs only.  Forms (per unit mu: A_mu = A0 + mu P):
  GS   : A_mu on Z_h = ker B                                    (the paper's mu*, Table tab:thresh)
  EFF  : A_mu - C^T M^{-1} C on Z_h   (= N_h - m_h, the effective-penalty form of notes/v7/robust, Lemma mh)
  EFFm : A_mu - Cm^T M^{-1} C (sym. part) on Z_h, Cm = mean-corrected C (production CNS* coupling)
  CNS* : the true saddle system [A_mu, -(B - Cm)^T; B, 0].  On Z_h it is a Petrov-Galerkin problem
         (trial Z_h = ker B, test Z* = ker(B - Cm)), so "coercivity" is not defined; the relevant
         quantity is the set of real mu at which the square system is SINGULAR.
Outputs per N:
  mu*_GS, mu*_EFF, mu*_EFFm  by inertia bisection (unpivoted symmetric LU of K + r B^T M^{-1} B, fresh per call)
  eigen-mu of the pencils  (M0 + mu M1) x = 0  near mu = sigma (shift-invert Arnoldi), for GS, EFF, CNS*:
     mu_c = largest real mu with a singular system.  For GS this must reproduce mu*_GS.
  at mu = 100: number of negative pivots (= negative eigenvalues on Z_h) of GS, EFF, EFFm.
usage: python cns_crit.py N1,N2,...  [sigma=100]
"""
import os, sys, json, time
for v_ in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v_, "2")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
import warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
import svn

zf = lambda X, Y: np.zeros(np.shape(X) + (2,))


class Prob:
    def __init__(self, N):
        v, t, c, o = svn.make_mesh(N)
        S = svn.Space(v, t, c, o); self.S = S
        s0 = svn.assemble(S, 0.0, zf, 0); s1 = svn.assemble(S, 1.0, zf, 0, want=("A",))
        dD = svn.dofs(S.dnodes).ravel(); fr = np.setdiff1d(np.arange(2 * S.nn), dD)
        self.A0 = s0["A"][fr][:, fr].tocsc(); self.P = (s1["A"] - s0["A"])[fr][:, fr].tocsc()
        B = s0["B"]; C = s0["C"]; Minv = s0["Minv"]
        e = np.zeros(B.shape[0]); e[0::10] = 1.0; phi = C.T @ e
        Cm = (C - sp.csr_matrix(s0["gGamma"][:, None]) @ sp.csr_matrix(phi[None, :])).tocsr()
        self.B = B[:, fr].tocsr(); self.C = C[:, fr].tocsr(); self.Cm = Cm[:, fr].tocsr(); self.Minv = Minv
        self.Pen = (self.B.T @ Minv @ self.B).tocsc()
        self.Mw = (self.C.T @ Minv @ self.C).tocsc()                      # m_h
        Mwm = (self.Cm.T @ Minv @ self.C).tocsc(); self.Mwm = ((Mwm + Mwm.T) / 2).tocsc()
        self.rho_G = S.h / S.hGamma                                          # gamma = mu / rho_G
        # first-layer heights H_e and chord lengths
        He, le = [], []
        for (tt, u, w) in S.bedges:
            a = S.verts[S.tris[tt, u]]; b = S.verts[S.tris[tt, w]]; L = np.linalg.norm(b - a)
            He.append(abs(S.detJ[tt]) / L); le.append(L)
        self.He, self.le = np.array(He), np.array(le)
        self.nf = len(fr); self.npr = B.shape[0]

    def K(self, mu, which):
        K = self.A0 + mu * self.P
        if which == "EFF":
            K = K - self.Mw
        elif which == "EFFm":
            K = K - self.Mwm
        return K

    def negpiv(self, mu, which, r=1e6):
        K = (self.K(mu, which) + r * self.Pen).tocsc(); K = ((K + K.T) / 2).tocsc()
        lu = spl.splu(K, permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0, options=dict(SymmetricMode=True))
        return int((lu.U.diagonal() < 0).sum())

    def threshold(self, which, lo=1.0, hi=1e4, tol=2e-4, r=1e6):
        if self.negpiv(hi, which, r) > 0:
            return np.inf
        if self.negpiv(lo, which, r) == 0:
            return lo
        while hi / lo > 1 + tol:
            mid = np.sqrt(lo * hi)
            if self.negpiv(mid, which, r) > 0:
                lo = mid
            else:
                hi = mid
        return hi

    def saddle(self, mu, which):
        Bt = {"GS": self.B, "EFF": self.B, "CNS": (self.B - self.Cm).tocsr(), "CNSnomean": (self.B - self.C).tocsr()}[which]
        KK = self.K(mu, "EFF" if which == "EFF" else "GS")
        return sp.bmat([[KK, -Bt.T], [self.B, None]], format="csc")

    def eigmu(self, which, sigma=100.0, nev=16):
        """real mu with (M0 + mu M1) singular, nearest sigma: (M0+sigma M1)^{-1} M1 x = x/(sigma - mu)."""
        Ms = self.saddle(sigma, which)
        lu = spl.splu(Ms, permc_spec="COLAMD", diag_pivot_thresh=1.0)
        n = Ms.shape[0]; P = self.P

        def mv(x):
            y = np.zeros(n, dtype=np.result_type(x, float)); y[:self.nf] = P @ x[:self.nf]
            return lu.solve(y)
        op = spl.LinearOperator((n, n), matvec=mv, dtype=float)
        nu = spl.eigs(op, k=nev, which="LM", return_eigenvectors=False, tol=1e-10, maxiter=5000)
        mus = sigma - 1.0 / nu
        return list(mus)


def main():
    Ns = [int(s) for s in sys.argv[1].split(",")]
    sigma = float(sys.argv[2]) if len(sys.argv) > 2 else 100.0
    for N in Ns:
        t0 = time.time(); pb = Prob(N)
        rec = dict(N=N, nf=pb.nf, npr=pb.npr, hG=pb.S.hGamma, h=pb.S.h, h_over_hG=pb.rho_G,
                   H_over_e=(float((pb.He / pb.le).min()), float((pb.He / pb.le).max())),
                   shift20_hG_over_He=(float((20 * pb.S.hGamma / pb.He).min()), float((20 * pb.S.hGamma / pb.He).max())))
        # local effective penalty at mu=100, in units of 1/hG:  gamma_loc_e = (mu/h - 20/H_e) * hG
        rec["gamma100"] = 100 / pb.rho_G
        rec["gamma_eff100_range"] = (float((100 / pb.S.h - 20 / pb.He).min() * pb.S.hGamma),
                                     float((100 / pb.S.h - 20 / pb.He).max() * pb.S.hGamma))
        for which in ("GS", "EFF", "EFFm"):
            mu = pb.threshold(which)
            rec[f"mustar_{which}"] = mu; rec[f"gammastar_{which}"] = mu / pb.rho_G
            rec[f"negpiv100_{which}"] = pb.negpiv(100.0, which)
            print(f"  N={N} {which}: mu*={mu:.4f} negpiv(100)={rec[f'negpiv100_{which}']}  [{time.time()-t0:.0f}s]", flush=True)
        # r-independence check of EFF threshold at small N
        if N <= 32:
            rec["mustar_EFF_r1e5"] = pb.threshold("EFF", r=1e5); rec["mustar_EFF_r1e7"] = pb.threshold("EFF", r=1e7)
        for which in ("GS", "EFF", "CNS"):
            ms = pb.eigmu(which, sigma)
            real = [float(z.real) for z in ms if abs(z.imag) < 1e-6 * max(1, abs(z))]
            cplx = [(float(z.real), float(z.imag)) for z in ms if abs(z.imag) >= 1e-6 * max(1, abs(z))]
            rec[f"eigmu_{which}_real"] = sorted(real, reverse=True)[:8]
            rec[f"eigmu_{which}_cplx"] = cplx[:6]
            rec[f"muc_{which}"] = max(real) if real else None
            print(f"  N={N} eig {which}: real {np.round(sorted(real, reverse=True)[:6], 3)} cplx {np.round(cplx[:3], 3)} [{time.time()-t0:.0f}s]", flush=True)
        # far check: eigen-mu of CNS* nearest mu = 400 (any real singular mu >= 100 would show here)
        ms = pb.eigmu("CNS", 400.0, nev=8)
        rec["eigmu_CNS_near400"] = [(float(z.real), float(z.imag)) for z in ms]
        print(f"  N={N} eig CNS near 400: {np.round(ms, 3)}", flush=True)
        rec["secs"] = round(time.time() - t0, 1)
        print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
