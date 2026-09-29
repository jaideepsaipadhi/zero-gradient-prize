"""notes/v7/mu100: what does the critical mode look like?  For GS, EFF (= N_h - m_h) and CNS*, take the singular
vector of the saddle pencil at the largest real critical mu (shift-invert, sigma = 100) and report its wall trace:
  fn  = ||u.n||^2_Gamma / ||u||^2_Gamma        (normal fraction)
  fnm = sum_e |e| (mean_e u.n)^2 / ||u.n||^2   (edge-mean share of the normal trace; m_h's sharp constant acts on it)
  m_h(u,u) / ((mu/h) ||u.n||^2)  at mu = mu_c   (how much of the normal penalty the layer removes on this mode)
usage: python crit_modes.py N1,N2
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, scipy.sparse.linalg as spl
from cns_crit import Prob, svn


def traces(pb, u):
    S = pb.S; full = np.zeros(2 * S.nn); dD = svn.dofs(S.dnodes).ravel()
    fr = np.setdiff1d(np.arange(2 * S.nn), dD); full[fr] = u
    nn = tt = mean2 = 0.0
    for (t, ref, Xe, we, n) in S.edge_data():
        ph = svn.basis(ref[:, 0], ref[:, 1]); d = svn.dofs(S.ids[t]); U = ph @ full[d]
        un = U @ n; ut = U @ np.array([-n[1], n[0]])
        nn += we @ un ** 2; tt += we @ ut ** 2; mean2 += (we @ un) ** 2 / we.sum()
    return nn, tt, mean2


def main():
    for N in [int(s) for s in sys.argv[1].split(",")]:
        pb = Prob(N); out = dict(N=N)
        for which in ("GS", "EFF", "CNS"):
            Ms = pb.saddle(100.0, which); lu = spl.splu(Ms, permc_spec="COLAMD", diag_pivot_thresh=1.0)
            n = Ms.shape[0]

            def mv(x):
                y = np.zeros(n); y[:pb.nf] = pb.P @ x[:pb.nf]; return lu.solve(y)
            nu, V = spl.eigs(spl.LinearOperator((n, n), matvec=mv, dtype=float), k=6, which="LM", tol=1e-10)
            mus = 100 - 1 / nu; real = [i for i in range(len(mus)) if abs(mus[i].imag) < 1e-6 * abs(mus[i])]
            i = max(real, key=lambda j: mus[j].real); mu = mus[i].real
            u = V[:pb.nf, i]; u = (u * np.exp(-1j * np.angle(u[np.argmax(abs(u))]))).real
            nn, tt, m2 = traces(pb, u)
            mh = float(u @ (pb.Mw @ u))
            out[which] = dict(mu_c=mu, fn=nn / (nn + tt), fnm=m2 / nn, mh_over_normal_pen=mh / (mu / pb.S.h * nn),
                              mh_over_full_pen=mh / (mu / pb.S.h * (nn + tt)))
        print(json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
