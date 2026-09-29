"""notes/v7/robust: NUMERICAL test of (S) on a genuinely LP^+ family with SMOOTHLY VARYING first-layer heights
(not rotation invariant, so not reducible to one Floquet cell).  Mesh of robust5/run_S.py with first-layer ratio
rho_i = 1 + A cos(3 theta_i) at wall vertex i ('smooth'); load a = cos 2 theta; k = 4; CNS* and GS (delta_a) solved.
Reports S = ||grad W^a|| gamma / h_Gamma and E = ||grad(W^a - delta_a)|| gamma / h_Gamma; (S) <=> both bounded in N.
usage: python smooth_S.py N1,N2 mu1,mu2 [A=0.35]"""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../../code")); sys.path.insert(0, os.path.join(HERE, "../../v6/robust2"))
sys.path.insert(0, os.path.join(HERE, "../../v6/robust5"))
import numpy as np
import svn_k as K
import leak as LK
from run_S import mesh, minangle


def run(N, mus, A):
    th = 2 * np.pi * np.arange(N) / N
    rhos = 1 + A * np.cos(3 * th)
    v, t, c, o = mesh(N, rhos)
    S = K.Space(v, t, c, o, K.Ref(4)); R = S.R
    for mu in mus:
        t0 = time.time()
        sysm = K.assemble(S, mu, lambda X, Y: np.zeros(np.shape(X) + (2,)))
        F = np.zeros(2 * S.nn)
        for (tt, ref, Xe, we, n) in S.edge_data():
            ph = R.basis(ref[:, 0], ref[:, 1]); d = K.dofs(S.ids[tt])
            a = np.cos(2 * np.arctan2(Xe[:, 1], Xe[:, 0]))
            for cc in range(2):
                F[d[:, cc]] += ph.T @ (we * a) * n[cc]
        uW, _, rel1 = LK.solve(S, sysm["A"], K.coupling(S, sysm, 1), sysm["B"], F)
        uD, _, rel2 = LK.solve(S, sysm["A"], K.coupling(S, sysm, 0), sysm["B"], F)
        gW = K.errors(S, uW, LK.ZERO, mu)["H1"]; gD = K.errors(S, uD, LK.ZERO, mu)["H1"]
        gE = K.errors(S, uW - uD, LK.ZERO, mu)["H1"]
        gam = mu * S.hGamma / S.h; hG = S.hGamma
        print(json.dumps(dict(N=N, A=A, mu=mu, gamma=gam, hG=hG, minangle=minangle(v, t), S=gW * gam / hG, SD=gD * gam / hG,
                              E=gE * gam / hG, E_layer=gE * gam**1.5 / hG**0.5, relres=max(rel1, rel2), time=round(time.time() - t0, 1))), flush=True)


if __name__ == "__main__":
    Ns = [int(x) for x in sys.argv[1].split(",")]; mus = [float(x) for x in sys.argv[2].split(",")]
    A = float(sys.argv[3]) if len(sys.argv) > 3 else 0.35
    for N in Ns:
        run(N, mus, A)
