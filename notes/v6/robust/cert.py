"""
Lower-bound certificate for the CNS* hydrostatic response (REPORT.md, Theorem R2.3).
For each boundary triangle T_e (Gamma_h edge e with end-vertices a,b) and r in P_{k-3}(T),
    v_r = curl(lam_a^2 lam_b^2 r)  in Z_h  (supported in T_e, zero on the two interior edges).
With eta = phi - pi_T phi, the local dual value
    d_e = sup_r <eta, v_r.n>_e / |||v_r|||_T,   |||v|||_T^2 = |v|_{1,T}^2 + (mu/h)|v|_e^2 + |e| |dn v|_e^2,
and D_h = (sum_e d_e^2)^{1/2} = sup over span of the local fields (disjoint supports).
Theorem R2.3:  nu * (M + C_xi) * |||u_h^{CNS*}||| >= D_h.
usage: python cert.py N1,N2,... [k] [mu]
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../code"))
import numpy as np
import svn_k as K
from hydro import PHIS


def cert(N, k=4, mu=100.0, phiname="sin"):
    phi, gphi = PHIS[phiname]
    S = K.build(N, k); R = S.R
    Ml_all = None
    rexp = [(a, b) for a in range(k - 2) for b in range(k - 2 - a)]       # P_{k-3} in (lam_a, lam_b)
    Dsq = 0.0; ds = []
    for (t, ref, Xe, we, n), (tt, u, w) in zip(S.edge_data(), S.bedges):
        assert t == tt
        V = S.verts[S.tris[t]]
        # barycentric gradients
        Jm = np.array([[V[0, 0], V[1, 0], V[2, 0]], [V[0, 1], V[1, 1], V[2, 1]], [1, 1, 1.0]])
        Ji = np.linalg.inv(Jm)                     # lam = Ji @ [x, y, 1]
        gl = Ji[:, :2]                             # grad lam_i (rows)
        lam = lambda X: (Ji @ np.vstack([X[:, 0], X[:, 1], np.ones(len(X))])).T
        nodes = S.P0[t][None, :] + (S.J[t] @ R.REFNODES.T).T
        def curlphi(X, a, b):
            L = lam(X); la, lb = L[:, u], L[:, w]
            r = la**a * lb**b
            gr = (a * la**max(a - 1, 0) * lb**b)[:, None] * gl[u] if a else 0 * X
            gr = gr + ((b * la**a * lb**max(b - 1, 0))[:, None] * gl[w] if b else 0 * X)
            g = (2 * la * lb**2 * r)[:, None] * gl[u] + (2 * la**2 * lb * r)[:, None] * gl[w] + (la**2 * lb**2)[:, None] * gr
            return np.stack([g[:, 1], -g[:, 0]], 1)          # curl = (d_y, -d_x)
        Vn = [curlphi(nodes, a, b) for (a, b) in rexp]      # nodal values -> P_k field (exact)
        # volume gradient Gram
        gx, gy = S.grads(R.QX, R.QY); gx, gy = gx[t], gy[t]
        wq = R.QW * abs(S.detJ[t])
        Gm = np.zeros((len(rexp),) * 2); ell = np.zeros(len(rexp))
        ph = R.basis(ref[:, 0], ref[:, 1]); dX, dY = R.dbasis(ref[:, 0], ref[:, 1])
        Jt = S.Jinv[t]; gxe = Jt[0, 0] * dX + Jt[1, 0] * dY; gye = Jt[0, 1] * dX + Jt[1, 1] * dY
        dn = n[0] * gxe + n[1] * gye
        # eta = phi - pi_T phi
        Xq = S.phys(R.QX, R.QY)[t]
        M = np.einsum('q,qi,qj->ij', wq, R.QP_Q, R.QP_Q)
        cpi = np.linalg.solve(M, R.QP_Q.T @ (wq * phi(Xq[:, 0], Xq[:, 1])))
        qe = R.pbasis(ref[:, 0], ref[:, 1])
        eta = phi(Xe[:, 0], Xe[:, 1]) - qe @ cpi
        le = we.sum()
        for i, Ui in enumerate(Vn):
            vi_e = ph @ Ui; dni = dn @ Ui
            ell[i] = we @ (eta * (vi_e @ n))
            for j, Uj in enumerate(Vn):
                vol = wq @ ((gx @ Ui) * (gx @ Uj) + (gy @ Ui) * (gy @ Uj)).sum(1)
                bd = (mu / S.h) * we @ ((ph @ Ui) * (ph @ Uj)).sum(1) + le * we @ ((dn @ Ui) * (dn @ Uj)).sum(1)
                Gm[i, j] = vol + bd
        d2 = ell @ np.linalg.solve(Gm, ell)
        ds.append(np.sqrt(d2)); Dsq += d2
    ds = np.array(ds)
    return dict(N=N, k=k, mu=mu, phi=phiname, hG=S.hGamma, h=S.h, D=np.sqrt(Dsq),
                D_over_hG45=np.sqrt(Dsq) / S.hGamma**(k + 0.5),
                frac_edges_de_gt_half_mean=float(np.mean(ds > 0.5 * ds.mean())))


if __name__ == "__main__":
    Ns = [int(s) for s in sys.argv[1].split(",")]
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    mu = float(sys.argv[3]) if len(sys.argv) > 3 else 100.0
    for N in Ns:
        print(json.dumps(cert(N, k, mu)), flush=True)
