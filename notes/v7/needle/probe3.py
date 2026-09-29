"""Theorem S check in the slit lab: (i) gamma_off with mean-zero q vanishing on the needles (preimage in V_h);
(ii) inf-sup of the SLIT-DOMAIN pair: velocities vanishing on the needles N,N' (all their P4 nodes), pressures on the
non-needle interior triangles, mean zero.  Usage: python3 probe3.py n eps..."""
import sys, json, time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, scipy.linalg as sla
from slit_lab import *
from infsup import laplace

def run(n, eps):
    t0 = time.time()
    S, g = build(n, eps)
    sysm = svn.assemble(S, 1.0, lambda X, Y: np.zeros(X.shape + (2,)), None)
    A = laplace(S); B = sysm["B"]
    M = sp.block_diag(list(sysm["Ml"]), format="csr").toarray()
    tN = tri_of(S, g["z"], g["dn"], g["up"]); tNp = tri_of(S, g["f"], g["up"], g["dn"])
    E = interior_basis(S)
    keep = [j for j in range(E.shape[1]) if E[:, j].nonzero()[0][0] // 10 not in (tN, tNp)]
    Eo = E[:, keep]
    # mean-zero within the kept space
    onek = Eo.T @ (M @ np.where(np.arange(M.shape[0]) % 10 == 0, 1.0, 0.0))
    Q, _ = np.linalg.qr(onek[:, None], mode="complete"); Eo0 = Eo @ Q[:, 1:]
    out = dict(n=n, eps=eps)
    for label, extra in (("full_V", []), ("slit_V", list(S.ids[tN]) + list(S.ids[tNp]))):
        fixed = np.union1d(svn.dofs(S.dnodes).ravel(), svn.dofs(np.array(extra, dtype=int)).ravel()) if extra else svn.dofs(S.dnodes).ravel()
        free = np.setdiff1d(np.arange(A.shape[0]), fixed)
        Af = A[free][:, free].tocsc(); Bf = B[:, free].tocsr()
        X = spl.splu(Af).solve(Bf.T.toarray()); Sd = Bf @ X; Sd = .5 * (Sd + Sd.T)
        if label == "slit_V":
            b, _ = beta_sub(Sd, M, Eo0, k=1); out["beta_slit"] = float(f"{b[0]:.4g}")
        else:
            # primal constant on Q_off,0 with preimages in the full V_h: 1/sqrt(max (Mq)^T Sd^+ (Mq) / q^T M q)
            ev, U = np.linalg.eigh(Sd); tol = 1e-12 * ev.max()
            P = U[:, ev > tol]; Sinv = P @ np.diag(1 / ev[ev > tol]) @ P.T
            G = M @ Eo0
            resid = np.linalg.norm(G - P @ (P.T @ G)) / np.linalg.norm(G)     # is M q in range(Sd)?
            K = G.T @ Sinv @ G; K = .5 * (K + K.T)
            evp = sla.eigh(K, Eo0.T @ M @ Eo0, eigvals_only=True)
            out["gamma_off0"] = float(f"{1/np.sqrt(evp[-1]):.4g}"); out["range_resid"] = float(f"{resid:.2g}")
            out["min_nonzero_eig_Sd"] = float(f"{ev[ev > tol].min():.3g}")
    out["secs"] = round(time.time() - t0, 1)
    print(json.dumps(out), flush=True)

if __name__ == "__main__":
    n = int(sys.argv[1])
    for e in sys.argv[2:]:
        run(n, float(e))
