"""Critical Nitsche penalty: mu_c = max_v [2<d_n v, v> - a_h(v,v)] / ((1/h)|v|^2_Gamma) over V_R,
so that N_h(v,v) >= (1 - mu_c/mu) (mu/h)|v|^2_Gamma - ... ; N_h is indefinite for mu < mu_c.
Computed exactly via the Schur complement onto the Gamma_h trace dofs.  gamma_c = mu_c hG/h.
usage: python mucrit.py N1,N2,..."""
import sys, json
sys.path.insert(0, "../../../code")
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, scipy.linalg as sl
import svn_k as K
for N in [int(s) for s in sys.argv[1].split(",")]:
    S = K.build(N, 4)
    sysm = K.assemble(S, 1.0, lambda X, Y: np.zeros(np.shape(X) + (2,)))
    Pen = sysm["Pen"].tocsr(); Q = -(sysm["A"] - 1.0 * Pen)        # Q = Dn + Dn^T - Avol
    ndof = Q.shape[0]
    dD = K.dofs(S.dnodes).ravel()
    bmask = np.abs(Pen).sum(1).A.ravel() > 0
    b = np.where(bmask)[0]; i = np.setdiff1d(np.setdiff1d(np.arange(ndof), dD), b)
    Qbb = Q[b][:, b].toarray(); Qbi = Q[b][:, i]; Aii = (-Q[i][:, i]).tocsc()
    X = spl.splu(Aii).solve(Qbi.T.toarray())
    Sch = Qbb + Qbi @ X
    Mb = Pen[b][:, b].toarray()
    lam = sl.eigh(0.5 * (Sch + Sch.T), Mb, eigvals_only=True)
    muc = lam.max()
    print(json.dumps(dict(N=N, h=S.h, hG=S.hGamma, mu_c=muc, gamma_c=muc * S.hGamma / S.h,
                          gamma_at_mu100=100 * S.hGamma / S.h, top5=list(np.sort(lam)[-5:]))), flush=True)
