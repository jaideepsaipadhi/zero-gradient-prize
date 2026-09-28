"""Referee5: V_R vs Z_h critical penalty on the same meshes (k=4).
V_R: as mucrit.py.  Z_h: add -r B^T Minv B (divergence penalty) before Schur elimination, r -> large."""
import sys, json
sys.path.insert(0, "/home/claude/zero-gradient-prize/code")
import numpy as np, scipy.sparse.linalg as spl, scipy.linalg as sl
import svn_k as K
for N in [int(s) for s in sys.argv[1].split(",")]:
    S = K.build(N, 4)
    sm = K.assemble(S, 1.0, lambda X, Y: np.zeros(np.shape(X) + (2,)))
    Pen = sm["Pen"].tocsr(); Q0 = -(sm["A"] - Pen)
    DtD = (sm["B"].T @ sm["Minv"] @ sm["B"]).tocsr()
    ndof = Q0.shape[0]; dD = K.dofs(S.dnodes).ravel()
    bm = np.abs(Pen).sum(1).A.ravel() > 0
    b = np.where(bm)[0]; i = np.setdiff1d(np.setdiff1d(np.arange(ndof), dD), b)
    out = dict(N=N, hG_over_h=S.hGamma / S.h)
    for r in [0.0, 1e3, 1e5, 1e7]:
        Q = (Q0 - r * DtD).tocsr()
        Qbb = Q[b][:, b].toarray(); Qbi = Q[b][:, i]; Aii = (-Q[i][:, i]).tocsc()
        X = spl.splu(Aii).solve(Qbi.T.toarray()); Sch = Qbb + Qbi @ X
        lam = sl.eigh(0.5 * (Sch + Sch.T), Pen[b][:, b].toarray(), eigvals_only=True)
        out["mu_c_r%g" % r] = lam.max()
    print(json.dumps(out), flush=True)
