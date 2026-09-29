"""Family II (interior slit lab): is the bad space spanned by the 8 needle 'vertex+mean' representers?
B = span{ 1_N, 1_N', R_N(z), R_N(up), R_N(dn), R_N'(f), R_N'(up), R_N'(dn) }  (R_T(x): L2(T) Riesz rep of q -> q|_T(x)).
Prints restricted inf-sup values on B (mean-zero part) /eps and the eigenvectors as functionals on the data
  d(q) = (mean_N q, mean_N' q, q_N(z), q_N(up), q_N(dn), q_N'(f), q_N'(up), q_N'(dn)).
Usage: python3 probe6.py n eps..."""
import sys, json, numpy as np, scipy.linalg as sla
from slit_lab import *

def run(n, eps):
    S, g = build(n, eps)
    Sd, M, one, sysm, _ = schur(S)
    tN = tri_of(S, g["z"], g["dn"], g["up"]); tNp = tri_of(S, g["f"], g["up"], g["dn"])
    V = S.verts
    reps = []
    for t in (tN, tNp):
        r = np.zeros(M.shape[0]); r[10 * t] = 1.0          # constant 1 on t (monomial 1)
        reps.append(r / (r @ M @ r))                        # (rep, q)_M = mean_t q ... up to |t|: (1_t,q)/|t|
    for t, xs in ((tN, ("z", "up", "dn")), (tNp, ("f", "up", "dn"))):
        for x in xs:
            reps.append(point_rep(S, sysm, t, V[g[x]]))
    R = np.array(reps).T                                    # columns
    # mean-zero combos: (1, R c)_M = 0
    c1 = one @ M @ R
    Q, _ = np.linalg.qr(c1[:, None], mode="complete"); Z = Q[:, 1:]
    E = R @ Z
    A = E.T @ Sd @ E; Bm = E.T @ M @ E
    ev, W = sla.eigh(.5 * (A + A.T), .5 * (Bm + Bm.T))
    b = np.sqrt(np.maximum(ev, 0))
    coef = Z @ W                                            # coefficients alpha on the 8 representers
    # (rho, q)_M = sum alpha_i d_i(q) with d_1 = mean? rep_1 = 1_t/|t| -> (rep_1,q)=mean_t q
    out = dict(n=n, eps=eps, betaB_over_eps=[float(f"{x/eps:.4g}") for x in b], betaB_over_sqrteps=[float(f"{x/np.sqrt(eps):.4g}") for x in b])
    for j in range(3):
        a = coef[:, j]; a = a / np.abs(a).max()
        out[f"f{j}"] = [float(f"{x:.4f}") for x in a]
    print(json.dumps(out), flush=True)

if __name__ == "__main__":
    n = int(sys.argv[1])
    for e in sys.argv[2:]:
        run(n, float(e))
