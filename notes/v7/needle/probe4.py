"""Structure of the 3 eps-modes and the sqrt(eps)-mode (interior slit, n=8): symmetry classes under the two
reflections of the kite (y -> -y : up<->dn ; x -> -x : N<->N') and vertex values.  Usage: python3 probe4.py eps"""
import sys, json, numpy as np
from slit_lab import *
eps = float(sys.argv[1]); n = 8
S, g = build(n, eps)
Sd, M, one, sysm, _ = schur(S)
E = interior_basis(S)
b, V = beta_sub(Sd, M, E, k=4)
tN = tri_of(S, g["z"], g["dn"], g["up"]); tNp = tri_of(S, g["f"], g["up"], g["dn"])
Z, F, U, D = (S.verts[g[k]] for k in ("z", "f", "up", "dn")); W0 = (U + D) / 2
for j in range(4):
    q = V[:, j]
    def ev(t, x): return eval_p(S, q, t, x)
    # sample on each needle at scaled points (s along axis from apex, tau in [-1,1] across)
    def samp(t, apex):
        pts = []
        for s in (0.25, 0.5, 0.75, 1.0):
            for tau in (-1, 0, 1):
                c = apex + s * (W0 - apex); half = s * (U - W0)
                pts.append(ev(t, c + tau * half))
        return np.array(pts)
    a = samp(tN, Z); c = samp(tNp, F)
    sc = np.abs(np.concatenate([a, c])).max()
    a /= sc; c /= sc
    oddy = np.allclose(a.reshape(4, 3)[:, 0], -a.reshape(4, 3)[:, 2], atol=2e-2)
    eveny = np.allclose(a.reshape(4, 3)[:, 0], a.reshape(4, 3)[:, 2], atol=2e-2)
    symx = np.allclose(a, c, atol=2e-2); asymx = np.allclose(a, -c, atol=2e-2)
    print(json.dumps(dict(mode=j, beta=float(f"{b[j]:.4g}"), y_odd=bool(oddy), y_even=bool(eveny), NNp_sym=bool(symx), NNp_antisym=bool(asymx),
                          N_rows_s_by_tau=np.round(a.reshape(4, 3), 3).tolist(), Np_rows=np.round(c.reshape(4, 3), 3).tolist())))
