import sys, json, time, numpy as np
from slit_lab import *
n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
for eps in [0.1, 0.03, 0.01, 0.003]:
    t0 = time.time()
    S, g = build(n, eps)
    Sd, M, one, sysm, _ = schur(S)
    E = interior_basis(S); b, V = beta_sub(Sd, M, E, k=3)
    q = V[:, 0]
    mass = np.array([q[10*t:10*t+10] @ sysm["Ml"][t] @ q[10*t:10*t+10] for t in range(len(S.tris))]); mass /= mass.sum()
    tN = tri_of(S, g["z"], g["dn"], g["up"]); tNp = tri_of(S, g["f"], g["up"], g["dn"])
    print(json.dumps(dict(n=n, eps=eps, beta=list(np.round(b, 6)), b_over_eps=list(np.round(b / eps, 4)),
          massN=round(mass[tN], 4), massNp=round(mass[tNp], 4), secs=round(time.time()-t0, 1))), flush=True)
    # values of the mode on the needles
    Z, F, U, D = (S.verts[g[k]] for k in ("z", "f", "up", "dn"))
    for nm, t, pts in (("N", tN, {"z": Z, "up": U, "dn": D, "mid_e": (U + D) / 2, "mid_long_up": (Z + U) / 2, "mid_long_dn": (Z + D) / 2, "c": (Z + U + D) / 3}),
                       ("N'", tNp, {"f": F, "up": U, "dn": D, "mid_e": (U + D) / 2, "mid_long_up": (F + U) / 2, "mid_long_dn": (F + D) / 2, "c": (F + U + D) / 3})):
        sc = np.sqrt(sum(mass)) 
        print("   ", nm, {k: round(eval_p(S, q, t, x), 4) for k, x in pts.items()})
