"""Bend angle scan: straight slit of slit_lab, then rotate the far apex f about w by angle g (so the two needle axes
meet at w at angle pi - g).  Prediction: beta ~ C (eps + g) (cross-difference functional, Prop. B) up to O(1).
Usage: python3 probe9.py n"""
import sys, json, time, numpy as np
from slit_lab import *

def run(n, eps, g):
    t0 = time.time()
    S0, gg = build(n, eps)
    V = S0.verts.copy(); w = (V[gg["up"]] + V[gg["dn"]]) / 2
    f = gg["f"]; r = V[f] - w; c, s = np.cos(g), np.sin(g)
    V[f] = w + np.array([c * r[0] - s * r[1], s * r[0] + c * r[1]])
    T = orient(V, S0.tris)
    S = svn.Space(V, T, {0, n + 1}, None)
    Sd, M, one, sysm, _ = schur(S)
    E = interior_basis(S)
    b, _ = beta_sub(Sd, M, E, k=2)
    tN = tri_of(S, gg["z"], gg["dn"], gg["up"]); tNp = tri_of(S, f, gg["up"], gg["dn"])
    rr = (point_rep(S, sysm, tN, V[gg["up"]]) - point_rep(S, sysm, tN, V[gg["dn"]])
          - point_rep(S, sysm, tNp, V[gg["up"]]) + point_rep(S, sysm, tNp, V[gg["dn"]]))
    print(json.dumps(dict(n=n, eps=eps, g=g, beta=[float(f"{x:.4g}") for x in b], beta_over_eps_plus_g=[float(f"{x/(eps+g):.4g}") for x in b],
                          cross_ratio=float(f"{dual_norm_ratio(Sd, M, rr):.4g}"), cross_over_eps_plus_g=float(f"{dual_norm_ratio(Sd, M, rr)/(eps+g):.4g}"),
                          secs=round(time.time() - t0, 1))), flush=True)

if __name__ == "__main__":
    n = int(sys.argv[1])
    for eps in (1e-3, 1e-4):
        for g in (0.0, 1e-4, 1e-3, 1e-2, 3e-2, 0.1, 0.3):
            run(n, eps, g)
