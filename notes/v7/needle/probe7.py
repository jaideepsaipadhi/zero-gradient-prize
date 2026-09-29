"""Oblique slit: split w along a direction tau_g at angle (pi/2 - g) to the axis z->w->f (g = 0: the symmetric slit of
slit_lab).  Prediction (REPORT, Prop. B): beta <= C (eps + g) from the cross-difference functional; the eps-modes should
disappear once g >> eps.  Also prints the cross-difference functional ratio (its sup_v / norm).  Usage: python3 probe7.py n"""
import sys, json, time, numpy as np
from slit_lab import *

def split_obl(verts, tris, w, z, f, eps, g):
    d = verts[f] - verts[z]; d = d / np.linalg.norm(d); nrm = np.array([-d[1], d[0]])
    tg = np.cos(g) * nrm + np.sin(g) * d
    r = np.linalg.norm(verts[w] - verts[z])
    delta = 2 * r * np.tan(eps / 2) / np.cos(g)
    up = len(verts); dn = w
    V = np.vstack([verts, verts[w]])
    newt = []
    for tri in tris:
        tri = list(tri)
        if w in tri:
            c = verts[tri].mean(0) - verts[w]
            if np.dot(c, nrm) > 0:
                tri = [up if v == w else v for v in tri]
        newt.append(tri)
    V[up] = verts[w] + 0.5 * delta * tg
    V[dn] = verts[w] - 0.5 * delta * tg
    newt += [[z, dn, up], [f, up, dn]]
    return V, orient(V, newt), up, dn

def run(n, eps, g):
    t0 = time.time()
    verts, tris, vid = square_mesh(n)
    i0 = n // 2
    w = vid[i0, n // 2]; z = vid[i0 - 1, n // 2]; f = vid[i0 + 1, n // 2]
    V, T, up, dn = split_obl(verts, tris, w, z, f, eps, g)
    S = svn.Space(V, T, {int(vid[0, 0]), int(vid[1, 0])}, None)
    Sd, M, one, sysm, _ = schur(S)
    E = interior_basis(S)
    b, _ = beta_sub(Sd, M, E, k=3)
    tN = tri_of(S, z, dn, up); tNp = tri_of(S, f, up, dn)
    r = (point_rep(S, sysm, tN, V[up]) - point_rep(S, sysm, tN, V[dn])
         - point_rep(S, sysm, tNp, V[up]) + point_rep(S, sysm, tNp, V[dn]))
    ang = lambda a, b, c: np.degrees(np.arccos(np.dot(b - a, c - a) / np.linalg.norm(b - a) / np.linalg.norm(c - a)))
    print(json.dumps(dict(n=n, eps=eps, g=g, apex_deg=round(ang(V[z], V[up], V[dn]), 4), beta=[float(f"{x:.4g}") for x in b],
                          beta_over_eps=[float(f"{x/eps:.4g}") for x in b],
                          cross_fn_ratio=float(f"{dual_norm_ratio(Sd, M, r):.4g}"), secs=round(time.time() - t0, 1))), flush=True)

if __name__ == "__main__":
    n = int(sys.argv[1])
    for eps in (1e-2, 1e-3):
        for g in (0.0, 0.03, 0.1, 0.3, 0.6):
            run(n, eps, g)
