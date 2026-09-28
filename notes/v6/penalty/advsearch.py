"""Adversarial float search: minimise lambda*_star over fans with straight Nitsche boundary,
subject to min triangle angle >= thmin (degrees).  Usage: python advsearch.py m thmin nstarts seed"""
import sys, numpy as np
from scipy.optimize import minimize
from patch import fan_polar

m = int(sys.argv[1]); thmin = float(sys.argv[2]) * np.pi / 180
ns = int(sys.argv[3]); seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
rng = np.random.default_rng(seed)


def decode(x):
    w = np.exp(x[:m]); al = np.pi * w / w.sum()
    ang = np.concatenate([[0], np.cumsum(al)]); ang[-1] = np.pi
    r = np.exp(x[m:]); r = r / r[0]
    return ang, r, al


def minangle(ang, r):
    out = []
    for i in range(m):
        O = np.zeros(2); A = r[i] * np.array([np.cos(ang[i]), np.sin(ang[i])])
        B = r[i + 1] * np.array([np.cos(ang[i + 1]), np.sin(ang[i + 1])])
        for p, q, s in ((O, A, B), (A, B, O), (B, O, A)):
            u, v = q - p, s - p
            out.append(np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1)))
    return min(out)


def obj(x):
    ang, r, al = decode(x)
    ma = minangle(ang, r)
    pen = 0.0
    if ma < thmin:
        pen = 1e3 * (thmin - ma) + 50
        return pen + 50
    lam, d = fan_polar(ang, r)
    return lam


best = (np.inf, None)
for s in range(ns):
    x0 = np.concatenate([rng.normal(0, .5, m), rng.normal(0, .5, m + 1)])
    if obj(x0) > 90:
        continue
    res = minimize(obj, x0, method="Nelder-Mead", options=dict(maxiter=3000, xatol=1e-6, fatol=1e-8))
    ang, r, al = decode(res.x)
    print(f"start {s}: lam={res.fun:.5f} alphas(deg)={np.round(al*180/np.pi,2)} radii={np.round(r,3)} minang={minangle(ang,r)*180/np.pi:.2f}", flush=True)
    if res.fun < best[0]:
        best = (res.fun, res.x)
ang, r, al = decode(best[1])
print("BEST", best[0], "alphas", np.round(al * 180 / np.pi, 3), "radii", np.round(r, 4), "dim", fan_polar(ang, r)[1])
