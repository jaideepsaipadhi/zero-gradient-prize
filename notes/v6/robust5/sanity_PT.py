"""Soundness check (floating point): the certified lower bound (cert_PT.kappa2) never exceeds the optimal numerical
kappa_T of robust4 (stars.py, full star, exact dual norms) on random admissible 3- and 4-triangle stars."""
import sys, math, json
import numpy as np
sys.path.insert(0, '../robust4'); sys.path.insert(0, '../robust3')
import stars as S
import cert_PT as C

def ang(P):
    out = []
    for i in range(3):
        u = P[(i + 1) % 3] - P[i]; w = P[(i + 2) % 3] - P[i]
        out.append(math.degrees(math.acos(np.clip(u @ w / np.linalg.norm(u) / np.linalg.norm(w), -1, 1))))
    return out

rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
amin = float(sys.argv[1]) if len(sys.argv) > 1 else 20
ratios = []; n = 0
while n < 120:
    m = 3 if rng.random() < 0.6 else 4
    d = rng.uniform(-0.3, 0.3)
    g2 = -rng.uniform(0.6, 1.6) * np.array([math.cos(d), math.sin(d)])
    angs = np.sort(rng.uniform(0, 1, m - 1)); th0 = 0; thE = math.atan2(g2[1], g2[0]) % (2 * math.pi)
    rays = [thE * a for a in angs]
    pts = [np.array([1.0, 0.0])] + [rng.uniform(0.6, 1.6) * np.array([math.cos(t), math.sin(t)]) for t in rays] + [g2]
    P = np.array([[0.0, 0.0]] + [list(p) for p in pts]); fan = [[0, j, j + 1] for j in range(1, m + 1)]
    fan[-1] = [len(P) - 1, 0, len(P) - 2]
    tris = [P[f] for f in fan]
    if not all(np.cross(t[1] - t[0], t[2] - t[0]) > 0 for t in tris) or min(min(ang(t)) for t in tris) < amin:
        continue
    n += 1
    S._ns["GAMMA"] = (1, len(P) - 1)
    M, info = S._star_forms_M(P, fan, 4)
    kn = S.kappa(M, S.KerT(P[fan[0]]))
    # cot parameters of T1 = (z, a, c), T2 = (z, c, b)
    z, a, c, b = P[0], P[1], P[2], P[3]
    cot = lambda A, B, Cc: np.dot(B - A, Cc - A) / abs(np.cross(B - A, Cc - A))
    x1, y1 = cot(z, a, c), cot(a, c, z); x2, y2 = cot(z, c, b), cot(c, b, z)
    kb = C.kappa2(x1, y1, C.Ffun(x2, y2), math.pi, max, math.sqrt) ** 0.5
    ratios.append(kn / kb)
    assert kb <= kn * (1 + 1e-6), (kb, kn)
print(json.dumps(dict(amin=amin, n=n, ratio_numeric_over_bound_min=float(min(ratios)),
                      median=float(np.median(ratios)), max=float(max(ratios)))))
