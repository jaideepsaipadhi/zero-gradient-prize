"""Soundness check (floating point) of the CORRECTED (P^T) bound: it must never exceed the optimal numerical kappa_T
of robust4 (stars.py: full star, exact dual norms, robust3 pk.py assembler, validated there against solves) on random
admissible 3/4-triangle stars (same sampler as robust5/sanity_PT.py).  usage: python sanity_fix.py AMIN SEED N"""
import sys, math, json
import numpy as np
sys.path.insert(0, '../robust4'); sys.path.insert(0, '../robust3'); sys.path.insert(0, '../robust5')
import stars as S
import cert_PT as C
import cert_fix as CF2
import sympy as sp
X, Y = CF2.X, CF2.Y
Fs = sp.Symbol('F')
q0 = {k: 13 * (CF2.Q1[k][0] + CF2.gC[k[0]] * CF2.gC[k[1]] * CF2.XX * Fs) + 60 * CF2.S1 * CF2.Q1[k][1] for k in CF2.Q1}
q1 = {k: CF2.Q1[k][1] + 3 * CF2.Q1[k][2] for k in CF2.Q1}
f0 = sp.lambdify((X, Y, Fs), CF2.comb(q0) / CF2.DL); f1 = sp.lambdify((X, Y), CF2.comb(q1) / CF2.DL)
def kbound(x1, y1, F):
    ct2 = CF2.CT2_of(x1 * x1 + 1, y1 * y1 + 1, x1 + y1, math.pi, max, math.sqrt, 1.0)
    return 1 / math.sqrt(f0(x1, y1, F) + ct2 * f1(x1, y1))
def ang(P):
    out = []
    for i in range(3):
        u = P[(i + 1) % 3] - P[i]; w = P[(i + 2) % 3] - P[i]
        out.append(math.degrees(math.acos(np.clip(u @ w / np.linalg.norm(u) / np.linalg.norm(w), -1, 1))))
    return out
amin = float(sys.argv[1]); rng = np.random.default_rng(int(sys.argv[2])); N = int(sys.argv[3])
ratios = []; n = 0
while n < N:
    m = 3 if rng.random() < 0.6 else 4
    d = rng.uniform(-0.3, 0.3)
    g2 = -rng.uniform(0.6, 1.6) * np.array([math.cos(d), math.sin(d)])
    angs = np.sort(rng.uniform(0, 1, m - 1)); thE = math.atan2(g2[1], g2[0]) % (2 * math.pi)
    pts = [np.array([1.0, 0.0])] + [rng.uniform(0.6, 1.6) * np.array([math.cos(thE * a), math.sin(thE * a)]) for a in angs] + [g2]
    P = np.array([[0.0, 0.0]] + [list(p) for p in pts]); fan = [[0, j, j + 1] for j in range(1, m + 1)]
    fan[-1] = [len(P) - 1, 0, len(P) - 2]
    tris = [P[f] for f in fan]
    if not all(np.cross(t[1] - t[0], t[2] - t[0]) > 0 for t in tris) or min(min(ang(t)) for t in tris) < amin:
        continue
    n += 1
    S._ns["GAMMA"] = (1, len(P) - 1)
    M, info = S._star_forms_M(P, fan, 4)
    kn = S.kappa(M, S.KerT(P[fan[0]]))
    z, a, c, b = P[0], P[1], P[2], P[3]
    cot = lambda A, B, Cc: np.dot(B - A, Cc - A) / abs(np.cross(B - A, Cc - A))
    x1, y1 = cot(z, a, c), cot(a, c, z); x2, y2 = cot(z, c, b), cot(c, b, z)
    kb = kbound(x1, y1, C.Ffun(x2, y2)); ratios.append(kn / kb)
print(json.dumps(dict(amin=amin, n=n, ratio_min=float(min(ratios)), median=float(np.median(ratios)), max=float(max(ratios)),
                      violations=int(sum(r < 1 - 1e-6 for r in ratios)))))
