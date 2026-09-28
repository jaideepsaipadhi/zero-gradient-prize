"""Theorem T1 made explicit on mesh (a) (strong_bc.make_mesh_alt) and mesh (s) (svn.make_mesh), test A.

For every wall vertex z whose star forces grad v|_T(z) = 0 (A_z = {0}: m <= 2, three distinct lines),
    ||grad(u~ - v)||_{L2(T)} >= gamma_k |T|^{1/2} |grad u~(z)|_F - ||grad u~ - grad u~(z)||_{L2(T)},
gamma_k = 2/(k(k+1)) (local_jets.py (D)).  Summing over the (disjoint) locked stars gives the rigorous bound
    E >= LB := gamma_k * S^{1/2} - R,   S = sum_{z locked} sum_{T in star z} |T| |grad u~(z)|^2,
    R = ( sum ||grad u~ - grad u~(z)||^2_{L2(T)} )^{1/2}   (computed by 12-point-rule quadrature on T).
Output: LB, S^{1/2}, R and the measured strong error from results/v5/strong.
"""
import os, sys, json, numpy as np
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "..", "..", "..", "code"))
os.environ.setdefault("OMP_NUM_THREADS", "1")
import svn, strong_bc

k = 4; gam = 2.0 / (k * (k + 1))
ex = svn.make_exact("A")
# degree-6 symmetric rule on the reference triangle via a tensor Duffy rule (plenty for smooth integrands)
g, w = np.polynomial.legendre.leggauss(8)
g = (g + 1) / 2; w = w / 2
A, B = np.meshgrid(g, g, indexing="ij"); WA, WB = np.meshgrid(w, w, indexing="ij")
QX = (A * (1 - B)).ravel(); QY = B.ravel(); QW = (WA * WB * (1 - B)).ravel()   # sum QW = 1/2

res = []
for mesh in ("alt", "std"):
    for N in (16, 24, 32, 48, 64, 96, 128, 192, 256):
        verts, tris, circ, outer = (strong_bc.make_mesh_alt(N) if mesh == "alt" else svn.make_mesh(N))
        oncirc = np.abs(np.linalg.norm(verts, axis=1) - 1) < 1e-12
        star = {}
        for t, T in enumerate(tris):
            for v in T:
                if oncirc[v]:
                    star.setdefault(v, []).append(t)
        locked = [v for v, st in star.items() if len(st) <= 2]
        mult = {}
        for v in locked:
            for t in star[v]:
                mult[t] = mult.get(t, 0) + 1
        assert all(c == 1 for c in mult.values()), "locked stars overlap"
        S = 0.0; R2 = 0.0
        for v in locked:
            Gz = ex["gu"](verts[v:v + 1, 0], verts[v:v + 1, 1])[0]
            for t in star[v]:
                P = verts[tris[t]]
                J = np.array([P[1] - P[0], P[2] - P[0]]).T; det = abs(np.linalg.det(J))
                S += 0.5 * det * np.sum(Gz ** 2)
                X = P[0][0] + J[0, 0] * QX + J[0, 1] * QY; Y = P[0][1] + J[1, 0] * QX + J[1, 1] * QY
                Gq = ex["gu"](X, Y)
                R2 += det * np.sum(QW[:, None, None] * (Gq - Gz) ** 2)
        hG = max(np.linalg.norm(verts[a] - verts[b]) for T in tris for a in T for b in T
                 if oncirc[a] and oncirc[b] and a != b)
        LB = gam * np.sqrt(S) - np.sqrt(R2)
        fn = os.path.join(ROOT, "..", "..", "..", "results", "v5", "strong",
                          f"N{N}_A_zero" + ("_alt" if mesh == "alt" else "") + ".jsonl")
        E = json.loads(open(fn).readline())["H1"] if os.path.exists(fn) else float("nan")
        r = dict(mesh=mesh, N=N, hG=hG, n_locked=len(locked), n_wall=len(star), sqrtS=np.sqrt(S),
                 gamma_sqrtS=gam * np.sqrt(S), R=np.sqrt(R2), LB=LB, E_strong=E)
        res.append(r)
        print(json.dumps({a: (round(b, 6) if isinstance(b, float) else b) for a, b in r.items()}))
with open(os.path.join(ROOT, "lock_bound.log"), "w") as f:
    for r in res:
        f.write(json.dumps(r) + "\n")
