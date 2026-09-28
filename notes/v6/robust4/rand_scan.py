"""Random admissible 3- and 4-triangle boundary stars (flat and tilted): scale-invariant kappa of sigma_z relative to
Ker_{T_e} (3-dim, Thm U) and to K_{T_e} (2-dim).  Minimum angle of every triangle >= amin degrees.
python rand_scan.py n amin seed"""
import sys, json
import numpy as np
import stars as S


def minang(P):
    a = []
    for i in range(3):
        u = P[(i + 1) % 3] - P[i]; w = P[(i + 2) % 3] - P[i]
        a.append(np.degrees(np.arccos(np.clip(u @ w / np.linalg.norm(u) / np.linalg.norm(w), -1, 1))))
    return min(a)


def sample(rng, amin):
    while True:
        l = rng.uniform(0.6, 1.6); d = rng.uniform(-0.25, 0.25) * (rng.random() < 0.7)
        g1 = np.array([np.cos(d / 2), -np.sin(d / 2)]); g2 = -l * np.array([np.cos(d / 2), np.sin(d / 2)])
        nt = 3 if rng.random() < 0.6 else 4
        if nt == 3:
            w1 = np.array([rng.uniform(0.2, 1.5), rng.uniform(0.4, 1.8)]); w0 = np.array([rng.uniform(-1.2, 0.5), rng.uniform(0.4, 1.8)])
            P = np.array([[0, 0], g1, w1, w0, g2]); fan = [[0, 1, 2], [0, 2, 3], [4, 0, 3]]
        else:
            w1 = np.array([rng.uniform(0.4, 1.5), rng.uniform(0.4, 1.8)]); wm = np.array([rng.uniform(-0.4, 0.4), rng.uniform(0.6, 1.8)])
            w0 = np.array([rng.uniform(-1.5, -0.4), rng.uniform(0.4, 1.8)])
            P = np.array([[0, 0], g1, w1, wm, w0, g2]); fan = [[0, 1, 2], [0, 2, 3], [0, 3, 4], [len(P) - 1, 0, 4]]
        ok = all(np.cross(P[b] - P[a], P[c] - P[a]) > 0 for a, b, c in fan) and all(minang(P[f]) >= amin for f in fan)
        if ok:
            return P, fan


def run(n, amin, seed):
    rng = np.random.default_rng(seed); worst = None; ks = []
    for it in range(n):
        P, fan = sample(rng, amin)
        S._ns["GAMMA"] = (1, len(P) - 1)
        M, info = S._star_forms_M(P, fan, 4)
        T = P[fan[0]]; Tp = P[fan[-1]]
        kKer = S.kappa(M, S.KerT(T)); kKerp = S.kappa(M, S.KerT(Tp)); kK = S.kappa(M, S.KT(T))
        ks.append((kKer, kKerp, kK))
        if worst is None or kKer < worst[0]:
            worst = (kKer, P.tolist(), len(fan), info["dimV"])
    ks = np.array(ks)
    print(json.dumps(dict(n=n, amin=amin, seed=seed, kKer_min=float(ks[:, 0].min()), kKer_q05=float(np.quantile(ks[:, 0], .05)),
                          kKerp_min=float(ks[:, 1].min()), kK_min=float(ks[:, 2].min()), worst=worst)), flush=True)


if __name__ == "__main__":
    run(int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]))
