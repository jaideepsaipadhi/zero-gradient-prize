"""Local (single-star) analysis of nearly singular wall stars, strong no-slip Scott-Vogelius.

Star at z = (0,0) on the circle of radius 1 centred at (0,-1) (Omega_h = exterior of the disc, above).
Wall neighbours z_+- at polar angle +-2t on that circle: chords from z point at angles -t (to z_+)
and pi+t (to z_-); turning angle phi = 2t; chord length h = 2 sin t.
Rays from z (ccw): e_first (angle -t), interior edges s_1..s_{m-1} (angles gam_j), e_last (pi+t).
Triangle T_i between consecutive rays, ray lengths r (default all = h).

Div-free vertex-gradient tuples (strong/REPORT.md, L1) in Hessian form H_i = J^{-1} G_i (symmetric):
   H_1 = alpha P(n_first),  H_{i+1} = H_i + beta_i P(s_i^perp),  H_m t_last = 0.
The exact gradient c = grad u(z) = d_nu u (x) nu corresponds to H = a P(nu), nu = (0,1).  a = 1.

Computed per star:
  dinf  = min over A_z of max_i |H_i - a P(nu)|_F                (SOCP via cvx-free LP on a polygonal norm)
  W     = Fekete wall oscillation  sum_{k in F} |l_k(n_last) - l_k(n_first)|   (Theorem N1 predicts dinf ~ W)
  dw    = min over A_z of (sum_i |T_i| |H_i - aP(nu)|_F^2)^{1/2} / h   (weighted, = d_z/h of S3; Thm N2)
Families (m = 3):
  I   : needle at the wall   : gam = (g0, pi+t-eps)            (angle eps between s_2 and e_last)
  II  : split interior edge  : gam = (g0-eps/2, g0+eps/2)      (angle eps between s_1 and s_2)
  III : paper-Theta small    : gam = (t+b, pi-t-b)             (s_1 || e_last, s_2 || e_first as b -> 0)
Predictions printed: pred_inf = min(1, phi/eps) (I, II; Thm N1 + Fekete asymptotics), min(1, phi/(phi+b)) (III);
  pred_w = min(1, phi/sqrt(eps)) (I, II; Thm N2), min(1, phi/sqrt(phi+b)) (III; observation only).
Usage: python3 local_ns.py
"""
import itertools
import numpy as np
from scipy.optimize import linprog


def P(u):
    return np.outer(u, u)


def vec(M):                       # isometry Sym2 -> R^3
    return np.array([M[0, 0], M[1, 1], np.sqrt(2) * M[0, 1]])


def unit(a):
    return np.array([np.cos(a), np.sin(a)])


def star(t, gam, r=None):
    angs = np.concatenate([[-t], np.asarray(gam, float), [np.pi + t]])
    m = len(angs) - 1
    h = 2 * np.sin(t)
    r = np.full(m + 1, h) if r is None else np.asarray(r, float)
    th = np.diff(angs)
    area = 0.5 * r[:-1] * r[1:] * np.sin(th)
    return angs, area, h, th


def constraint_system(angs):
    """unknowns x = (alpha, beta_1..beta_{m-1}); returns list of affine maps x -> vec(H_i) (as 3 x n matrices)
    and the closure matrix C (2 x n) with C x = 0."""
    m = len(angs) - 1
    n = m
    nrm = [unit(a + np.pi / 2) for a in angs]            # unit normals of the rays
    Hs = []
    cur = np.zeros((3, n)); cur[:, 0] = vec(P(nrm[0]))
    Hs.append(cur.copy())
    for i in range(1, m):
        cur = cur.copy(); cur[:, i] = vec(P(nrm[i]))
        Hs.append(cur)
    tl = unit(angs[-1])
    # H_m t_last = 0 : H_m = sum_j x_j Q_j,  (Q_j t)
    Qs = [P(nrm[0])] + [P(nrm[i]) for i in range(1, m)]
    C = np.array([[(Q @ tl)[k] for Q in Qs] for k in range(2)])
    return Hs, C


def target():
    return vec(P(np.array([0.0, 1.0])))


def d_weighted(angs, area, h):
    Hs, C = constraint_system(angs)
    N = null_space(C)
    g = target()
    A = np.vstack([np.sqrt(ar) * (H @ N) for H, ar in zip(Hs, area)])
    b = np.concatenate([np.sqrt(ar) * g for ar in area])
    y, *_ = np.linalg.lstsq(A, b, rcond=None)
    return np.linalg.norm(A @ y - b) / h


def null_space(C, tol=1e-13):
    u, s, vt = np.linalg.svd(C)
    rank = int((s > tol * max(1, s[0])).sum())
    return vt[rank:].T


def d_inf(angs):
    """min_x max_i ||H_i - aP(nu)||_inf (entrywise, on R^3 coordinates) s.t. C x = 0 : LP"""
    Hs, C = constraint_system(angs)
    n = C.shape[1]; g = target()
    # variables (x, s); minimise s
    c = np.zeros(n + 1); c[-1] = 1
    Aub = []; bub = []
    for H in Hs:
        for k in range(3):
            row = np.zeros(n + 1); row[:n] = H[k]; row[-1] = -1; Aub.append(row); bub.append(g[k])
            row = np.zeros(n + 1); row[:n] = -H[k]; row[-1] = -1; Aub.append(row); bub.append(-g[k])
    Aeq = np.hstack([C, np.zeros((2, 1))])
    res = linprog(c, A_ub=np.array(Aub), b_ub=np.array(bub), A_eq=Aeq, b_eq=np.zeros(2),
                  bounds=[(None, None)] * (n + 1), method="highs")
    return res.fun


def fekete_W(angs):
    nrm = [unit(a + np.pi / 2) for a in angs]
    V = [vec(P(u)) for u in nrm]
    best = None
    for tri in itertools.combinations(range(len(V)), 3):
        d = abs(np.linalg.det(np.array([V[i] for i in tri])))
        if best is None or d > best[0]:
            best = (d, tri)
    d, F = best
    M = np.array([V[i] for i in F]).T                   # columns
    # Lagrange functions: l_k(x) = coefficient of V[k] in the expansion of the dual ... use Cramer:
    def ell(k, x):
        Mk = M.copy(); Mk[:, k] = x
        return np.linalg.det(Mk) / np.linalg.det(M)
    W = sum(abs(ell(k, V[-1]) - ell(k, V[0])) for k in range(3))
    return W, d, F


def fam(name, t, eps, g0=np.pi / 3):
    if name == "I":
        gam = (g0, np.pi + t - eps)
    elif name == "II":
        gam = (g0 - eps / 2, g0 + eps / 2)
    elif name == "III":
        gam = (t + eps, np.pi - t - eps)
    else:
        raise ValueError
    return star(t, gam)


def theta_paper(th):
    return max(abs(np.sin(th[i] + th[i + 1])) for i in range(len(th) - 1))


if __name__ == "__main__":
    np.set_printoptions(precision=4)
    print("Check: m=2 lock (gam = pi/2): dinf, W, dw =",
          d_inf(star(0.02, [np.pi / 2])[0]), fekete_W(star(0.02, [np.pi / 2])[0])[0],
          d_weighted(*star(0.02, [np.pi / 2])[:3]))
    for name in ("I", "II", "III"):
        print(f"\nFamily {name}")
        print(" phi      eps       Theta_paper  dinf     W       dinf/W   pred_inf     dw      pred_w    dw/pred")
        for t in (0.04, 0.01, 0.0025):
            phi = 2 * t
            for eps in (0.5, 0.2, 0.05, phi * 4, phi, phi / 4, phi**2 * 4, phi**2, phi**2 / 4):
                if name == "III" and eps > 0.5:
                    continue
                angs, area, h, th = fam(name, t, eps)
                if np.any(th <= 0):
                    continue
                di = d_inf(angs); W, V, F = fekete_W(angs); dw = d_weighted(angs, area, h)
                if name == "I":
                    pinf = min(1, phi / eps); pw = min(1, phi / np.sqrt(eps))
                elif name == "II":
                    pinf = min(1, phi / eps); pw = min(1, phi / np.sqrt(eps))
                else:
                    pinf = min(1, phi / (phi + eps)); pw = min(1, phi / np.sqrt(phi + eps))
                print(f" {phi:.4f}  {eps:.3e}  {theta_paper(th):.3e}   {di:.3e}  {W:.3e}  {di / W:6.3f}  {pinf:.3e}  "
                      f"{dw:.3e}  {pw:.3e}  {dw / pw:6.3f}")
