"""NUMERICAL: the actual quotient Y(2B-A)/C of v = curl(I_Argyris F) on half-plane meshes of size kappa*Y
(Y=1, X=2, F = phi(x/2) g(y)), to compare the rigorous kappa_0 with what is really needed.
Meshes: squares of side kappa split by a diagonal (min angle 45 deg), optionally with random vertex perturbation
(boundary vertices move along y=0), reporting the min angle."""
import numpy as np, sympy as sp, math, sys
from validate_interp import argyris, mon_d, der

# degree-8 symmetric triangle rule via collapsed Gauss (tensor Gauss-Jacobi, exact for high degree)
gx, gw = np.polynomial.legendre.leggauss(8)


def tri_rule(P):
    a = (gx + 1) / 2; wa = gw / 2
    pts, wts = [], []
    for i in range(8):
        for j in range(8):
            u = a[i]; v = a[j] * (1 - u)
            pts.append(P[0] + u * (P[1] - P[0]) + v * (P[2] - P[0])); wts.append(wa[i] * wa[j] * (1 - u))
    u_, w_ = P[1] - P[0], P[2] - P[0]; area = abs(u_[0] * w_[1] - u_[1] * w_[0])
    return np.array(pts), np.array(wts) * area


def run(kappa, pert=0.0, seed=0):
    rng = np.random.default_rng(seed)
    nx = int(round(4.4 / kappa)); ny = int(round(1.2 / kappa))
    xs = -2.2 + kappa * np.arange(nx + 1); ys = kappa * np.arange(ny + 1)
    V = np.array([[x, y] for y in ys for x in xs])
    if pert:
        d = rng.uniform(-pert, pert, V.shape) * kappa
        d[V[:, 1] == 0, 1] = 0
        V = V + d
    idx = lambda i, j: j * (nx + 1) + i
    T = []
    for j in range(ny):
        for i in range(nx):
            a, b, c, dd = idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)
            T += [(a, b, c), (a, c, dd)]
    A = B = C = 0.0; minang = 180
    for tri in T:
        P = V[list(tri)]
        if P[:, 1].min() > 1.0 or P[:, 0].min() > 2.0 or P[:, 0].max() < -2.0:
            continue  # all DOFs of F vanish (F=0 with all derivatives there)
        for i in range(3):
            u = P[(i + 1) % 3] - P[i]; w = P[(i + 2) % 3] - P[i]
            minang = min(minang, math.degrees(math.acos(u @ w / np.linalg.norm(u) / np.linalg.norm(w))))
        c, coef = argyris(P, kappa)
        Q, W = tri_rule(P)
        H = [mon_d(c, Q, 2, 0, kappa) @ coef, mon_d(c, Q, 1, 1, kappa) @ coef, mon_d(c, Q, 0, 2, kappa) @ coef]
        # v = (psi_y, -psi_x): grad v = [[psi_xy, psi_yy],[-psi_xx,-psi_xy]]
        v1x, v1y, v2x, v2y = H[1], H[2], -H[0], -H[1]
        A += np.sum(W * 0.5 * ((2 * v1x)**2 + 2 * (v1y + v2x)**2 + (2 * v2y)**2))
        # boundary edge on y = 0
        on = [k for k in range(3) if abs(P[k, 1]) < 1e-14]
        if len(on) == 2:
            p, q = P[on[0]], P[on[1]]
            s = (gx + 1) / 2; Qb = p[None] + s[:, None] * (q - p)[None]; Wb = gw / 2 * np.linalg.norm(q - p)
            px = mon_d(c, Qb, 1, 0, kappa) @ coef; py = mon_d(c, Qb, 0, 1, kappa) @ coef
            hb = [mon_d(c, Qb, 2, 0, kappa) @ coef, mon_d(c, Qb, 1, 1, kappa) @ coef, mon_d(c, Qb, 0, 2, kappa) @ coef]
            v1, v2 = py, -px
            dn1, dn2 = -hb[2], hb[1]           # d_n = -d_y : -(psi_yy), -(-psi_xy)
            B += np.sum(Wb * (dn1 * v1 + dn2 * v2)); C += np.sum(Wb * (v1**2 + v2**2))
    return (2 * B - A) / C, minang


if __name__ == "__main__":
    for pert in (0.0, 0.25):
        for kappa in (1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05):
            q, ma = run(kappa, pert)
            print(f"pert={pert} kappa={kappa:5.3f} diam/Y={kappa*math.sqrt(2)*(1+pert):5.3f} min angle={ma:5.1f}  Y(2B-A)/C = {q:.4f}", flush=True)
