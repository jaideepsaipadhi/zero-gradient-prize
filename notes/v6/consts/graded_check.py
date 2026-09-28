"""NUMERICAL: the continuum-profile witness v = curl(I_Argyris F_Y) on a DYADICALLY GRADED half-plane mesh
(boundary edges of length 1; the mesh size doubles every layer: row 0 = unit squares split by the diagonal,
row j>=1 = squares of side 2^j on top of pairs of row-(j-1) squares, each cut into 3 triangles from the hanging
bottom midpoint; min angle atan(1/2)=26.57 deg). Every ball B(x0,3Y), Y>=1, contains a triangle of diameter >= Y,
so hypothesis (LQ_kappa) fails for every kappa<1 at every scale.  We report Y(2B-A)/C for several Y.
Also reports the same quotient on the UNIFORM unit mesh for comparison."""
import numpy as np, math
from validate_interp import mon_d
from realistic_kappa import tri_rule
import sympy as sp

gx, gw = np.polynomial.legendre.leggauss(8)


def dyadic_mesh(J, W):
    T = []
    for i in range(-W, W):                 # row 0
        a, b, c, d = (i, 0), (i + 1, 0), (i + 1, 1), (i, 1)
        T += [(a, b, c), (a, c, d)]
    for j in range(1, J + 1):
        s = 2**j; y0 = 2**j - 1
        for i in range(-W // s, W // s):
            x0 = i * s
            bl, br, tr, tl, mid = (x0, y0), (x0 + s, y0), (x0 + s, y0 + s), (x0, y0 + s), (x0 + s / 2, y0)
            T += [(bl, mid, tl), (mid, br, tr), (mid, tr, tl)]
    return [np.array(t, float) for t in T]


def uniform_mesh(H, W):
    T = []
    for j in range(H):
        for i in range(-W, W):
            a, b, c, d = (i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)
            T += [(a, b, c), (a, c, d)]
    return [np.array(t, float) for t in T]


def make_F(Y):
    x, y = sp.symbols('x y')
    F = (1 - (x / (2 * Y))**2)**6 * Y * (y / Y) * (1 - y / Y)**6
    der = {}
    for a in range(3):
        for b in range(3 - a):
            e = sp.diff(F, x, a, y, b) if a + b else F
            der[(a, b)] = sp.lambdify((x, y), e, 'numpy')
    return der


def argyris_F(P, der, s):
    c = P.mean(0)
    rows, rhs = [], []
    for Pt in P:
        for (a, b) in [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]:
            rows.append(mon_d(c, Pt[None], a, b, s)[0]); rhs.append(der[(a, b)](*Pt) if abs(Pt[0]) < 2 * Yg and Pt[1] < Yg else 0.0)
    for i in range(3):
        p, q = P[(i + 1) % 3], P[(i + 2) % 3]
        mid = (p + q) / 2; tv = q - p; n = np.array([tv[1], -tv[0]]) / np.linalg.norm(tv)
        rows.append(n[0] * mon_d(c, mid[None], 1, 0, s)[0] + n[1] * mon_d(c, mid[None], 0, 1, s)[0])
        inside = abs(mid[0]) < 2 * Yg and mid[1] < Yg
        rhs.append((n[0] * der[(1, 0)](*mid) + n[1] * der[(0, 1)](*mid)) if inside else 0.0)
    return c, np.linalg.solve(np.array(rows), np.array(rhs))


def quotient(T, Y):
    global Yg
    Yg = Y
    der = make_F(Y)
    A = B = C = 0.0; kap = 0.0
    for P in T:
        if P[:, 1].min() >= Y or P[:, 0].min() >= 2 * Y or P[:, 0].max() <= -2 * Y:
            continue
        s = max(np.linalg.norm(P[i] - P[j]) for i in range(3) for j in range(3)); kap = max(kap, s / Y)
        c, coef = argyris_F(P, der, s)
        Q, Wt = tri_rule(P)
        H = [mon_d(c, Q, 2, 0, s) @ coef, mon_d(c, Q, 1, 1, s) @ coef, mon_d(c, Q, 0, 2, s) @ coef]
        v1x, v1y, v2x, v2y = H[1], H[2], -H[0], -H[1]
        A += np.sum(Wt * 0.5 * ((2 * v1x)**2 + 2 * (v1y + v2x)**2 + (2 * v2y)**2))
        on = [k for k in range(3) if P[k, 1] == 0]
        if len(on) == 2:
            p, q = P[on[0]], P[on[1]]
            u = (gx + 1) / 2; Qb = p[None] + u[:, None] * (q - p)[None]; Wb = gw / 2 * np.linalg.norm(q - p)
            px = mon_d(c, Qb, 1, 0, s) @ coef; py = mon_d(c, Qb, 0, 1, s) @ coef
            hyy = mon_d(c, Qb, 0, 2, s) @ coef; hxy = mon_d(c, Qb, 1, 1, s) @ coef
            B += np.sum(Wb * (-hyy * py + hxy * (-px))); C += np.sum(Wb * (py**2 + px**2))
    return Y * (2 * B - A) / C, kap


if __name__ == "__main__":
    Td = dyadic_mesh(7, 256); Tu = uniform_mesh(40, 256)
    print("Y   dyadic: Y(2B-A)/C  (max diam/Y over support)   uniform unit mesh: Y(2B-A)/C")
    for Y in (1, 1.5, 2, 3, 4, 6, 8, 16, 32):
        qd, kd = quotient(Td, Y); qu, ku = quotient(Tu, Y)
        print(f"{Y:5.1f}  {qd:10.4f}  ({kd:5.2f})      {qu:10.4f}  ({ku:5.2f})", flush=True)
