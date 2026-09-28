"""NUMERICAL sanity check of the explicit Argyris bounds |D^j(F - I F)| <= E_j(theta) M6 h^{6-j}, j=1,2.
Physical Argyris interpolation (DOFs: value, gradient, Hessian at vertices; physical unit-normal derivative at edge
midpoints), random triangles with min angle >= theta, F = phi(x/2) g(y) (Y=1), placed at random in the support,
various sizes h. Reports max ratio actual/bound (must be <= 1)."""
import numpy as np, sympy as sp, math, sys, json, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from assemble import E, m

xs, ys = sp.symbols('x y')
Fs = (1 - (xs / 2)**2)**6 * ys * (1 - ys)**6
der = {}
for a in range(4):
    for b in range(4 - a):
        der[(a, b)] = sp.lambdify((xs, ys), sp.diff(Fs, xs, a, ys, b) if a + b else Fs, 'numpy')
M6 = m[6]   # Y = 1
mons = [(i, d - i) for d in range(6) for i in range(d + 1)]


def mon_d(c, p, a, b, s):  # d^a_x d^b_y of ((x-c)/s)^i ((y-c)/s)^j at points p
    out = np.zeros((len(p), len(mons)))
    for k, (i, j) in enumerate(mons):
        if i < a or j < b:
            continue
        ci = math.factorial(i) / math.factorial(i - a); cj = math.factorial(j) / math.factorial(j - b)
        out[:, k] = ci * cj * ((p[:, 0] - c[0]) / s)**(i - a) * ((p[:, 1] - c[1]) / s)**(j - b) / s**(a + b)
    return out


def argyris(Pts, s):
    c = Pts.mean(0)
    rows, rhs = [], []
    for P in Pts:
        P1 = P[None]
        for (a, b) in [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]:
            rows.append(mon_d(c, P1, a, b, s)[0]); rhs.append(der[(a, b)](*P))
    for i in range(3):
        p, q = Pts[(i + 1) % 3], Pts[(i + 2) % 3]
        mid = (p + q) / 2; tv = q - p; n = np.array([tv[1], -tv[0]]) / np.linalg.norm(tv)
        r = n[0] * mon_d(c, mid[None], 1, 0, s)[0] + n[1] * mon_d(c, mid[None], 0, 1, s)[0]
        rows.append(r); rhs.append(n[0] * der[(1, 0)](*mid) + n[1] * der[(0, 1)](*mid))
    coef = np.linalg.solve(np.array(rows), np.array(rhs))
    return c, coef


if __name__ == '__main__':
    rng = np.random.default_rng(1)
    worst = {1: 0, 2: 0}
    for trial in range(600):
        th = math.radians(rng.choice([10, 15, 20, 30, 45]))
        while True:
            P = rng.random((3, 2))
            ang = []
            for i in range(3):
                u = P[(i + 1) % 3] - P[i]; v = P[(i + 2) % 3] - P[i]
                ang.append(math.acos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1)))
            if min(ang) >= th and min(ang) < th + math.radians(8):
                break
        h = 10**rng.uniform(-2, -0.7)
        diam = max(np.linalg.norm(P[i] - P[j]) for i in range(3) for j in range(3))
        P = P / diam * h + np.array([rng.uniform(-1.8, 1.8), rng.uniform(0.02, 0.95)])
        c, coef = argyris(P, h)
        # sample points
        lam = rng.dirichlet([1, 1, 1], 400); Q = lam @ P
        E1, E2 = E(min(ang))
        # j=1
        g_err = np.stack([mon_d(c, Q, 1, 0, h) @ coef - der[(1, 0)](*Q.T), mon_d(c, Q, 0, 1, h) @ coef - der[(0, 1)](*Q.T)], 1)
        e1 = np.max(np.linalg.norm(g_err, axis=1))
        H = np.stack([mon_d(c, Q, 2, 0, h) @ coef - der[(2, 0)](*Q.T), mon_d(c, Q, 1, 1, h) @ coef - der[(1, 1)](*Q.T),
                      mon_d(c, Q, 0, 2, h) @ coef - der[(0, 2)](*Q.T)], 1)
        e2 = max(np.max(np.abs(np.linalg.eigvalsh(np.array([[a, b], [b, d]])))) for a, b, d in H)
        worst[1] = max(worst[1], e1 / (E1 * M6 * h**5)); worst[2] = max(worst[2], e2 / (E2 * M6 * h**4))
    print("max actual/bound over 600 triangles: j=1:", worst[1], " j=2:", worst[2])
