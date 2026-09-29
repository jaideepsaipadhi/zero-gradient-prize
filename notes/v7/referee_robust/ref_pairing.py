"""Referee: independent recomputation of c_f = lim P_f / l^{k+2} for phi=(r-1)^k on UP (k=4), rho = 1 and 2.
Exact UP geometry (z=(1,0), a=(cos al, sin al), c=(1+rho l) a, l = 2 sin(pi/N)); exact phi; projection onto P_{k-1}(T_e)
in PHYSICAL centred/scaled coordinates via high-order Gauss quadrature at 50 digits; Richardson in l.
t_f = B_1 - B_{f+1} (Bernstein on e, parameter from z to a)."""
import sys
import mpmath as mp
mp.mp.dps = 50
k = int(sys.argv[1]) if len(sys.argv) > 1 else 4
NQ = 30
def gl(n):
    xs, ws = [], []
    # Golub-Welsch-free: use mpmath legendre roots
    for i in range(1, n + 1):
        x0 = mp.cos(mp.pi * (i - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
        for _ in range(100):
            p0, p1 = mp.mpf(1), x0
            for j in range(2, n + 1):
                p0, p1 = p1, ((2 * j - 1) * x0 * p1 - (j - 1) * p0) / j
            dp = n * (x0 * p1 - p0) / (x0**2 - 1)
            dx = p1 / dp; x0 -= dx
            if abs(dx) < mp.mpf(10)**(-45): break
        xs.append((x0 + 1) / 2); ws.append(1 / ((1 - x0**2) * dp**2))   # weight 2/((1-x^2)P'^2) /2
    return xs, ws
GX, GW = gl(NQ)
def binom(n, r): return mp.binomial(n, r)
def pairing(N, rho):
    al = 2 * mp.pi / N; l = 2 * mp.sin(mp.pi / N)
    z = (mp.mpf(1), mp.mpf(0)); a = (mp.cos(al), mp.sin(al)); c = ((1 + rho * l) * a[0], (1 + rho * l) * a[1])
    cx = (z[0] + a[0] + c[0]) / 3; cy = (z[1] + a[1] + c[1]) / 3
    phi = lambda X, Y: ((X * X + Y * Y - 1) / (mp.sqrt(X * X + Y * Y) + 1))**k
    mons = [(i, j) for i in range(k) for j in range(k - i)]
    bas = lambda X, Y: [((X - cx) / l)**i * ((Y - cy) / l)**j for i, j in mons]
    m = len(mons)
    G = mp.matrix(m, m); b = mp.matrix(m, 1)
    for s, ws in zip(GX, GW):
        for t, wt in zip(GX, GW):
            u = s * (1 - t); w = t; W = ws * wt * (1 - t)
            X = z[0] + (a[0] - z[0]) * u + (c[0] - z[0]) * w
            Y = z[1] + (a[1] - z[1]) * u + (c[1] - z[1]) * w
            B = bas(X, Y); f = phi(X, Y)
            for i in range(m):
                b[i] += W * B[i] * f
                for j in range(m):
                    G[i, j] += W * B[i] * B[j]
    co = mp.lu_solve(G, b)
    P = [mp.mpf(0)] * (k - 2)
    for u, wu in zip(GX, GW):
        X = z[0] + (a[0] - z[0]) * u; Y = z[1] + (a[1] - z[1]) * u
        B = bas(X, Y)
        eta = phi(X, Y) - sum(co[i] * B[i] for i in range(m))
        Bs = [binom(k, i) * u**i * (1 - u)**(k - i) for i in range(k + 1)]
        for f in range(1, k - 1):
            P[f - 1] += wu * eta * (Bs[1] - Bs[f + 1]) * l
    return l, [p / l**(k + 2) for p in P]
for rho in (1, 2):
    rows = []
    for N in (256, 512, 1024, 2048, 4096):
        l, sc = pairing(N, mp.mpf(rho))
        rows.append(sc)
        rich = [2 * sc[i] - rows[-2][i] for i in range(len(sc))] if len(rows) > 1 else None
        print('rho', rho, 'N', N, [mp.nstr(v, 12) for v in sc], 'rich', [mp.nstr(v, 12) for v in rich] if rich else '', flush=True)
    print('claimed k=4: rho^3/900 =', mp.nstr(mp.mpf(rho)**3 / 900, 12), ' rho^3/525 =', mp.nstr(mp.mpf(rho)**3 / 525, 12))
