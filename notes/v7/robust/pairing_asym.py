"""notes/v7/robust: the order-(k+3/2) lower bound on UP for phi = (r-1)^k, k = 4 (REPORT.tex, Sec. 3.3).
On UP every boundary star is congruent (rotation), so Thm rb:cert (first form, NO remainder) with robust6's corrected
fields v_1, v_2 (traces B_1-B_2, B_1-B_3 on e; their stream parts have zero normal trace on e and do not enter)
gives ||grad W|| >= K^{-1/2} N^{1/2} max_f |P_f| / Phi*(v_f), Phi* the same for every star and convergent as N -> oo.
So the order is (at most) k+3/2 from below as soon as  c_f := lim P_f / l^{k+2} != 0, where
     P_f(N) = int_e (I - pi_{T_e}) phi . (B_1 - B_{f+1}) ds,   l = 2 sin(pi/N) = |e|,
T_e = (z, a, c), z = (1,0), a = (cos al, sin al), c = (1 + rho l)(cos al, sin al) (forward diagonal, as on UP).
Double precision, Duffy/Gauss quadrature of order 24 in reference coordinates, (r-1) evaluated as (r^2-1)/(r+1)
(no cancellation); the projection is solved in the reference monomial basis. Loss of accuracy ~ l^{-2} (P ~ l^2 |eta| |e|)
is harmless (1e-9 relative at N=4096); cross-checked by the two orders of quadrature.
usage: python pairing_asym.py [rho=1] [phi: r (default) | rc]"""
import sys
import numpy as np
k = 4
rho = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
which = sys.argv[2] if len(sys.argv) > 2 else "r"


def phi(X, Y):
    r = np.hypot(X, Y); d = (X * X + Y * Y - 1) / (r + 1)
    if which == "r":
        return d**k
    return d**k * np.cos(4 * np.arctan2(Y, X))       # (r-1)^k cos 4 theta (A not constant)


def quad(n):
    g, wg = np.polynomial.legendre.leggauss(n); g = (g + 1) / 2; wg = wg / 2
    U, V = np.meshgrid(g, g, indexing="ij"); W = np.outer(wg, wg)
    return (U * (1 - V)).ravel(), V.ravel(), (W * (1 - V)).ravel()      # Duffy: u = s(1-t), w = t


def pairing(N, nq=24):
    al = 2 * np.pi / N; l = 2 * np.sin(np.pi / N)
    z = np.array([1.0, 0.0]); a = np.array([np.cos(al), np.sin(al)]); c = (1 + rho * l) * a
    J = np.stack([a - z, c - z], 1)
    uq, wq, Wq = quad(nq)
    X = z[0] + J[0, 0] * uq + J[0, 1] * wq; Y = z[1] + J[1, 0] * uq + J[1, 1] * wq
    f = phi(X, Y)
    mons = [(i, j) for i in range(k) for j in range(k - i)]
    Mq = np.stack([uq**i * wq**j for i, j in mons], 1)
    G = Mq.T @ (Wq[:, None] * Mq); b = Mq.T @ (Wq * f); cf = np.linalg.solve(G, b)
    g, wg = np.polynomial.legendre.leggauss(nq); g = (g + 1) / 2; wg = wg / 2
    Xe = z[0] + (a - z)[0] * g; Ye = z[1] + (a - z)[1] * g
    eta = phi(Xe, Ye) - np.stack([g**i * 0**j if j else g**i for i, j in mons], 1) @ cf
    B = [np.math.comb(k, i) * g**i * (1 - g)**(k - i) if hasattr(np, "math") else None for i in range(k + 1)]
    from math import comb
    B = [comb(k, i) * g**i * (1 - g)**(k - i) for i in range(k + 1)]
    return l, [np.sum(wg * eta * (B[1] - B[f])) * l for f in (2, 3)]


if __name__ == "__main__":
    print("phi =", "(r-1)^4" if which == "r" else "(r-1)^4 cos4theta", " rho =", rho)
    prev = None
    for N in (64, 128, 256, 512, 1024, 2048, 4096, 8192):
        l, P = pairing(N); l2, P2 = pairing(N, 32)
        sc = [p / l**(k + 2) for p in P]
        line = "N=%5d P1/l^6=% .10e P2/l^6=% .10e  quad-check %.1e" % (N, sc[0], sc[1], abs(P2[0] - P[0]) / abs(P[0]))
        if prev is not None:
            line += "  Richardson(O(l)): % .8e % .8e" % (2 * sc[0] - prev[0], 2 * sc[1] - prev[1])
        print(line, flush=True); prev = sc
