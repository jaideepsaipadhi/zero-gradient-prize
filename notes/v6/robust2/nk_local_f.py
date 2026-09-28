"""
Float version of nk_local.py (fast; any k).  For T = conv{(0,0),(1,0),(xi,zeta)}, e = [0,1]x{0}:
  Lam_T : Hom_k -> P_e,  H |-> Pi_{P_e} ((I - pi_T) H)|_e    (P_e = traces of the local div-free fields X_e)
reports singular values, the kernel, |Lam_T(x^k)|, the component of y^k (pure normal) and the L_k pairing.
Shape sweep: xi in [-0.5, 1.5], zeta in [0.3, 1.5] (covers the production meshes' boundary triangles).
usage: python nk_local_f.py k
"""
import sys
import numpy as np
from numpy.polynomial import legendre as Lg


def quad_tri(n):
    g, w = Lg.leggauss(n); g = (g + 1) / 2; w = w / 2
    U, V = np.meshgrid(g, g, indexing="ij"); WU, WV = np.meshgrid(w, w, indexing="ij")
    return U.ravel(), (V * (1 - U)).ravel(), (WU * WV * (1 - U)).ravel()


def lam_T(k, xi, ze, nq=None):
    nq = nq or k + 6
    s, r, w = quad_tri(nq)
    X = s + r * xi; Y = r * ze; W = w * ze
    Pk1 = [(a, b) for a in range(k) for b in range(k - a)]
    B = np.stack([X**a * Y**b for (a, b) in Pk1], 1)
    G = B.T @ (W[:, None] * B)
    gl, gw = Lg.leggauss(k + 6); tt = (gl + 1) / 2; tw = gw / 2
    Be = np.stack([tt**a * 0.0**b if b else tt**a for (a, b) in Pk1], 1)
    R = []
    for j in range(k + 1):
        H = X**(k - j) * Y**j
        c = np.linalg.solve(G, B.T @ (W * H))
        He = tt**(k - j) if j == 0 else 0 * tt
        R.append(He - Be @ c)                                   # (I - pi_T)H on e at GL nodes
    R = np.array(R).T                                            # (nq, k+1)
    # P_e basis: d/dt[t^2(1-t)^2 t^j]
    Pe = np.stack([(2 * tt * (1 - tt)**2 - 2 * tt**2 * (1 - tt)) * tt**j + tt**2 * (1 - tt)**2 * j * tt**max(j - 1, 0)
                   for j in range(k - 2)], 1)
    Ge = Pe.T @ (tw[:, None] * Pe)
    Lam = Pe.T @ (tw[:, None] * R)
    ev, V = np.linalg.eigh(Ge); M = (V @ np.diag(ev**-0.5) @ V.T) @ Lam   # orthonormalised
    Lk = Lg.legval(2 * tt - 1, [0] * k + [1])
    lk = Lk @ (tw[:, None] * R)
    flux = tw @ R
    return M, lk, flux


def main():
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    from math import factorial
    kap = factorial(k)**2 / factorial(2 * k + 1)
    print(f"k={k}  kappa_k = k!^2/(2k+1)! = {kap:.6e}")
    worst_tan = np.inf; worst_normal = 0; worst_lk = 0; smin = np.inf
    for xi in np.linspace(-0.5, 1.5, 21):
        for ze in np.linspace(0.3, 1.5, 13):
            M, lk, flux = lam_T(k, xi, ze)
            sv = np.linalg.svd(M, compute_uv=False)
            worst_tan = min(worst_tan, np.linalg.norm(M[:, 0]))
            worst_normal = max(worst_normal, np.linalg.norm(M[:, k]) / max(np.linalg.norm(M), 1e-300))
            worst_lk = max(worst_lk, abs(lk[0] - kap) / kap + np.abs(lk[1:]).max() / kap)
            smin = min(smin, sv[-1])
    print(f"  min over shapes |Lam_T(x^k)|              = {worst_tan:.3e}   (tangential k-th derivative always visible to X_e)")
    print(f"  max over shapes |Lam_T(y^k)|/|Lam_T|       = {worst_normal:.3e}  (pure normal k-th derivative invisible to X_e)")
    print(f"  max rel. dev. of <(I-pi)H,L_k> from kappa*delta_(j,0) = {worst_lk:.3e}")
    print(f"  min over shapes smallest sing. value of Lam_T (rank k-2 map) = {smin:.3e}")
    M, lk, flux = lam_T(k, 0.5, np.sqrt(3) / 2)
    _, s, Vt = np.linalg.svd(M)
    print("  equilateral: sing. values", np.array2string(s, precision=3), " kernel dim", k + 1 - np.sum(s > 1e-10 * s[0]))
    print("  equilateral: flux functional int_e (I-pi)x^(k-j)y^j =", np.array2string(flux, precision=3))


if __name__ == "__main__":
    main()
