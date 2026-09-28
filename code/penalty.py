"""Coercivity threshold of the Nitsche form on the exactly-divergence-free space Z_h.
A_mu is positive definite on ker B  <=>  K_r = A_mu + r B^T M^{-1} B is PD for r large.
PD is tested by the sign of the pivots of an unpivoted symmetric LU (inertia)."""
import sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, warnings
warnings.filterwarnings("ignore")
import svn

zero_f = lambda X, Y: np.zeros(np.shape(X) + (2,))


def build(N, L=2.5, layers=None):
    v, t, c, o = svn.make_mesh(N, L=L, layers=layers)
    return svn.Space(v, t, c, o)


def negpivots(S, mu, r=1e6, cache={}):
    key = id(S)
    if key not in cache:
        s0 = svn.assemble(S, 0.0, zero_f, None)
        s1 = svn.assemble(S, 1.0, zero_f, None)
        dD = svn.dofs(S.dnodes).ravel(); free = np.setdiff1d(np.arange(2 * S.nn), dD)
        A0 = s0["A"][free][:, free]; P = (s1["A"] - s0["A"])[free][:, free]   # penalty part per unit mu
        Bf = s0["B"][:, free]
        Pen = (Bf.T @ s0["Minv"] @ Bf).tocsc()
        cache.clear(); cache[key] = (A0, P, Pen)
    A0, P, Pen = cache[key]
    K = (A0 + mu * P + r * Pen).tocsc()
    K = (K + K.T) / 2
    lu = spl.splu(K, permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0,
                  options=dict(SymmetricMode=True))
    return int((lu.U.diagonal() < 0).sum())


def threshold(S, lo=1.0, hi=1e5, tol=1e-4, r=1e6):
    if negpivots(S, hi, r) > 0:
        return np.inf
    if negpivots(S, lo, r) == 0:
        return lo
    while hi / lo > 1 + tol:
        mid = np.sqrt(lo * hi)
        if negpivots(S, mid, r) > 0:
            lo = mid
        else:
            hi = mid
    return hi


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "validate"
    if mode == "validate":
        for N in (8, 12, 16, 20):
            S = build(N)
            for r in (1e5, 1e6, 1e7):
                mu = threshold(S, r=r)
                print(f"N={N} r={r:.0e}  mu*={mu:.5f}  gamma*={mu*S.hGamma/S.h:.4f}  rho={S.h/S.hGamma:.3f}", flush=True)
    elif mode == "rho":
        N = int(sys.argv[2]) if len(sys.argv) > 2 else 16
        for L in (2.5, 5.0, 10.0, 20.0, 40.0):
            m = max(N // 4, int(round(N / 4 * (L - 1) / 1.5)))
            S = build(N, L=L, layers=m)
            mu = threshold(S)
            print(f"N={N} L={L:5.1f} layers={m:3d} ntri={len(S.tris):5d} hG={S.hGamma:.4f} hOmega={S.h:.4f} "
                  f"rho={S.h/S.hGamma:7.3f}  mu*={mu:10.4f}  gamma*=mu*/rho={mu*S.hGamma/S.h:.4f}", flush=True)
    elif mode == "refine":
        for N in (8, 12, 16, 20, 24, 32, 48, 64):
            S = build(N)
            mu = threshold(S)
            print(f"N={N} ntri={len(S.tris)} rho={S.h/S.hGamma:.3f} mu*={mu:.4f} gamma*={mu*S.hGamma/S.h:.4f}", flush=True)
