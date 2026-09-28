"""Coercivity threshold of the Nitsche form on the exactly-divergence-free space Z_h, degree-general
(adapted from penalty.py to svn_k).  A_mu is PD on ker B  <=>  K_r = A_mu + r B^T M^{-1} B is PD for
r large; PD is tested by the sign of the pivots of an unpivoted symmetric LU (inertia).

  python3 penalty_k.py validate K MESH   -> checks the k=4 std numbers against penalty.py / r-independence
"""
import sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, warnings
warnings.filterwarnings("ignore")
import svn_k

zero_f = lambda X, Y: np.zeros(np.shape(X) + (2,))
_cache = {}


def negpivots(S, mu, r=1e6):
    key = id(S)
    if key not in _cache:
        s0 = svn_k.assemble(S, 0.0, zero_f)
        dD = svn_k.dofs(S.dnodes).ravel(); free = np.setdiff1d(np.arange(2 * S.nn), dD)
        A0 = s0["A"][free][:, free]; P = s0["Pen"][free][:, free]
        Bf = s0["B"][:, free]
        Pen = (Bf.T @ s0["Minv"] @ Bf).tocsc()
        _cache.clear(); _cache[key] = (A0, P, Pen)
    A0, P, Pen = _cache[key]
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
    k = int(sys.argv[2]); mesh = sys.argv[3]
    for N in (8, 12, 16):
        S = svn_k.build(N, k, mesh)
        for r in (1e5, 1e6, 1e7):
            mu = threshold(S, r=r)
            print(f"k={k} {mesh} N={N} r={r:.0e}  mu*={mu:.5f}  gamma*={mu*S.hGamma/S.h:.4f}  rho={S.h/S.hGamma:.3f}", flush=True)
