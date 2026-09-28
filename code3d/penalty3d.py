"""Coercivity threshold of the 3D Nitsche form on the exactly-divergence-free space Z_h (box nodes
Dirichlet, Gamma_h Nitsche, flux 'grad'), as code/penalty_k.py:  A_mu PD on ker B  <=>  K_r = A_mu + r B^T M^-1 B
PD for r large.  Inertia from MKL PARDISO symmetric-indefinite LDL^T (iparm 22/23).
  python3 penalty3d.py K N [N ...]"""
import os, sys, numpy as np, scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svn3d

zero = lambda a, b, c: np.zeros(np.shape(a) + (3,))


def pieces(S):
    s0 = svn3d.assemble(S, 0.0, zero, [(S.fGamma, None, "grad")], [])
    s1 = svn3d.assemble(S, 1.0, zero, [(S.fGamma, None, "grad")], [])
    dD = (3 * S.boxnodes[:, None] + np.arange(3)).ravel()
    free = np.setdiff1d(np.arange(3 * S.nn), dD)
    A0 = s0["A"][free][:, free]; P = (s1["A"] - s0["A"])[free][:, free]
    Bf = s0["B"][:, free]
    Pen = (Bf.T @ s0["Minv"] @ Bf)
    return A0, P, Pen


def inertia(K):
    import pypardiso
    K = sp.triu((K + K.T) / 2, format="csr"); K.sort_indices()
    slv = pypardiso.PyPardisoSolver(mtype=-2)
    slv.factorize(K)
    npos, nneg = int(slv.get_iparm(22)), int(slv.get_iparm(23))
    slv.free_memory(everything=True)
    return npos, nneg


def threshold(S, lo=1.0, hi=3e3, tol=1e-3, r=1e6, pc=None):
    A0, P, Pen = pc or pieces(S)
    neg = lambda mu: inertia(A0 + mu * P + r * Pen)[1]
    if neg(hi) > 0:
        return np.inf
    if neg(lo) == 0:
        return lo
    while hi / lo > 1 + tol:
        mid = np.sqrt(lo * hi)
        if neg(mid) > 0:
            lo = mid
        else:
            hi = mid
    return hi


if __name__ == "__main__":
    k = int(sys.argv[1])
    for N in map(int, sys.argv[2:]):
        S = svn3d.build("sph", N, k, L=2.0)
        pc = pieces(S)
        for r in (1e5, 1e6, 1e7):
            mu = threshold(S, r=r, pc=pc)
            print(f"k={k} N={N} m={S.info['m']} r={r:.0e} muZ*={mu:.3f} gamma*={mu*S.hGamma/S.h:.4f} "
                  f"h/hG={S.h/S.hGamma:.3f} localV*={svn3d.local_threshold(S, S.fGamma).max():.2f}", flush=True)
