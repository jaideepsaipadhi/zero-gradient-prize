"""Referee: global threshold on the polar column mesh (H=5, N=128, L=2.5) -- sensitivity to the divergence penalty r,
comparison Z_h vs full V_h (r=0), fresh-process (no cache) value, and an eigen-check at mu slightly below/above mu*."""
import sys, os, numpy as np, scipy.sparse.linalg as spl
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "penalty")); sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
import global_mesh as gm, svn, penalty

H = float(sys.argv[1]) if len(sys.argv) > 1 else 5.0
V, T, circ, outer, layer = gm.build_mesh(128, H, 2, 2.5)
S = svn.Space(V, T, circ, outer)
print(f"H={H} ntri={len(T)} h={S.h:.4f} hG={S.hGamma:.5f}")
for r in (1e4, 1e5, 1e6, 1e7, 1e8):
    penalty.negpivots.__defaults__[-1].clear()
    mu = penalty.threshold(S, lo=1e-2, hi=1e6, tol=1e-5, r=r)
    print(f"  r={r:.0e}: mu*={mu:.5f} gamma*={mu*S.hGamma/S.h:.5f}", flush=True)
penalty.negpivots.__defaults__[-1].clear()
mu0 = penalty.threshold(S, lo=1e-2, hi=1e6, tol=1e-5, r=0.0)
print(f"  r=0 (full V_h, NOT Z_h): mu*={mu0:.5f} gamma*={mu0*S.hGamma/S.h:.5f}")

# independent: null-space basis of B (dense SVD) and smallest eigenvalue of A_mu restricted to ker B
s0 = svn.assemble(S, 0.0, gm.penalty.zero_f, None); s1 = svn.assemble(S, 1.0, gm.penalty.zero_f, None)
dD = svn.dofs(S.dnodes).ravel(); free = np.setdiff1d(np.arange(2 * S.nn), dD)
A0 = s0["A"][free][:, free].toarray(); P = (s1["A"] - s0["A"])[free][:, free].toarray(); B = s0["B"][:, free].toarray()
print("  dense sizes", A0.shape, B.shape, flush=True)
U, sv, Vt = np.linalg.svd(B, full_matrices=True)
rank = int((sv > 1e-10 * sv[0]).sum()); Z = Vt[rank:].T
print(f"  rank B={rank}, dim Z_h={Z.shape[1]}", flush=True)
Az = Z.T @ A0 @ Z; Pz = Z.T @ P @ Z
mu_star = None
import scipy.linalg as sl
for mu in (0.95, 0.999, 1.001, 1.05):
    m = mu * 65.936 * (H == 5.0) + mu * (1 - (H == 5.0))
    ev = sl.eigvalsh((Az + m * Pz + (Az + m * Pz).T) / 2)
    print(f"  mu={m:.4f}: min eig of N_h on Z_h = {ev[0]:.3e}", flush=True)
