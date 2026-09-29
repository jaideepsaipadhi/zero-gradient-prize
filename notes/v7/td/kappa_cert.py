"""Certified reference constants on the Alfeld-split reference tetrahedron (exact rational arithmetic).

For k given on the command line:
  Y^_k   = div-free continuous piecewise P_k on the Alfeld split of T^, zero on dT^;
  Y^1_k  = {v in Y^_k : int_F g = 0},  g = d_z v on F^  (z = inward normal direction nu^);
  E^     = (int mu1mu2 g, int mu1mu3 g, int mu2mu3 g) in R^6           [quadratic edge moments, Prop. sh:joint]
  E^1    = (int mu2 g, int mu3 g) in R^4  (= the first-moment matrix M^1 = int g (x)(x - x_c) on Y^1)
Both on (Y^1, |grad .|_{L2(T^)}).  Right-inverse norms ||R^|| = lambda_min(E G^{-1} E^T)^{-1/2}.
The smallest eigenvalue is bracketed EXACTLY by LDL^T tests of E G^{-1} E^T - t I (Sylvester).
Also: the local Stokes right-inverse constant C^_loc = 1/sigma_min(div : V0(T^A) -> P_{k-1}^0) (float).
"""
import sys, time
from fractions import Fraction as Fr
import numpy as np
import refalfeld as ra
import ratlin as rl

k = int(sys.argv[1]) if len(sys.argv) > 1 else 4
t0 = time.time()
A = ra.Alfeld(k, 'all')
rows = A.div_rows()
rk, N = ra.nullspace(rows, A.nunk)
print(f"k={k}: domain pts={A.nd} unknowns={A.nunk} div rows={len(rows)} rank(div)={rk} dim Y={len(N)}  ({time.time()-t0:.1f}s)", flush=True)
G = A.grad_gram()
GY = ra.congruence(G, N)
W = [('m', (0, 0, 0)), ('m1', (1, 0, 0)), ('m2', (0, 1, 0)), ('m3', (0, 0, 1)),
     ('E12', (1, 1, 0)), ('E13', (1, 0, 1)), ('E23', (0, 1, 1))]
FR = ra.face_moment_rows(A, W)
L = {key: ra.apply_rows([r], N)[0] for key, r in FR.items()}
mean = [L[('m', 0)], L[('m', 1)]]
print("rank of mean on Y:", rl.rank(mean))
K1 = rl.nullspace(mean)                      # coefficient vectors (in Y basis) spanning Y^1
d1 = len(K1)
G1 = rl.matmul(rl.matmul(K1, GY), rl.transpose(K1))
restr = lambda key: [sum(a * b for a, b in zip(L[key], v)) for v in K1]
Eq = [restr((n, c)) for n in ('E12', 'E13', 'E23') for c in (0, 1)]
E1 = [restr((n, c)) for n in ('m2', 'm3') for c in (0, 1)]
E1all = [restr((n, c)) for n in ('m1', 'm2', 'm3') for c in (0, 1)]
print(f"dim Y^1={d1}; rank quadratic moments on Y^1 = {rl.rank(Eq)}; rank first moments on Y^1 = {rl.rank(E1)} "
      f"(with mu1 too: {rl.rank(E1all)})", flush=True)


def rinv(E, name):
    if rl.rank(E) < len(E):
        print(f"  {name}: not onto (rank {rl.rank(E)} < {len(E)})")
        Ef = np.array([[float(x) for x in r] for r in E]); Gf = np.array([[float(x) for x in r] for r in G1])
        M = Ef @ np.linalg.solve(Gf, Ef.T)
        print(f"  {name}: eigenvalues of E G^-1 E^T (float): {np.linalg.eigvalsh(M)}")
        return None
    X = rl.solve(G1, rl.transpose(E))            # G^{-1} E^T
    M = rl.matmul(E, X)                          # E G^{-1} E^T (exact)
    Mf = np.array([[float(x) for x in r] for r in M])
    ev = np.linalg.eigvalsh(Mf)
    lo, hi = Fr(ev[0] * (1 - 1e-6)).limit_denominator(10 ** 15), Fr(ev[0] * (1 + 1e-6)).limit_denominator(10 ** 15)
    lo, hi = rl.lmin_bracket(M, lo, hi, iters=30)
    Rlo, Rhi = float(hi) ** -0.5, float(lo) ** -0.5
    print(f"  {name}: singular values (float) {np.sqrt(ev[::-1])}")
    print(f"  {name}: CERTIFIED lambda_min(E G^-1 E^T) in [{float(lo):.12e}, {float(hi):.12e}]  (exact LDL^T)")
    print(f"  {name}: CERTIFIED ||R^|| in [{Rlo:.9f}, {Rhi:.9f}]", flush=True)
    return lo, hi


rinv(Eq, "quadratic edge moments E^ (kappa_0 of Thm sh:alfeld)")
rinv(E1, "first moments E^1 (leak witness)")

# ---- reference local Stokes constant (float): min over q in P_{k-1}^0 of sup (div v,q)/(|grad v| |q|)
Gf = np.zeros((A.nunk, A.nunk))
for u, row in G.items():
    for w, c in row.items():
        Gf[u, w] = float(c)
Df = np.zeros((len(rows), A.nunk))
for i, r in enumerate(rows):
    for u, c in r.items():
        Df[i, u] = float(c)
# pressure mass in the degree-(k-1) Bernstein basis, block diagonal over the 4 sub-tets
nb = len(A.bk1)
Mp = np.zeros((4 * nb, 4 * nb))
for s in range(4):
    M = A.mass_k1(s)
    for i, a in enumerate(A.bk1):
        for j, b in enumerate(A.bk1):
            Mp[s * nb + i, s * nb + j] = float(M[(a, b)])
# (div v, q) = sum_s int (k sum_beta B_beta divcoef_beta) q ; q in Bernstein basis => Bmat = Mp @ (k*Df)
Bm = Mp @ (k * Df)
S = Bm @ np.linalg.solve(Gf, Bm.T)
import scipy.linalg as sl
ev = sl.eigh(S, Mp, eigvals_only=True)
print(f"local Stokes (V0(T^A) x P_(k-1)^disc): smallest gen. eigenvalues {ev[:3]} -> beta^_loc = {np.sqrt(ev[1]):.6f}, "
      f"C^_loc = 1/beta^_loc = {1/np.sqrt(ev[1]):.4f}  (float)", flush=True)
print(f"done ({time.time()-t0:.1f}s)")
