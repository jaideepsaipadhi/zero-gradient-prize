"""Check: ker Lam_T = span{ (nu_a.x)^k, (nu_b.x)^k, (nu_c.x)^k }, nu_w = normal of the edge opposite vertex w
(k-th powers of the three barycentric coordinates, modulo P_{k-1}).  Reports, over a shape sweep,
max |Lam_T(ell_w^k)| / (|Lam_T| |ell_w^k|coeff)  and  rank of Lam_T."""
import sys, numpy as np
from math import comb
from nk_local_f import lam_T
k = int(sys.argv[1]) if len(sys.argv) > 1 else 4
worst = 0; ranks = set()
for xi in np.linspace(-0.6, 1.6, 23):
    for ze in np.linspace(0.25, 1.6, 10):
        M, _, _ = lam_T(k, xi, ze)
        V = np.array([[0, 0], [1, 0], [xi, ze]])
        s = np.linalg.svd(M, compute_uv=False); ranks.add(int(np.sum(s > 1e-9 * s[0])))
        for w in range(3):
            P, Q = V[(w + 1) % 3], V[(w + 2) % 3]
            d = Q - P; nu = np.array([d[1], -d[0]]) / np.linalg.norm(d)
            # (nu.x)^k = sum_j C(k,j) nu0^(k-j) nu1^j x^(k-j) y^j
            c = np.array([comb(k, j) * nu[0]**(k - j) * nu[1]**j for j in range(k + 1)])
            worst = max(worst, np.linalg.norm(M @ c) / (np.linalg.norm(M, 2) * np.linalg.norm(c)))
print(f"k={k}: max relative |Lam_T((nu_w.x)^k)| = {worst:.2e}; ranks of Lam_T seen: {sorted(ranks)} (dim Hom_k = {k+1}, dim P_e = {k-2})")
