"""lambda_T(k) along apex height (apex above the circumcentre-ish point (1/2,7/24) of the near-equilateral base)."""
import sys, alfeld_tet as a
from fractions import Fraction as Fr
k = int(sys.argv[1]); hs = [Fr(s) for s in sys.argv[2:]]
for h in hs:
    T = [(Fr(1, 2), Fr(7, 24), h), (0, 0, 0), (1, 0, 0), (Fr(1, 2), Fr(7, 8), 0)]
    d, nr, A, B, C, dF = a.build(T, k)
    L, vec = a.lam_float(A, B, C)
    print(f"k={k} h={h} ({float(h):.4f}) dim={d} lambda_T={L*dF:.5f}", flush=True)
