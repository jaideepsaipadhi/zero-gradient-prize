"""Where is the pure tangential tensor x^k invisible to the local fields X_e?  |Lam_T(x^k)| over (xi, zeta)."""
import sys, numpy as np
from nk_local_f import lam_T
k = int(sys.argv[1]) if len(sys.argv) > 1 else 4
for ze in [0.3, 0.5, 0.866, 1.2, 1.5]:
    row = []
    for xi in np.linspace(-0.5, 1.5, 9):
        M, lk, fl = lam_T(k, xi, ze)
        row.append(np.linalg.norm(M[:, 0]))
    print(f"zeta={ze:5.3f}: " + " ".join(f"{v:8.2e}" for v in row))
print("xi grid:", np.linspace(-0.5, 1.5, 9))
