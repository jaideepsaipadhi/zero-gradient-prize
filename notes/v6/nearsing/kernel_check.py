"""Check of the explicit kernel-vector asymptotics used in the proof of Theorem N2 (m = 3, one needle).
Lines (angles of the rays mod pi, measured from e_first): Family II: (0, g, g+eps, phi); Family I: (0, g, phi-eps, phi).
A_z is spanned by x = kernel of [P_1 P_2 P_3 P_4] (P_j = n_j n_j^T); H_1 = x_1 P_1, H_2 = H_1 + x_2 P_2, H_3 = -x_4 P_4.
Normalise x_1 = 1 and print:  II: |x_2| * eps/phi  (-> const = 1),  |x_4 + 1| (-> O(phi));
                               I : |(-x_4) - 1| * eps/phi  (-> 1),  |x_2| (-> O(phi)).
"""
import numpy as np
from local_ns import P, vec, unit
def kern(psis):
    M = np.array([vec(P(unit(p + np.pi / 2))) for p in psis]).T
    x = np.linalg.svd(M)[2][-1]; return x / x[0]
g = np.pi / 3
print(" phi     eps     II:|x2|eps/phi  II:|x4+1|/phi   I:|-x4-1|eps/phi  I:|x2|/phi")
for phi in (0.04, 0.01):
    for eps in (0.3, 0.1, 0.03, 0.01, 0.003):
        x = kern([0, g, g + eps, phi]); y = kern([0, g, phi - eps, phi])
        print(f" {phi:.3f}  {eps:.3f}   {abs(x[1]) * eps / phi:8.4f}       {abs(x[3] + 1) / phi:8.4f}      "
              f"{abs(-y[3] - 1) * eps / phi:8.4f}        {abs(y[1]) / phi:8.4f}")
