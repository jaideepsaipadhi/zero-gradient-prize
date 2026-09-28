# witness definition (extracted from verify_witness.py); needs sympy as sp and symbols a1..b3,x,y
from sympy import Rational as Q
P = [sp.Matrix([1, 0]), sp.Matrix([a1, b1]), sp.Matrix([a2, b2]), sp.Matrix([a3, b3])]
cr = lambda u, v: u[0] * v[1] - u[1] * v[0]
D = {(i, j): cr(P[i], P[j]) for i in range(4) for j in range(4)}
X = sp.Matrix([x, y])
Z = sp.Matrix([0, 0])

def bary(A_, B_, C_):
    """barycentric coords (lam_A, lam_B, lam_C) of X in triangle A_,B_,C_ (rational in coords)"""
    area2 = cr(B_ - A_, C_ - A_)
    lA = cr(B_ - X, C_ - X) / area2
    lB = cr(C_ - X, A_ - X) / area2
    lC = cr(A_ - X, B_ - X) / area2
    return lA, lB, lC

A = D[1, 3] * D[2, 3]
A3 = D[0, 1] * D[0, 2]
# T1=(0,p0,p1), T2=(0,p1,p2), T3=(0,p2,p3); lambda0 at z, lambda1 at p_{j-1}, lambda2 at p_j
l = [bary(Z, P[j - 1], P[j]) for j in (1, 2, 3)]
L0, L1, L2 = l[0]
phi1 = A * D[0, 1]**2 * L0**2 * L2**2 * (L0 - 3 * L1)
L0, L1, L2 = l[2]
phi3 = A3 * D[2, 3]**2 * L0**2 * L1**2 * (L0 - 3 * L2)
L0, L1, L2 = l[1]
phi2 = L0**2 * (A * D[0, 1]**2 * L1**2 * (L0 + L2) + A3 * D[2, 3]**2 * L2**2 * (L0 + L1)
                 + 2 * A * D[0, 1] * D[0, 2] * L1 * L2 * (L0 + L1 + L2)
                 + A * D[0, 1] * (6 * D[1, 2] - 5 * D[0, 2] + 2 * D[0, 1]) * L1**2 * L2
                 + A3 * D[2, 3] * (6 * D[1, 2] - 5 * D[1, 3] + 2 * D[2, 3]) * L1 * L2**2)
PHIS = [phi1, phi2, phi3]
