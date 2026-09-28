"""Exact negative certificate: on a tall Alfeld tetrahedron the single-tet space W_3(T) has NO witness,
i.e. A-2B-cC is positive definite on W_3(T) (exact LDL^T over Q, all pivots > 0)."""
import sys, alfeld_tet as a
from fractions import Fraction as Fr
def posdef(Q):
    n = len(Q); Q = [row[:] for row in Q]
    for i in range(n):
        if Q[i][i] <= 0: return False
        for j in range(i + 1, n):
            f = Q[j][i] / Q[i][i]
            for l in range(i, n): Q[j][l] -= f * Q[i][l]
    return True
for nm, c in [(n, Fr(int(v))) for n, v in zip(sys.argv[1::2], sys.argv[2::2])] or [("tall(h=2)", Fr(10)), ("tall(h=3)", Fr(57))]:
    d, nr, A, B, C, dF = a.build(a.TETS[nm], 3)
    Q = [[A[i][j] - 2 * B[i][j] - c * C[i][j] for j in range(d)] for i in range(d)]
    Qc = [[A[i][j] - 2 * B[i][j] - (c + 1) * C[i][j] for j in range(d)] for i in range(d)]
    print(f"{nm}: dim={d}; A-2B-{c}C > 0 exactly: {posdef(Q)}  (control A-2B-{c+1}C > 0: {posdef(Qc)})", flush=True)
