"""Small exact linear algebra over Fractions: solve, inverse, LDL^T positivity test, rank."""
from fractions import Fraction as Fr


def solve(A, B):
    """A X = B, A square invertible (lists of Fractions), B matrix (list of rows)."""
    n = len(A); m = len(B[0])
    M = [list(A[i]) + list(B[i]) for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [row[n:] for row in M]


def matmul(A, B):
    Bt = list(zip(*B))
    return [[sum(a * b for a, b in zip(r, c)) for c in Bt] for r in A]


def transpose(A):
    return [list(r) for r in zip(*A)]


def is_pd(S):
    """exact test that the symmetric matrix S is positive definite (LDL^T without pivoting)."""
    n = len(S); A = [list(r) for r in S]
    for c in range(n):
        if A[c][c] <= 0:
            return False
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            for j in range(c, n):
                A[r][j] -= f * A[c][j]
    return True


def rank(A):
    M = [list(r) for r in A]; rk = 0; n = len(M); m = len(M[0]) if n else 0
    for c in range(m):
        p = next((r for r in range(rk, n) if M[r][c] != 0), None)
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        for r in range(n):
            if r != rk and M[r][c] != 0:
                f = M[r][c] / M[rk][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[rk])]
        rk += 1
    return rk


def nullspace(A):
    """basis of {x: A x = 0} (rows of A), exact."""
    M = [list(r) for r in A]; n = len(M); m = len(M[0]); piv = []; rk = 0
    for c in range(m):
        p = next((r for r in range(rk, n) if M[r][c] != 0), None)
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        pv = M[rk][c]; M[rk] = [x / pv for x in M[rk]]
        for r in range(n):
            if r != rk and M[r][c] != 0:
                f = M[r][c]; M[r] = [x - f * y for x, y in zip(M[r], M[rk])]
        piv.append(c); rk += 1
    free = [c for c in range(m) if c not in piv]
    basis = []
    for f in free:
        x = [Fr(0)] * m; x[f] = Fr(1)
        for i, c in enumerate(piv):
            x[c] = -M[i][f]
        basis.append(x)
    return basis


def lmin_bracket(S, lo, hi, iters=40):
    """exact bisection bracket for the smallest eigenvalue of symmetric S: returns (a,b) with
    S - a I > 0 and S - b I not > 0.  Requires S - lo I > 0 and not (S - hi I > 0)."""
    n = len(S)
    shift = lambda t: [[S[i][j] - (t if i == j else 0) for j in range(n)] for i in range(n)]
    assert is_pd(shift(lo)) and not is_pd(shift(hi))
    for _ in range(iters):
        mid = (lo + hi) / 2
        mid = Fr(mid).limit_denominator(10 ** 12)
        if is_pd(shift(mid)):
            lo = mid
        else:
            hi = mid
    return lo, hi
