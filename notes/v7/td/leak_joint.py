"""Leak witness with vanishing quadratic edge moments (Alfeld and barycentric WF, exact).
On Y^1 (zero face mean(s) of g), the map v -> (first moments E^1 in R^4, quadratic edge moments E^ in R^6) of g = d_z v.
If it has rank 10, the leak witness can be chosen with E^(v)=0; its right inverse restricted to {E^=0} is certified here:
   ||R^_1,0|| := sup{ |grad v| : E^1(v)=m, E^(v)=0, v minimal } over |m|=1   (exact LDL^T bracket).
usage: python3 leak_joint.py alfeld|wf k"""
import sys
from fractions import Fraction as Fr
import numpy as np
import refalfeld as ra, ratlin as rl
kind = sys.argv[1]; k = int(sys.argv[2])
if kind == 'alfeld':
    A = ra.Alfeld(k, 'all'); rk, N = ra.nullspace(A.div_rows(), A.nunk)
    W = [('m', (0, 0, 0)), ('m2', (0, 1, 0)), ('m3', (0, 0, 1)), ('E12', (1, 1, 0)), ('E13', (1, 0, 1)), ('E23', (0, 1, 1))]
    FR = ra.face_moment_rows(A, W)
    L = {key: ra.apply_rows([r], N)[0] for key, r in FR.items()}
    get = lambda n, c: L[(n, c)]
    means = [get('m', 0), get('m', 1)]
else:
    sys.argv = ['wf_witness.py', str(k)]
    src = open('wf_witness.py').read().split("k = int(sys.argv[1])")[0]
    exec(src)
    A = WF(k); rk, N = ra.nullspace(A.div_rows(), A.nunk)
    Wd = {'m': {(0, 0, 0): Fr(1)}, 'm2': {(0, 1, 0): Fr(1)}, 'm3': {(0, 0, 1): Fr(1)},
          'E12': {(1, 1, 0): Fr(1)}, 'E13': {(1, 0, 1): Fr(1)}, 'E23': {(0, 1, 1): Fr(1)}}
    FR = A.face_rows(Wd); subs = [s for s, _ in A.fsub]
    L = {key: ra.apply_rows([r], N)[0] for key, r in FR.items()}
    get = lambda n, c: [sum(L[(n, s, c)][j] for s in subs) for j in range(len(N))]
    means = [L[('m', s, c)] for s in subs for c in (0, 1)]
K1 = rl.nullspace(means)
restr = lambda vec: [sum(a * b for a, b in zip(vec, v)) for v in K1]
E1 = [restr(get(n, c)) for n in ('m2', 'm3') for c in (0, 1)]
Eq = [restr(get(n, c)) for n in ('E12', 'E13', 'E23') for c in (0, 1)]
print(f"{kind} k={k}: dim Y={len(N)}, dim Y^1={len(K1)}, rank E^1={rl.rank(E1)}, rank E^={rl.rank(Eq)}, joint rank={rl.rank(E1+Eq)} (of 10)")
if rl.rank(E1 + Eq) == 10:
    K2 = rl.nullspace(Eq)                                   # coefficients (in K1 basis) with E^ = 0
    G = A.grad_gram(); GY = ra.congruence(G, N)
    G1 = rl.matmul(rl.matmul(K1, GY), rl.transpose(K1))
    G2 = rl.matmul(rl.matmul(K2, G1), rl.transpose(K2))
    E12 = [[sum(a * b for a, b in zip(r, v)) for v in K2] for r in E1]
    X = rl.solve(G2, rl.transpose(E12)); M = rl.matmul(E12, X)
    ev = np.linalg.eigvalsh(np.array([[float(x) for x in r] for r in M]))
    lo, hi = Fr(ev[0] * (1 - 1e-6)).limit_denominator(10 ** 15), Fr(ev[0] * (1 + 1e-6)).limit_denominator(10 ** 15)
    lo, hi = rl.lmin_bracket(M, lo, hi, iters=30)
    print(f"  CERTIFIED ||R^_1,0|| (first moments prescribed, quadratic edge moments zero) in [{float(hi)**-0.5:.9f}, {float(lo)**-0.5:.9f}]")
