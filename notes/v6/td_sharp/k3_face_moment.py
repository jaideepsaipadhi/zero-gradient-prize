import sympy as sp
import alfeld_patch as ap
from k3_trace import traces
s, t = ap.s_, ap.t_
mac = [[(0,0,2),(0,0,1),r1,r2] for r1, r2 in zip([(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)],[(0,1,0),(-1,0,0),(0,-1,0),(1,0,0)])]
G = traces(mac, 0, 3)
mu = [1 - s - t, s, t]
Q = mu[0]*mu[1] + mu[0]*mu[2] + mu[1]*mu[2]   # equal weights (equilateral)
M0 = sp.Matrix([[ap.ref_int(g[c]) for c in range(3)] for g in G]).T
M2 = sp.Matrix([[ap.ref_int(Q * g[c]) for c in range(3)] for g in G]).T
print('rank mean on F1 =', M0.rank())
K = M0.nullspace()
Km = sp.Matrix.hstack(*K)
print('dim {v: mean_F1 = 0} =', Km.shape[1], ' rank of F1-moment there =', (M2 * Km).rank())
