import sympy as sp
import alfeld_patch as ap
from k3_trace import traces
s, t = ap.s_, ap.t_
ring = [(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)]
mac = [[(0,0,2),(0,0,1),ring[j],ring[(j+1)%4]] for j in range(4)]
mu = [1 - s - t, s, t]
Gs = [traces(mac, j, 3) for j in range(4)]   # same Y basis each time (deterministic construction)
def total(weights):
    M0 = []; M2 = sp.zeros(3, len(Gs[0]))
    for j in range(4):
        w12, w13, w23 = weights[j]
        Q = w12*mu[0]*mu[1] + w13*mu[0]*mu[2] + w23*mu[1]*mu[2]
        M0.append(sp.Matrix([[ap.ref_int(g[c]) for c in range(3)] for g in Gs[j]]).T)
        M2 += sp.Matrix([[ap.ref_int(Q * g[c]) for c in range(3)] for g in Gs[j]]).T
    C = sp.Matrix.vstack(*M0)
    K = sp.Matrix.hstack(*C.nullspace())
    return (M2 * K).rank(), (M2[0:2, :] * K).rank()
print('equal weights  : rank (3d, xy) =', total([(1,1,1)]*4))
print('weights (1,2,3): rank (3d, xy) =', total([(1,2,3)]*4))
print('face-1 only weighted:', total([(1,1,1),(0,0,0),(0,0,0),(0,0,0)]))
