"""k=3: can larger Alfeld patches supply both tangential moment directions?  Exact rational."""
import sympy as sp, sys
from alfeld_patch import analyse
R = sp.Rational

def ring(m, rad=R(1), curv=R(1, 4), jitter=None):
    # rational approximations of a ring of m points around z=0 on the paraboloid zz = -curv*(x^2+y^2)
    pts = []
    for j in range(m):
        th = 2 * sp.pi * j / m
        c = sp.nsimplify(round(float(sp.cos(th)), 3)); s = sp.nsimplify(round(float(sp.sin(th)), 3))
        if jitter: c += jitter[j][0]; s += jitter[j][1]
        X, Y = rad * c, rad * s
        pts.append((X, Y, -curv * (X * X + Y * Y)))
    return pts

def fan(m, apex_h=R(1), **kw):
    a = ring(m, **kw); zz = (0, 0, 0); qq = (0, 0, apex_h)
    macro = [[qq, zz, a[j], a[(j + 1) % m]] for j in range(m)]
    gf = [(j, (1, 2, 3)) for j in range(m)]
    return macro, gf

k = int(sys.argv[1]) if len(sys.argv) > 1 else 3
tang = [(1, 0, 0), (0, 1, 0)]
for m in (4, 5, 6):
    print(f'vertex fan, m={m} boundary faces (regular ring):')
    analyse(*fan(m), k, tang)
print('vertex fan, m=5, irregular ring:')
jit = [(R(1,10),0),(0,R(1,7)),(R(-1,9),R(1,11)),(0,0),(R(1,13),R(-1,8))]
analyse(*fan(5, jitter=jit), k, tang)
