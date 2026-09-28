import sympy as sp
import alfeld_patch as ap
from k3_trace import traces
s, t = ap.s_, ap.t_
mac = [[(0,0,2),(0,0,1),r1,r2] for r1, r2 in zip([(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)],[(0,1,0),(-1,0,0),(0,-1,0),(1,0,0)])]
G = traces(mac, 0, 3)
# F1 = (e3, e1, e2): P1=e3, P2=e1, P3=e2 ; frame: a1 = P2-P1, a2 = P3-P1 ; express g = c1 a1 + c2 a2
P1, P2, P3 = sp.Matrix([0,0,1]), sp.Matrix([1,0,0]), sp.Matrix([0,1,0])
A = sp.Matrix.hstack(P2 - P1, P3 - P1)
pinv = (A.T * A).inv() * A.T
mu = [1 - s - t, s, t]
# Bernstein quadratic basis on F
bern = {}
for a in [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)]:
    c = sp.factorial(2) / sp.prod([sp.factorial(ai) for ai in a])
    bern[a] = sp.expand(c * mu[0]**a[0] * mu[1]**a[1] * mu[2]**a[2])
keys = list(bern)
# solve for Bernstein coefficients of each tangential component
mons = [s**a * t**b for a in range(3) for b in range(3) if a + b <= 2]
Bmat = sp.Matrix([[sp.Poly(bern[kk], s, t).coeff_monomial(m) for kk in keys] for m in mons])
rows = []
for g in G:
    gv = sp.Matrix(g)
    comps = pinv * gv          # 2 components in (a1,a2) frame
    # check tangential
    assert sp.simplify((A * comps - gv).applyfunc(sp.expand)) == sp.zeros(3, 1)
    row = []
    for c in range(2):
        rhs = sp.Matrix([sp.Poly(sp.expand(comps[c]), s, t).coeff_monomial(m) for m in mons])
        row += list(Bmat.solve(rhs))
    rows.append(row)
T = sp.Matrix(rows)
print('trace space rank', T.rank())
print('labels:', [(c, kk) for c in ('a1','a2') for kk in keys])
cons = T.nullspace()
print('constraints (annihilators of the trace space), in Bernstein coefficients (a1-comp then a2-comp):')
for v in cons:
    print(list(v.T))
