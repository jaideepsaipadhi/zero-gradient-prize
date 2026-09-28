"""k=3 Alfeld: structure of the face trace g = d_nu v on a Gamma-face, over Y(patch)."""
import sympy as sp, sys
from sympy.polys.matrices import DomainMatrix
from sympy import QQ
import alfeld_patch as ap
x, y, z = ap.X
def traces(macro, face_tet, k):
    # replicate analyse() internals to get g on face (tet face_tet, local (1,2,3)) for a Y basis
    import itertools
    M = [[sp.Matrix(v) for v in T] for T in macro]
    cnt = {}
    fk = lambda pts: frozenset(tuple(p) for p in pts)
    for T in M:
        for i in range(4):
            f = fk([T[j] for j in range(4) if j != i]); cnt[f] = cnt.get(f, 0) + 1
    bfaces = [f for f, c in cnt.items() if c == 1]
    subs = []
    for ti, T in enumerate(M):
        B = (T[0] + T[1] + T[2] + T[3]) / 4
        for i in range(4): subs.append(([B] + [T[j] for j in range(4) if j != i], ti, i))
    def onb(pt):
        for f in bfaces:
            P = [sp.Matrix(p) for p in f]; n = (P[1] - P[0]).cross(P[2] - P[0])
            if n.dot(pt - P[0]) != 0: continue
            A = sp.Matrix.hstack(P[1] - P[0], P[2] - P[0]); a1, a2 = (A.T * A).inv() * A.T * (pt - P[0])
            if a1 >= 0 and a2 >= 0 and a1 + a2 <= 1: return True
        return False
    gidx = {}; polys = []; bc = {}
    for (K, ti, i) in subs:
        lam = ap.bary(K); bl = []
        for a in ap.mi(k):
            pt = sum((sp.Rational(a[j], k) * K[j] for j in range(4)), sp.zeros(3, 1)); key = tuple(pt)
            if key not in bc: bc[key] = onb(pt)
            if bc[key]: continue
            if key not in gidx: gidx[key] = len(gidx)
            c = sp.factorial(k) / sp.prod([sp.factorial(ai) for ai in a])
            bl.append((gidx[key], sp.Poly(sp.expand(c * sp.prod([lam[j] ** a[j] for j in range(4)])), *ap.X)))
        polys.append(bl)
    nunk = 3 * len(gidx); rows = []
    for bl in polys:
        co = {}
        for gi, p in bl:
            for c in range(3):
                for mon, cv in p.diff(ap.X[c]).terms():
                    co.setdefault(mon, {}); u = 3 * gi + c; co[mon][u] = co[mon].get(u, 0) + cv
        for mon, d in co.items():
            r = [QQ(0)] * nunk
            for u, cv in d.items(): r[u] = ap.q(cv)
            rows.append(r)
    N = DomainMatrix(rows, (len(rows), nunk), QQ).nullspace().to_Matrix()
    T = M[face_tet]; si = [i for i, (K, t, o) in enumerate(subs) if t == face_tet and o == 0][0]
    P1, P2, P3 = T[1], T[2], T[3]; NF = (P2 - P1).cross(P3 - P1)
    if NF.dot(T[0] - P1) < 0: NF = -NF
    s, t = ap.s_, ap.t_
    sub = {ap.X[c]: P1[c] + s * (P2 - P1)[c] + t * (P3 - P1)[c] for c in range(3)}
    gb = []
    for gi, p in polys[si]:
        gb.append((gi, sp.expand(sum(p.diff(ap.X[c]).as_expr() * NF[c] for c in range(3)).subs(sub))))
    G = []
    for r in range(N.shape[0]):
        g = [sp.expand(sum(N[r, 3 * gi + c] * e for gi, e in gb)) for c in range(3)]
        G.append(g)
    return G
if __name__ == '__main__':
    k = 3
    cases = {
     'single tet (general)': ([[(0, 0, 3), (0, 0, 0), (5, 1, 0), (1, 4, 0)]], 0),
     'octa fan': ([[(0,0,2),(0,0,1),r1,r2] for r1, r2 in zip([(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)],[(0,1,0),(-1,0,0),(0,-1,0),(1,0,0)])], 0),
    }
    s, t = ap.s_, ap.t_
    for name, (mac, ft) in cases.items():
        G = traces(mac, ft, k)
        # coefficient vectors of g (3 comps, quadratic in s,t)
        mons = [s**a * t**b for a in range(3) for b in range(3) if a + b <= 2]
        Mat = sp.Matrix([[sp.Poly(g[c], s, t).coeff_monomial(m) for c in range(3) for m in mons] for g in G])
        print(name, ': dim Y =', len(G), ' rank of trace map g =', Mat.rank())
        vals = {'vertex1': {s: 0, t: 0}, 'vertex2': {s: 1, t: 0}, 'vertex3': {s: 0, t: 1}, 'mid12': {s: sp.Rational(1,2), t: 0}, 'centroid': {s: sp.Rational(1,3), t: sp.Rational(1,3)}}
        for vn, vv in vals.items():
            V = sp.Matrix([[g[c].subs(vv) for c in range(3)] for g in G])
            print('   rank of g at', vn, '=', V.rank())
        # tangential curl / divergence structure: express g in tangential frame and print a basis of traces
        print('   basis of traces (components along NF-frame omitted):')
        B = Mat.T.columnspace()
        print('   ', len(B))
