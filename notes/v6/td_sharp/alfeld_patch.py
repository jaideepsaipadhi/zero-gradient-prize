"""Exact rational computation of divergence-free Alfeld fields on a PATCH of macro-tetrahedra.

patch: list of macro-tets (4 rational points each). Gamma-faces: list of (tet index, local face
vertex triple).  Fields: C^0 piecewise P_k on the Alfeld refinement of the patch, zero on the
boundary of the patch, div-free.  Face functionals (exact, area-weighted normal N_F = cross
product oriented into the tet):
    mean_F(v) = int_F d_nu v dA  = int_ref (grad v . N_F)
    M_F(v)    = 1/2 sum_{i<j} l_ij^2 int_F mu_i mu_j d_nu v dA
Returns ranks of the constraint map and of the tangential moment map on Y^1.
"""
import itertools
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy import QQ

x, y, z = sp.symbols('x y z'); X = (x, y, z)
s_, t_ = sp.symbols('s t')

def q(v):
    v = sp.Rational(v); return QQ(int(v.p), int(v.q))

def bary(verts):
    M = sp.Matrix([[*v, 1] for v in verts]).T
    lam = M.inv() * sp.Matrix([x, y, z, 1])
    return [sp.expand(l) for l in lam]

def mi(k, n=4):
    for a in itertools.product(range(k + 1), repeat=n):
        if sum(a) == k: yield a

def ref_int(expr):
    p = sp.Poly(sp.expand(expr), s_, t_); r = sp.Rational(0)
    for (a, b), c in p.terms():
        r += c * sp.factorial(a) * sp.factorial(b) / sp.factorial(a + b + 2)
    return r

def analyse(macro, gfaces, k, tang, verbose=True):
    macro = [[sp.Matrix(v) for v in T] for T in macro]
    # boundary faces of the patch (macro faces not shared)
    def fkey(pts): return frozenset(tuple(p) for p in pts)
    cnt = {}
    for T in macro:
        for i in range(4):
            f = fkey([T[j] for j in range(4) if j != i]); cnt[f] = cnt.get(f, 0) + 1
    bfaces = [f for f, c in cnt.items() if c == 1]
    subs = []   # (verts, parent, opposite index)
    for ti, T in enumerate(macro):
        B = (T[0] + T[1] + T[2] + T[3]) / 4
        for i in range(4):
            subs.append(([B] + [T[j] for j in range(4) if j != i], ti, i))
    def on_bface(pt):
        # a domain point lies on the patch boundary iff it lies in a boundary macro face
        for f in bfaces:
            P = [sp.Matrix(p) for p in f]
            nrm = (P[1] - P[0]).cross(P[2] - P[0])
            if nrm.dot(pt - P[0]) != 0: continue
            # inside triangle? barycentric via solving
            A = sp.Matrix.hstack(P[1] - P[0], P[2] - P[0])
            sol = (A.T * A).inv() * A.T * (pt - P[0])
            a1, a2 = sol
            if a1 >= 0 and a2 >= 0 and a1 + a2 <= 1: return True
        return False
    gidx = {}; polys = []
    bcache = {}
    for (K, ti, i) in subs:
        lam = bary(K); bl = []
        for a in mi(k):
            pt = sum((sp.Rational(a[j], k) * K[j] for j in range(4)), sp.zeros(3, 1))
            key = tuple(pt)
            if key not in bcache: bcache[key] = on_bface(pt)
            if bcache[key]: continue
            if key not in gidx: gidx[key] = len(gidx)
            c = sp.factorial(k) / sp.prod([sp.factorial(ai) for ai in a])
            bl.append((gidx[key], sp.Poly(sp.expand(c * sp.prod([lam[j] ** a[j] for j in range(4)])), *X)))
        polys.append(bl)
    nunk = 3 * len(gidx)
    rows = []
    for bl in polys:
        co = {}
        for gi, p in bl:
            for c in range(3):
                for mon, cv in p.diff(X[c]).terms():
                    co.setdefault(mon, {}); u = 3 * gi + c
                    co[mon][u] = co[mon].get(u, 0) + cv
        for mon, d in co.items():
            r = [QQ(0)] * nunk
            for u, cv in d.items(): r[u] = q(cv)
            rows.append(r)
    A = DomainMatrix(rows, (len(rows), nunk), QQ)
    N = A.nullspace().to_Matrix()   # basis rows
    if verbose: print('  nullspace done', flush=True)
    dimY = N.shape[0]
    # face functionals
    cons = []; mom = sp.zeros(3, nunk)
    for (ti, loc) in gfaces:
        T = macro[ti]
        opp = [j for j in range(4) if j not in loc][0]
        P1, P2, P3 = [T[j] for j in loc]
        NF = (P2 - P1).cross(P3 - P1)
        if NF.dot(T[opp] - P1) < 0: NF = -NF
        # sub-tet containing the face: the one of parent ti opposite 'opp'
        si = [idx for idx, (K, pt, i) in enumerate(subs) if pt == ti and i == opp][0]
        l2 = {(1, 2): (P2 - P1).dot(P2 - P1), (1, 3): (P3 - P1).dot(P3 - P1), (2, 3): (P3 - P2).dot(P3 - P2)}
        mu = {1: 1 - s_ - t_, 2: s_, 3: t_}
        W = sp.Rational(1, 2) * sum(l2[(i, j)] * mu[i] * mu[j] for (i, j) in l2)
        sub = {x: P1[0] + s_ * (P2 - P1)[0] + t_ * (P3 - P1)[0],
               y: P1[1] + s_ * (P2 - P1)[1] + t_ * (P3 - P1)[1],
               z: P1[2] + s_ * (P2 - P1)[2] + t_ * (P3 - P1)[2]}
        crow = sp.zeros(3, nunk); mrow = sp.zeros(3, nunk)
        for gi, p in polys[si]:
            dnu = sum(p.diff(X[c]).as_expr() * NF[c] for c in range(3))
            dn = sp.expand(dnu.subs(sub))
            m0 = ref_int(dn); m1 = ref_int(W * dn)
            for c in range(3):
                crow[c, 3 * gi + c] = m0; mrow[c, 3 * gi + c] = m1
        cons.append(crow); mom += mrow
    Nd = DomainMatrix([[q(N[i, j]) for j in range(nunk)] for i in range(dimY)], (dimY, nunk), QQ).transpose()
    def dm(Mx):
        return DomainMatrix([[q(Mx[i, j]) for j in range(Mx.shape[1])] for i in range(Mx.shape[0])], Mx.shape, QQ)
    C = dm(sp.Matrix.vstack(*cons)) * Nd
    T2 = dm(sp.Matrix([list(tv) for tv in tang]))
    Mt = T2 * dm(mom) * Nd
    Kc = C.nullspace()             # rows
    dimY1 = Kc.shape[0]
    rk = (Mt * Kc.transpose()).rank() if dimY1 else 0
    if verbose:
        print(f'  k={k}: unknowns={nunk}, dim Y={dimY}, rank constraints={C.rank()}, dim Y1={dimY1}, rank tangential moment on Y1={rk}')
    return dimY, dimY1, rk

if __name__ == '__main__':
    import sys
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    # sanity: single macro-tet, random rational shape
    T = [(0, 0, 3), (0, 0, 0), (5, 1, 0), (1, 4, 0)]
    print('single macro-tet, general shape:')
    analyse([T], [(0, (1, 2, 3))], k, [(1, 0, 0), (0, 1, 0)])
