"""Exact (rational) computation of the local divergence-free space on the Alfeld split of a
reference tetrahedron with one face on Gamma_h, and of the face functionals used in (M3^3).

Reference macro-tet  T = conv(P0,P1,P2,P3), P0=(0,0,1) apex, face F=conv(P1,P2,P3) in z=0,
P1=(0,0,0), P2=(1,0,0), P3=(0,1,0).  Inward normal of T on F is e_z.
Alfeld split: barycentre B, sub-tets K_i = conv(B, face opposite P_i); K_0 contains F.

Space  Y  = { v in C^0 piecewise P_k (on the split), v = 0 on dT, div v = 0 }.
Face data (from K_0):  g = d_z v on F  (tangential: g_z = 0 automatically).
Functionals:  mean  m(v)   = int_F g           (2 comps)
              E_ij(v)      = int_F mu_i mu_j g   (2 comps each), mu1=1-x-y, mu2=x, mu3=y.
Usage: python3 alfeld_witness.py k
"""
import sys, itertools
from fractions import Fraction as Fr
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy import QQ

x, y, z = sp.symbols('x y z')
def q(v):
    v = sp.Rational(v)
    return QQ(int(v.p), int(v.q))
X = (x, y, z)

def bary_funcs(verts):
    # affine barycentric coordinates of tet with vertices verts (list of 4 rational 3-tuples)
    M = sp.Matrix([[*v, 1] for v in verts]).T  # 4x4 columns (v,1)
    Minv = M.inv()
    vec = sp.Matrix([x, y, z, 1])
    lam = Minv * vec
    return [sp.expand(l) for l in lam]

def multi_indices(k, n=4):
    for a in itertools.product(range(k + 1), repeat=n):
        if sum(a) == k:
            yield a

def build(k):
    P = [sp.Matrix([0, 0, 1]), sp.Matrix([0, 0, 0]), sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 0])]
    B = (P[0] + P[1] + P[2] + P[3]) / 4
    subs = []
    for i in range(4):
        face = [P[j] for j in range(4) if j != i]
        subs.append([B] + face)          # vertex 0 of each sub-tet is B
    # global interior domain points
    gidx = {}
    local = []   # per sub-tet: list of (alpha, global index or None)
    for K in subs:
        lst = []
        for a in multi_indices(k):
            pt = sum((sp.Rational(a[j], k) * K[j] for j in range(4)), sp.zeros(3, 1))
            key = tuple(pt)
            if a[0] == 0:
                lst.append((a, None))    # on dT: coefficient zero
            else:
                if key not in gidx:
                    gidx[key] = len(gidx)
                lst.append((a, gidx[key]))
        local.append(lst)
    nd = len(gidx)
    nunk = 3 * nd
    # Bernstein polys per sub-tet
    polys = []
    for K, lst in zip(subs, local):
        lam = bary_funcs(K)
        bl = []
        for a, gi in lst:
            if gi is None:
                continue
            c = sp.factorial(k) / sp.prod([sp.factorial(ai) for ai in a])
            bl.append((gi, sp.Poly(sp.expand(c * sp.prod([lam[j] ** a[j] for j in range(4)])), *X)))
        polys.append(bl)
    # divergence rows
    rows = []
    for bl in polys:
        coeffs = {}  # monomial -> {unknown: value}
        for gi, p in bl:
            for comp in range(3):
                dp = p.diff(X[comp])
                for mon, cval in dp.terms():
                    coeffs.setdefault(mon, {})
                    u = 3 * gi + comp
                    coeffs[mon][u] = coeffs[mon].get(u, 0) + cval
        for mon, d in coeffs.items():
            r = [QQ(0)] * nunk
            for u, cval in d.items():
                r[u] = q(cval)
            rows.append(r)
    A = DomainMatrix(rows, (len(rows), nunk), QQ)
    return subs, polys, nd, nunk, A

def tri_int(expr):
    # integral over reference triangle {x,y>=0, x+y<=1} of polynomial in x,y
    p = sp.Poly(sp.expand(expr), x, y)
    s = sp.Rational(0)
    for (a, b), c in p.terms():
        s += c * sp.factorial(a) * sp.factorial(b) / sp.factorial(a + b + 2)
    return s

def face_functionals(polys, nunk):
    # g = d_z v on z=0 from sub-tet K_0 ; returns dict name -> row vector (len nunk)
    mu = {1: 1 - x - y, 2: x, 3: y}
    weights = {'m': 1, 'E12': mu[1] * mu[2], 'E13': mu[1] * mu[3], 'E23': mu[2] * mu[3]}
    out = {}
    for name, wgt in weights.items():
        for comp in range(2):   # tangential components x,y
            r = [QQ(0)] * nunk
            for gi, p in polys[0]:
                g = p.diff(z).as_expr().subs(z, 0)
                val = tri_int(wgt * g)
                r[3 * gi + comp] = q(val)
            out[(name, comp)] = r
    return out

if __name__ == '__main__':
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    subs, polys, nd, nunk, A = build(k)
    rk = A.rank()
    print(f'k={k}: interior domain pts={nd}, unknowns={nunk}, div rows={A.shape[0]}, rank(div)={rk}')
    N = A.nullspace()               # rows span the kernel
    Nm = N.to_Matrix()
    dimY = Nm.shape[0]
    print(f'dim Y (local div-free, zero on dT) = {dimY}')
    # check divergence surjectivity onto P_{k-1}^disc mean-zero: dim = 4*C(k+2,3)-1
    target = 4 * sp.binomial(k + 2, 3) - 1
    print(f'dim P_(k-1)^disc(T^A) - 1 = {target}; rank(div) = {rk}  -> onto mean-zero: {rk == target}')
    F = face_functionals(polys, nunk)
    order = [('m', 0), ('m', 1), ('E12', 0), ('E12', 1), ('E13', 0), ('E13', 1), ('E23', 0), ('E23', 1)]
    L = sp.Matrix([[sum(sp.Rational(F[o][j].numerator, F[o][j].denominator) * Nm[i, j] for j in range(nunk)) for i in range(dimY)] for o in order])
    # L: 8 x dimY  (functionals on the basis of Y)
    print('rank of [mean; E12; E13; E23] on Y =', L.rank())
    print('rank of [mean] on Y =', L[0:2, :].rank())
    import pickle
    with open(f'alfeld_k{k}.pkl', 'wb') as fh:
        pickle.dump({'L': L, 'N': Nm, 'nunk': nunk}, fh)
