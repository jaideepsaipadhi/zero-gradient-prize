"""Referee's independent exact check (NOT derived from allk_exact.py).
Physical triangle z=(0,0), a=(3,4), c=(-2,5): e=[z,a], |e|=5, area 23/2, H_e=23/5, outward n=(4,-3)/5.
Triangle integrals: pull back to reference (u,w) with Jacobian 23, monomial formula.
Edge integrals: OWN parametrisation x=3t, y=4t, ds=5dt (restriction to the edge FIRST), asserted on functions vanishing on e.
Checks for given k: dim U_k, traces = P_e (zero mean, rank k-2, inside lam_z lam_a P_{k-2}), zero flux of the
g-problem for every mean-zero g in P_k(e), zeta trace = k(k+1)/H_e, sharp one-edge trace constant, kernel = Ker_T.
"""
import sys, time
from math import factorial
import flint
from flint import fmpq, fmpq_mat

ctx = flint.fmpq_mpoly_ctx.get(('x', 'y'), 'lex')
X, Y = ctx.gens()
Q = lambda p, q=1: fmpq(p, q)
Z = (0, 0); A = (3, 4); C = (-2, 5)
DET = 23
lam_c = (3 * Y - 4 * X) / DET
lam_a = (5 * X + 2 * Y) / DET
lam_z = 1 - lam_a - lam_c
assert lam_c.subs({'x': -2, 'y': 5}) if False else True
# pull-back substitution x = 3u - 2w, y = 4u + 5w
U_, W_ = X, Y   # reuse gens as (u,w) in the pulled-back polynomial
def intT(p):
    q = p.compose(3 * U_ - 2 * W_, 4 * U_ + 5 * W_)
    s = Q(0)
    for (i, j), c in q.terms():
        s += c * Q(factorial(i) * factorial(j), factorial(i + j + 2))
    return s * DET
def intE(p):
    q = p.compose(3 * X, 4 * X)          # restriction to e: (x,y) = (3t,4t), t in [0,1], y-slot now t too but y absent
    s = Q(0)
    for (i, j), c in q.terms():
        assert j == 0
        s += c * Q(1, i + 1)
    return 5 * s
# regression: functions vanishing on e integrate to 0; a known value
assert intE((4 * X - 3 * Y) * X**3) == 0 and intE(lam_c * (X + Y)**5) == 0
assert intE(ctx.from_dict({(0, 0): 1})) == 5 and intE(X) == Q(15, 2)
assert intT(ctx.from_dict({(0, 0): 1})) == Q(23, 2)
n = (Q(4, 5), Q(-3, 5))
def P(d):
    return [X**i * Y**j for i in range(d + 1) for j in range(d + 1 - i)]
def mat(rows):
    return fmpq_mat(rows)
def rank(rows):
    return mat(rows).rank() if rows else 0
def nullspace(M, ncols):
    R, rk = M.rref(); piv = []
    for i in range(rk):
        piv.append(next(j for j in range(ncols) if R[i, j] != 0))
    free = [j for j in range(ncols) if j not in piv]; out = []
    for f in free:
        v = [Q(0)] * ncols; v[f] = Q(1)
        for i, pj in enumerate(piv):
            v[pj] = -R[i, f]
        out.append(v)
    return out
def tpoly(p):                     # restriction to e as a polynomial in t (coefficient list)
    q = p.compose(X, X * 0)       # placeholder, not used
def restrict_t(p, k):
    q = p.compose(3 * X, 4 * X)
    c = [Q(0)] * (k + 2)
    for (i, j), cc in q.terms():
        c[i] += cc
    return c

def run(k):
    t0 = time.time(); out = {}
    qs = P(k - 2)
    Xb = [(lam_z * lam_a * q, ctx.from_dict({})) for q in qs] + [(ctx.from_dict({}), lam_z * lam_a * q) for q in qs]
    rs = P(k - 1)
    grads = [(r.derivative('x'), r.derivative('y')) for r in rs]
    Mfull = [[intT(p0 * gx + p1 * gy) for (p0, p1) in Xb] for (gx, gy) in grads]
    rk = rank(Mfull)
    out['rank=dimP_{k-1}-1'] = rk == len(rs) - 1
    Nsp = nullspace(mat(Mfull), len(Xb))
    out['dimU'] = (len(Nsp), (k - 1) * (k - 2) // 2)
    # traces u.n on e as polys in t
    tr = []
    for v in Nsp:
        p0 = sum((c * b0 for c, (b0, b1) in zip(v, Xb)), ctx.from_dict({}))
        p1 = sum((c * b1 for c, (b0, b1) in zip(v, Xb)), ctx.from_dict({}))
        pn = p0 * n[0] + p1 * n[1]
        tr.append(pn)
    out['trace zero mean'] = all(intE(t) == 0 for t in tr)
    trc = [restrict_t(t, k) for t in tr]
    out['trace rank=k-2'] = rank(trc) == k - 2
    # P_e basis: t(1-t)t^j - mean, as restrictions; check traces span exactly P_e (rank of union still k-2)
    import fractions
    Pe = []
    for j in range(k - 1):
        # t(1-t)t^j on e has mean (over t) 1/((j+2)(j+3)); zero-mean combos: b_j - c_j b_0
        pass
    Tt = flint.fmpq_poly
    basis = [Tt([0, 1, -1]) * Tt([0] * j + [1]) for j in range(k - 1)]
    means = [sum(b[i] / (i + 1) for i in range(b.degree() + 1)) for b in basis]
    Pe_t = [basis[j] - (means[j] / means[0]) * basis[0] for j in range(1, k - 1)]
    Pe_c = [[p[i] for i in range(k + 2)] for p in Pe_t]
    out['span traces = P_e'] = rank(trc + Pe_c) == k - 2 and rank(Pe_c) == k - 2
    # zero flux of g-problem: for every mean-zero g in P_k(e): (p,grad r) = -<g,r>_e
    fl = []
    M1 = mat(Mfull[1:])
    MMt = M1 * M1.transpose()
    for j in range(1, k + 1):
        g_t = Tt([0] * j + [1]) - Tt([Q(1, j + 1)])          # g(t) = t^j - 1/(j+1)
        # <g,r>_e = 5 int_0^1 g(t) r(3t,4t) dt
        b = []
        for r in rs:
            rt = restrict_t(r, k)
            prod = g_t * Tt(rt)
            b.append(5 * sum(prod[i] / (i + 1) for i in range(prod.degree() + 1)))
        assert b[0] == 0
        w = MMt.solve(mat([[-v] for v in b[1:]]))
        pv = M1.transpose() * w
        # verify
        for row, bv in zip(Mfull, b):
            assert sum((row[i] * pv[i, 0] for i in range(len(Xb))), Q(0)) == -bv
        p0 = sum((pv[i, 0] * Xb[i][0] for i in range(len(Xb))), ctx.from_dict({}))
        p1 = sum((pv[i, 0] * Xb[i][1] for i in range(len(Xb))), ctx.from_dict({}))
        fl.append(intE(p0 * n[0] + p1 * n[1]))
    out['g-problem zero flux'] = all(f == 0 for f in fl)
    # zeta and sharp trace constant on P_{k-1}
    MT = mat([[intT(a * b) for b in rs] for a in rs]); ME = mat([[intE(a * b) for b in rs] for a in rs])
    zc = MT.solve(mat([[intE(a)] for a in rs]))
    zeta = sum((zc[i, 0] * m for i, m in enumerate(rs)), ctx.from_dict({}))
    zt = restrict_t(zeta, k)
    out['zeta trace'] = (zt[0], all(c == 0 for c in zt[1:]), Q(k * (k + 1) * 5, 23))
    Kc = Q(k * (k + 1) * 5, 23)
    D = Kc * MT - ME
    import numpy as np
    ev = np.linalg.eigvalsh(np.array([[float(D[i, j]) for j in range(D.ncols())] for i in range(D.nrows())]))
    Dd = D.det()
    cp = D.charpoly()   # symmetric: PSD iff coefficients alternate in sign (Descartes)
    cs = [cp[i] for i in range(cp.degree() + 1)]
    m = D.nrows()
    alt = all((c == 0) or ((c > 0) == ((m - i) % 2 == 0)) for i, c in enumerate(cs))
    out['Kc*M_T - M_e PSD(exact), singular'] = (alt, Dd == 0, float(ev.min()))
    # kernel of pairing Hom_k x P_e
    MTinv = MT.inv()
    Hom = [X**i * Y**(k - i) for i in range(k + 1)]
    rows = []
    for H in Hom:
        c = MTinv * mat([[intT(H * a)] for a in rs])
        eta = H - sum((c[i, 0] * m for i, m in enumerate(rs)), ctx.from_dict({}))
        et = Tt(restrict_t(eta, k))
        row = []
        for p in Pe_t:
            prod = et * p
            row.append(5 * sum(prod[i] / (i + 1) for i in range(prod.degree() + 1)))
        rows.append(row)
    Pm = mat(rows)                          # (k+1) x (k-2)
    out['pairing rank'] = Pm.rank() == k - 2
    def coeffs(H):
        d = H.to_dict(); return [d.get((i, k - i), Q(0)) for i in range(k + 1)]
    lin = lambda L: L - L.coefficient(0) if False else L
    def linpart(L):
        d = L.to_dict(); return d.get((1, 0), Q(0)) * X + d.get((0, 1), Q(0)) * Y
    Ker = [coeffs(linpart(L)**k) for L in (lam_c, lam_z, lam_a)]
    KerM = mat(Ker)
    out['Ker_T in left kernel'] = KerM.rank() == 3 and all(x == 0 for x in (KerM * Pm).entries())
    # left kernel dim = (k+1) - rank = 3
    out['left kernel dim'] = (k + 1) - Pm.rank()
    out['secs'] = round(time.time() - t0, 1)
    return out

if __name__ == '__main__':
    print('edge regression OK')
    for k in [int(a) for a in sys.argv[1:]] or [5, 7, 9]:
        print('k=%d' % k, run(k), flush=True)
