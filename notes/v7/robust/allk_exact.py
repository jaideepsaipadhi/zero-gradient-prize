"""notes/v7/robust: exact (rational) verification, k = 4..KMAX, of the structural facts behind R-final (b),(c) for all k.
Reference triangle T^ = (z,a,c) = (0,0),(1,0),(0,1); Gamma edge e = [z,a] = {y=0}; outward normal n = (0,-1);
lam_z = 1-x-y, lam_a = x, lam_c = y.  Triangle integrals by the exact monomial formula int x^i y^j = i! j!/(i+j+2)!;
EDGE INTEGRALS BY RESTRICTION y=0 THEN int_0^1 dx  (so a polynomial with a factor y integrates to 0 on e:
this is exactly the robust5 bug, and it is regression-tested below).  Linear algebra over Q with python-flint.

Checks (numbering as in REPORT.tex, Sec. 1):
 F1 constraint map X -> (P_{k-1}/R)^*, p -> ((p,grad r))_r, has full rank dim P_{k-1}-1;  dim U_k = (k-1)(k-2)/2
 F2 traces of U_k: zero mean, rank k-2 (= dim P_e)
 F3 flux certificate: the 1D function f(y) (f' = g in P_{k-2}) with (p, grad f(lam_c)) = int_e p.n for ALL p in X exists
    (overdetermined system solved exactly)  ==> Lemma flux0 for this k; plus the direct check (flux of the solution
    for every mean-zero g in P_k(e) is 0)
 F4 zeta: the Riesz representative in P_{k-1}(T^) of w -> int_e w is a function of y alone, trace k(k+1)
 F5 pairing Hom_k x P_e, (H,p) -> <(I-pi_T)H, p>_e has rank k-2; kernel = span{l_c^k, l_z^k, l_a^k}
 F6 restriction of O_k = (I-pi_T)Hom_k to e is onto P_k(e) (rank k+1)
 F7 (I-pi_T) lam_a^k = J_k(lam_a) (J_k monic, orthogonal to P_{k-1} for weight 1-s), and beta_k = int J_k = (-1)^k J_k(0)/(k+1)
 F8 (k>=5) psi_z = lam_z^{k-2} lam_a^2 lam_c, psi_c = lam_z^{k-3} lam_a^2 lam_c^2: C^1 across the two interior edges
    (grad = 0 there), psi = 0 on e, grad psi(z) = 0, M_t(curl psi_c) = 0, M_t(curl psi_z) != 0, D_t(curl psi_c) != 0.
usage: python allk_exact.py [KMAX=12]
"""
import sys, time
from math import factorial
import sympy as sp
import flint

x, y = sp.symbols('x y')
Q = sp.Rational


def intT(p):
    p = sp.Poly(sp.expand(p), x, y)
    return sum(c * Q(factorial(i) * factorial(j), factorial(i + j + 2)) for (i, j), c in p.terms()) if not p.is_zero else Q(0)


def intE(p):                     # int over e = {y=0, 0<=x<=1}; restriction FIRST
    q = sp.Poly(sp.expand(sp.sympify(p).subs(y, 0)), x)
    return sum(c * Q(1, i + 1) for (i,), c in q.terms()) if not q.is_zero else Q(0)


def fm(rows):
    return flint.fmpq_mat([[flint.fmpq(int(sp.numer(v)), int(sp.denom(v))) for v in r] for r in rows])


def rank(rows):
    return fm(rows).rank() if rows and rows[0] else 0


def nullspace(rows, ncols):
    R, rk = fm(rows).rref()
    piv = []
    for i in range(rk):
        j = next(j for j in range(ncols) if R[i, j] != 0); piv.append(j)
    free = [j for j in range(ncols) if j not in piv]
    out = []
    for f in free:
        v = [sp.Integer(0)] * ncols; v[f] = sp.Integer(1)
        for i, pj in enumerate(piv):
            v[pj] = -sp.Rational(int(R[i, f].p), int(R[i, f].q))
        out.append(v)
    return out


def solve_exact(rows, rhs):
    """least-squares-free exact solve of overdetermined consistent system; returns solution or None"""
    A = fm(rows); b = fm([[v] for v in rhs])
    r = A.rank()
    if fm([list(rw) + [bv] for rw, bv in zip(rows, rhs)]).rank() != r:
        return None
    # solve via normal equations (A^T A full rank if columns independent)
    AtA = A.transpose() * A; Atb = A.transpose() * b
    s = AtA.solve(Atb)
    return [sp.Rational(int(s[i, 0].p), int(s[i, 0].q)) for i in range(s.nrows())]


lz, la, lc = 1 - x - y, x, y
P = lambda d: [x**i * y**j for i in range(d + 1) for j in range(d + 1 - i)]


def run(k):
    t0 = time.time(); out = {}
    mons = P(k - 2)
    Xb = [(lz * la * m, 0) for m in mons] + [(0, lz * la * m) for m in mons]
    rs = P(k - 1)[1:]
    A = [[intT(u0 * sp.diff(r, x) + u1 * sp.diff(r, y)) for (u0, u1) in Xb] for r in rs]
    rkA = rank(A)
    out['F1'] = (rkA == len(rs), len(Xb) - rkA == (k - 1) * (k - 2) // 2)
    N = nullspace(A, len(Xb))
    tr = [sp.expand(-sum(c * u1 for c, (u0, u1) in zip(v, Xb)).subs(y, 0)) for v in N]
    trm = [[sp.Poly(t, x).coeff_monomial(x**j) for j in range(k + 1)] for t in tr]
    out['F2'] = (all(intE(t) == 0 for t in tr), rank(trm) == k - 2)
    # F3: f' = g(y) in P_{k-2}; (p, grad f(y)) = int_e p.n  for all p in X  <=> for the e_y-part:
    #     int_T lz la q g(y) = - int_e lz la q  for all q in P_{k-2}   (e_x-part: both sides 0)
    gb = [y**j for j in range(k - 1)]
    rows = [[intT(lz * la * q * gj) for gj in gb] for q in mons]
    rhs = [-intE(lz * la * q) for q in mons]
    gsol = solve_exact(rows, rhs)
    ok3 = gsol is not None
    if ok3:
        g = sum(c * b for c, b in zip(gsol, gb)); f = sp.integrate(g, y)
        ok3 = all(intT(u0 * sp.diff(f, x) + u1 * sp.diff(f, y)) == intE(-u1) for (u0, u1) in Xb)
    # direct check: particular solution for each mean-zero g_e, then flux
    rows_all = [[intT(u0 * sp.diff(r, x) + u1 * sp.diff(r, y)) for (u0, u1) in Xb] for r in P(k - 1)]
    fl = []
    for j in range(1, k + 1):
        ge = x**j - Q(1, j + 1)
        b = [-intE(ge * r) for r in P(k - 1)]
        s = solve_exact(rows_all + [[1 if i == jj else 0 for i in range(len(Xb))] for jj in []], b) if False else None
        # solve with free variables fixed: use flint solve on a maximal independent column subset via lstsq-free route
        Mf = fm(rows_all); bf = fm([[v] for v in b])
        aug = fm([list(rw) + [bv] for rw, bv in zip(rows_all, b)])
        assert aug.rank() == Mf.rank()
        # particular solution: minimum-norm-free: solve (M M^T) w = b, p = M^T w  (M has full row rank after dropping r=1)
        M1 = fm(rows_all[1:]); b1 = fm([[v] for v in b[1:]])
        w = (M1 * M1.transpose()).solve(b1); pvec = M1.transpose() * w
        pv = [sp.Rational(int(pvec[i, 0].p), int(pvec[i, 0].q)) for i in range(pvec.nrows())]
        assert all(sum(c * a for c, a in zip(pv, rw)) == bv for rw, bv in zip(rows_all, b))
        fl.append(intE(-sum(c * u1 for c, (u0, u1) in zip(pv, Xb))))
    out['F3'] = (ok3, all(v == 0 for v in fl))
    # F4 zeta
    Pk1 = P(k - 1)
    G = [[intT(a * b) for b in Pk1] for a in Pk1]
    zc = fm(G).solve(fm([[intE(a)] for a in Pk1]))
    zeta = sum(sp.Rational(int(zc[i, 0].p), int(zc[i, 0].q)) * m for i, m in enumerate(Pk1))
    ztr = sp.expand(zeta.subs(y, 0))
    out['F4'] = (sp.expand(sp.diff(zeta, x)) == 0, ztr == k * (k + 1))
    # F5/F6
    Gi = fm(G).inv()
    def proj_perp(H):
        b = fm([[intT(H * a)] for a in Pk1]); c = Gi * b
        return sp.expand(H - sum(sp.Rational(int(c[i, 0].p), int(c[i, 0].q)) * m for i, m in enumerate(Pk1)))
    Hom = [x**i * y**(k - i) for i in range(k + 1)]
    Ob = [proj_perp(H) for H in Hom]
    Pe = [sp.expand(la * lz * x**j - (Q(1, 6) if j == 0 else 0)) for j in range(k - 1)]
    Pe = [sp.expand(la * lz * (x**j - Q(6, (j + 2) * (j + 3)))) for j in range(1, k - 1)]   # zero mean on e
    assert all(intE(p) == 0 for p in Pe) and len(Pe) == k - 2
    Pm = [[intE(o * p) for p in Pe] for o in Ob]
    out['F5_rank'] = rank(Pm) == k - 2
    def coeffs(H):
        pH = sp.Poly(sp.expand(H), x, y); return [pH.coeff_monomial(x**i * y**(k - i)) for i in range(k + 1)]
    kerT = [coeffs(l**k) for l in (y, -(x + y), x)]           # l_c, l_z, l_a (linear parts)
    ok5 = rank(kerT) == 3 and all(all(sum(c * Pm[i][j] for i, c in enumerate(h)) == 0 for j in range(k - 2)) for h in kerT)
    out['F5_kernel=Ker_T'] = ok5
    Or = [[sp.Poly(sp.expand(o.subs(y, 0)), x).coeff_monomial(x**j) for j in range(k + 1)] for o in Ob]
    out['F6'] = rank(Or) == k + 1
    # F7
    s = sp.symbols('s'); cs = sp.symbols('c0:%d' % k)
    J = s**k + sum(c * s**i for i, c in enumerate(cs))
    J = J.subs(sp.solve([sp.integrate(J * (1 - s) * s**i, (s, 0, 1)) for i in range(k)], cs))
    beta = sp.integrate(J, (s, 0, 1))
    out['F7'] = (sp.expand(proj_perp(x**k) - J.subs(s, x)) == 0, beta == (-1)**k * J.subs(s, 0) / (k + 1))
    out['beta_k'] = beta
    # F8
    if k >= 5:
        res = []
        for psi in (lz**(k - 2) * la**2 * lc, lz**(k - 3) * la**2 * lc**2):
            gx, gy = sp.diff(psi, x), sp.diff(psi, y)
            c1 = all(sp.expand(g.subs(x, 0)) == 0 for g in (gx, gy)) and all(sp.expand(g.subs(x, 1 - y)) == 0 for g in (gx, gy))
            on_e = sp.expand(psi.subs(y, 0)) == 0
            atz = gx.subs({x: 0, y: 0}) == 0 and gy.subs({x: 0, y: 0}) == 0
            Mt = intE(gy)                         # curl psi = (psi_y, -psi_x); t = (1,0)
            Dt = intE(-sp.diff(gy, y))            # d_n = -d_y
            res.append((c1, on_e, atz, Mt, Dt))
        (a1, b1, c1_, Mz, Dz), (a2, b2, c2, Mc, Dc) = res
        out['F8'] = (a1 and b1 and c1_ and a2 and b2 and c2, Mc == 0, Mz != 0, Dc != 0)
    out['secs'] = round(time.time() - t0, 1)
    return out


if __name__ == "__main__":
    # regression test of the edge integral (robust5 bug): functions with a factor lam_c = y integrate to 0 on e
    assert intE(y * x**3) == 0 and intE(lz * la * y) == 0 and intE(x**2) == Q(1, 3)
    print("edge-integral regression: OK")
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    for k in range(4, KMAX + 1):
        print("k=%d" % k, run(k), flush=True)
