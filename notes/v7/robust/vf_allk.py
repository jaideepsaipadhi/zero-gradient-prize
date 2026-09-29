"""notes/v7/robust: exact test (sympy rationals) of the all-k vertex field (REPORT.tex, Prop. 1.6 / Thm 1.7) for k >= 5.
Same star and step 1 as notes/v6/robust6/vertex_field.py (the stream function Psi is the SAME C^1 piecewise QUINTIC for
every k >= 4, since curl P_5 subset P_4 subset P_k), then
  step 2 (k-dependent): p on T1 in lam_z lam_a P_{k-2}^2 with (p, grad r)_T1 = -<curl Psi.n_e, r>_e for r in P_{k-1};
         likewise on T3.  (Lemma flux0 for this k says: flux 0.)
  step 3 (k >= 5, single-triangle, no (G) needed): curls of psi_z = lam_z^{k-2} lam_a^2 lam_c and
         psi_c = lam_z^{k-3} lam_a^2 lam_c^2 on T1, and the analogues on T3; solve combined (M) = (D) = 0 (4 x 4).
Every constraint of V#(omega_z) at degree k is re-checked independently at the end (moments against ALL of P_{k-1},
edge integrals by parametrisation of the segment, so functions vanishing on the edge integrate to 0).
usage: python vf_allk.py k [star index 0..3]"""
import sys, time
import sympy as sp
X, Y, t = sp.symbols('X Y t')
R = sp.Rational


def build(k, tilt, hc, hb, Gvec=(1, 0)):
    z = sp.Matrix([0, 0]); a = sp.Matrix([1, 0]); c = sp.Matrix([R(1, 2), hc]); b = sp.Matrix([R(-1, 2), hb])
    g = sp.Matrix([-1, -tilt])
    mons5 = [X**i * Y**j for i in range(6) for j in range(6 - i)]
    seg = lambda e, P, Q: sp.expand(sp.sympify(e).subs({X: P[0] + t * (Q[0] - P[0]), Y: P[1] + t * (Q[1] - P[1])}, simultaneous=True))
    ct = lambda e: sp.Poly(e, t).coeffs() if e != 0 else []
    grad = lambda p: [sp.diff(p, X), sp.diff(p, Y)]
    curl = lambda p: [sp.diff(p, Y), -sp.diff(p, X)]
    # ---- step 1: Psi
    polys, unk = [], []
    for i in range(3):
        cs = sp.symbols('p%d_0:%d' % (i, len(mons5))); polys.append(sum(ci * m for ci, m in zip(cs, mons5))); unk += cs
    eqs = []
    for (i, j, P, Q) in ((0, 1, z, c), (1, 2, z, b)):
        for f1, f2 in zip([polys[i]] + grad(polys[i]), [polys[j]] + grad(polys[j])):
            eqs += ct(seg(f1 - f2, P, Q))
    for i, (P, Q) in enumerate(((a, c), (c, b), (b, g))):
        for f in [polys[i]] + grad(polys[i]):
            eqs += ct(seg(f, P, Q))
    eqs += [polys[0].subs({X: 0, Y: 0}), sp.diff(polys[0], X).subs({X: 0, Y: 0}) - Gvec[0],
            sp.diff(polys[0], Y).subs({X: 0, Y: 0}) - Gvec[1]]
    sol = list(sp.linsolve(eqs, unk)); assert sol
    sol = dict(zip(unk, sol[0])); free = set().union(*[sp.sympify(v).free_symbols for v in sol.values()]) & set(unk)
    Psi = [sp.expand(p.subs({kk: sp.sympify(v).subs({f: 0 for f in free}) for kk, v in sol.items()})) for p in polys]
    v = [curl(p) for p in Psi]

    def bary(P0, P1, P2):
        M = sp.Matrix([[P0[0], P1[0], P2[0]], [P0[1], P1[1], P2[1]], [1, 1, 1]])
        return list(M.inv() * sp.Matrix([X, Y, 1]))

    def intT(expr, P0, P1, P2):
        u, w = sp.symbols('u w'); J = sp.Matrix.hstack(P1 - P0, P2 - P0)
        e = sp.expand(sp.sympify(expr).subs({X: P0[0] + J[0, 0] * u + J[0, 1] * w, Y: P0[1] + J[1, 0] * u + J[1, 1] * w}, simultaneous=True))
        return sp.integrate(sp.integrate(e, (u, 0, 1 - w)), (w, 0, 1)) * abs(J.det())

    def intE(expr, P, Q):
        return sp.integrate(seg(expr, P, Q), (t, 0, 1)) * sp.sqrt((Q[0] - P[0])**2 + (Q[1] - P[1])**2)

    def outward(P, Q, opp):
        tt = Q - P; n = sp.Matrix([tt[1], -tt[0]]) / sp.sqrt(tt.dot(tt))
        return n if n.dot(opp - P) < 0 else -n
    ne = outward(z, a, c); nep = outward(g, z, b)
    rs = [X**i * Y**j for i in range(k) for j in range(k - i)]
    # regression: edge integral of a function vanishing on the edge
    Lz1 = bary(z, a, c)
    assert intE(Lz1[2] * X**3, z, a) == 0

    def bcorr(ti, P0, P1, opp, n, name):
        L = bary(P0, P1, opp); m2 = [X**i * Y**j for i in range(k - 1) for j in range(k - 1 - i)]
        cs = sp.symbols(name + '0:%d' % (2 * len(m2)))
        p = [sum(cs[kk] * L[0] * L[1] * m for kk, m in enumerate(m2)),
             sum(cs[len(m2) + kk] * L[0] * L[1] * m for kk, m in enumerate(m2))]
        gtr = v[ti][0] * n[0] + v[ti][1] * n[1]
        E = [intT(p[0] * sp.diff(r, X) + p[1] * sp.diff(r, Y), P0, P1, opp) + intE(gtr * r, P0, P1) for r in rs]
        s = list(sp.linsolve(E, cs)); assert s, "boundary correction not solvable"
        s = dict(zip(cs, s[0])); fr = set().union(*[sp.sympify(q).free_symbols for q in s.values()]) & set(cs)
        s = {kk: sp.sympify(q).subs({f: 0 for f in fr}) for kk, q in s.items()}
        return [sp.expand(q.subs(s)) for q in p]
    p1 = bcorr(0, z, a, c, ne, 'q'); p3 = bcorr(2, g, z, b, nep, 'r')
    v[0] = [v[0][i] + p1[i] for i in range(2)]; v[2] = [v[2][i] + p3[i] for i in range(2)]
    flux_before = (sp.simplify(intE(v[0][0] * ne[0] + v[0][1] * ne[1], z, a)), sp.simplify(intE(v[2][0] * nep[0] + v[2][1] * nep[1], g, z)))
    # ---- step 3: single-triangle stream corrections
    L1 = bary(z, a, c)          # lam_z, lam_a, lam_c on T1
    L3 = bary(g, z, b)          # lam_g, lam_z, lam_b on T3  (Gamma edge [g,z], apex b)
    psis = [(0, L1[0]**(k - 2) * L1[1]**2 * L1[2]), (0, L1[0]**(k - 3) * L1[1]**2 * L1[2]**2),
            (2, L3[1]**(k - 2) * L3[0]**2 * L3[2]), (2, L3[1]**(k - 3) * L3[0]**2 * L3[2]**2)]
    basis = []
    for ti, ps in psis:
        w = [[0, 0], [0, 0], [0, 0]]; w[ti] = [sp.expand(q) for q in curl(ps)]; basis.append(w)

    def MD(vv):
        out = []
        for comp in range(2):
            out.append(intE(vv[0][comp], z, a) + intE(vv[2][comp], g, z))
            dn1 = sp.diff(vv[0][comp], X) * ne[0] + sp.diff(vv[0][comp], Y) * ne[1]
            dn3 = sp.diff(vv[2][comp], X) * nep[0] + sp.diff(vv[2][comp], Y) * nep[1]
            out.append(intE(dn1, z, a) + intE(dn3, g, z))
        return [sp.nsimplify(sp.simplify(q)) for q in out]
    base = MD(v); Am = sp.Matrix([MD(w) for w in basis]).T
    ks = sp.symbols('k0:4'); ksol = list(sp.linsolve((Am, -sp.Matrix(base)), ks)); assert ksol, "(M),(D) not solvable"
    kv = [sp.sympify(q).subs({kk: 0 for kk in ks}) for q in ksol[0]]
    for kk, w in zip(kv, basis):
        v = [[sp.expand(v[i][cc] + kk * w[i][cc]) for cc in range(2)] for i in range(3)]
    # ---- independent checks
    chk = {}
    chk["C0"] = all(seg(v[i][cc] - v[j][cc], P, Q) == 0 for (i, j, P, Q) in ((0, 1, z, c), (1, 2, z, b)) for cc in range(2))
    chk["zero outer"] = all(seg(v[i][cc], P, Q) == 0 for i, (P, Q) in enumerate(((a, c), (c, b), (b, g))) for cc in range(2))
    chk["div T2"] = sp.expand(sp.diff(v[1][0], X) + sp.diff(v[1][1], Y)) == 0
    def moments(i, P0, P1, opp, n):
        dv = sp.diff(v[i][0], X) + sp.diff(v[i][1], Y); gtr = v[i][0] * n[0] + v[i][1] * n[1]
        return all(sp.simplify(intT(dv * r, P0, P1, opp) - intE(gtr * r, P0, P1)) == 0 for r in rs)
    chk["moments T1,T3 (P_%d)" % (k - 1)] = (moments(0, z, a, c, ne), moments(2, g, z, b, nep))
    chk["flux e,e' (before step 3)"] = flux_before
    chk["flux e,e'"] = (sp.simplify(intE(v[0][0] * ne[0] + v[0][1] * ne[1], z, a)), sp.simplify(intE(v[2][0] * nep[0] + v[2][1] * nep[1], g, z)))
    chk["(M),(D)"] = all(q == 0 for q in MD(v))
    vz = [sp.simplify(v[0][cc].subs({X: 0, Y: 0})) for cc in range(2)]
    chk["v(z)"] = vz; chk["v(z).n_e"] = sp.simplify(vz[0] * ne[0] + vz[1] * ne[1]); chk["rank(M,D)"] = Am.rank()
    chk["deg<=k"] = all(sp.Poly(v[i][cc], X, Y).total_degree() <= k for i in range(3) for cc in range(2))
    return chk


STARS = ((R(0), R(1), R(1)), (R(0), R(1), R(3, 2)), (R(1, 5), R(1), R(1)), (R(-1, 3), R(4, 5), R(6, 5)))
if __name__ == "__main__":
    k = int(sys.argv[1]); idx = [int(s) for s in sys.argv[2].split(",")] if len(sys.argv) > 2 else range(4)
    for i in idx:
        t0 = time.time(); tilt, hc, hb = STARS[i]
        print("k", k, "tilt", tilt, "hc", hc, "hb", hb, ":", build(k, tilt, hc, hb), " secs", round(time.time() - t0, 1), flush=True)
