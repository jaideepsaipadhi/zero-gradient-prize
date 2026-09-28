"""Exact test (sympy rationals) of the vertex-field construction for (V) (robust6 REPORT, Sec. 3), k = 4, on a
3-triangle boundary star z; T1 = (z,a,c) with Gamma edge e = [z,a]; T2 = (z,c,b) interior; T3 = (z,b,g) with Gamma
edge e' = [g,z].
  step 1  Psi: a C^1 piecewise quintic on the star with Psi = grad Psi = 0 on the outer edges [a,c],[c,b],[b,g],
          Psi(z) = 0, grad Psi(z) = G   (an Argyris-type function; any particular solution);  v0 = curl Psi.
  step 2  p on T1 in lam_z lam_a P_2^2 with (p, grad r)_T1 = -<v0.n_e, r>_e for r in P_3; p' likewise on T3.
  step 3  corrections curl(chi), chi in S = {C^1 quintics, chi = grad chi = 0 on outer edges, chi = 0 on e, e',
          grad chi(z) = 0}; solve (M) = (D) = 0 (vectors, combined over e, e').
Checks every constraint of V#(omega_z) independently at the end, prints v(z), per-edge fluxes and v(z).n_e.
usage: python vertex_field.py tilt hdiff     (rationals; g = (-1 - hdiff?, ...) see code)"""
import sys
import sympy as sp
X, Y, t = sp.symbols('X Y t')
R = sp.Rational
def build(tilt, hc, hb, Gvec):
    z = sp.Matrix([0, 0]); a = sp.Matrix([1, 0]); c = sp.Matrix([R(1, 2), hc]); b = sp.Matrix([R(-1, 2), hb])
    g = sp.Matrix([-1, -tilt])                  # e' = [g, z]; tilt = 0 is flat
    tris = [(z, a, c), (z, c, b), (z, b, g)]
    mons5 = [X**i * Y**j for i in range(6) for j in range(6 - i)]
    def gen_poly(name, mons):
        cs = sp.symbols(name + '0:%d' % len(mons)); return sum(ci * m for ci, m in zip(cs, mons)), list(cs)
    def seg(expr, P, Q): return sp.expand(expr.subs({X: P[0] + t * (Q[0] - P[0]), Y: P[1] + t * (Q[1] - P[1])}, simultaneous=True))
    def coeffs_t(e): return sp.Poly(e, t).coeffs() if e != 0 else []
    def c1_system(prefix):
        polys, unk = [], []
        for i in range(3):
            p, cs = gen_poly(prefix + str(i) + '_', mons5); polys.append(p); unk += cs
        eqs = []
        grad = lambda p: [sp.diff(p, X), sp.diff(p, Y)]
        for (i, j, P, Q) in ((0, 1, z, c), (1, 2, z, b)):              # C^1 across spokes
            for f1, f2 in zip([polys[i]] + grad(polys[i]), [polys[j]] + grad(polys[j])):
                eqs += coeffs_t(seg(f1 - f2, P, Q))
        for i, (P, Q) in enumerate(((a, c), (c, b), (b, g))):          # Psi = grad Psi = 0 on outer edges
            for f in [polys[i]] + grad(polys[i]):
                eqs += coeffs_t(seg(f, P, Q))
        return polys, unk, eqs, grad
    # step 1
    Ps, unk, eqs, grad = c1_system('p')
    eqs += [Ps[0].subs({X: 0, Y: 0}), sp.diff(Ps[0], X).subs({X: 0, Y: 0}) - Gvec[0], sp.diff(Ps[0], Y).subs({X: 0, Y: 0}) - Gvec[1]]
    sol = sp.linsolve(eqs, unk); sol = list(sol)
    assert sol, "no Argyris-type Psi"
    sol = dict(zip(unk, sol[0])); free = set().union(*[sp.sympify(v).free_symbols for v in sol.values()]) & set(unk)
    sol = {k: sp.sympify(v).subs({f: 0 for f in free}) for k, v in sol.items()}
    Psi = [sp.expand(p.subs(sol)) for p in Ps]
    curl = lambda p: [sp.diff(p, Y), -sp.diff(p, X)]
    v = [curl(p) for p in Psi]
    # step 2: boundary corrections on T1 (Gamma edge e=[z,a]) and T3 (Gamma edge e'=[g,z])
    def bary(P0, P1, P2):
        M = sp.Matrix([[P0[0], P1[0], P2[0]], [P0[1], P1[1], P2[1]], [1, 1, 1]])
        return list(M.inv() * sp.Matrix([X, Y, 1]))
    def intT(expr, P0, P1, P2):
        u, w = sp.symbols('u w')
        J = sp.Matrix.hstack(P1 - P0, P2 - P0)
        e = sp.expand(expr.subs({X: P0[0] + J[0, 0] * u + J[0, 1] * w, Y: P0[1] + J[1, 0] * u + J[1, 1] * w}, simultaneous=True))
        return sp.integrate(sp.integrate(e, (u, 0, 1 - w)), (w, 0, 1)) * abs(J.det())
    def intE(expr, P, Q): return sp.integrate(seg(expr, P, Q), (t, 0, 1)) * sp.sqrt((Q[0] - P[0])**2 + (Q[1] - P[1])**2)
    def outward(P, Q, opp):
        tt = Q - P; n = sp.Matrix([tt[1], -tt[0]]) / sp.sqrt(tt.dot(tt))
        return n if n.dot(opp - P) < 0 else -n
    ne = outward(z, a, c); nep = outward(g, z, b)
    rs = [X**i * Y**j for i in range(4) for j in range(4 - i)]
    def bcorr(ti, P0, P1, opp, n, name):   # Gamma edge [P0,P1] of triangle (P0,P1,opp); p in lam_P0 lam_P1 P_2^2
        L = bary(P0, P1, opp); mons2 = [X**i * Y**j for i in range(3) for j in range(3 - i)]
        cs = sp.symbols(name + '0:12')
        p = [sum(cs[k] * L[0] * L[1] * m for k, m in enumerate(mons2)), sum(cs[6 + k] * L[0] * L[1] * m for k, m in enumerate(mons2))]
        gtr = v[ti][0] * n[0] + v[ti][1] * n[1]
        E = [intT(p[0] * sp.diff(r, X) + p[1] * sp.diff(r, Y), P0, P1, opp) + intE(gtr * r, P0, P1) for r in rs]
        s = list(sp.linsolve(E, cs)); assert s, "boundary correction not solvable"
        s = dict(zip(cs, s[0])); fr = set().union(*[sp.sympify(q).free_symbols for q in s.values()]) & set(cs)
        s = {k: sp.sympify(q).subs({f: 0 for f in fr}) for k, q in s.items()}
        return [sp.expand(q.subs(s)) for q in p]
    p1 = bcorr(0, z, a, c, ne, 'q'); p3 = bcorr(2, g, z, b, nep, 'r')
    v[0] = [v[0][i] + p1[i] for i in range(2)]; v[2] = [v[2][i] + p3[i] for i in range(2)]
    # step 3: stream-function corrections
    Cs, cunk, ceqs, _ = c1_system('s')
    for P, Q in ((z, a),): ceqs += coeffs_t(seg(Cs[0], P, Q))
    ceqs += coeffs_t(seg(Cs[2], g, z))
    ceqs += [sp.diff(Cs[0], X).subs({X: 0, Y: 0}), sp.diff(Cs[0], Y).subs({X: 0, Y: 0})]
    csol = list(sp.linsolve(ceqs, cunk))[0]
    csol = dict(zip(cunk, csol)); cfree = sorted(set().union(*[sp.sympify(q).free_symbols for q in csol.values()]) & set(cunk), key=str)
    basis = []
    for f in cfree:
        sub = {q: (1 if q == f else 0) for q in cfree}
        basis.append([curl(sp.expand(Cs[i].subs(csol).subs(sub))) for i in range(3)])
    def MD(vv):   # combined (M), (D) over e and e' (vectors); d_n = outward normal derivative
        out = []
        for comp in range(2):
            out.append(intE(vv[0][comp], z, a) + intE(vv[2][comp], g, z))
            dn1 = sp.diff(vv[0][comp], X) * ne[0] + sp.diff(vv[0][comp], Y) * ne[1]
            dn3 = sp.diff(vv[2][comp], X) * nep[0] + sp.diff(vv[2][comp], Y) * nep[1]
            out.append(intE(dn1, z, a) + intE(dn3, g, z))
        return [sp.nsimplify(sp.simplify(q)) for q in out]
    base = MD(v); A = sp.Matrix([MD(w) for w in basis]).T
    ks = sp.symbols('k0:%d' % len(basis))
    ksol = list(sp.linsolve((A, -sp.Matrix(base)), ks))
    assert ksol, "(M),(D) not solvable"
    kv = [sp.sympify(q).subs({k: 0 for k in ks}) for q in ksol[0]]
    for kk, w in zip(kv, basis):
        v = [[sp.expand(v[i][cc] + kk * w[i][cc]) for cc in range(2)] for i in range(3)]
    # independent final checks
    chk = {}
    chk["cont"] = all(seg(v[i][cc] - v[j][cc], P, Q) == 0 for (i, j, P, Q) in ((0, 1, z, c), (1, 2, z, b)) for cc in range(2))
    chk["zero on outer edges"] = all(seg(v[i][cc], P, Q) == 0 for i, (P, Q) in enumerate(((a, c), (c, b), (b, g))) for cc in range(2))
    chk["div=0 on T2"] = sp.expand(sp.diff(v[1][0], X) + sp.diff(v[1][1], Y)) == 0
    def moments(i, P0, P1, opp, n):
        dv = sp.diff(v[i][0], X) + sp.diff(v[i][1], Y); gtr = v[i][0] * n[0] + v[i][1] * n[1]
        return all(sp.simplify(intT(dv * r, P0, P1, opp) - intE(gtr * r, P0, P1)) == 0 for r in rs)
    chk["moments T1"] = moments(0, z, a, c, ne); chk["moments T3"] = moments(2, g, z, b, nep)
    Fe = sp.simplify(intE(v[0][0] * ne[0] + v[0][1] * ne[1], z, a)); Fep = sp.simplify(intE(v[2][0] * nep[0] + v[2][1] * nep[1], g, z))
    chk["flux e, e'"] = (Fe, Fep)
    chk["(M),(D) combined"] = all(q == 0 for q in MD(v))
    vz = [v[0][cc].subs({X: 0, Y: 0}) for cc in range(2)]
    chk["v(z)"] = vz; chk["v(z).n_e"] = sp.simplify(vz[0] * ne[0] + vz[1] * ne[1])
    chk["dim S (stream corrections)"] = len(basis); chk["rank of corrections on (M,D)"] = A.rank()
    return chk
if __name__ == "__main__":
    for (tilt, hc, hb) in ((R(0), R(1), R(1)), (R(0), R(1), R(3, 2)), (R(1, 5), R(1), R(1)), (R(-1, 3), R(4, 5), R(6, 5))):
        print("tilt", tilt, "hc", hc, "hb", hb, ":", build(tilt, hc, hb, (1, 0)), flush=True)
