"""Independent exact check that the pairing <(I - pi_T1) m_j, v_f.n>_e is the constant matrix of ref_pairing.py on a
non-reference T1 = (0,0),(1,0),(c1,c2) (sympy, exact), with m_j = l_z^j l_a^(4-j) and the traces B1-B2, B1-B3."""
import sympy as sp
x, y, s = sp.symbols('x y s')
for c1, c2 in [(sp.Rational(3, 10), sp.Rational(7, 5)), (sp.Rational(-1, 4), sp.Rational(2, 3)), (sp.Rational(13, 10), sp.Rational(1, 2))]:
    lz = 1 - x + (c1 - 1) / c2 * y; la = x - c1 / c2 * y
    mons = [x**i * y**j for i in range(4) for j in range(4 - i)]
    def intT(f):
        u, w = sp.symbols('u w')
        g = sp.expand(f.subs({x: u + w * c1, y: w * c2}, simultaneous=True)) * c2
        return sp.integrate(sp.integrate(g, (u, 0, 1 - w)), (w, 0, 1))
    Gm = sp.Matrix(len(mons), len(mons), lambda i, j: intT(mons[i] * mons[j]))
    def eta(H):
        cc = Gm.LUsolve(sp.Matrix([intT(H * m) for m in mons]))
        return sp.expand((H - sum(cc[i] * mons[i] for i in range(len(mons)))).subs(y, 0).subs(x, s))
    lzl = -x - (1 - c1) / c2 * y; lal = x - c1 / c2 * y          # linear parts
    B = [sp.binomial(4, i) * s**i * (1 - s)**(4 - i) for i in range(5)]
    g = [B[1] - B[2], B[1] - B[3]]
    P = sp.Matrix(2, 5, lambda f, j: sp.integrate(eta(sp.expand(lzl**j * lal**(4 - j))) * g[f], (s, 0, 1)))
    print((c1, c2), P)
