"""PROVED (exact rational arithmetic): on the stacked column half-plane mesh (cells [i,i+1] x [jH,(j+1)H] split by
the forward diagonal; every boundary star is the counterexample star S^(H)), the separable stream function
    F(x,y) = Phi(x) G(y),
Phi = uniform cubic B-spline with knot spacing d (integer), support [-2d, 2d]  (C^2, piecewise cubic, knots on
      vertical mesh lines),
G   = C^1 piecewise quadratic with knots at multiples of H:
      G = y - y^2/(2H) on [0,H]; then G'' = -c on [H,(1+M)H], +c on [(1+M)H,(1+2M)H], G = 0 after, c = 1/(2 M^2 H),
is C^1 and piecewise P5 on the mesh, so v = curl F lies in Z_h for every k >= 4 (continuous piecewise P4, div v = 0,
compact support).  Because F is in the discrete space there is NO interpolation error, and with h_G = 1
    Q := (2B - A)/C = 2/H - int G''^2 - 2 R1 int G'^2 - R2 int G^2,  R1 = int Phi'^2/int Phi^2, R2 = int Phi''^2/int Phi^2
(A = 2 int Phi'^2 int G'^2 + int Phi^2 int G''^2 + int Phi''^2 int G^2, B = (1/H) int Phi^2, C = int Phi^2).
Q > 0 certifies that N_h is indefinite for mu/h < Q/h_G on any mesh containing this patch (support) exactly."""
import sympy as sp

x, y = sp.symbols('x y', real=True)


def bspline(d):
    """uniform cubic B-spline on knots -2d,-d,0,d,2d as a list of (a,b,poly)"""
    t = sp.Symbol('t')
    pieces_t = [(-2, -1, (2 + t)**3 / 6), (-1, 0, (4 - 6 * t**2 - 3 * t**3) / 6),
                (0, 1, (4 - 6 * t**2 + 3 * t**3) / 6), (1, 2, (2 - t)**3 / 6)]
    return [(a * d, b * d, sp.expand(p.subs(t, x / d))) for a, b, p in pieces_t]


def G_pieces(H, M):
    c = sp.Rational(1, 2) / (M**2 * H)
    P1 = y - y**2 / (2 * H)
    P2 = H / 2 - c * (y - H)**2 / 2
    y2 = (1 + M) * H
    P3 = c * (y - (1 + 2 * M) * H)**2 / 2
    assert sp.simplify(P2.subs(y, y2) - P3.subs(y, y2)) == 0
    assert sp.simplify(sp.diff(P2, y).subs(y, y2) - sp.diff(P3, y).subs(y, y2)) == 0
    return [(0, H, P1), (H, y2, P2), (y2, (1 + 2 * M) * H, P3)]


def I(pieces, var, f):
    return sum(sp.integrate(f(p), (var, a, b)) for a, b, p in pieces)


def Q(H, M, d):
    H = sp.nsimplify(H); Ph = bspline(d); Gp = G_pieces(H, M)
    i0 = I(Ph, x, lambda p: p**2); i1 = I(Ph, x, lambda p: sp.diff(p, x)**2); i2 = I(Ph, x, lambda p: sp.diff(p, x, 2)**2)
    g0 = I(Gp, y, lambda p: p**2); g1 = I(Gp, y, lambda p: sp.diff(p, y)**2); g2 = I(Gp, y, lambda p: sp.diff(p, y, 2)**2)
    A = 2 * i1 * g1 + i0 * g2 + i2 * g0
    B = i0 / H; C = i0
    return sp.nsimplify(sp.simplify((2 * B - A) / C), rational=True)


if __name__ == "__main__":
    print("exact Q = h_G (2B-A)/C of the separable witness on stacked column meshes (h_G = 1)")
    for H in (1, 5, 8, 20):
        for (M, d) in ((2, 20), (2, 60), (3, 200)):
            q = Q(H, M, d)
            print(f"H={H:3d} M={M} knot spacing d={d:4d} (support width {4*d}):  Q = {q}  = {float(q):.6f}"
                  f"   (limit (1 - 1/(2M^3))/H = {(1 - 1/(2*M**3))/H:.6f})", flush=True)
