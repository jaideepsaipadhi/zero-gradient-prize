"""Exact (rational) checks of the local ingredients of Theorem L / L' (REPORT.tex).

(A) Locked star on the unit circle: z = (1,0), z_-, z_+ rational points at parameter +-t (chord ~ 2t),
    apex w = (1 + a t, 0) outside the circle (Omega_h is exterior).  T1 = (z_-, z, w), T2 = (z, z_+, w).
    For v = c * lam_w * lam_z^2 (c in R^2):  grad v|_{Ti}(z) = c (x) grad lam_w|_{Ti}, so the vertex
    divergences are D c with D = [grad lam_w|_{T1}; grad lam_w|_{T2}].  We report, with h = |e1|:
      h*sigma_max(D), h*sigma_min(D)/phi   (phi = turning angle)   -> both O(1)
      |c|/h for data (1,1) [the Q_lock direction]          -> bounded   (Theorem L')
      |c|/h * phi for data (1,0) [the Q_lock^perp direction]-> bounded, i.e. |c|/h ~ 1/phi
    and exactness facts: v vanishes on e1, e2 and on the outer edges, grad v = 0 at z_-, z_+, w.
(B) Bubble lemma: div maps (P_k cap H^1_0(T))^2 ONTO {q in P_{k-1}: int q = 0, q(vertices) = 0}; exact rank, k=3..8.
(C) The flux bubble b = lam_z^2 lam_w^2: grad b = 0 at all vertices, int_s b > 0.
Usage: python3 local_L.py
"""
import sympy as sy
from sympy import Rational as R

x, y = sy.symbols("x y")


def bary(P):
    (x1, y1), (x2, y2), (x3, y3) = P
    M = sy.Matrix([[x1, x2, x3], [y1, y2, y3], [1, 1, 1]])
    Mi = M.inv()
    return [sy.expand(Mi[i, 0] * x + Mi[i, 1] * y + Mi[i, 2]) for i in range(3)]


def circ(t):
    return (sy.nsimplify((1 - t**2) / (1 + t**2)), sy.nsimplify(2 * t / (1 + t**2)))


def partA(a=R(3, 2), b=0, ts=(R(1, 10), R(1, 20), R(1, 40), R(1, 80), R(1, 160))):
    print("(A) locked star, apex w=(1+a t, b t), a=%s b=%s" % (a, b))
    print("   t      h*smax   h*smin/phi   |c|/h (1,1)   |c|phi/h (1,0)")
    for t in ts:
        z = (1, 0); zm = circ(-t); zp = circ(t); w = (1 + a * t, b * t)
        T1 = [zm, z, w]; T2 = [z, zp, w]
        l1 = bary(T1); l2 = bary(T2)          # order: (z_-, z, w) and (z, z_+, w)
        lw1, lz1 = l1[2], l1[1]; lw2, lz2 = l2[2], l2[0]
        g1 = sy.Matrix([sy.diff(lw1, x), sy.diff(lw1, y)]).T; g2 = sy.Matrix([sy.diff(lw2, x), sy.diff(lw2, y)]).T
        D = sy.Matrix.vstack(g1, g2)
        h = sy.sqrt((zm[0] - 1)**2 + zm[1]**2)
        phi = 2 * sy.atan(t)   # z_+- at polar angles +-2atan(t): turning angle of the inscribed polygon at z
        import numpy as np
        sv = np.linalg.svd(np.array(D.evalf(30).tolist(), dtype=float), compute_uv=False)   # D exact, svd in floats
        smin, smax = sv[-1], sv[0]
        c11 = D.solve(sy.Matrix([1, 1])); c10 = D.solve(sy.Matrix([1, 0]))
        assert sy.simplify(c11 - sy.Matrix([w[0] - 1, w[1]])) == sy.zeros(2, 1)   # c = w - z exactly (Lemma K(c))
        nrm = lambda c: sy.sqrt(c[0]**2 + c[1]**2)
        print("  1/%-4d %.4f   %.4f       %.4f        %.4f" % (1 / t, float(h * smax), float(h * smin / phi),
              float(nrm(c11) / h), float(nrm(c10) * phi / h)))
        # exactness checks on the field f = lam_w lam_z^2 (scalar factor of v)
        f1 = lw1 * lz1**2; f2 = lw2 * lz2**2
        ok = True
        # continuity along s = [z,w]
        s = lambda r: (1 + r * (w[0] - 1), r * w[1])
        ok &= sy.simplify((f1 - f2).subs({x: s(sy.Symbol('r'))[0], y: s(sy.Symbol('r'))[1]})) == 0
        # vanishing on e1, e2 (lam_w = 0 there) and on outer edges (lam_z = 0): by construction; check grads at vertices
        for f, P in ((f1, T1), (f2, T2)):
            for V in P:
                gv = [sy.diff(f, x).subs({x: V[0], y: V[1]}), sy.diff(f, y).subs({x: V[0], y: V[1]})]
                if V != z:
                    ok &= all(sy.simplify(g) == 0 for g in gv)
        print("         continuity on s, grad=0 at z_-, z_+, w:", ok)


def partB():
    print("(B) bubble lemma: rank of div on (P_k cap H^1_0(T))^2 vs dim{q in P_{k-1}: mean 0, q(vertices)=0}")
    b = x * y * (1 - x - y)
    for k in range(3, 9):
        mons = [x**i * y**j for i in range(k - 2) for j in range(k - 2 - i)]    # b * P_{k-3}
        fields = [(b * m, 0) for m in mons] + [(0, b * m) for m in mons]
        qmons = [x**i * y**j for i in range(k) for j in range(k - i)]
        rows = []
        for (u, v) in fields:
            d = sy.Poly(sy.expand(sy.diff(u, x) + sy.diff(v, y)), x, y)
            rows.append([d.coeff_monomial(m) for m in qmons])
        rank = sy.Matrix(rows).rank()
        target = len(qmons) - 1 - 3
        print("   k=%d  dim bubbles=%d  rank(div)=%d  target dim=%d  onto=%s" % (k, len(fields), rank, target, rank == target))


def partC():
    print("(C) flux bubble lam_z^2 lam_w^2 on T = (z, z', w) reference: grads at vertices and edge integral")
    lz, lzp, lw = 1 - x - y, x, y
    b = lz**2 * lw**2
    g = [sy.diff(b, x), sy.diff(b, y)]
    print("   grads at vertices:", [[gg.subs({x: X, y: Y}) for gg in g] for (X, Y) in [(0, 0), (1, 0), (0, 1)]])
    r = sy.Symbol('r')
    print("   int over edge z-w (x=0, y=r) of b:", sy.integrate(b.subs({x: 0, y: r}), (r, 0, 1)))


if __name__ == "__main__":
    partA(); partA(R(3, 1))
    partA(R(3, 2), R(2, 3), (R(1, 10), R(1, 20), R(1, 40)))   # off-axis apex: sympy simplify slow beyond 1/40
    partB(); partC()
