"""Construct test potentials phi = Re sum_j c_j (r-1)^j e^{i m theta} whose 4th-derivative tensor on the unit circle lies,
at every point, in a prescribed subspace of Hom_4 written in the local frame x = e_theta (ccw tangent), y = e_r
(outward radial = inward normal of Omega_h).   H_* = x^4 - (x - y/rho)^4   (rho = height/edge of the limit cell,
apex radially above the ccw-forward endpoint), K = span{y^4, H_*}.
Modes: 'B'  : grad^4 phi/4! in span{H_*}          (pure H_* direction)
       'K'  : grad^4 phi/4! in span{y^4, H_*}     (whole invisible set; choose the solution with nonzero H_* part)
       'Q'  : grad^4 phi/4! in span{x^4}  (visible control: tangential power, not in K)
Prints complex coefficients c_j (j=0..4) for mode m.   python kphi.py m rho
"""
import sys, json
import sympy as sp
X, Y = sp.symbols('X Y', real=True)
xl, yl = sp.symbols('xl yl')


def tensor_poly(f):
    """H(xl,yl) = sum_{|a|=4} d^a f(1,0) / a! * (point displacement)^a, with displacement = yl*e_r + xl*e_theta
    at the point (1,0): e_r = (1,0), e_theta = (0,1)  =>  dX = yl, dY = xl."""
    H = 0
    for i in range(5):
        d = sp.diff(f, X, i, Y, 4 - i).subs({X: 1, Y: 0})
        H += sp.simplify(d) / (sp.factorial(i) * sp.factorial(4 - i)) * yl**i * xl**(4 - i)
    return sp.expand(H)


def coeffs(H):
    P = sp.Poly(H, xl, yl)
    return [P.coeff_monomial(xl**j * yl**(4 - j)) for j in range(5)]


def solve(m, rho, mode):
    r = sp.sqrt(X**2 + Y**2)
    e = ((X + sp.I * Y) / r)**m
    Hs = [coeffs(tensor_poly((r - 1)**j * e)) for j in range(5)]
    Hstar = sp.expand(xl**4 - (xl - yl / rho)**4)
    if mode == 'B':
        basis = [Hstar]
    elif mode == 'K':
        basis = [yl**4, Hstar]
    else:
        basis = [xl**4]
    c = sp.symbols('c0:5'); t = sp.symbols('t0:%d' % len(basis))
    expr = [sum(c[j] * Hs[j][i] for j in range(5)) - sum(t[b] * coeffs(basis[b])[i] for b in range(len(basis)))
            for i in range(5)]
    M = sp.Matrix([[sp.diff(q, v) for v in list(c) + list(t)] for q in expr])
    ns = M.nullspace()
    return ns, Hs


if __name__ == "__main__":
    m = int(sys.argv[1]); rho = sp.Rational(sys.argv[2])
    for mode in ('B', 'K', 'Q'):
        ns, Hs = solve(m, rho, mode)
        out = []
        for v in ns:
            v = [sp.nsimplify(sp.simplify(z)) for z in v]
            out.append([str(z) for z in v])
        print(mode, len(ns), json.dumps(out))
