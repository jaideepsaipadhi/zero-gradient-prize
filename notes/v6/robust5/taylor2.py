"""Sharp ball enclosures of bivariate polynomials with rational coefficients over boxes (Taylor form at the centre).
p(x0+u, y0+w) = sum c'_ij u^i w^j, computed EXACTLY (flint fmpq_mpoly.compose, x0, y0 rational);
then p(box) is contained in [c'_00 - rho, c'_00 + rho], rho = sum_{(i,j) != 0} |c'_ij| rx^i ry^j (exact rational)."""
import sympy as sp
import flint
from flint import arb, fmpq

CTX = flint.fmpq_mpoly_ctx.get(('u', 'w'), 'lex')
U, W = CTX.gens()


def to_arb(q):
    return arb(q.p) / arb(q.q)


class P2:
    def __init__(self, expr, x, y):
        pol = sp.Poly(sp.expand(expr), x, y)
        self.terms = [(a, b, sp.Rational(c)) for (a, b), c in pol.terms()]
        self.fterms = [(a, b, float(c)) for a, b, c in self.terms]
        self.mp = CTX.from_dict({(a, b): fmpq(int(c.p), int(c.q)) for a, b, c in self.terms})

    def f(self, x, y):
        return sum(c * x**a * y**b for a, b, c in self.fterms)

    def exact_range(self, x0, y0, rx, ry):
        """x0, y0, rx, ry: fmpq.  Returns (centre value, rho) as fmpq."""
        sh = self.mp.compose(U + x0, W + y0)
        c00 = fmpq(0); rho = fmpq(0)
        for (i, j), c in sh.to_dict().items():
            if i == 0 and j == 0:
                c00 = c
            else:
                rho += abs(c) * rx**i * ry**j
        return c00, rho

    def enclose(self, x0, y0, rx, ry):
        c00, rho = self.exact_range(x0, y0, rx, ry)
        return to_arb(c00) + arb(0, 1) * to_arb(rho)
