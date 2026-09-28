"""Corrected reference fields for robust5's (P^T) construction (k=4): u^_1, u^_2 in U(T^) with normal traces
B_1-B_2 and B_1-B_3, returned as BB coefficient dicts {(comp, alpha): Fraction} in robust5's pt_fields format
(alpha = (alpha_z, alpha_a, alpha_c), T^ = z(0,0), a(1,0), c(0,1)).  The free curl direction is fixed by setting
the free nullspace coordinate to zero; the final v_f do not depend on this choice (w_A spans that direction)."""
import sympy as sp
from fractions import Fraction as Fr
from math import factorial as f
from uref_correct import run
x, y = sp.symbols('x y')
L = [1 - x - y, x, y]
I4 = [(a, b, 4 - a - b) for a in range(4, -1, -1) for b in range(4 - a, -1, -1)]
def B(a, n=4): return sp.Rational(f(n), f(a[0]) * f(a[1]) * f(a[2])) * L[0]**a[0] * L[1]**a[1] * L[2]**a[2]
def to_bb(p):
    cs = sp.symbols('c0:15')
    expr = sp.expand(sum(c * B(a) for c, a in zip(cs, I4)) - p)
    sol = sp.solve(sp.Poly(expr, x, y).coeffs(), cs, dict=True)[0]
    return {a: sol[c] for c, a in zip(cs, I4)}
def get():
    fields, traces = run(4)
    s = x
    Bs = [sp.binomial(4, j) * s**j * (1 - s)**(4 - j) for j in range(5)]
    targets = [sp.expand(Bs[1] - Bs[2]), sp.expand(Bs[1] - Bs[3])]
    cs = sp.symbols('k0:3')
    out = []
    for tg in targets:
        comb0 = sum(c * fl[0] for c, fl in zip(cs, fields)); comb1 = sum(c * fl[1] for c, fl in zip(cs, fields))
        eqs = sp.Poly(sp.expand(-comb1.subs(y, 0) - tg), x).coeffs()
        sol = sp.solve(eqs, cs, dict=True)[0]
        free = [c for c in cs if c not in sol]
        sol = {c: sp.sympify(v).subs({q: 0 for q in free}) for c, v in sol.items()}
        sol.update({q: 0 for q in free})
        u0 = sp.expand(comb0.subs(sol)); u1 = sp.expand(comb1.subs(sol))
        d = {}
        for comp, u in ((0, u0), (1, u1)):
            for a, v in to_bb(u).items():
                if v != 0:
                    d[(comp, a)] = Fr(int(sp.fraction(v)[0]), int(sp.fraction(v)[1]))
        out.append(d)
    return out
if __name__ == "__main__":
    for d in get():
        print(d)
