"""Independent check (sympy, monomials, no BB integration formulas) of robust5's reference fields UREF on
T^=(0,0),(1,0),(0,1) with (z,a,c) = vertices 0,1,2.  Correct pressure-blind condition for u vanishing on the two
edges other than e={y=0}:  (u, grad r)_T^ = 0 for all r in P3."""
import sys; sys.path.insert(0, '../robust5')
import sympy as sp
from math import factorial as f
import pt_fields as PF
x, y = sp.symbols('x y')
L = [1 - x - y, x, y]
def bb(coef, n):
    return sum(sp.Rational(c) * sp.Rational(f(n), f(a[0]) * f(a[1]) * f(a[2])) * L[0]**a[0] * L[1]**a[1] * L[2]**a[2]
               for a, c in coef.items())
def intT(p):
    return sp.integrate(sp.integrate(sp.expand(p), (y, 0, 1 - x)), (x, 0, 1))
for i, ur in enumerate(PF.UREF):
    u = [bb({a: c for (cp, a), c in ur.items() if cp == comp}, 4) for comp in range(2)]
    u = [sp.expand(w) for w in u]
    # boundary: vanish on x=0 and x+y=1
    bc = [sp.expand(w.subs(x, 0)) for w in u] + [sp.expand(w.subs(y, 1 - x)) for w in u]
    res = [intT(u[0] * sp.diff(r, x) + u[1] * sp.diff(r, y)) for r in [x**p * y**q for p in range(4) for q in range(4 - p)]]
    tr = sp.expand(-u[1].subs(y, 0))  # v.n with n=(0,-1)
    # Bernstein on e in s = x (lam_a): B_j = C(4,j) s^j (1-s)^(4-j)
    s = x
    B = [sp.binomial(4, j) * s**j * (1 - s)**(4 - j) for j in range(5)]
    print("field", i + 1, "BC zero:", all(b == 0 for b in bc), " (u,grad r) residuals:", res)
    print("   trace v.n =", sp.factor(tr), " B1-B2:", sp.expand(tr - (B[1] - B[2])) == 0, " B1-B3:", sp.expand(tr - (B[1] - B[3])) == 0)
    # the buggy constraint: <v.n, r~>_e with r~ = B^3_r restricted ignoring lambda_c
    print("   int_e tr * s^j, j=0..3:", [sp.integrate(tr * s**j, (x, 0, 1)) for j in range(4)])
