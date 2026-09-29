"""notes/v7/robust: EXACT leading coefficient of the UP pairing for phi = (r-1)^4 (companion of pairing_asym.py).
Blow-up at z = (1,0): X = 1 + l xi, Y = l ups.  To first order in l:  (r-1)^4 / l^4 = xi^4 + 2 l xi^3 ups^2 + O(l^2),
T_e(l) = {(0,0), (-l/2, 1), (rho - l/2, 1 + rho l)} + O(l^2) (exact expansion of the UP vertices).
P_f / l^5 = int_0^1 [(I - pi_{T(l)}) phihat](z + u(a-z)) (B_1 - B_{f+1})(u) du =: p_f(l),  p_f(0) = 0 (xi^4 in Ker_T), and
c_f = lim P_f/l^6 = p_f'(0), computed exactly in Q(rho)."""
import sympy as sp
from math import factorial
l, rho, u, w, xi, up = sp.symbols('l rho u w xi ups')
import sys
k = int(sys.argv[1]) if len(sys.argv) > 1 else 4
z = sp.Matrix([0, 0]); a = sp.Matrix([-l / 2, 1]); c = sp.Matrix([rho - l / 2, 1 + rho * l])
J = sp.Matrix.hstack(a - z, c - z)
phihat = xi**k + sp.Rational(k, 2) * l * xi**(k - 1) * up**2      # (r-1)^k / l^k to first order
ph = sp.expand(phihat.subs({xi: J[0, 0] * u + J[0, 1] * w, up: J[1, 0] * u + J[1, 1] * w}, simultaneous=True))
P3 = [u**i * w**j for i in range(k) for j in range(k - i)]
def I(p):
    p = sp.Poly(sp.expand(p), u, w)
    return sum(cf * sp.Rational(factorial(i) * factorial(j), factorial(i + j + 2)) for (i, j), cf in p.terms())
G = sp.Matrix([[I(p * q) for q in P3] for p in P3])          # reference Gram (the Jacobian cancels)
b = sp.Matrix([I(ph * p) for p in P3])
cf = G.LUsolve(b)
eta_e = sp.expand((ph - sum(ci * p for ci, p in zip(cf, P3))).subs(w, 0))
B = [sp.binomial(k, i) * u**i * (1 - u)**(k - i) for i in range(k + 1)]
tests = [B[1] - B[f] for f in range(2, k)]          # a basis of P_e (zero ends, zero mean), k-2 elements
print("k =", k)
for f, tf in enumerate(tests, 1):
    pf = sp.integrate(sp.expand(eta_e * tf), (u, 0, 1))
    print("f=%d: p_f(0) =" % f, sp.simplify(pf.subs(l, 0)), "  c_f = p_f'(0) =", sp.factor(sp.diff(pf, l).subs(l, 0)))
