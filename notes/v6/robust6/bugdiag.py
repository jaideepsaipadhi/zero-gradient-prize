"""Diagnose: do robust5's UREF fields satisfy the BUGGY constraint of fan.py / sym_verify.py
(edge term computed with int_e(al, r) even when r has lambda_c-exponent > 0, where B^3_r vanishes on e)?"""
import sys; sys.path.insert(0, '../robust5')
import sympy as sp
from math import factorial as f
from fractions import Fraction as Fr
import pt_fields as PF
x, y = sp.symbols('x y')
L = [1 - x - y, x, y]
def B(a, n): return sp.Rational(f(n), f(a[0]) * f(a[1]) * f(a[2])) * L[0]**a[0] * L[1]**a[1] * L[2]**a[2]
def intT(p): return sp.integrate(sp.integrate(sp.expand(p), (y, 0, 1 - x)), (x, 0, 1))
for i, ur in enumerate(PF.UREF):
    cf = [{a: c for (cp, a), c in ur.items() if cp == comp} for comp in range(2)]
    u = [sum(sp.Rational(c) * B(a, 4) for a, c in cf[comp].items()) for comp in range(2)]
    div = sp.diff(u[0], x) + sp.diff(u[1], y)
    good, bad = [], []
    for r in PF.I3:
        dterm = intT(div * B(r, 3))
        corr_edge = sp.integrate((-u[1] * B(r, 3)).subs(y, 0), (x, 0, 1))
        bug_edge = sum(-sp.Rational(cf[1].get(al, 0)) * sp.Rational(PF.int_e(al, r, 4, 3)) for al in PF.I4 if al[2] == 0)
        good.append(dterm - corr_edge); bad.append(dterm - bug_edge)
    print("field", i + 1, "\n  correct residuals:", good, "\n  buggy residuals  :", bad)
