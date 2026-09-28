"""Covariance matrix of the edge bubbles mu_i mu_j under the weight B=(mu1 mu2 mu3)^2 on a triangle
(area-normalised: integrals divided by |F|).  Used in the k>=9 single-tetrahedron witness."""
import sympy as sp
from math import factorial
def I(a):  # (1/|F|) int_F mu^a
    return sp.Rational(2 * factorial(a[0]) * factorial(a[1]) * factorial(a[2]), factorial(sum(a) + 2))
e = [(1, 1, 0), (1, 0, 1), (0, 1, 1)]
add = lambda *v: tuple(sum(t) for t in zip(*v))
B = (2, 2, 2)
m0 = I(B)
V = sp.Matrix(3, 3, lambda i, j: I(add(B, e[i], e[j])) - I(add(B, e[i])) * I(add(B, e[j])) / m0)
print('int B / |F| =', m0)
print('V =', V)
print('eigenvalues:', V.eigenvals())
