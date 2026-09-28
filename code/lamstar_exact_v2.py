import sympy as sy
from lamstar_exact import run
lam0 = sy.Rational(9, 1)
run([(1, 0), (1, 1), (0, sy.Rational(6, 5)), (-1, 1), (-1, 0)], lam0)   # S4': apex moved off y=1, nonsingular
