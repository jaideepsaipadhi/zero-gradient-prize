"""Rigorous per-edge inverse trace constant for P_3 on triangles:  ||u||_e^2 <= c3 (|e|/|T|) ||u||_T^2  for all u in P_3(T).
The ratio is affine invariant, so c3 = lambda_max(M_e, M_T) * |T^|/|e^| on the reference (exact charpoly + rigorous roots)."""
import flint
from fractions import Fraction as Fr
import pt_fields as PF
I3 = PF.I3
MT = flint.fmpq_mat(10, 10); Me = flint.fmpq_mat(10, 10)
for i, a in enumerate(I3):
    for j, b in enumerate(I3):
        v = PF.int_T(1, a, b, 3, 3)          # times area2 = 1 for T^
        MT[i, j] = flint.fmpq(v.numerator, v.denominator)
        if a[2] == 0 and b[2] == 0:
            w = PF.int_e(a, b, 3, 3); Me[i, j] = flint.fmpq(w.numerator, w.denominator)
A = MT.inv() * Me
cp = A.charpoly()
num, den = cp.numer_denom() if hasattr(cp, 'numer_denom') else (cp, 1)
zp = flint.fmpz_poly([int(x) for x in (cp * cp.denom() if hasattr(cp,'denom') else cp).coeffs()]) if hasattr(cp,'coeffs') else None
roots = zp.complex_roots()
lam = max((r[0].real for r in roots), key=lambda x: float(x.mid()))
print("charpoly:", cp)
print("lambda_max =", lam, " c3 = lambda_max * |T^|/|e^| =", lam / 2)
