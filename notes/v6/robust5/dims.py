"""Exact dimension counts (robust5 REPORT, Sec. 2.1): fields on the two-triangle patch Q = T1 u T2 and on a 3-triangle
fan, with per-edge (M), (D); normal-trace rank on e.  Random rational shapes, k = 4, 5."""
import random
from fractions import Fraction as Fr
import fan
random.seed(1)
for trial in range(4):
    c = (Fr(random.randint(2, 8), 10), Fr(random.randint(6, 12), 10))
    p2 = (Fr(-random.randint(0, 6), 10), Fr(random.randint(8, 12), 10))
    p3 = (Fr(-random.randint(8, 12), 10), Fr(random.randint(-1, 1), 10))
    for P, lg, name in [([(0, 0), (1, 0), c, p2], False, "Q = T1 u T2"), ([(0, 0), (1, 0), c, p2, p3], True, "3-star, v=0 on e'")]:
        for k in (4, 5):
            F = fan.Fan(P, k=k, last_is_gamma=lg); B, r = F.nullspace()
            print("trial %d  %-18s k=%d  dim=%d  trace rank=%d  (dim P_e=%d)" % (trial, name, k, len(B), fan.trmat(F, B).rank(), k - 2))
for k in (4, 5, 6):
    F, B = fan.U_space([(0, 0), (1, 0), (0, 1)], k)
    print("reference U(T^), k=%d: dim %d, trace rank %d" % (k, len(B), fan.trmat(F, B).rank()))
