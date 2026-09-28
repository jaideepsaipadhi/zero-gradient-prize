"""referee11: independent exact check of (i) the T1 Gram closed forms (robust6 closed_forms.pkl), (ii) the T2 factorisation
gC_f gC_g (x1^2+1) F(x2,y2) with robust5's F, (iii) the pairing matrix P_W and annihilation of Ker_T, using my own
integrators from verify_PT.py and my own L2 projection onto P3(T1)."""
import sys, pickle
sys.path.insert(0, '../robust5'); sys.path.insert(0, '../robust6')
from fractions import Fraction as Fr
import sympy as sp
from verify_PT import *
import pt_fields as PF, uref_fix, cert_PT as C
CF = pickle.load(open('../robust6/closed_forms.pkl', 'rb'))
def solve(A, b):
    n = len(b); M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for i in range(n):
        p = next(r for r in range(i, n) if M[r][i] != 0); M[i], M[p] = M[p], M[i]
        for r in range(n):
            if r != i and M[r][i] != 0:
                f = M[r][i] / M[i][i]; M[r] = [u - f * w for u, w in zip(M[r], M[i])]
    return [M[i][n] / M[i][i] for i in range(n)]
def run(x1, y1, x2, y2):
    PF.UREF = uref_fix.get()
    G, V, info = PF.star_fields(x1, y1, x2, y2, PF.NumExact)
    s1, s2 = x1 + y1, x2 + y2
    z = (Fr(0), Fr(0)); a = (Fr(1), Fr(0)); c = (x1 / s1, 1 / s1); b = ((x2 * c[0] - c[1]) / s2, (x2 * c[1] + c[0]) / s2)
    L1 = bary(z, a, c); L2 = bary(z, c, b)
    F1 = [[bb_to_poly(F[1][d], L1) for d in range(2)] for F in V]
    F2 = [[bb_to_poly(F[2][d], L2) for d in range(2)] for F in V]
    ok = True
    sub = {sp.Symbol('x1'): sp.Rational(x1.numerator, x1.denominator), sp.Symbol('y1'): sp.Rational(y1.numerator, y1.denominator)}
    for f in range(2):
        for g in range(f, 2):
            gT1 = sum(int_tri(padd(pmul(dx(p), dx(q)), pmul(dy(p), dy(q))), z, a, c) for p, q in zip(F1[f], F1[g]))
            gT2 = sum(int_tri(padd(pmul(dx(p), dx(q)), pmul(dy(p), dy(q))), z, c, b) for p, q in zip(F2[f], F2[g]))
            ee = sum(int_seg_param(pmul(p, q), z, a) for p, q in zip(F1[f], F1[g]))
            dd = sum(int_seg_param(pmul(dy(p), dy(q)), z, a) for p, q in zip(F1[f], F1[g]))
            cf = [sp.Rational(q.subs({sy: (sub[sp.Symbol("x1")] if sy.name == "x1" else sub[sp.Symbol("y1")]) for sy in q.free_symbols})) for q in CF["Q1"][(f, g)]]
            gc = [Fr(1, 5), Fr(2, 5)]
            e2 = gc[f] * gc[g] * (x1 * x1 + 1) * C.Ffun(x2, y2)
            r = [Fr(str(cf[0])) == gT1, Fr(str(cf[1])) == ee, Fr(str(cf[2])) == dd, e2 == gT2]
            print(f"  ({f},{g}) T1grad/edge/dn closed forms, T2 factorisation:", r); ok &= all(r)
    # L2 projection onto P3(T1) and pairing
    mons = [{(p, q): Fr(1)} for p in range(4) for q in range(4 - p)]
    A = [[int_tri(pmul(m, n), z, a, c) for n in mons] for m in mons]
    def proj(H):
        co = solve(A, [int_tri(pmul(H, m), z, a, c) for m in mons]); o = {}
        for cc, m in zip(co, mons): o = padd(o, pscale(m, cc))
        return o
    la = {(1, 0): Fr(1), (0, 1): -x1}; lz = {(1, 0): Fr(-1), (0, 1): -y1}; lc = pscale(padd(la, lz), -1)
    Ms = [pmul(lz, ppow(la, 3)), pmul(ppow(lz, 2), ppow(la, 2))]; Ks = [ppow(lz, 4), ppow(la, 4), ppow(lc, 4)]
    def pair(H, f): return int_seg_param(pmul(padd(H, pscale(proj(H), -1)), pscale(F1[f][1], -1)), z, a)
    P = [[pair(m, f) for m in Ms] for f in range(2)]
    Kp = [pair(k, f) for k in Ks for f in range(2)]
    r = [P == [[0, Fr(-1, 3780)], [Fr(-1, 2520), 0]], all(q == 0 for q in Kp)]
    print("  P_W =", P, " Ker_T pairings all 0:", r[1]); ok &= all(r)
    print("shape", x1, y1, x2, y2, "GRAM/PAIRING OK" if ok else "FAIL", flush=True)
if __name__ == "__main__":
    run(*[Fr(q) for q in sys.argv[1:5]])
