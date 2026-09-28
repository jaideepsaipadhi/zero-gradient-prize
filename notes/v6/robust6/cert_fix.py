"""Rigorous interval certificate for (P^T), k = 4, with the CORRECTED fields (robust6 REPORT, Sec. 1).
Same method and code path as robust5/cert_PT.py (python-flint arb balls; exact Taylor-form enclosures of the
polynomial pieces over exact rational boxes; covering branch and bound), but with
  * the corrected T1 Gram closed forms (closed_forms.pkl, checked independently in check_forms.py), and
  * the corrected gamma_C = (1/5, 2/5) (robust5 hard-coded (-1/20, -1/10) in three places).
F(x2, y2) (the T2 part) and the Frobenius-distance Gram are unchanged (they do not involve the reference fields).
usage: python cert_fix.py THETA_DEG KAPPA0"""
import sys, math, time, pickle
sys.path.insert(0, '../robust5')
from fractions import Fraction as Fr
import sympy as sp
import flint
from flint import arb, fmpq
import cert_PT as C
from taylor2 import P2
CF = pickle.load(open('closed_forms.pkl', 'rb'))
X, Y = sp.symbols('x1 y1')
def _ren(e): return e.subs({s: (X if s.name == 'x1' else Y) for s in e.free_symbols})
Q1 = {k: tuple(_ren(q) for q in v) for k, v in CF["Q1"].items()}
gC = [sp.nsimplify(g) for g in CF["gC"]]
assert gC == [sp.Rational(1, 5), sp.Rational(2, 5)], gC
S1 = X + Y; XX = X**2 + 1
DL = (12*X**4*Y**4 + 8*X**4*Y**2 + 2*X**4 - 16*X**3*Y**3 - 4*X**3*Y + 8*X**2*Y**4 + 18*X**2*Y**2
      + 3*X**2 - 4*X*Y**3 - 6*X*Y + 2*Y**4 + 3*Y**2 + 2)
a_, b_, c_ = C.GWperp_cf(X, Y)
A_, B_, C_ = (sp.cancel(q * DL) for q in (a_, b_, c_))
def comb(q):   # tr(P_W^{-1} Q P_W^{-T} Gperp) * DL
    return 2520**2 * q[(1, 1)] * A_ + 2 * 2520 * 3780 * q[(0, 1)] * B_ + 3780**2 * q[(0, 0)] * C_
def build(Fmax):
    Fq = sp.Rational(str(Fmax))
    q0 = {k: 13 * (Q1[k][0] + gC[k[0]] * gC[k[1]] * XX * Fq) + 60 * S1 * Q1[k][1] for k in Q1}   # 3 C_I^2 = 60 s1
    q1 = {k: Q1[k][1] + 3 * Q1[k][2] for k in Q1}                                               # coefficient of C_T^2
    N0, D0 = sp.fraction(sp.cancel(sp.together(comb(q0))))
    N1, D1 = sp.fraction(sp.cancel(sp.together(comb(q1))))
    print("  denominators:", sp.factor(D0), "|", sp.factor(D1), flush=True)
    return dict(N0=P2(N0, X, Y), D0=P2(D0, X, Y), N1=P2(N1, X, Y), D1=P2(D1, X, Y), DL=P2(DL, X, Y), S1=P2(S1, X, Y),
                XX=P2(XX, X, Y), YY=P2(Y**2 + 1, X, Y)), (N0, D0, N1, D1)
def CT2_of(xx, yy, s1, PI, amax, asqrt, one):
    R2 = amax(xx, yy) / (s1 * s1); D2 = amax(R2, one)
    return 2 * s1 * (D2 / (PI * PI) + asqrt(R2 * D2) / PI)
def kappa2_box(P, x0, y0, rx, ry):
    E = {k: p.enclose(x0, y0, rx, ry) for k, p in P.items()}
    if E["S1"].lower() <= 0 or E["DL"].lower() <= 0 or E["D0"].lower() <= 0 or E["D1"].lower() <= 0:
        return None
    CT2 = CT2_of(E["XX"], E["YY"], E["S1"], arb.pi(), lambda p, q: p.max(q), lambda p: p.sqrt(), arb(1))
    return 1 / ((E["N0"] / E["D0"] + CT2 * E["N1"] / E["D1"]) / E["DL"])
def bb2(theta_deg, P, kappa0, maxboxes=3_000_000, minwidth=Fr(1, 2 ** 24)):
    th = arb.pi() * theta_deg / 180
    Cc = 1 / th.tan(); C2 = 1 / (2 * th).tan(); SIN = th.sin()
    lo = -Fr(math.ceil(float(C2.mid()) * 10**6) + 10, 10**6); hi = Fr(math.ceil(float(Cc.mid()) * 10**6) + 10, 10**6)
    k2 = (arb(kappa0) ** 2).upper()
    stack = [C.Box([lo, lo], [hi, hi])]; st = dict(accepted=0, discarded=0, split=0, min_lower=None); t0 = time.time()
    q = lambda f: fmpq(f.numerator, f.denominator)
    while stack:
        B = stack.pop()
        x = C.ball_of(B.lo[0], B.hi[0]); y = C.ball_of(B.lo[1], B.hi[1])
        if C.infeasible(x, y, Cc, SIN):
            st["discarded"] += 1; continue
        v = kappa2_box(P, q((B.lo[0] + B.hi[0]) / 2), q((B.lo[1] + B.hi[1]) / 2), q((B.hi[0] - B.lo[0]) / 2), q((B.hi[1] - B.lo[1]) / 2))
        if v is not None and v.is_finite() and v.lower() >= k2:
            st["accepted"] += 1
            if st["min_lower"] is None or v.lower() < st["min_lower"]: st["min_lower"] = v.lower()
            continue
        if B.hi[0] - B.lo[0] < minwidth and B.hi[1] - B.lo[1] < minwidth:
            raise RuntimeError("REFUTED/unresolved box %s %s value %s" % (B.lo, B.hi, v))
        st["split"] += 1
        if st["split"] > maxboxes: raise RuntimeError("too many boxes")
        stack.extend(B.split())
    st["time_s"] = round(time.time() - t0, 1); st["min_lower"] = str(st["min_lower"])
    return st
if __name__ == "__main__":
    th = float(sys.argv[1]); kappa0 = sys.argv[2]
    Cc = 1 / math.tan(math.radians(th)); C2 = 1 / math.tan(math.radians(2 * th))
    n = 400; best = 0
    pts = [(-C2 + (Cc + C2) * i / n, -C2 + (Cc + C2) * j / n) for i in range(n + 1) for j in range(n + 1)]
    adm = [(p, r) for p, r in pts if p + r > 0 and 1 - p * r <= Cc * (p + r) + 1e-12 and p <= Cc and r <= Cc]
    best = max(C.Ffun(p, r) for p, r in adm)
    Fmax = math.ceil(best * 1.02)
    print("theta = %g deg: float max of F = %.6g; certify F <= %d" % (th, best, Fmax), flush=True)
    print("  F certificate:", C.bb(th, C.Ffun, lambda v: v.upper() <= Fmax), flush=True)
    P, (N0, D0, N1, D1) = build(Fmax)
    fN0, fD0, fN1, fD1, fDL = (sp.lambdify((X, Y), e) for e in (N0, D0, N1, D1, DL))
    def kf(p, r):
        ct2 = CT2_of(p * p + 1, r * r + 1, p + r, math.pi, max, math.sqrt, 1.0)
        return (fDL(p, r) / (fN0(p, r) / fD0(p, r) + ct2 * fN1(p, r) / fD1(p, r))) ** 0.5
    vals = [(kf(p, r), p, r) for p, r in adm]
    print("  float min of the kappa bound on a 401^2 grid: %.4g at (x1,y1)=(%.4f,%.4f)" % min(vals), flush=True)
    st = bb2(th, P, kappa0)
    print("  CERTIFIED: kappa_T >= %s on all admissible (T1, T2) with min angle >= %g deg:" % (kappa0, th), st, flush=True)
