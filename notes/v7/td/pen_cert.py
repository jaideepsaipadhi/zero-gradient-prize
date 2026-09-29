"""Single-cell penalty necessity (Theorem pen:single) over a SHAPE CLASS, by interval arithmetic (python-flint arb).

Shape space (mod similarity): p1=(0,0,0), p2=(1,0,0), p3=(a,b,0), apex p0=(x,y,z); F=conv(p1,p2,p3) on Gamma_h.
For a fixed rational witness v^ in W_3(T^) (Alfeld, zero on the 3 faces through the apex) the pulled-back forms are
   A_T = (1/(2 det J)) sum K_{ab,cd}(J) SG[ab,cd],   B_T = -(1/det J) sum (J^T J)_{ac} w_d SB[a,cd],  w = J^{-1}J^{-T}e3,
   C_T = (|J^{-T}e3|/det J) sum (J^T J)_{ab} SC[a,b],
with SG = int_T^ G(x)G, SB = int_F^ v^(x)G, SC = int_F^ v^(x)v^  (exact rationals, G = grad^ v^).
lambda_T >= diam F * (2B_T - A_T)/C_T at the witness.  A box of shapes is CERTIFIED if the arb enclosure of that
quotient is >= lambda0; otherwise it is bisected.  Each box gets its own witness (float optimum at the box centre,
rounded to a rational vector; the enclosure is then rigorous for that rational witness).
usage: python3 pen_cert.py lambda0 delta [maxboxes]
"""
import sys, time
from fractions import Fraction as Fr
import numpy as np
import flint
import refalfeld as ra

flint.ctx.prec = 80
k = 3
t0 = time.time()
Al = ra.Alfeld(k, 'apex')
rk, N = ra.nullspace(Al.div_rows(), Al.nunk)
d = len(N)
print(f"W_3(T^): unknowns={Al.nunk}, rank(div)={rk}, dim W={d}", flush=True)

# ---------------- exact tensors per basis pair
def field_D(v):
    """per sub s, beta -> 3x3 Fraction matrix D (grad v = k sum B_beta D_beta)."""
    out = []
    for s in range(4):
        dd = {}
        for beta in Al.bk1:
            Dl = Al.Dlin(s, beta)
            dd[beta] = [[sum(c * v[u] for u, c in Dl.get((a, b), {}).items()) for b in range(3)] for a in range(3)]
        out.append(dd)
    return out

Ds = [field_D(v) for v in N]
Ms = [Al.mass_k1(s) for s in range(4)]
tr = Al.face_trace()
from itertools import product
def face_vals(v):
    return {ap: [v[3 * n + a] for a in range(3)] for ap, n in tr.items()}
FV = [face_vals(v) for v in N]
bk1f = ra.mindex(k - 1, 3)
TG = np.empty((d, d), object); TB = np.empty((d, d), object); TC = np.empty((d, d), object)
for i in range(d):
    for j in range(d):
        SG = [[Fr(0)] * 9 for _ in range(9)]
        for s in range(4):
            M = Ms[s]
            for b1 in Al.bk1:
                D1 = Ds[i][s][b1]
                for b2 in Al.bk1:
                    m = M[(b1, b2)] * k * k
                    if m == 0: continue
                    D2 = Ds[j][s][b2]
                    for p in range(9):
                        x1 = D1[p // 3][p % 3]
                        if x1 == 0: continue
                        for q in range(9):
                            x2 = D2[q // 3][q % 3]
                            if x2: SG[p][q] += m * x1 * x2
        SB = [[Fr(0)] * 9 for _ in range(3)]
        for ap, vv in FV[i].items():
            for bp in bk1f:
                Dj = Ds[j][0][(0,) + bp]
                c = k * ra.face_int(k, ap, k - 1, bp)
                for a in range(3):
                    if vv[a] == 0: continue
                    for q in range(9):
                        if Dj[q // 3][q % 3]: SB[a][q] += c * vv[a] * Dj[q // 3][q % 3]
        SC = [[Fr(0)] * 3 for _ in range(3)]
        for ap, vi in FV[i].items():
            for bp, vj in FV[j].items():
                c = ra.face_int(k, ap, k, bp)
                for a in range(3):
                    for b in range(3):
                        if vi[a] and vj[b]: SC[a][b] += c * vi[a] * vj[b]
        TG[i, j], TB[i, j], TC[i, j] = SG, SB, SC
print(f"tensors done ({time.time()-t0:.1f}s)", flush=True)
TGf = np.array([[np.array(TG[i, j], float) for j in range(d)] for i in range(d)])   # d,d,9,9
TBf = np.array([[np.array(TB[i, j], float) for j in range(d)] for i in range(d)])   # d,d,3,9
TCf = np.array([[np.array(TC[i, j], float) for j in range(d)] for i in range(d)])   # d,d,3,3


def geo(J, sqrt):
    """J as 3x3 nested list (any ring). Returns K (9x9), JtJ, w, |J^-T e3|^2 pieces for the formulas."""
    # inverse of upper-triangular-like J = [[1,a,x],[0,b,y],[0,0,z]] computed generally via adjugate
    a, b_, c = J[0]; dd, e, f = J[1]; g, h, i = J[2]
    det = a * (e * i - f * h) - b_ * (dd * i - f * g) + c * (dd * h - e * g)
    adj = [[e * i - f * h, c * h - b_ * i, b_ * f - c * e],
           [f * g - dd * i, a * i - c * g, c * dd - a * f],
           [dd * h - e * g, b_ * g - a * h, a * e - b_ * dd]]
    Ji = [[adj[r][s] / det for s in range(3)] for r in range(3)]
    return det, Ji


def forms_float(p):
    a, b, x, y, z = p
    J = np.array([[1, a, x], [0, b, y], [0, 0, z]], float)
    Ji = np.linalg.inv(J); det = np.linalg.det(J)
    T = np.einsum('ia,bj->ijab', J, Ji); T = T + T.transpose(1, 0, 2, 3)          # T^{ij}_{ab}
    K = np.einsum('ijab,ijcd->abcd', T, T).reshape(9, 9)
    A = np.einsum('pq,xypq->xy', K, TGf) / (2 * det)
    JtJ = J.T @ J; w = Ji @ Ji.T @ np.array([0, 0, 1.0])
    Bm = -np.einsum('ac,d,xyacd->xy', JtJ, w, TBf.reshape(d, d, 3, 3, 3)) / det
    B = (Bm + Bm.T) / 2
    nrm = np.linalg.norm(Ji.T @ np.array([0, 0, 1.0]))
    C = np.einsum('ab,xyab->xy', JtJ, TCf) * nrm / det
    diamF = max(1.0, np.hypot(a, b), np.hypot(a - 1, b))
    return A, B, C, diamF


def best_witness(p):
    A, B, C, dF = forms_float(p)
    Q = 2 * B - A
    lo, hi = -1e3, 1e4
    f = lambda t: np.linalg.eigvalsh(Q - t * C)[-1]
    for _ in range(100):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    w, V = np.linalg.eigh(Q - (lo - 1e-7) * C)
    return lo * dF, V[:, -1]


def contract(qv):
    q = [Fr(float(t)).limit_denominator(10 ** 6) for t in qv]
    SG = [[sum(q[i] * q[j] * TG[i, j][r][s] for i in range(d) for j in range(d)) for s in range(9)] for r in range(9)]
    SB = [[sum(q[i] * q[j] * TB[i, j][r][s] for i in range(d) for j in range(d)) for s in range(9)] for r in range(3)]
    SC = [[sum(q[i] * q[j] * TC[i, j][r][s] for i in range(d) for j in range(d)) for s in range(3)] for r in range(3)]
    A = lambda M: [[flint.arb(Fr(v).numerator) / Fr(v).denominator for v in row] for row in M]
    return A(SG), A(SB), A(SC)


class Du:
    """forward-mode dual number over arb with n partials (for the mean-value form)."""
    n = 5
    def __init__(self, v, g=None):
        self.v = v if isinstance(v, flint.arb) else flint.arb(v)
        self.g = g if g is not None else [flint.arb(0)] * Du.n
    def _c(o):
        return o if isinstance(o, Du) else Du(o)
    def __add__(s, o):
        o = Du._c(o); return Du(s.v + o.v, [a + b for a, b in zip(s.g, o.g)])
    __radd__ = __add__
    def __sub__(s, o):
        o = Du._c(o); return Du(s.v - o.v, [a - b for a, b in zip(s.g, o.g)])
    def __rsub__(s, o):
        return Du._c(o) - s
    def __neg__(s):
        return Du(-s.v, [-a for a in s.g])
    def __mul__(s, o):
        o = Du._c(o); return Du(s.v * o.v, [s.v * b + o.v * a for a, b in zip(s.g, o.g)])
    __rmul__ = __mul__
    def __truediv__(s, o):
        o = Du._c(o); q = s.v / o.v
        return Du(q, [(a - q * b) / o.v for a, b in zip(s.g, o.g)])
    def __rtruediv__(s, o):
        return Du._c(o) / s
    def __pow__(s, m):
        r = Du(1)
        for _ in range(m): r = r * s
        return r
    def sqrt(s):
        r = s.v.sqrt(); return Du(r, [a / (2 * r) for a in s.g])


def quotient_expr(prm, S):
    """the witness quotient (2B-A)/C (without diam F) as a function of the 5 shape parameters (arb or Du)."""
    SG, SB, SC = S
    a, b, x, y, z = prm
    one, zero = 1, 0
    J = [[one, a, x], [zero, b, y], [zero, zero, z]]
    det, Ji = geo(J, None)
    A2 = 0
    for i in range(3):
        for j in range(3):
            T = [J[i][p // 3] * Ji[p % 3][j] + J[j][p // 3] * Ji[p % 3][i] for p in range(9)]
            for p in range(9):
                u = 0
                for q in range(9):
                    if SG[p][q] != 0: u = u + T[q] * SG[p][q]
                A2 = A2 + T[p] * u
    JtJ = [[J[0][aa] * J[0][cc] + J[1][aa] * J[1][cc] + J[2][aa] * J[2][cc] for cc in range(3)] for aa in range(3)]
    e3t = [Ji[2][0], Ji[2][1], Ji[2][2]]
    w = [Ji[r][0] * e3t[0] + Ji[r][1] * e3t[1] + Ji[r][2] * e3t[2] for r in range(3)]
    Bd = 0
    for aa in range(3):
        for cc in range(3):
            for dd in range(3):
                if SB[aa][3 * cc + dd] != 0: Bd = Bd - JtJ[aa][cc] * w[dd] * SB[aa][3 * cc + dd]
    nrm = (e3t[0] * e3t[0] + e3t[1] * e3t[1] + e3t[2] * e3t[2]).sqrt()
    Cd = 0
    for aa in range(3):
        for bb in range(3):
            if SC[aa][bb] != 0: Cd = Cd + JtJ[aa][bb] * SC[aa][bb]
    Cd = Cd * nrm
    return (2 * Bd - A2 / 2) / Cd


def quotient_arb(box, S):
    """rigorous lower bound of diamF*(2B-A)/C over the box, mean-value form."""
    ctr = [flint.arb(m) for (m, r) in box]
    fc = quotient_expr(ctr, S)                                     # value at the centre (tight ball)
    prm = [Du(flint.arb(m, r), [flint.arb(1 if i == j else 0) for i in range(5)]) for j, (m, r) in enumerate(box)]
    fb = quotient_expr(prm, S)                                     # gradient enclosure over the box
    enc = fc
    for gi, (m, r) in zip(fb.g, box):
        enc = enc + gi * flint.arb(0, r)
    enc = flint.arb(enc)
    a = flint.arb(box[0][0], box[0][1]); b = flint.arb(box[1][0], box[1][1])
    d2 = (a * a + b * b).sqrt(); d3 = ((a - 1) ** 2 + b * b).sqrt()
    lo = float(enc.lower())
    if lo > 0:
        return lo * max(1.0, float(d2.lower()), float(d3.lower())), enc
    return lo * max(1.0, float(d2.upper()), float(d3.upper())), enc


if __name__ == "__main__":
    lam0 = float(sys.argv[1]) if len(sys.argv) > 1 else 20.0
    delta = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
    maxb = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
    c0 = [0.5, np.sqrt(3) / 2, 0.5, np.sqrt(3) / 6, np.sqrt(2 / 3)]        # regular tetrahedron
    # sanity check against the tabulated near-regular value lambda_T = 32.50
    lt, _ = best_witness([0.5, 7 / 8, 0.5, 7 / 24, 13 / 16])
    print(f"check: near-regular lambda_T (float) = {lt:.4f}  (table: 32.50)")
    lt, _ = best_witness(c0)
    print(f"regular tetrahedron lambda_T (float) = {lt:.4f}", flush=True)
    stack = [[(c, delta) for c in c0]]
    done = 0; worst = 1e9; nb = 0; wit = {}
    while stack and nb < maxb:
        box = stack.pop(); nb += 1
        centre = [m for m, r in box]
        lt, qv = best_witness(centre)
        S = contract(qv)
        lb, _ = quotient_arb(box, S)
        if lb >= lam0:
            done += 1; worst = min(worst, lb)
            continue
        if lt < lam0:
            print(f"FAIL: float lambda_T={lt:.3f} < lambda0 at {centre}"); break
        _, enc = quotient_arb(box, S)
        prm = [Du(flint.arb(m, r), [flint.arb(1 if i == jj else 0) for i in range(5)]) for jj, (m, r) in enumerate(box)]
        gb = quotient_expr(prm, S).g
        j = int(np.argmax([float(abs(gi).upper()) * r for gi, (m, r) in zip(gb, box)]))
        m, r = box[j]
        for s in (-1, 1):
            nbx = list(box); nbx[j] = (m + s * r / 2, r / 2); stack.append(nbx)
    print(f"lambda0={lam0} delta={delta}: boxes processed={nb}, certified={done}, remaining={len(stack)}, "
          f"min certified lower bound={worst:.4f}  ({time.time()-t0:.1f}s)")
    if not stack:
        print(f"CERTIFIED: lambda_T >= {lam0} for every tetrahedron with p1=0,p2=e1,p3=(a,b,0),p0=(x,y,z) and "
              f"|(a,b,x,y,z)-regular|_inf <= {delta}")
