"""Exact (rational) Bernstein--Bezier toolkit for ONE Alfeld-split reference tetrahedron.

T^ = conv(P0,P1,P2,P3), P0=(0,0,1) apex, P1=0, P2=e1, P3=e2; Gamma-face F^ = conv(P1,P2,P3) in z=0.
Alfeld split at the barycentre B: K_i = conv(B, face opposite P_i); K_0 carries F^.
A field is given by vector Bernstein coefficients at the domain points (shared points = C^0).

zero='all'  : v = 0 on all of dT        (the spaces Y^_k of the sharpness section / the local Stokes space)
zero='apex' : v = 0 on the three faces through the apex (the space W_k(T) of the penalty section)

Everything here is exact (fractions.Fraction / sympy QQ).  Used by kappa_cert.py, leak_witness.py, pen_cert.py.
"""
from fractions import Fraction as Fr
from itertools import product
from math import factorial
from sympy.polys.matrices import DomainMatrix
from sympy import QQ

P = [(Fr(0), Fr(0), Fr(1)), (Fr(0), Fr(0), Fr(0)), (Fr(1), Fr(0), Fr(0)), (Fr(0), Fr(1), Fr(0))]


def mindex(k, n):
    return [a for a in product(range(k + 1), repeat=n) if sum(a) == k]


def mfact(a):
    r = 1
    for x in a:
        r *= factorial(x)
    return r


def inv3(M):
    a, b, c = M[0]; d, e, f = M[1]; g, h, i = M[2]
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    adj = [[e * i - f * h, c * h - b * i, b * f - c * e],
           [f * g - d * i, a * i - c * g, c * d - a * f],
           [d * h - e * g, b * g - a * h, a * e - b * d]]
    return [[x / det for x in row] for row in adj], det


def bary_grads(V):
    """gradients of the 4 barycentric coordinates of tet V (list of 4 points)."""
    E = [[V[j][c] - V[0][c] for c in range(3)] for j in (1, 2, 3)]      # rows = edge vectors
    Mt = [[E[j][c] for j in range(3)] for c in range(3)]               # columns = edges
    Mi, det = inv3(Mt)                                                  # rows of Mi = grad lambda_1..3
    g = [Mi[0], Mi[1], Mi[2]]
    g0 = [-(g[0][c] + g[1][c] + g[2][c]) for c in range(3)]
    return [g0] + g, abs(det) / 6


class Alfeld:
    def __init__(self, k, zero='all'):
        self.k = k
        B = tuple(sum(p[c] for p in P) / 4 for c in range(3))
        self.B = B
        self.subs = [[B] + [P[j] for j in range(4) if j != i] for i in range(4)]
        # zero faces
        if zero == 'all':
            zf = [[P[j] for j in range(4) if j != i] for i in range(4)]
        else:
            zf = [[P[j] for j in range(4) if j != i] for i in (1, 2, 3)]
        zpts = set()
        for tri in zf:
            for b in mindex(k, 3):
                zpts.add(tuple(sum(Fr(b[j], k) * tri[j][c] for j in range(3)) for c in range(3)))
        self.key = {}
        self.loc = []
        for V in self.subs:
            d = {}
            for a in mindex(k, 4):
                pt = tuple(sum(Fr(a[j], k) * V[j][c] for j in range(4)) for c in range(3))
                if pt in zpts:
                    d[a] = None
                else:
                    if pt not in self.key:
                        self.key[pt] = len(self.key)
                    d[a] = self.key[pt]
            self.loc.append(d)
        self.nd = len(self.key)
        self.nunk = 3 * self.nd
        self.gr = []
        self.vol = []
        for V in self.subs:
            g, vol = bary_grads(V)
            self.gr.append(g); self.vol.append(vol)
        self.bk1 = mindex(k - 1, 4)

    # --- D[s][beta] : 3x3 matrix (as dict of unknown -> coefficient) giving grad v on sub s in the
    #     degree-(k-1) Bernstein basis: grad v = k sum_beta B_beta D_beta,  (D_beta)_{a,b} = sum_l c^a_{beta+e_l} (grad lam_l)_b
    def Dlin(self, s, beta):
        """returns dict (a,b) -> {unknown: coef} (linear in the unknowns), without the factor k."""
        out = {}
        d = self.loc[s]; g = self.gr[s]
        for l in range(4):
            al = tuple(beta[m] + (1 if m == l else 0) for m in range(4))
            n = d[al]
            if n is None:
                continue
            for a in range(3):
                for b in range(3):
                    if g[l][b] != 0:
                        dd = out.setdefault((a, b), {})
                        dd[3 * n + a] = dd.get(3 * n + a, 0) + g[l][b]
        return out

    def div_rows(self):
        rows = []
        for s in range(len(self.subs)):
            for beta in self.bk1:
                D = self.Dlin(s, beta)
                r = {}
                for a in range(3):
                    for u, c in D.get((a, a), {}).items():
                        r[u] = r.get(u, 0) + c
                rows.append(r)
        return rows

    def mass_k1(self, s):
        """exact L2 Gram of the degree-(k-1) Bernstein basis on sub s."""
        m = self.k - 1; vol = self.vol[s]
        C = lambda a: Fr(factorial(m), mfact(a))
        out = {}
        denom_all = Fr(factorial(2 * m) * 6 * factorial(2 * m), factorial(2 * m + 3))  # unused helper
        for a in self.bk1:
            for b in self.bk1:
                ab = tuple(x + y for x, y in zip(a, b))
                # int B_a B_b = vol * C(a)C(b) * 3! (a+b)! / (2m+3)!
                out[(a, b)] = vol * C(a) * C(b) * 6 * mfact(ab) / factorial(2 * m + 3)
        return out

    def grad_gram(self, sym=False):
        """exact matrix (dict of dicts) of int grad v : grad w  (or 1/2 int Dv:Dw if sym) over T^, in unknowns."""
        k = self.k
        G = {}
        for s in range(len(self.subs)):
            M = self.mass_k1(s)
            Ds = {beta: self.Dlin(s, beta) for beta in self.bk1}
            for b1 in self.bk1:
                D1 = Ds[b1]
                for b2 in self.bk1:
                    m = M[(b1, b2)] * k * k
                    D2 = Ds[b2]
                    pairs = []
                    for a in range(3):
                        for b in range(3):
                            if sym:
                                # 1/2 (G+G^T):(H+H^T) = G:H + G:H^T
                                pairs.append(((a, b), (a, b)))
                                pairs.append(((a, b), (b, a)))
                            else:
                                pairs.append(((a, b), (a, b)))
                    for p1, p2 in pairs:
                        x1 = D1.get(p1); x2 = D2.get(p2)
                        if not x1 or not x2:
                            continue
                        for u, c1 in x1.items():
                            row = G.setdefault(u, {})
                            for w, c2 in x2.items():
                                row[w] = row.get(w, 0) + m * c1 * c2
        return G

    # ---------------- face F^ (in K_0, B is vertex 0 of K_0, face vertices P1,P2,P3 -> mu1,mu2,mu3)
    def face_trace(self):
        """v on F^: dict alpha' (deg k, 3 indices) -> {unknown(comp a): coef} with v_a = sum c B^k_alpha'(mu)."""
        d = self.loc[0]; out = {}
        for ap in mindex(self.k, 3):
            n = d[(0,) + ap]
            if n is not None:
                out[ap] = n
        return out

    def face_grad(self):
        """G|_F = k sum_{beta'} B^{k-1}_{beta'}(mu) D_{(0,beta')}: dict beta' -> Dlin."""
        return {bp: self.Dlin(0, (0,) + bp) for bp in mindex(self.k - 1, 3)}


def face_int(m, a, n, b, extra=(0, 0, 0)):
    """int_{F^} B^m_a B^n_b mu^extra dA over the reference triangle (area 1/2)."""
    g = tuple(x + y + z for x, y, z in zip(a, b, extra))
    return Fr(factorial(m), mfact(a)) * Fr(factorial(n), mfact(b)) * Fr(mfact(g), factorial(sum(g) + 2))


def face_moment_rows(A, weights):
    """rows (dict unknown->coef) for int_F w(mu) g_c, g = d_z v = G e3, c=0,1 (tangential comps).
    weights: list of (name, exponent-tuple e) meaning w = mu^e (mu1,mu2,mu3)."""
    FG = A.face_grad(); k = A.k
    rows = {}
    for name, e in weights:
        for c in (0, 1):
            r = {}
            for bp, D in FG.items():
                coef = k * Fr(factorial(k - 1), mfact(bp)) * Fr(mfact(tuple(x + y for x, y in zip(bp, e))), factorial(k - 1 + sum(e) + 2))
                for u, cc in D.get((c, 2), {}).items():
                    r[u] = r.get(u, 0) + coef * cc
            rows[(name, c)] = r
    return rows


def nullspace(rows, nunk):
    M = DomainMatrix([[QQ(r.get(j, 0).numerator, r.get(j, 0).denominator) if j in r else QQ(0) for j in range(nunk)]
                      for r in rows], (len(rows), nunk), QQ)
    rk = M.rank()
    N = M.nullspace()   # rows span kernel
    Nl = [[Fr(int(N[i, j].element.numerator), int(N[i, j].element.denominator)) for j in range(nunk)] for i in range(N.shape[0])]
    return rk, Nl


def apply_rows(rows, N):
    """rows (dict) applied to basis vectors N (list of lists) -> matrix len(rows) x len(N)."""
    return [[sum(c * v[u] for u, c in r.items()) for v in N] for r in rows]


def congruence(G, N):
    """N G N^T for sparse dict-of-dict G."""
    d = len(N)
    GN = []
    for v in N:
        w = {}
        for u, row in G.items():
            if v[u] == 0:
                continue
            for x, g in row.items():
                w[x] = w.get(x, 0) + v[u] * g
        GN.append(w)
    return [[sum(c * N[j][x] for x, c in GN[i].items()) for j in range(d)] for i in range(d)]
