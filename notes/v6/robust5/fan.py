"""
Exact (rational) assembly of continuous P_k vector fields on a vertex fan, in the Bernstein basis.

Fan: z = P[0], outer points q_0..q_m = P[1..m+1]; triangles T_j = (z, q_{j-1}, q_j), j = 1..m (counterclockwise).
e = [z, q_0] is the Gamma_h edge of T_1 (interior of the fan on its left when walking z -> q_0 ... we only need
the OUTWARD normal of Omega_h on e, which points away from T_1).
Fields v are allowed to be nonzero only in the open edge e and the interior of the fan: v = 0 on every other boundary
edge of the fan, in particular on [z, q_m] (for m = 3 stars this is e'; for m >= 4 it is an interior edge of the star).

Constraint set  F(fan) = { v :
   continuity across [z, q_j], j = 1..m-1;
   v = 0 on [q_{j-1}, q_j] (all j), on [z, q_m];
   pressure-blind:  (div v, r)_{T_1} = <v.n, r>_e, (div v, r)_{T_j} = 0 (j >= 2), r in P_{k-1};
   (M)  int_e v = 0 (vector);  (D) int_e d_n v = 0 (vector, from T_1);  (D') int_{[z,q_m]} d_n v = 0 (vector, from T_m) }
Every v in F(fan), extended by zero, lies in V#(omega_z) of robust3 (Lemma char + (M),(D)), and v.n|_e lies in
P_e = { p in P_k(e): p(z) = p(q_0) = 0, int_e p = 0 }.

Everything is exact over QQ given rational vertices (python Fractions; nullspaces via flint.fmpq_mat.rref).
"""
from fractions import Fraction as Fr
from math import factorial as fact, comb
import itertools
import flint

K_DEFAULT = 4


def midx(n):
    """multi-indices (a0,a1,a2) with sum n, fixed order."""
    return [(a, b, n - a - b) for a in range(n, -1, -1) for b in range(n - a, -1, -1)]


def mfact(a):
    r = 1
    for x in a:
        r *= fact(x)
    return r


class Tri:
    """triangle with rational vertices V0,V1,V2; barycentric gradients, area."""
    def __init__(self, V):
        self.V = [tuple(Fr(c) for c in p) for p in V]
        (x0, y0), (x1, y1), (x2, y2) = self.V
        det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
        assert det != 0
        self.det = det                     # 2*signed area
        self.area2 = abs(det)
        # grad lambda_i = rot(opposite edge)/det
        self.g = [((y1 - y2) / det, (x2 - x1) / det), ((y2 - y0) / det, (x0 - x2) / det), ((y0 - y1) / det, (x1 - x0) / det)]


def int_T(T, a, b, m, n):
    """int_T B^m_a B^n_b"""
    ab = tuple(x + y for x, y in zip(a, b))
    return Fr(T.area2 * fact(m) * fact(n) * mfact(ab), mfact(a) * mfact(b) * fact(m + n + 2))


def int_edge(len_, a, b, m, n, i, j):
    """int over the edge {lambda_l = 0} (l the third index) of B^m_a B^n_b; i, j the two live indices.
    len_ = edge length (may be a Fraction or a symbol-free float; caller decides)."""
    ai, aj, bi, bj = a[i], a[j], b[i], b[j]
    return len_ * Fr(fact(m) * fact(n) * fact(ai + bi) * fact(aj + bj),
                     fact(ai) * fact(aj) * fact(bi) * fact(bj) * fact(m + n + 1))


class Fan:
    def __init__(self, P, k=K_DEFAULT, last_is_gamma=True):
        self.k = k
        self.P = [tuple(Fr(c) for c in p) for p in P]
        self.m = len(P) - 2
        self.tris = [Tri([self.P[0], self.P[j], self.P[j + 1]]) for j in range(1, self.m + 1)]
        self.I4 = midx(k); self.I3 = midx(k - 1)
        self.nl = len(self.I4)
        self.nu = 2 * self.nl * self.m                    # unknowns: tri j, comp c, alpha
        self.last_is_gamma = last_is_gamma
        z, q0 = self.P[0], self.P[1]
        # e = [z, q0]; require |e| = 1 rational normalisation: caller puts z = 0, q0 = (1,0)
        t = (q0[0] - z[0], q0[1] - z[1])
        assert t[0] ** 2 + t[1] ** 2 == 1
        # outward normal: away from T_1's apex q1
        n = (t[1], -t[0]); a = self.P[2]
        if n[0] * (a[0] - z[0]) + n[1] * (a[1] - z[1]) > 0:
            n = (-n[0], -n[1])
        self.n_e = n
        # direction normal (unnormalised) of [z, q_m], pointing away from T_m
        qm = self.P[-1]; te = (qm[0] - z[0], qm[1] - z[1]); ne = (te[1], -te[0]); a = self.P[-2]
        if ne[0] * (a[0] - z[0]) + ne[1] * (a[1] - z[1]) > 0:
            ne = (-ne[0], -ne[1])
        self.N_last = ne

    def u(self, j, c, a):
        """unknown index; j = 1..m"""
        return ((j - 1) * 2 + c) * self.nl + self.I4.index(a)

    # ---- linear functionals as dict {unknown: coeff}
    def dx_coeffs(self, j, d, c, beta):
        """coefficient functional of  (d/dx_d v_c) in Bernstein B^{k-1}_beta on T_j"""
        T = self.tris[j - 1]; out = {}
        for i in range(3):
            a = list(beta); a[i] += 1; a = tuple(a)
            out[self.u(j, c, a)] = out.get(self.u(j, c, a), 0) + self.k * T.g[i][d]
        return out

    def div_poly(self, j):
        """div v on T_j as list over beta in I3 of functionals"""
        res = []
        for beta in self.I3:
            f = {}
            for c in range(2):
                for key, val in self.dx_coeffs(j, c, c, beta).items():
                    f[key] = f.get(key, 0) + val
            res.append(f)
        return res

    def constraints(self):
        rows = []
        k, m = self.k, self.m
        # continuity across [z, q_j]: T_j alpha (a,0,b) == T_{j+1} alpha (a,b,0)
        for j in range(1, m):
            for a in range(k + 1):
                b = k - a
                for c in range(2):
                    rows.append({self.u(j, c, (a, 0, b)): 1, self.u(j + 1, c, (a, b, 0)): -1})
        # zero on outer edges (alpha0 = 0) of all T_j, and on [z, q_m] (alpha1 = 0 in T_m)
        for j in range(1, m + 1):
            for al in self.I4:
                if al[0] == 0 or (j == m and al[1] == 0):
                    for c in range(2):
                        rows.append({self.u(j, c, al): 1})
        # pressure-blind
        for j in range(1, m + 1):
            T = self.tris[j - 1]; dv = self.div_poly(j)
            for r in self.I3:
                f = {}
                for beta, fb in zip(self.I3, dv):
                    w = int_T(T, beta, r, k - 1, k - 1)
                    for key, val in fb.items():
                        f[key] = f.get(key, 0) + w * val
                if j == 1:  # - <v.n, r>_e ; e = {lambda2 = 0} of T_1, |e| = 1
                    for al in self.I4:
                        if al[2] == 0:
                            w = int_edge(1, al, r, k, k - 1, 0, 1)
                            for c in range(2):
                                key = self.u(1, c, al); f[key] = f.get(key, 0) - w * self.n_e[c]
                rows.append(f)
        # (M): int_e v_c = 0 ; int_e B^k_al = |e|/(k+1)
        for c in range(2):
            rows.append({self.u(1, c, al): Fr(1, k + 1) for al in self.I4 if al[2] == 0})
        # (D): int_e d_n v_c  (T_1)
        rows += self.dn_int_rows(1, self.n_e, edge_third=2)
        if self.last_is_gamma:
            rows += self.dn_int_rows(m, self.N_last, edge_third=1)
        return rows

    def dn_int_rows(self, j, n, edge_third):
        """int over edge {lambda_{edge_third} = 0} of T_j of (n . grad) v_c, in units of the edge length."""
        out = []
        live = [i for i in range(3) if i != edge_third]
        for c in range(2):
            f = {}
            for beta in self.I3:
                if beta[edge_third] != 0:
                    continue
                wint = Fr(1, self.k)           # int_0^1 B^{k-1}_beta dt
                for d in range(2):
                    for key, val in self.dx_coeffs(j, d, c, beta).items():
                        f[key] = f.get(key, 0) + wint * n[d] * val
            out.append(f)
        return out

    def matrix(self, rows):
        M = flint.fmpq_mat(len(rows), self.nu)
        for i, r in enumerate(rows):
            for key, val in r.items():
                if val != 0:
                    M[i, key] = flint.fmpq(val.numerator, val.denominator) if isinstance(val, Fr) else flint.fmpq(val)
        return M

    def nullspace(self):
        M = self.matrix(self.constraints())
        R, rank = M.rref()
        piv = []; r = 0
        for col in range(self.nu):
            if r < rank and R[r, col] != 0:
                piv.append(col); r += 1
        free = [c for c in range(self.nu) if c not in set(piv)]
        basis = []
        for fcol in free:
            v = [flint.fmpq(0)] * self.nu; v[fcol] = flint.fmpq(1)
            for i, pc in enumerate(piv):
                v[pc] = -R[i, fcol]
            basis.append(v)
        return basis, rank

    # ---- traces on e
    def trace_e(self, v):
        """v.n on e as Bernstein coefficients in t (t = 0 at z, 1 at q0): alpha = (k-i, i, 0)."""
        return [v[self.u(1, 0, (k_ - i, i, 0))] * self.n_e[0] + v[self.u(1, 1, (k_ - i, i, 0))] * self.n_e[1]
                for k_ in [self.k] for i in range(self.k + 1)]


# ---------------------------------------------------------------- single-triangle reference space U(T) and traces
def U_space(P, k=4):
    """U(T) for T = (P[0], P[1], P[2]): fields vanishing on the two edges other than e = [P0, P1], pressure-blind on T
    (no (M), (D)); returns (Fan, exact basis)."""
    F = Fan(P, k=k, last_is_gamma=False)
    rows = F.constraints()[:-4]                    # drop (M) and (D)
    M = F.matrix(rows); R, rank = M.rref()
    piv = []; r = 0
    for col in range(F.nu):
        if r < rank and R[r, col] != 0:
            piv.append(col); r += 1
    free = [c for c in range(F.nu) if c not in set(piv)]
    B = []
    for fc in free:
        v = [flint.fmpq(0)] * F.nu; v[fc] = flint.fmpq(1)
        for i, pc in enumerate(piv):
            v[pc] = -R[i, fc]
        B.append(v)
    return F, B


def trmat(F, basis):
    """normal traces v.n on e (Bernstein coefficients in s, s = 0 at z) of the basis fields, as an fmpq_mat"""
    ne = [flint.fmpq(x.numerator, x.denominator) for x in F.n_e]
    return flint.fmpq_mat([[v[F.u(1, 0, (F.k - i, i, 0))] * ne[0] + v[F.u(1, 1, (F.k - i, i, 0))] * ne[1]
                            for i in range(F.k + 1)] for v in basis])
