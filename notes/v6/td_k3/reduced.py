"""EXACT reduced model of the joint trace space of Y(R), R = ring (star) of an interior edge [z,q]
(Alfeld split, k=3).  See REPORT.tex, Section 3 (Ring Trace Theorem).

Ring: link cycle x_0..x_{m-1}; macro T_j = [q,z,x_j,x_{j+1}]; interior faces f_j = [q,z,x_j].
Gamma faces: a subset of the link faces F_j = [z,x_j,x_{j+1}] (all of them are on the ring boundary).

Unknowns:
  ghat_i  : P2 tangential field on F_i (derivative of v along B_j - x, x in F_i), 6 coeffs x 3 comps
  d       : d_{q-z} v (z)               (3)
  dp      : the coefficient p_q of v|_[z,q] = lam_q lam_z (p_q lam_q + d lam_z)   (3)
  p_j     : face coefficient of v|_{f_j} = lam_q lam_z (dp lam_q + d lam_z + p_j lam_x)  (3 per j)
  Phi     : common flux (1)
  d1_i    : d_{B-z} v (z) in macro of F_i (3 per Gamma face)
Relations (necessary for v in Y(R), and sufficient by the construction in the report):
  U_i:  ghat_i . N_F = 0 (all coeffs),  ghat_i(x_j)=ghat_i(x_{j+1})=0,  ghat_i(M_ab).n_{q x_j x_{j+1}} = 0
  jet:  d1.grad lamB^{s1}=0, d.grad lamq^{s2} + d1.grad lamB^{s2}=0, d.grad lamq^{s3}+d1.grad lamB^{s3}=0,
        ghat_i(z) = d1_i
  edge: ghat_i(M_{z x_j}).grad lamB^{s3} + (d+p_j)/4 . grad lamq^{s3} = 0     (s3 = [B,q,z,x_j])
        ghat_i(M_{z x_{j+1}}).grad lamB^{s2} + (d+p_{j+1})/4 . grad lamq^{s2} = 0   (s2 = [B,q,z,x_{j+1}])
  flux: A_j . ((d+dp)/30 + p_j/60) = Phi,   A_j = (q-z) x (x_j - z) / 2
Then  J = projection of the solution space onto (ghat_i).  Total moment: sum_i c_i int_{F_i} Q_i ghat_i,
c_i = |grad lamB^{s1}| |F_i| / |F_i| ... (we use the exact physical scaling:  M_i = fac_i * 2 int_ref Q ghat,
fac_i = |F_i|^2/(3 vol(s1)) ), mean_i = fac_i * 2 int_ref ghat.
"""
from fractions import Fraction as Fr
import itertools
import sympy as sp
from patchlib import F3, sub, add, scal, dot, cross, bary_grads, tri_mono_int

MONS = [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)]  # Bernstein-like monomials mu^alpha

def evalmon(al, mu):
    r = Fr(1)
    for a, m in zip(al, mu): r *= Fr(m) ** a
    return r

class Reduced:
    def __init__(self, link, q, z, gam, weights=None):
        self.link = [F3(x) for x in link]; self.q = F3(q); self.z = F3(z); m = len(link); self.m = m
        self.gam = list(gam)
        # variable layout
        names = []
        for i in self.gam:
            for c in range(3):
                for al in MONS: names.append(('g', i, c, al))
        names += [('d', c) for c in range(3)] + [('dp', c) for c in range(3)]
        for j in range(m): names += [('p', j, c) for c in range(3)]
        names += [('Phi',)]
        for i in self.gam: names += [('d1', i, c) for c in range(3)]
        self.names = names; self.ix = {n: k for k, n in enumerate(names)}
        self.eqs = []
        q_, z_ = self.q, self.z
        self.fac = {}; self.Q = {}
        for i in self.gam:
            a, b = self.link[i], self.link[(i + 1) % m]
            T = [q_, z_, a, b]
            B = scal(Fr(1, 4), add(add(q_, z_), add(a, b)))
            s1 = [B, z_, a, b]; s2 = [B, q_, z_, b]; s3 = [B, q_, z_, a]
            G1, D1 = bary_grads(s1); G2, _ = bary_grads(s2); G3, _ = bary_grads(s3)
            NF = cross(sub(a, z_), sub(b, z_))
            nqab = cross(sub(a, q_), sub(b, q_))
            # ghat as function of face barycentrics (mu_z, mu_a, mu_b)
            def gval(mu, c):
                return {('g', i, c, al): evalmon(al, mu) for al in MONS}
            # tangential: ghat . NF = 0 identically -> for each monomial
            for al in MONS:
                self.eqs.append({('g', i, c, al): NF[c] for c in range(3)})
            for mu in [(0, 1, 0), (0, 0, 1)]:
                for c in range(3): self.eqs.append(gval(mu, c))
            Mab = (0, Fr(1, 2), Fr(1, 2))
            e = {}
            for c in range(3):
                for k_, v in gval(Mab, c).items(): e[k_] = e.get(k_, 0) + v * nqab[c]
            self.eqs.append(e)
            # jet
            self.eqs.append({('d1', i, c): G1[0][c] for c in range(3)})
            for (G, qidx) in [(G2, 1), (G3, 1)]:
                e = {('d1', i, c): G[0][c] for c in range(3)}
                for c in range(3): e[('d', c)] = e.get(('d', c), 0) + G[qidx][c]
                self.eqs.append(e)
            for c in range(3):
                e = gval((1, 0, 0), c); e[('d1', i, c)] = -1; self.eqs.append(e)
            # edge conditions
            for (G, jj, mu) in [(G3, i, (Fr(1, 2), Fr(1, 2), 0)), (G2, (i + 1) % m, (Fr(1, 2), 0, Fr(1, 2)))]:
                e = {}
                for c in range(3):
                    for k_, v in gval(mu, c).items(): e[k_] = e.get(k_, 0) + v * G[0][c]
                    e[('d', c)] = e.get(('d', c), 0) + Fr(1, 4) * G[1][c]
                    e[('p', jj, c)] = e.get(('p', jj, c), 0) + Fr(1, 4) * G[1][c]
                self.eqs.append(e)
            area2 = dot(NF, NF)
            self.fac[i] = (area2 / 4) / (abs(D1) / 2)
            l2 = [dot(sub(z_, a), sub(z_, a)), dot(sub(z_, b), sub(z_, b)), dot(sub(a, b), sub(a, b))]
            w = weights[i] if weights and i in weights else l2   # (w_za, w_zb, w_ab)
            self.Q[i] = w
        for j in range(m):
            A = scal(Fr(1, 2), cross(sub(q_, z_), sub(self.link[j], z_)))
            e = {('Phi',): -1}
            for c in range(3):
                e[('d', c)] = e.get(('d', c), 0) + A[c] / 30
                e[('dp', c)] = e.get(('dp', c), 0) + A[c] / 30
                e[('p', j, c)] = e.get(('p', j, c), 0) + A[c] / 60
            self.eqs.append(e)

    def mat(self, rows):
        return sp.Matrix([[sp.Rational(r.get(n, 0).numerator, r.get(n, 0).denominator) if n in r else 0 for n in self.names] for r in rows])

    def face_int_rows(self, i, weight_fn):
        """rows (3 comps) of fac_i * 2 * int_ref weight * ghat_i ; weight_fn: dict alpha->coef (poly in mu)."""
        rows = []
        for c in range(3):
            r = {}
            for al in MONS:
                s = 0
                for wa, wc in weight_fn.items():
                    s += wc * tri_mono_int(tuple(x + y for x, y in zip(al, wa)))
                r[('g', i, c, al)] = self.fac[i] * 2 * s
            rows.append(r)
        return rows

    def mean_rows(self, i): return self.face_int_rows(i, {(0, 0, 0): 1})
    def mom_rows_face(self, i):
        wza, wzb, wab = self.Q[i]
        return self.face_int_rows(i, {(1, 1, 0): Fr(1, 2) * wza, (1, 0, 1): Fr(1, 2) * wzb, (0, 1, 1): Fr(1, 2) * wab})

    def analyse(self, frame):
        Eq = self.mat(self.eqs)
        N = Eq.nullspace()
        Nm = sp.Matrix.hstack(*N)  # columns: basis of solution space
        gidx = [self.ix[n] for n in self.names if n[0] == 'g']
        Gm = Nm.extract(gidx, list(range(Nm.shape[1])))
        jdim = Gm.rank()
        Mean = self.mat(sum([self.mean_rows(i) for i in self.gam], [])) * Nm
        K = Mean.nullspace()
        Kn = sp.Matrix.hstack(*K) if K else sp.zeros(Nm.shape[1], 0)
        tot = [dict() for _ in frame]
        for i in self.gam:
            for c, r in enumerate(self.mom_rows_face(i)):
                for fi, t in enumerate(frame):
                    for k_, v in r.items(): tot[fi][k_] = tot[fi].get(k_, 0) + Fr(t[c]) * v
        Mt = self.mat(tot) * Nm * Kn
        return jdim, Mt.rank(), Mt
