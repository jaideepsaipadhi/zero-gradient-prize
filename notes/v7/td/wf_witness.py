"""Worsey--Farin (WF) macro-element: exact rank / right-inverse computations for the k=3 witnesses (sharpness (M3^3)
and leak lower bound), analogous to kappa_cert.py for Alfeld.

Reference T^ = conv(P0=e3 apex, P1=0, P2=e1, P3=e2), Gamma-face F^ = conv(P1,P2,P3).
WF split with the interior point z at the barycentre and each face point at the face barycentre (12 sub-tetrahedra
conv(z, f_face, edge)).  Fields: continuous piecewise P_k, div-free, zero on all of dT^.  The boundary faces of the
mesh on F^ are the three sub-triangles conv(f, Pi, Pj); Y^1 imposes zero mean of g = d_z v on EACH of them.
Moments (as in kappa_cert.py) use the barycentrics mu1,mu2,mu3 of the MACRO face F^:
   E^  = (int_F^ mu_a mu_b g)_{a<b}  in R^6,     E^1 = (int_F^ mu2 g, int_F^ mu3 g) in R^4.
usage: python3 wf_witness.py k
"""
import sys, time
from fractions import Fraction as Fr
import numpy as np
import refalfeld as ra
import ratlin as rl

P = ra.P


class WF(ra.Alfeld):
    def __init__(self, k):
        self.k = k
        z = tuple(sum(p[c] for p in P) / 4 for c in range(3))
        faces = [[j for j in range(4) if j != i] for i in range(4)]
        self.subs = []
        self.fsub = []          # (sub index, (a,b)) for sub-tets with a face on F^ (face opposite P0 = [1,2,3])
        for i, fc in enumerate(faces):
            f = tuple(sum(P[j][c] for j in fc) / 3 for c in range(3))
            for a in range(3):
                for b in range(a + 1, 3):
                    ia, ib = fc[a], fc[b]
                    self.subs.append([z, f, P[ia], P[ib]])
                    if i == 0:
                        self.fsub.append((len(self.subs) - 1, (ia, ib)))
        self.B = z
        onb = lambda pt: pt[0] == 0 or pt[1] == 0 or pt[2] == 0 or pt[0] + pt[1] + pt[2] == 1
        self.key = {}; self.loc = []
        for V in self.subs:
            d = {}
            for al in ra.mindex(k, 4):
                pt = tuple(sum(Fr(al[j], k) * V[j][c] for j in range(4)) for c in range(3))
                if onb(pt):
                    d[al] = None
                else:
                    if pt not in self.key:
                        self.key[pt] = len(self.key)
                    d[al] = self.key[pt]
            self.loc.append(d)
        self.nd = len(self.key); self.nunk = 3 * self.nd
        self.gr = []; self.vol = []
        for V in self.subs:
            g, vol = ra.bary_grads(V); self.gr.append(g); self.vol.append(vol)
        self.bk1 = ra.mindex(k - 1, 4)

    def face_rows(self, weights):
        """weights: dict name -> polynomial in macro barycentrics (mu1,mu2,mu3) given as {exponent(3): coef}.
        Returns rows[(name, subface or 'all', comp)] for int over the sub-face(s) of w * g_comp."""
        k = self.k; rows = {}
        for s, (ia, ib) in self.fsub:
            V = self.subs[s]                      # [z, f, Pa, Pb]; sub-face = conv(f, Pa, Pb): sub-barycentrics nu
            # macro barycentric mu_m (m=1,2,3 <-> P1,P2,P3) of sub-face vertices f, Pa, Pb
            vert = [(Fr(1, 3), Fr(1, 3), Fr(1, 3)), tuple(Fr(1 if m == ia else 0) for m in (1, 2, 3)),
                    tuple(Fr(1 if m == ib else 0) for m in (1, 2, 3))]
            area2 = Fr(1, 3)                      # 2|subface| / (2|F^|) ... face_int integrates over area 1/2; subface area 1/6
            FG = {bp: self.Dlin(s, (0,) + bp) for bp in ra.mindex(k - 1, 3)}
            for name, poly in weights.items():
                # expand poly(mu) with mu_m = sum_j vert[j][m] nu_j  into monomials in nu
                expo = {}
                for e, c in poly.items():
                    terms = {(0, 0, 0): Fr(1)}
                    for m in range(3):
                        for _ in range(e[m]):
                            new = {}
                            for tt, cc in terms.items():
                                for j in range(3):
                                    if vert[j][m] != 0:
                                        t2 = tuple(tt[q] + (1 if q == j else 0) for q in range(3))
                                        new[t2] = new.get(t2, 0) + cc * vert[j][m]
                            terms = new
                    for tt, cc in terms.items():
                        expo[tt] = expo.get(tt, 0) + c * cc
                for comp in (0, 1):
                    r = {}
                    for bp, D in FG.items():
                        for e, c in expo.items():
                            coef = area2 * c * k * Fr(ra.factorial(k - 1), ra.mfact(bp)) * \
                                Fr(ra.mfact(tuple(x + y for x, y in zip(bp, e))), ra.factorial(k - 1 + sum(e) + 2))
                            for u, cc in D.get((comp, 2), {}).items():
                                r[u] = r.get(u, 0) + coef * cc
                    rows[(name, s, comp)] = r
        return rows


k = int(sys.argv[1]) if len(sys.argv) > 1 else 3
t0 = time.time()
A = WF(k)
rows = A.div_rows()
rk, N = ra.nullspace(rows, A.nunk)
print(f"WF k={k}: unknowns={A.nunk}, div rows={len(rows)}, rank(div)={rk}, dim Y={len(N)} ({time.time()-t0:.1f}s)", flush=True)
W = {'m': {(0, 0, 0): Fr(1)}, 'm2': {(0, 1, 0): Fr(1)}, 'm3': {(0, 0, 1): Fr(1)},
     'E12': {(1, 1, 0): Fr(1)}, 'E13': {(1, 0, 1): Fr(1)}, 'E23': {(0, 1, 1): Fr(1)}}
FR = A.face_rows(W)
subs = [s for s, _ in A.fsub]
L = {key: ra.apply_rows([r], N)[0] for key, r in FR.items()}
tot = lambda name, c: [sum(L[(name, s, c)][j] for s in subs) for j in range(len(N))]
means = [L[('m', s, c)] for s in subs for c in (0, 1)]
print("rank of the 6 sub-face means on Y:", rl.rank(means))
K1 = rl.nullspace(means) if rl.rank(means) else [list(r) for r in __import__('numpy').eye(len(N), dtype=object)]; d1 = len(K1)
if d1 == 0:
    print('dim Y^1 = 0: no witness'); sys.exit(0)
G = A.grad_gram(); GY = ra.congruence(G, N)
G1 = rl.matmul(rl.matmul(K1, GY), rl.transpose(K1))
restr = lambda vec: [sum(a * b for a, b in zip(vec, v)) for v in K1]
Eq = [restr(tot(n, c)) for n in ('E12', 'E13', 'E23') for c in (0, 1)]
E1 = [restr(tot(n, c)) for n in ('m2', 'm3') for c in (0, 1)]
print(f"dim Y^1={d1}; rank quadratic edge moments on Y^1 = {rl.rank(Eq)} (of 6); rank first moments = {rl.rank(E1)} (of 4)", flush=True)
for E, name in ((Eq, 'quadratic edge moments (sharpness)'), (E1, 'first moments (leak witness)')):
    Ef = np.array([[float(x) for x in r] for r in E]); Gf = np.array([[float(x) for x in r] for r in G1])
    M = Ef @ np.linalg.solve(Gf, Ef.T); ev = np.linalg.eigvalsh(M)
    if rl.rank(E) == len(E):
        X = rl.solve(G1, rl.transpose(E)); Mx = rl.matmul(E, X)
        lo, hi = Fr(ev[0] * (1 - 1e-6)).limit_denominator(10 ** 15), Fr(ev[0] * (1 + 1e-6)).limit_denominator(10 ** 15)
        lo, hi = rl.lmin_bracket(Mx, lo, hi, iters=30)
        print(f"  {name}: CERTIFIED ||R^|| in [{float(hi)**-0.5:.9f}, {float(lo)**-0.5:.9f}]  (exact LDL^T)")
    else:
        print(f"  {name}: NOT onto; eigenvalues of E G^-1 E^T: {ev}")
        if name.startswith('first'):
            # image direction(s) as 2x2 matrices [[g1.x1,g1.x2],[g2.x1,g2.x2]]
            U, S_, Vt = np.linalg.svd(Ef); print("   image basis (columns):", np.round(U[:, :rl.rank(E)].T, 6))
print(f"done ({time.time()-t0:.1f}s)")
