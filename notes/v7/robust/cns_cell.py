"""notes/v7/robust: EXACT Floquet reduction of the CNS* smooth-data problem (S) on rotation-invariant meshes, k = 4.

Omega_h = (disc of radius ROUT, polygonal) minus (regular N-gon inscribed in the unit circle); the mesh is invariant
under the rotation by p*alpha, alpha = 2 pi/N (p = number of wall edges per period cell).  First layer: vertex i of
the wall has its first-layer neighbour at radius 1 + rho_{i mod p} * l (l = 2 sin(pi/N) = h_Gamma); further layers
geometric with ratio (1 + a l).  Quads split on the forward diagonal (pattern of code/svn.py and unifmu2/floquet_cell.py).
  p = 1, rho = (1,)          : the UP family (LP, constant heights)                        -> 'u'
  p = 2, rho = (0.7, 1.4)    : robust5's 'alt' (first-layer height alternating; NOT LP^+)   -> 'alt'
Load v -> <a, v.n>_{Gamma_h}, a = e^{i m theta} (pure Floquet mode m != 0, so the Gamma-mean terms of B* vanish and
the flux constraint is automatic).  Solved on ONE period cell with the quasi-periodic condition
u(R x) = e^{i m p alpha} R u(x), complex saddle system, all 10 pressure moments per triangle:
  CNS*: N_h(W,v) + B*(xi,v) = <a, v.n>, (q, div W) = 0;  B*(q,v) = -(q,div v) + <q, v.n>   (mode m != 0)
  GS  : same with B(q,v) = -(q, div v)                   (delta_a of robust5 Thm S1)
Penalty per unit length lam = gamma/l (so gamma = mu h_Gamma/h with h := h_Gamma).
Output per (N, gamma): ||grad W||, ||grad delta||, ||grad(W - delta)||, S = ||grad W|| gamma / h_Gamma,
E = ||grad(W-delta)|| gamma / h_Gamma, and the fitted h-exponent across N is computed by the caller.
usage: python cns_cell.py kind N1,N2,... gamma1,gamma2,... [m=2] [a=1]      kind in {u, alt, alt2:r1:r2, per:r1:...}
"""
import os, sys, json, time
for v_ in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v_, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
import svn

ROUT = 2.5
R3 = np.array([[0, 0], [1, 0], [0, 1]], float)


def cell_mesh(N, rhos, a=1.0, rout=ROUT):
    p = len(rhos); al = 2 * np.pi / N; ell = 2 * np.sin(np.pi / N)
    r1 = [1 + rhos[i % p] * ell for i in range(p + 1)]
    R1 = max(r1); radii = [R1 * (1 + a * ell)]
    while radii[-1] < rout:
        radii.append(radii[-1] * (1 + a * ell))
    radii[-1] = rout
    if len(radii) > 1 and radii[-1] - radii[-2] < 0.5 * a * ell * radii[-2]:
        radii.pop(-2)
    layers = [[1.0] * (p + 1), r1] + [[rr] * (p + 1) for rr in radii]
    verts = []
    for j, lay in enumerate(layers):
        for i in range(p + 1):
            verts.append([lay[i] * np.cos(i * al), lay[i] * np.sin(i * al)])
    verts = np.array(verts); vid = lambda j, i: (p + 1) * j + i
    tris = []
    for j in range(len(layers) - 1):
        for i in range(p):
            tris.append([vid(j, i), vid(j, i + 1), vid(j + 1, i + 1)])
            tris.append([vid(j, i), vid(j + 1, i + 1), vid(j + 1, i)])
    tris = np.array(tris)
    circ = {vid(0, i) for i in range(p + 1)}; J = len(layers) - 1
    outer = {vid(J, i) for i in range(p + 1)}
    return verts, tris, circ, outer, ell, al, p


class PCell:
    def __init__(self, N, rhos, a=1.0):
        verts, tris, circ, outer, ell, al, p = cell_mesh(N, rhos, a)
        S = svn.Space(verts, tris, circ, outer); S.h = 1.0
        self.S, self.N, self.ell, self.al, self.p = S, N, ell, al, p
        xy = S.xy; ang = np.arctan2(xy[:, 1], xy[:, 0]); rad = np.hypot(xy[:, 0], xy[:, 1])
        be = p * al
        left = np.where(np.abs(xy[:, 1]) < 1e-9 * np.maximum(1, rad))[0]
        right = np.where(np.abs(ang - be) < 1e-9)[0]
        Rm = np.array([[np.cos(be), -np.sin(be)], [np.sin(be), np.cos(be)]])
        pre = xy[right] @ Rm; match = {}
        for kk, node in enumerate(right):
            d = np.linalg.norm(xy[left] - pre[kk], axis=1); j = np.argmin(d); assert d[j] < 1e-9
            match[node] = left[j]
        rmax = rad.max()
        outer_nodes = set(np.where(np.abs(rad - rmax) < 1e-9)[0].tolist())
        self.Rm, self.match, self.outer_nodes = Rm, match, outer_nodes
        self.red = [n for n in range(S.nn) if n not in match and n not in outer_nodes]
        self.pos = {n: kk for kk, n in enumerate(self.red)}
        zf = lambda X, Y: np.zeros(np.shape(X) + (2,))
        s1 = svn.assemble(S, 1.0, zf, 0, want=("A", "B")); s2 = svn.assemble(S, 2.0, zf, 0, want=("A",))
        self.Avol = s1["Avol"].tocsr(); self.Mg = (s2["A"] - s1["A"]).tocsr()
        self.Cn = (s1["A"] - s1["Avol"] - self.Mg).tocsr()
        self.B = s1["B"].tocsr(); self.C = s1["C"].tocsr(); self.Minv = s1["Minv"].tocsr()
        assert len(S.bedges) == p

    def T(self, phase):
        S = self.S; rows, cols, vals = [], [], []; om = np.exp(1j * phase)
        for n in range(S.nn):
            if n in self.outer_nodes:
                continue
            if n in self.match:
                kk = self.pos[self.match[n]]
                for c in range(2):
                    for d in range(2):
                        rows.append(2 * n + c); cols.append(2 * kk + d); vals.append(om * self.Rm[c, d])
            else:
                kk = self.pos[n]
                for c in range(2):
                    rows.append(2 * n + c); cols.append(2 * kk + c); vals.append(1.0 + 0j)
        return sp.coo_matrix((vals, (rows, cols)), shape=(2 * S.nn, 2 * len(self.red))).tocsr()

    def load(self, m):
        S = self.S; F = np.zeros(2 * S.nn, complex)
        for (t, ref, Xe, we, n) in S.edge_data():
            ph = svn.basis(ref[:, 0], ref[:, 1]); d = svn.dofs(S.ids[t])
            aa = np.exp(1j * m * np.arctan2(Xe[:, 1], Xe[:, 0]))
            for c in range(2):
                np.add.at(F, d[:, c], ph.T @ (we * aa) * n[c])
        return F

    def gradnorm(self, vec):
        S = self.S; gx, gy = S.grads(svn.QX, svn.QY)
        wq = svn.QW[None, :] * np.abs(S.detJ)[:, None]; Vl = vec[svn.dofs(S.ids)]
        vx = np.einsum('tqa,tac->tqc', gx, Vl); vy = np.einsum('tqa,tac->tqc', gy, Vl)
        return np.sqrt((wq * (np.abs(vx) ** 2 + np.abs(vy) ** 2).sum(-1)).sum() * self.N / self.p)


def solve(cell, gam, m, T, TH, F):
    lam = gam / cell.ell
    Kr = (TH @ (cell.Avol + cell.Cn + lam * cell.Mg) @ T).tocsc()
    BT = (cell.B @ T).tocsr(); Fr = TH @ F
    out = {}
    Mw = (cell.C.T @ cell.Minv @ cell.C).tocsr()          # <xi_b[u.n], v.n>: the Nitsche pressure layer as a form
    Keff = (TH @ (cell.Avol + cell.Cn + lam * cell.Mg - Mw) @ T).tocsc()
    # equal-height shift: gamma' = gamma - k(k+1) h/H  (k = 4), meaningful when all wall triangles have the same height
    Hs = [abs(cell.S.detJ[t]) / np.linalg.norm(cell.S.verts[cell.S.tris[t, w]] - cell.S.verts[cell.S.tris[t, u]]) for (t, u, w) in cell.S.bedges]
    cell.Hmin, cell.Hmax = min(Hs), max(Hs)
    lamp = lam - 20.0 / np.mean(Hs)
    Kp = (TH @ (cell.Avol + cell.Cn + lamp * cell.Mg) @ T).tocsc()
    for name, Pm, KK in (("CNS", ((cell.B - cell.C) @ T).tocsr(), Kr), ("GS", BT, Kr), ("EFF", BT, Keff), ("GSP", BT, Kp)):
        M = sp.bmat([[KK, -Pm.conj().T], [BT, None]], format="csc")
        rhs = np.concatenate([Fr, np.zeros(BT.shape[0], complex)])
        lu = spl.splu(M, permc_spec="COLAMD", diag_pivot_thresh=1.0); x = lu.solve(rhs)
        for _ in range(2):
            x = x + lu.solve(rhs - M @ x)
        out[name] = (T @ x[:Kr.shape[0]], float(np.linalg.norm(M @ x - rhs) / np.linalg.norm(rhs)))
    return out


def rhos_of(kind):
    if kind == "u":
        return (1.0,)
    if kind == "alt":
        return (0.7, 1.4)
    if kind.startswith("per:") or kind.startswith("alt2:"):
        return tuple(float(s) for s in kind.split(":")[1:])
    raise ValueError(kind)


if __name__ == "__main__":
    kind = sys.argv[1]; Ns = [int(s) for s in sys.argv[2].split(",")]; gams = [float(s) for s in sys.argv[3].split(",")]
    m = int(sys.argv[4]) if len(sys.argv) > 4 else 2; a = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
    for N in Ns:
        t0 = time.time(); cell = PCell(N, rhos_of(kind), a)
        T = cell.T(m * cell.p * cell.al); TH = T.conj().T.tocsr(); F = cell.load(m)
        for gam in gams:
            o = solve(cell, gam, m, T, TH, F)
            W, rW = o["CNS"]; D, rD = o["GS"]; Df, rF = o["EFF"]
            gW, gD, gE = cell.gradnorm(W), cell.gradnorm(D), cell.gradnorm(W - D)
            gF, gR = cell.gradnorm(Df), cell.gradnorm(W - Df)
            Dp, rP = o["GSP"]; gP, gRP = cell.gradnorm(Dp), cell.gradnorm(W - Dp)
            gp = gam - 20.0 * cell.ell / (0.5 * (cell.Hmin + cell.Hmax))
            h = cell.ell
            rec = dict(kind=kind, N=N, m=m, a=a, gamma=gam, hG=h, gradW=gW, gradD=gD, gradE=gE,
                       S=gW * gam / h, SD=gD * gam / h, E=gE * gam / h, E_layer=gE * gam ** 1.5 / h ** 0.5, SF=gF * gam / h, R=gR * gam / h, R_h32=gR * gam / h ** 1.5, gamma_p=gp, SP=gP * gam / h, RP_h32=gRP * gam / h ** 1.5, H_over_h=(cell.Hmin / h, cell.Hmax / h),
                       relres=max(rW, rD, rF), ndof=int(T.shape[1]), secs=round(time.time() - t0, 1))
            print(json.dumps({kk: (float(f"{vv:.6g}") if isinstance(vv, float) else vv) for kk, vv in rec.items()}), flush=True)
