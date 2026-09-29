"""notes/v7/P: the discrete first-layer cell problem for the CNS* pressure moment (REPORT.tex, Section 2).

Cell frame (cell units, wall chord = [0,1] x {0}, fluid y > 0):
  tau = (1,0), e_in = (0,1) (into the fluid), outward normal of Omega_h n = (0,-1),  sigma = x - 1/2.
Lattice: rows j = 0..J, left vertex L_j, right vertex L_j + (1,0) (translation periodic, period 1);
row j -> j+1 split 'f' (forward: [L_j, R_j, R_{j+1}], [L_j, R_{j+1}, L_{j+1}]) or 'b' (backward).
Traction-free top (the far field is free: constant velocity, pressure -> 0), no Dirichlet anywhere.

Forms (code/svn.py, S.h = 1 so the Nitsche penalty is lam = lambda_hat = mu*l/h):
  N(U,v) = 1/2 (DU,Dv) - <d_n U, v> - <U, d_n v> + lam <U, v>,   b(P,v) = -(P, div v) + <P, v.n>
  N(U,v) + b(P,v) = F(v),   (q, div U) = 0 for all q in P3-disc.
Unit problems:
  'bump' : F(v) = <t1, v> - <g1, d_n v> + lam <g1, v>,  g1 = -beta(sigma) tau, beta = (1/4 - sigma^2)/2,
           t1 = sigma e_in                                   (chordal bump + transposed load, kappa W l^2 = 1)
  'trac' : F(v) = <n, v>                                    (uniform wall traction, global flux multiplier)
Outputs: A = int_0^1 P(sigma,0) dsigma (wall mean of the pressure, trace from the wall triangle),
  V = int_0^1 U.n, identity lam*V + A = int F-normal part (0 for bump, 1 for trac), ||P||_{L2(cell)},
  far-field velocity, per-row decay of |P| and |grad U|.
usage: python3 cell.py scan            (delta/aspect/lam tables)
"""
import os, sys, json
for v_ in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v_, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
import svn


def lattice_rows(a, delta, J):
    """uniform sheared lattice: L_j = (j*delta, j*a)."""
    return np.array([[j * delta, j * a] for j in range(J + 1)], float)


class CellProblem:
    def __init__(self, rows, splits, lam):
        rows = np.asarray(rows, float); J = len(rows) - 1
        verts = np.zeros((2 * (J + 1), 2))
        verts[0::2] = rows; verts[1::2] = rows + np.array([1.0, 0.0])
        vid = lambda j, i: 2 * j + i
        tris = []
        for j in range(J):
            s = splits if (isinstance(splits, str) and len(splits) == 1) else splits[j]
            if s == 'f':
                tris += [[vid(j, 0), vid(j, 1), vid(j + 1, 1)], [vid(j, 0), vid(j + 1, 1), vid(j + 1, 0)]]
            else:
                tris += [[vid(j, 0), vid(j, 1), vid(j + 1, 0)], [vid(j, 1), vid(j + 1, 1), vid(j + 1, 0)]]
        tris = np.array(tris)
        P = verts[tris]; ar = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
        assert np.all(np.abs(ar) > 1e-12)
        tris[ar < 0] = tris[ar < 0][:, [0, 2, 1]]
        S = svn.Space(verts, tris, {0, 1}, {2 * J, 2 * J + 1})
        S.h = 1.0
        assert len(S.bedges) == 1
        self.S, self.J, self.lam, self.rows, self.tris = S, J, lam, rows, tris
        zero = lambda X, Y: np.zeros(np.shape(X) + (2,))
        sysm = svn.assemble(S, lam, zero, 1)
        self.A, self.B, self.C = sysm["A"].tocsr(), sysm["B"].tocsr(), sysm["C"].tocsr()
        # periodic identification x -> x - 1
        xy = S.xy
        key = {(round(x, 9), round(y, 9)): n for n, (x, y) in enumerate(xy)}
        match = {}
        for n, (x, y) in enumerate(xy):
            m = key.get((round(x - 1.0, 9), round(y, 9)))
            if m is not None:
                match[n] = m
        red = [n for n in range(S.nn) if n not in match]
        pos = {n: k for k, n in enumerate(red)}
        r_, c_ = [], []
        for n in range(S.nn):
            k = pos[match[n]] if n in match else pos[n]
            for c in range(2):
                r_.append(2 * n + c); c_.append(2 * k + c)
        self.T = sp.coo_matrix((np.ones(len(r_)), (r_, c_)), shape=(2 * S.nn, 2 * len(red))).tocsr()
        T = self.T
        Bt = (self.B - self.C)
        self.K = sp.bmat([[T.T @ self.A @ T, -(Bt @ T).T], [self.B @ T, None]], format="csc")
        self.lu = spl.splu(self.K, permc_spec="COLAMD", diag_pivot_thresh=1.0)
        (t, ref, Xe, we, n) = S.edge_data()[0]
        self.wall = (t, ref, Xe, we, n)
        assert np.allclose(n, [0, -1])

    def rhs_edge(self, g, tload):
        """F(v) = <tload, v> - <g, d_n v> + lam <g, v> on the wall; g, tload: (ng,2) at the edge Gauss points."""
        S = self.S; (t, ref, Xe, we, n) = self.wall
        ph = svn.basis(ref[:, 0], ref[:, 1]); dX, dY = svn.dbasis(ref[:, 0], ref[:, 1]); Ji = S.Jinv[t]
        dn = n[0] * (Ji[0, 0] * dX + Ji[1, 0] * dY) + n[1] * (Ji[0, 1] * dX + Ji[1, 1] * dY)
        F = np.zeros(2 * S.nn); d = svn.dofs(S.ids[t])
        for c in range(2):
            np.add.at(F, d[:, c], np.einsum('g,g,ga->a', we, tload[:, c], ph)
                      + np.einsum('g,g,ga->a', we, g[:, c], -dn + self.lam * ph))
        return F

    def solve(self, F):
        rhs = np.concatenate([self.T.T @ F, np.zeros(self.B.shape[0])])
        x = self.lu.solve(rhs)
        for _ in range(2):
            x = x + self.lu.solve(rhs - self.K @ x)
        nr = self.T.shape[1]
        return self.T @ x[:nr], x[nr:], float(np.linalg.norm(self.K @ x - rhs) / max(np.linalg.norm(rhs), 1e-300))

    def unit_data(self, kind):
        (t, ref, Xe, we, n) = self.wall
        sig = Xe[:, 0] - 0.5 - self.rows[0, 0]
        if kind == "bump":
            g = np.stack([-(0.25 - sig ** 2) / 2, 0 * sig], -1)
            tl = np.stack([0 * sig, sig], -1)
        elif kind == "bumpg":          # bump only (no transposed load)
            g = np.stack([-(0.25 - sig ** 2) / 2, 0 * sig], -1); tl = 0 * g
        elif kind == "tl":             # transposed load only
            tl = np.stack([0 * sig, sig], -1); g = 0 * tl
        elif kind == "trac":
            g = 0 * np.stack([sig, sig], -1); tl = np.stack([0 * sig, -1 + 0 * sig], -1)   # n = (0,-1)
        return g, tl

    def analyse(self, U, P):
        S = self.S; (t, ref, Xe, we, n) = self.wall
        q3 = svn._mono(svn.MON3, ref[:, 0], ref[:, 1])
        Pw = q3 @ P[10 * t:10 * t + 10]
        ph = svn.basis(ref[:, 0], ref[:, 1])
        Uw = ph @ U[svn.dofs(S.ids[t])]
        A = float(we @ Pw); V = float(we @ (Uw @ n))
        sig = Xe[:, 0] - 0.5 - self.rows[0, 0]
        # volume norms per row
        wq = svn.QW[None, :] * np.abs(S.detJ)[:, None]
        Pq = np.einsum('qj,tj->tq', svn.Q3_Q, P.reshape(-1, 10))
        gx, gy = S.grads(svn.QX, svn.QY)
        Ul = U[svn.dofs(S.ids)]
        ux = np.einsum('tqa,tac->tqc', gx, Ul); uy = np.einsum('tqa,tac->tqc', gy, Ul)
        cy = S.verts[S.tris].mean(1)[:, 1]
        rowP, rowG = [], []
        a = self.rows[1, 1]
        for j in range(self.J):
            m = (cy > j * a) & (cy < (j + 1) * a)
            rowP.append(float(np.sqrt((wq[m] * Pq[m] ** 2).sum())))
            rowG.append(float(np.sqrt((wq[m] * (ux[m] ** 2 + uy[m] ** 2).sum(-1)).sum())))
        # half-cell norms (only the first ~5 rows matter)
        PL2 = float(np.sqrt((wq * Pq ** 2).sum())); GL2 = float(np.sqrt((wq * (ux ** 2 + uy ** 2).sum((-1))).sum()))
        # odd/even parts of the wall pressure about sigma = 0 (Gauss points are symmetric)
        Pw_odd = 0.5 * (Pw - Pw[::-1]); Pw_even = 0.5 * (Pw + Pw[::-1])
        # wall traction trace: T = (1/2)(DU)n - P n with (1/2)DU = sym grad? code form: a = 1/2 (DU,DU), D = grad+grad^T
        dX, dY = svn.dbasis(ref[:, 0], ref[:, 1]); Ji = S.Jinv[t]
        gxe = (Ji[0, 0] * dX + Ji[1, 0] * dY) @ U[svn.dofs(S.ids[t])]
        gye = (Ji[0, 1] * dX + Ji[1, 1] * dY) @ U[svn.dofs(S.ids[t])]
        G = np.stack([gxe, gye], -1)                 # G[g, c, d] = d_d U_c
        Dn = np.einsum('gcd,d->gc', G + np.transpose(G, (0, 2, 1)), n)
        Tr = Dn - Pw[:, None] * n[None, :]
        top = S.xy[2 * self.J]
        Utop = U[2 * (2 * self.J):2 * (2 * self.J) + 2]
        return dict(A=A, V=V, lamV_plus_A=self.lam * V + A, PL2=PL2, GL2=GL2,
                    Pw_odd=float(np.sqrt(we @ Pw_odd ** 2)), Pw_even_mean0=float(np.sqrt(we @ (Pw_even - A) ** 2)),
                    trac_L2=float(np.sqrt(we @ (Tr ** 2).sum(-1))), Utop=[float(Utop[0]), float(Utop[1])],
                    rowP=rowP, rowG=rowG)

    def run(self, kind):
        g, tl = self.unit_data(kind)
        U, P, rr = self.solve(self.rhs_edge(g, tl))
        out = self.analyse(U, P); out["relres"] = rr
        return out, U, P


def cell_from_lattice(a, delta, lam, J=None, split='f'):
    J = J or int(np.ceil(9.0 / a))
    return CellProblem(lattice_rows(a, delta, J), split, lam)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "scan"
    if cmd == "scan":
        # (1) truncation check
        for J in (6, 10, 14, 20):
            c = CellProblem(lattice_rows(1.0, -1.0, J), 'f', 30.0)
            o, _, _ = c.run("bump"); oc, _, _ = c.run("trac")
            print(json.dumps(dict(test="J", J=J, A=o["A"], lamV_A=o["lamV_plus_A"], Ac=oc["A"],
                                  lamV_A_c=oc["lamV_plus_A"], relres=o["relres"])), flush=True)
        # (2) delta scan at a = 1, lam = 20 : A odd about delta = -1/2
        for d in (-1.0, -0.875, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0.0):
            c = cell_from_lattice(1.0, d, 30.0)
            o, _, _ = c.run("bump"); oc, _, _ = c.run("trac")
            print(json.dumps(dict(test="delta", a=1.0, delta=d, apex_offset=0.5 + d, lam=30.0, A=o["A"], Ac=oc["A"],
                                  PL2=o["PL2"], Pw_odd=o["Pw_odd"], Pw_even_mean0=o["Pw_even_mean0"],
                                  trac=o["trac_L2"], lamV_A=o["lamV_plus_A"], Utop=o["Utop"])), flush=True)
        # (3) row decay (production-type cell)
        c = cell_from_lattice(1.0, -1.0, 30.0, J=16)
        o, _, _ = c.run("bump")
        print(json.dumps(dict(test="decay", rowP=o["rowP"], rowG=o["rowG"])), flush=True)
        # (4) lam scan (resonance at lam = k(k+1)/a = 20 for a = 1), rectangle lattice delta = -1
        for lam in (5, 10, 15, 19, 19.9, 20.1, 21, 25, 30, 50, 100, 300, 1e3, 1e4, 1e6):
            c = cell_from_lattice(1.0, -1.0, lam)
            o, _, _ = c.run("bump"); oc, _, _ = c.run("trac")
            print(json.dumps(dict(test="lam", a=1.0, delta=-1.0, lam=lam, A=o["A"], Ac=oc["A"], V=o["V"],
                                  lamA_scaled=(lam - 20) * o["A"], PL2=o["PL2"], lamV_A=o["lamV_plus_A"])), flush=True)
        # (5) bump-only vs transposed-load-only split (delta = -1, lam = 30)
        c = cell_from_lattice(1.0, -1.0, 30.0)
        for kd in ("bumpg", "tl", "bump"):
            o, _, _ = c.run(kd)
            print(json.dumps(dict(test="split", kind=kd, A=o["A"])), flush=True)
