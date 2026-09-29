"""
notes/v7/nsconst: isoparametric Taylor--Hood P_k / P_{k-1} (continuous pressure) on the EXACT domain
Omega = (-L,L)^2 minus the closed unit disk, on the topology of the paper's production meshes (svn.make_mesh).
Boundary triangles are curved: P_k-isoparametric interpolation of the Zlamal blending map of the radial projection.
Only used to evaluate the continuous constants C_F, beta_O, C_R of Theorem n2:L (REPORT.tex, Section 2).

Everything is assembled with the exact (curved) Jacobian at every quadrature point.
"""
import os, sys
for v_ in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v_, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
import numpy as np, scipy.sparse as sp
import svn, svn_k


class Iso:
    def __init__(self, N, k=4, L=2.5, nq=None):
        verts, tris, circ, outer = svn.make_mesh(N, L=L)
        self.N, self.k, self.L = N, k, L
        self.R = svn_k.Ref(k); self.Rp = svn_k.Ref(k - 1)
        if nq:
            self.R.QX, self.R.QY, self.R.QW = svn_k.tri_quad(nq)
        S = svn_k.Space(verts, tris, circ, outer, self.R)
        Sp = svn_k.Space(verts, tris, circ, outer, self.Rp)
        self.S, self.Sp = S, Sp
        xy = S.xy.copy()
        # ---- curve the boundary triangles (Zlamal blending of the radial projection, interpolated at the P_k nodes)
        lb = np.array(self.R.LBARY, float) / k          # (l1, l2, l0) per local node
        lam_of = lambda v: {0: lb[:, 2], 1: lb[:, 0], 2: lb[:, 1]}[v]
        self.curved = np.zeros(len(tris), bool)
        for (t, u, w) in S.bedges:
            A, B = verts[tris[t, u]], verts[tris[t, w]]
            la, lbb = lam_of(u), lam_of(w)
            s = la + lbb
            with np.errstate(invalid="ignore", divide="ignore"):
                tt = np.where(s > 1e-14, lbb / np.maximum(s, 1e-300), 0.0)
            P = A[None, :] + tt[:, None] * (B - A)[None, :]
            dlt = P / np.linalg.norm(P, axis=1)[:, None] - P
            D = s[:, None] * dlt
            base = S.P0[t][None, :] + (S.J[t] @ self.R.REFNODES.T).T
            xy[S.ids[t]] = base + D
            self.curved[t] = True
        self.xy = xy
        # ---- geometry at volume quadrature points
        R = self.R
        dX, dY = R.dbasis(R.QX, R.QY)                     # (nq, nl)
        Xn = xy[S.ids]                                     # (nt, nl, 2)
        J = np.stack([np.einsum('qa,tac->tqc', dX, Xn), np.einsum('qa,tac->tqc', dY, Xn)], -1)  # J[t,q,c,d] = dx_c/dxi_d
        self.detJ = J[..., 0, 0] * J[..., 1, 1] - J[..., 0, 1] * J[..., 1, 0]
        assert (self.detJ > 0).all() or (self.detJ < 0).all()
        Ji = np.linalg.inv(J)                              # Ji[t,q,d,c] = dxi_d/dx_c
        self.gx = np.einsum('tq,qa->tqa', Ji[..., 0, 0], dX) + np.einsum('tq,qa->tqa', Ji[..., 1, 0], dY)
        self.gy = np.einsum('tq,qa->tqa', Ji[..., 0, 1], dX) + np.einsum('tq,qa->tqa', Ji[..., 1, 1], dY)
        self.wq = R.QW[None, :] * np.abs(self.detJ)
        self.Xq = np.einsum('qa,tac->tqc', R.PHI_Q, Xn)    # physical quadrature points
        self.nt = len(tris); self.nu = S.nn; self.np_ = Sp.nn
        # pressure basis at the same points (reference P_{k-1} Lagrange), and its physical gradient
        Rp = self.Rp
        self.PP = Rp.basis(R.QX, R.QY)                     # (nq, npl)
        dpx, dpy = Rp.dbasis(R.QX, R.QY)
        self.pgx = np.einsum('tq,qa->tqa', Ji[..., 0, 0], dpx) + np.einsum('tq,qa->tqa', Ji[..., 1, 0], dpy)
        self.pgy = np.einsum('tq,qa->tqa', Ji[..., 0, 1], dpx) + np.einsum('tq,qa->tqa', Ji[..., 1, 1], dpy)
        # Dirichlet velocity nodes: square + circle (all of dOmega)
        cn = set()
        for (t, u, w) in S.bedges:
            for kk, (l1, l2, l0) in enumerate(R.LBARY):
                lam = {0: l0, 1: l1, 2: l2}
                if lam[3 - u - w] == 0:
                    cn.add(S.ids[t, kk])
        self.cnodes = np.array(sorted(cn))
        self.dnodes = np.union1d(S.dnodes, self.cnodes)
        self.free = np.setdiff1d(np.arange(self.nu), self.dnodes)
        self.vfree = np.sort(np.concatenate([2 * self.free, 2 * self.free + 1]))

    # ------------------------------------------------------------------ scalar / vector assembly helpers
    def _scal(self, loc, ids_r, ids_c, nr, nc):
        nt = self.nt
        rows = np.repeat(ids_r, ids_c.shape[1], axis=1).ravel(); cols = np.tile(ids_c, (1, ids_r.shape[1])).ravel()
        return sp.coo_matrix((loc.reshape(nt, -1).ravel(), (rows, cols)), shape=(nr, nc)).tocsr()

    def scalar_stiff(self):
        loc = np.einsum('tq,tqa,tqb->tab', self.wq, self.gx, self.gx) + np.einsum('tq,tqa,tqb->tab', self.wq, self.gy, self.gy)
        return self._scal(loc, self.S.ids, self.S.ids, self.nu, self.nu)

    def scalar_mass(self):
        loc = np.einsum('tq,qa,qb->tab', self.wq, self.R.PHI_Q, self.R.PHI_Q)
        return self._scal(loc, self.S.ids, self.S.ids, self.nu, self.nu)

    @staticmethod
    def vec(Ms):
        """block-diagonal vector version, dof ordering 2*node + comp"""
        Ms = Ms.tocoo()
        r = np.concatenate([2 * Ms.row, 2 * Ms.row + 1]); c = np.concatenate([2 * Ms.col, 2 * Ms.col + 1])
        return sp.coo_matrix((np.concatenate([Ms.data, Ms.data]), (r, c)), shape=(2 * Ms.shape[0], 2 * Ms.shape[1])).tocsr()

    def div(self):
        """B[q, v] = (q, div v), q continuous P_{k-1}"""
        Bx = np.einsum('tq,qi,tqa->tia', self.wq, self.PP, self.gx); By = np.einsum('tq,qi,tqa->tia', self.wq, self.PP, self.gy)
        nt, npl, nl = Bx.shape
        loc = np.stack([Bx, By], -1).reshape(nt, npl, 2 * nl)
        D = np.stack([2 * self.S.ids, 2 * self.S.ids + 1], -1).reshape(nt, 2 * nl)
        return self._scal(loc, self.Sp.ids, D, self.np_, 2 * self.nu)

    def pmass(self):
        loc = np.einsum('tq,qa,qb->tab', self.wq, self.PP, self.PP)
        return self._scal(loc, self.Sp.ids, self.Sp.ids, self.np_, self.np_)

    def pstiff(self):
        loc = np.einsum('tq,tqa,tqb->tab', self.wq, self.pgx, self.pgx) + np.einsum('tq,tqa,tqb->tab', self.wq, self.pgy, self.pgy)
        return self._scal(loc, self.Sp.ids, self.Sp.ids, self.np_, self.np_)

    def pmean(self):
        v = np.zeros(self.np_); np.add.at(v, self.Sp.ids.ravel(), np.einsum('tq,qa->ta', self.wq, self.PP).ravel())
        return v

    def oseen_conv(self, ufun, gufun):
        """matrix of (w.grad)u + (u.grad)w tested with v (standard form c1), u given analytically."""
        X = self.Xq
        U = ufun(X[..., 0], X[..., 1])                     # (nt,nq,2)
        G = gufun(X[..., 0], X[..., 1])                    # (nt,nq,2,2) [c,d] = d_d u_c
        PHI = self.R.PHI_Q
        adv = U[..., 0:1] * self.gx + U[..., 1:2] * self.gy        # u.grad phi_b
        K1 = np.einsum('tq,qa,tqb->tab', self.wq, PHI, adv)
        nt, nl = self.nt, self.R.nl
        loc = np.zeros((nt, nl, 2, nl, 2))
        loc[:, :, 0, :, 0] += K1; loc[:, :, 1, :, 1] += K1
        loc += np.einsum('tq,qa,qb,tqcd->tacbd', self.wq, PHI, PHI, G)   # (w.grad u)_c = G[c,d] w_d
        D = np.stack([2 * self.S.ids, 2 * self.S.ids + 1], -1).reshape(nt, 2 * nl)
        return self._scal(loc.reshape(nt, 2 * nl, 2 * nl), D, D, 2 * self.nu, 2 * self.nu)

    def load(self, gfun):
        """F[v] = (g, v) for vector g given on physical points; gfun(X, Y) -> (..., 2)"""
        X = self.Xq
        g = gfun(X[..., 0], X[..., 1])
        loc = np.einsum('tq,qa,tqc->tac', self.wq, self.R.PHI_Q, g)
        F = np.zeros(2 * self.nu)
        D = np.stack([2 * self.S.ids, 2 * self.S.ids + 1], -1)
        np.add.at(F, D.ravel(), loc.ravel())
        return F

    def circle_edges(self):
        """quadrature on the curved Gamma: for each boundary triangle, points (ref), weights (arc length), outward normal of Omega"""
        R = self.R; out = []
        Rr = np.array([[0, 0], [1, 0], [0, 1]], float)
        for (t, u, w) in self.S.bedges:
            ref = Rr[u][None, :] + R.GL[:, None] * (Rr[w] - Rr[u])[None, :]
            Xn = self.xy[self.S.ids[t]]
            dX, dY = R.dbasis(ref[:, 0], ref[:, 1])
            J = np.stack([dX @ Xn, dY @ Xn], -1)            # (ng, 2, 2)
            tang = J @ (Rr[w] - Rr[u])
            ds = np.linalg.norm(tang, axis=1)
            X = R.basis(ref[:, 0], ref[:, 1]) @ Xn
            n = -X / np.linalg.norm(X, axis=1)[:, None]      # outward normal of Omega on the circle: -e_r
            Ji = np.linalg.inv(J)
            gx = Ji[:, 0, 0][:, None] * dX + Ji[:, 1, 0][:, None] * dY
            gy = Ji[:, 0, 1][:, None] * dX + Ji[:, 1, 1][:, None] * dY
            out.append((t, X, R.GLW * ds, n, gx, gy))
        return out
