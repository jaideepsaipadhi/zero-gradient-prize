"""
3D Scott-Vogelius-Nitsche solver (NUMERICAL companion to paper/sections/12_three_d.tex).

  continuous P_k velocity (k >= 3), discontinuous P_{k-1} pressure (full P_{k-1} on every
  micro-tetrahedron) on an ALFELD split tetrahedral mesh.

Forms (same conventions as code/svn_k.py, nu = 1):
  a_h(u,v)  = 1/2 (Du, Dv),                Du = grad u + grad u^T
  N_h(u,v)  = a_h(u,v) - <flux(u), v> - <u, flux(v)> + (mu/h) <u, v>     on Nitsche faces
              flux(u) = d_n u       ("grad", the paper's Gjerde-Scott form)
              flux(u) = (Du) n      ("sym",  fully consistent; used only for the box verification)
  momentum  N_h(u,v) - (p, div v) + theta <p - pbar_Gamma(p), v.n>_{Gamma_h} = (f,v) + Nitsche data terms
  mass      (q, div u) = 0
  h = max tet diameter,  n = outward unit normal of the fluid domain (i.e. pointing INTO P_h on Gamma_h).
  GS  = theta 0,  CNS* = theta 1.

Pressure is eliminated by the iterated-penalty method (as code/svn.py "ipm"):
  (A + r Bt^T M^{-1} B) u^{n+1} = F + Bt^T p^n,   p^{n+1} = p^n - r M^{-1} B u^{n+1}
(fixed point: A u - Bt^T p = F, B u = 0).  Only the velocity matrix is factorised (MKL PARDISO).
The rank-one pbar_Gamma part of Bt for CNS* is handled by Sherman-Morrison.
"""
import os
os.environ.setdefault("PYPARDISO_MKL_RT", "/usr/local/lib/libmkl_rt.so.3")
import numpy as np
import scipy.sparse as sp

CH = 4000   # tets per chunk in vectorised loops


# ------------------------------------------------------------------------------ quadrature
def tet_quad(n):
    g, w = np.polynomial.legendre.leggauss(n)
    g = (g + 1) / 2; w = w / 2
    U, V, W = np.meshgrid(g, g, g, indexing="ij")
    WU, WV, WW = np.meshgrid(w, w, w, indexing="ij")
    x = U; y = V * (1 - U); z = W * (1 - U) * (1 - V)
    wt = WU * WV * WW * (1 - U) ** 2 * (1 - V)
    return np.stack([x.ravel(), y.ravel(), z.ravel()], 1), wt.ravel()


def tri_quad(n):
    g, w = np.polynomial.legendre.leggauss(n)
    g = (g + 1) / 2; w = w / 2
    U, V = np.meshgrid(g, g, indexing="ij"); WU, WV = np.meshgrid(w, w, indexing="ij")
    return np.stack([U.ravel(), (V * (1 - U)).ravel()], 1), (WU * WV * (1 - U)).ravel()   # weights sum 1/2


def _mono(M, X):
    return np.stack([X[:, 0] ** a * X[:, 1] ** b * X[:, 2] ** c for (a, b, c) in M], -1)


def _dmono(M, X):
    out = []
    for d in range(3):
        cols = []
        for m in M:
            if m[d] == 0:
                cols.append(0 * X[:, 0]); continue
            e = list(m); e[d] -= 1
            cols.append(m[d] * X[:, 0] ** e[0] * X[:, 1] ** e[1] * X[:, 2] ** e[2])
        out.append(np.stack(cols, -1))
    return out


class Ref:
    REFV = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)

    def __init__(self, k, nq=None, nqe=None):
        self.k = k
        self.MONV = [(a, b, c) for a in range(k + 1) for b in range(k + 1 - a) for c in range(k + 1 - a - b)]
        self.MONP = [(a, b, c) for a in range(k) for b in range(k - a) for c in range(k - a - b)]
        self.LBARY = [(k - a - b - c, a, b, c) for a in range(k + 1) for b in range(k + 1 - a) for c in range(k + 1 - a - b)]
        self.nl = len(self.MONV); self.npl = len(self.MONP)
        self.NODES = np.array([[l[1], l[2], l[3]] for l in self.LBARY], float) / k
        self.COEF = np.linalg.inv(_mono(self.MONV, self.NODES))
        # exact reference integrals (degree <= 2k)
        X, w = tet_quad(k + 3)
        D = self.dbasis(X); Q = self.pbasis(X); P = self.basis(X)
        self.K = np.einsum('q,iqa,jqb->ijab', w, np.array(D), np.array(D))
        self.Bq = np.einsum('q,qj,iqa->jia', w, Q, np.array(D))
        self.Mp = np.einsum('q,qi,qj->ij', w, Q, Q)
        self.Mpinv = np.linalg.inv(self.Mp)
        self.Mv = np.einsum('q,qa,qb->ab', w, P, P)
        # quadrature for loads / errors
        self.QX, self.QW = tet_quad(nq or (k + 3))
        self.FX, self.FW = tri_quad(nqe or (k + 4))

    def basis(self, X):
        return _mono(self.MONV, X) @ self.COEF

    def dbasis(self, X):
        return [d @ self.COEF for d in _dmono(self.MONV, X)]

    def pbasis(self, X):
        return _mono(self.MONP, X)


# ------------------------------------------------------------------------------ meshes
def _orient(V, T):
    P = V[T]
    det = np.einsum('ti,ti->t', P[:, 1] - P[:, 0], np.cross(P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]))
    T = T.copy(); T[det < 0] = T[det < 0][:, [0, 2, 1, 3]]
    return T


def alfeld(V, T):
    """Alfeld split: each tet -> 4 tets sharing its barycentre (boundary faces unchanged)."""
    nv = len(V)
    V2 = np.vstack([V, V[T].mean(1)])
    c = nv + np.arange(len(T))
    faces = [(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)]
    T2 = np.concatenate([np.stack([T[:, a], T[:, b], T[:, d], c], 1) for (a, b, d) in faces])
    return V2, _orient(V2, T2)


def _split_prisms(bot, top, tri):
    """prisms tri x [bot, top] -> 3 tets each; conforming via the min-index rule (sorted vertices)."""
    s = np.sort(tri, 1)
    a = [bot[s[:, i]] for i in range(3)]; b = [top[s[:, i]] for i in range(3)]
    return np.concatenate([np.stack([a[0], a[1], a[2], b[2]], 1),
                           np.stack([a[0], a[1], b[1], b[2]], 1),
                           np.stack([a[0], b[0], b[1], b[2]], 1)])


def cube_surface(N, L):
    """equiangular N x N grid on each face of the cube [-L,L]^3; returns unique points and triangles
    (diagonal of each quad chosen so that the radial projection to the unit sphere is CONVEX)."""
    t = L * np.tan(np.linspace(-np.pi / 4, np.pi / 4, N + 1))
    pts, quads = [], []
    for ax in range(3):
        for sgn in (-1, 1):
            U, W = np.meshgrid(t, t, indexing="ij")
            P = np.zeros((N + 1, N + 1, 3))
            o = [d for d in range(3) if d != ax]
            P[..., ax] = sgn * L; P[..., o[0]] = U; P[..., o[1]] = W
            base = len(pts) and sum(len(p) for p in pts)
            pts.append(P.reshape(-1, 3))
            idx = base + np.arange((N + 1) ** 2).reshape(N + 1, N + 1)
            q = np.stack([idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]], -1).reshape(-1, 4)
            quads.append(q)
    P = np.vstack(pts); Q = np.vstack(quads)
    key = np.round(P / L * 1e9).astype(np.int64)
    _, uidx, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    P = P[uidx]; Q = inv.ravel()[Q]
    S = P / np.linalg.norm(P, axis=1)[:, None]
    # choose diagonal: convex on the sphere <=> the 4th vertex is below the plane of each triangle
    A, B, C, D = (S[Q[:, i]] for i in range(4))
    nrm = np.cross(B - A, C - A); nrm *= np.sign(np.einsum('ij,ij->i', nrm, A))[:, None]
    diag_ac = np.einsum('ij,ij->i', nrm, D - A) <= 0      # split along A-C keeps D inside
    T1 = np.where(diag_ac[:, None], Q[:, [0, 1, 2]], Q[:, [0, 1, 3]])
    T2 = np.where(diag_ac[:, None], Q[:, [0, 2, 3]], Q[:, [1, 2, 3]])
    return P, np.vstack([T1, T2])


def sphere_box_mesh(N, L=2.0, m=None):
    """Omega_h = (-L,L)^3 minus the inscribed polyhedron P_h (vertices on the unit sphere).
    Layers are log-spaced along rays: r_j = |b|^(j/m); m defaults to N//2."""
    Pb, tri = cube_surface(N, L)
    rb = np.linalg.norm(Pb, axis=1); S = Pb / rb[:, None]
    if m is None:
        m = max(1, N // 2)                                    # self-similar family for even N
    ns = len(Pb)
    V = np.concatenate([S * (rb ** (j / m))[:, None] for j in range(m + 1)])
    V[m * ns:] = Pb                                           # exact box points
    T = np.concatenate([_split_prisms(j * ns + np.arange(ns), (j + 1) * ns + np.arange(ns), tri) for j in range(m)])
    T = _orient(V, T)
    return V, T, dict(N=N, L=L, m=m, hole=True, gamma_tris=tri)


def box_mesh(N, L=1.0):
    """structured (-L,L)^3, each cube -> 6 tets (Kuhn), no hole."""
    x = np.linspace(-L, L, N + 1)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    V = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)
    idx = np.arange((N + 1) ** 3).reshape(N + 1, N + 1, N + 1)
    c = [idx[i:N + i, j:N + j, l:N + l].ravel() for i in (0, 1) for j in (0, 1) for l in (0, 1)]
    # corner (i,j,l) -> c[4i+2j+l]
    perms = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
    T = []
    for p in perms:
        path = [np.zeros(3, int)]
        for d in p:
            nx = path[-1].copy(); nx[d] = 1; path.append(nx)
        T.append(np.stack([c[4 * q[0] + 2 * q[1] + q[2]] for q in path], 1))
    T = _orient(V, np.concatenate(T))
    return V, T, dict(N=N, L=L, hole=False)


# ------------------------------------------------------------------------------ FE space
class Space:
    def __init__(self, V, T, info, R):
        self.R, self.V, self.T, self.info = R, V, T, info
        k, nl, nt = R.k, R.nl, len(T)
        lam = np.array(R.LBARY)                                        # (nl,4)
        # global node key: sorted (vertex, lambda) pairs of the support
        gv = np.broadcast_to(T[:, None, :], (nt, nl, 4)).astype(np.int64)
        lm = np.broadcast_to(lam[None], (nt, nl, 4)).astype(np.int64)
        gvm = np.where(lm > 0, gv, np.int64(1) << 40)
        o = np.argsort(gvm, axis=2)
        gk = np.take_along_axis(gvm, o, 2); lk = np.take_along_axis(np.where(lm > 0, lm, 0), o, 2)
        key = np.concatenate([gk, lk], 2).reshape(-1, 8)
        _, inv = np.unique(key, axis=0, return_inverse=True)
        self.ids = inv.reshape(nt, nl); self.nn = self.ids.max() + 1
        P = V[T]
        self.P0 = P[:, 0]
        self.J = np.stack([P[:, 1] - P[:, 0], P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]], 2)
        self.det = np.linalg.det(self.J); assert (self.det > 0).all()
        self.Jinv = np.linalg.inv(self.J)
        xyz = np.zeros((self.nn, 3))
        xyz[self.ids.ravel()] = (self.P0[:, None, :] + np.einsum('tij,qj->tqi', self.J, R.NODES)).reshape(-1, 3)
        self.xyz = xyz
        E = np.stack([np.linalg.norm(P[:, i] - P[:, j], axis=1) for i in range(4) for j in range(i + 1, 4)], 1)
        self.h = E.max()
        # boundary faces (opposite local vertex i)
        FL = [(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)]
        allf = np.concatenate([np.sort(T[:, f], 1) for f in FL])
        _, fi, cnt = np.unique(allf, axis=0, return_inverse=True, return_counts=True)
        bnd = cnt[fi.ravel()] == 1
        tt = np.tile(np.arange(nt), 4); ii = np.repeat(np.arange(4), nt)
        bt, bi = tt[bnd], ii[bnd]
        L = info["L"]
        FLa = np.array(FL)
        FV = np.stack([V[T[bt, FLa[bi, j]]] for j in range(3)], 1)
        onbox = np.zeros(len(bt), bool)
        for d in range(3):
            onbox |= np.all(np.abs(np.abs(FV[:, :, d]) - L) < 1e-9, 1)
        onsph = np.all(np.abs(np.linalg.norm(FV, axis=2) - 1) < 1e-9, 1)
        assert np.all(onbox ^ onsph), "boundary face neither on box nor on sphere"
        self.fGamma = (bt[onsph], bi[onsph]); self.fBox = (bt[onbox], bi[onbox])
        if info.get("hole"):
            assert len(self.fGamma[0]) == 2 * len(info["gamma_tris"]) or len(self.fGamma[0]) == len(info["gamma_tris"])
            gF = FV[onsph]
            self.hGamma = max(np.linalg.norm(gF[:, i] - gF[:, j], axis=1).max() for i, j in ((0, 1), (1, 2), (0, 2)))
        else:
            self.hGamma = np.nan
        on = np.zeros(self.nn, bool)
        for d in range(3):
            on |= np.abs(np.abs(xyz[:, d]) - L) < 1e-9
        self.boxnodes = np.where(on)[0]
        self.dofs = 3 * self.ids[:, :, None] + np.arange(3)[None, None, :]    # (nt,nl,3)

    def phys(self, t, X):
        return self.P0[t][:, None, :] + np.einsum('tij,qj->tqi', self.J[t], X)

    def face_data(self, faces):
        """per boundary face: tet, reference points, physical points, weights, outward normal (all vectorised)."""
        R = self.R
        bt, bi = faces
        FL = [(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)]
        out = []
        for i in range(4):
            sel = np.where(bi == i)[0]
            if len(sel) == 0:
                continue
            t = bt[sel]
            a, b, c = (Ref.REFV[j] for j in FL[i])
            ref = a[None] + R.FX[:, :1] * (b - a)[None] + R.FX[:, 1:] * (c - a)[None]
            Vt = self.V[self.T[t]]
            A, B, C = Vt[:, FL[i][0]], Vt[:, FL[i][1]], Vt[:, FL[i][2]]
            nv = np.cross(B - A, C - A); area2 = np.linalg.norm(nv, axis=1); n = nv / area2[:, None]
            s = np.sign(np.einsum('ti,ti->t', n, Vt[:, i] - A)); n = -s[:, None] * n
            X = A[:, None] + R.FX[None, :, :1] * (B - A)[:, None] + R.FX[None, :, 1:] * (C - A)[:, None]
            w = R.FW[None, :] * area2[:, None]
            diam = np.max(np.stack([np.linalg.norm(B - A, axis=1), np.linalg.norm(C - B, axis=1),
                                    np.linalg.norm(A - C, axis=1)], 1), 1)
            out.append(dict(t=t, ref=ref, X=X, w=w, n=n, diam=diam))
        return out

    def pgrads(self, t, dref):
        """physical gradients: (nt_chunk, nq, nl, 3)"""
        D = np.stack(dref, -1)                                  # (nq,nl,3) reference
        return np.einsum('tim,qai->tqam', self.Jinv[t], D)


def build(kind, N, k, **kw):
    R = Ref(k)
    if kind == "box":
        V, T, info = box_mesh(N, **kw)
    else:
        V, T, info = sphere_box_mesh(N, **kw)
    V, T = alfeld(V, T)
    return Space(V, T, info, R)


# ------------------------------------------------------------------------------ assembly
def assemble(S, mu, f, nitsche, pcouple):
    """nitsche: list of (faces, g or None, flux 'grad'|'sym');  pcouple: list of (faces, mean_subtract)."""
    R = S.R; nl, npl = R.nl, R.npl; n3 = 3 * nl; nt = len(S.T); ndof = 3 * S.nn
    D = S.dofs.reshape(nt, n3)
    rows, cols, vals = [], [], []
    for s in range(0, nt, CH):
        t = np.arange(s, min(nt, s + CH))
        G = np.einsum('tim,tjn,ijab->tmnab', S.Jinv[t], S.Jinv[t], R.K) * S.det[t][:, None, None, None, None]
        loc = np.zeros((len(t), nl, 3, nl, 3))
        tr = G[:, 0, 0] + G[:, 1, 1] + G[:, 2, 2]
        for c in range(3):
            loc[:, :, c, :, c] += tr
            for d in range(3):
                loc[:, :, c, :, d] += G[:, d, c]
        Dt = D[t]
        rows.append(np.repeat(Dt, n3, 1).ravel()); cols.append(np.tile(Dt, (1, n3)).ravel()); vals.append(loc.ravel())
    A = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), (ndof, ndof)).tocsr()
    del rows, cols, vals
    # divergence  B[(t,j),(a,c)] = int q_j d_c phi_a
    Bl = np.einsum('tic,jia->tjac', S.Jinv, R.Bq) * S.det[:, None, None, None]
    prow = np.arange(nt)[:, None] * npl + np.arange(npl)[None, :]
    B = sp.coo_matrix((Bl.reshape(nt, npl, n3).ravel(), (np.repeat(prow, n3, 1).ravel(), np.tile(D, (1, npl)).ravel())),
                      (npl * nt, ndof)).tocsr()
    del Bl
    Minv = sp.block_diag([R.Mpinv / d for d in S.det], format="csr")
    # load
    F = np.zeros(ndof)
    PHI = R.basis(R.QX)
    for s in range(0, nt, CH):
        t = np.arange(s, min(nt, s + CH))
        X = S.phys(t, R.QX)
        fv = f(X[..., 0], X[..., 1], X[..., 2])
        Fl = np.einsum('q,t,qa,tqc->tac', R.QW, S.det[t], PHI, fv)
        np.add.at(F, S.dofs[t].reshape(-1), Fl.reshape(-1))
    # Nitsche faces
    rN, cN, vN = [], [], []
    for faces, g, flux in nitsche:
        for fd in S.face_data(faces):
            t, w, n = fd["t"], fd["w"], fd["n"]
            ph = R.basis(fd["ref"])                                  # (nq,nl)
            gr = S.pgrads(t, R.dbasis(fd["ref"]))                    # (nf,nq,nl,3)
            dn = np.einsum('tqam,tm->tqa', gr, n)
            M = np.einsum('tq,qa,qb->tab', w, ph, ph)
            Dn = np.einsum('tq,qa,tqb->tab', w, ph, dn)              # int phi_a d_n phi_b
            T4 = np.zeros((len(t), nl, 3, nl, 3))                    # int v.flux(u), v=phi_a e_c, u=phi_b e_d
            for c in range(3):
                T4[:, :, c, :, c] += Dn
            if flux == "sym":
                for c in range(3):
                    for d in range(3):
                        T4[:, :, c, :, d] += np.einsum('tq,qa,tqb->tab', w, ph, gr[..., c]) * n[:, d][:, None, None]
            loc = -T4 - np.transpose(T4, (0, 3, 4, 1, 2))
            for c in range(3):
                loc[:, :, c, :, c] += (mu / S.h) * M
            Dt = S.dofs[t].reshape(len(t), n3)
            rN.append(np.repeat(Dt, n3, 1).ravel()); cN.append(np.tile(Dt, (1, n3)).ravel()); vN.append(loc.ravel())
            if g is not None:
                X = fd["X"]; gv = g(X[..., 0], X[..., 1], X[..., 2])     # (nf,nq,3)
                Fl = (mu / S.h) * np.einsum('tq,qa,tqc->tac', w, ph, gv) - np.einsum('tq,tqa,tqc->tac', w, dn, gv)
                if flux == "sym":
                    gdot = np.einsum('tqam,tqm->tqa', gr, gv)
                    Fl -= np.einsum('tq,tqa,tc->tac', w, gdot, n)
                np.add.at(F, S.dofs[t].reshape(-1), Fl.reshape(-1))
    if rN:
        A = A + sp.coo_matrix((np.concatenate(vN), (np.concatenate(rN), np.concatenate(cN))), (ndof, ndof)).tocsr()
    # pressure traction coupling C[(t,j),(a,c)] = int q_j phi_a n_c
    Cs = []
    for faces, mean in pcouple:
        rC, cC, vC = [], [], []
        gG = np.zeros(npl * nt); area = 0.0
        for fd in S.face_data(faces):
            t, w, n = fd["t"], fd["w"], fd["n"]
            ph = R.basis(fd["ref"]); q = R.pbasis(fd["ref"])
            Cl = np.einsum('tq,qj,qa,tc->tjac', w, q, ph, n)
            rC.append(np.repeat(prow[t], n3, 1).ravel()); cC.append(np.tile(S.dofs[t].reshape(len(t), n3), (1, npl)).ravel())
            vC.append(Cl.ravel())
            np.add.at(gG, prow[t].ravel(), (w @ q).ravel()); area += w.sum()
        C = sp.coo_matrix((np.concatenate(vC), (np.concatenate(rC), np.concatenate(cC))), (npl * nt, ndof)).tocsr()
        Cs.append((C, gG / area if mean else None))
    return dict(A=A, B=B, Minv=Minv, F=F, C=Cs)


# ------------------------------------------------------------------------------ solve (iterated penalty)
class _SymSolver:
    """PARDISO symmetric indefinite (mtype -2) on the upper triangle; solve(K, b) with the full K."""
    def __init__(self, K):
        import pypardiso
        self.U = sp.triu(K, format="csr"); self.U.sort_indices()
        self.s = pypardiso.PyPardisoSolver(mtype=-2)
        self.s.factorize(self.U)

    def solve(self, K, b):
        return self.s.solve(self.U, b)

    def free_memory(self, everything=True):
        self.s.free_memory(everything=everything)


def _factor(K, sym):
    K = K.tocsr(); K.sort_indices()
    if sym:
        return _SymSolver(K), K
    import pypardiso
    slv = pypardiso.PyPardisoSolver(mtype=11)
    slv.set_iparm(11, 1); slv.set_iparm(13, 1)
    slv.factorize(K)
    return slv, K


def solve_ipm(S, sysm, dirichlet=None, r=1e4, tol=1e-13, rtol=1e-11, maxit=200, verbose=False, symK=True):
    """dirichlet: function g(x,y,z) imposed strongly at box nodes (None: no strong BC).
    symK=True: iteration matrix K = A + r B^T M^-1 B (symmetric; PARDISO LDL^T, half the memory); the
    CNS* traction coupling Bt enters only the right-hand side Bt^T p^n (same fixed point, since B u = 0 there).
    symK=False: K = A + r Bt^T M^-1 B (nonsymmetric LU; rank-one pbar part by Sherman-Morrison)."""
    A, B, Minv, F = sysm["A"], sysm["B"], sysm["Minv"], sysm["F"]
    ndof = A.shape[0]
    # Bt = B - sum C + sum gG (x) phi ; phi = C^T e (constant pressure modes)
    npl = S.R.npl
    Bt = B.copy()
    rank1 = []
    for C, gG in sysm["C"]:
        Bt = Bt - C
        if gG is not None:
            e = np.zeros(B.shape[0]); e[0::npl] = 1.0          # MONP[0] = (0,0,0) constant
            rank1.append((gG, C.T @ e))
    Bt = Bt.tocsr()
    if dirichlet is not None:
        dn = S.boxnodes
        dD = (3 * dn[:, None] + np.arange(3)).ravel()
        ug = np.zeros(ndof); ug[dD] = dirichlet(S.xyz[dn, 0], S.xyz[dn, 1], S.xyz[dn, 2]).ravel()
    else:
        dD = np.zeros(0, int); ug = np.zeros(ndof)
    free = np.setdiff1d(np.arange(ndof), dD)
    MB = (Minv @ B).tocsr()
    K0 = A + r * ((B if symK else Bt).T @ MB)
    if symK:
        K0 = (K0 + K0.T) / 2
    Kff = K0[free][:, free]
    b0 = F - K0[:, dD] @ ug[dD] if len(dD) else F.copy()
    # rank-one pieces: K = K0 + r * phi (gG^T M^-1 B)
    Z = [] if symK else [(r * ph, (gG @ MB)) for gG, ph in rank1]
    slv, Kc = _factor(Kff, sym=symK)

    def solve0(rhs):
        x = slv.solve(Kc, rhs)
        for _ in range(2):
            x = x + slv.solve(Kc, rhs - Kc @ x)
        return x
    corr = []
    for (ph, z) in Z:
        y = solve0(ph[free]); corr.append((y, z[free], ph[free]))
        b0 = b0 - ph * (z[dD] @ ug[dD]) if len(dD) else b0

    def solveK(rhs):
        x = solve0(rhs)
        for (y, z, ph) in corr:
            x = x - y * (z @ x) / (1 + z @ y)
        return x
    p = np.zeros(B.shape[0]); u = ug.copy(); hist = []; rhist = []

    def momres(u, p):
        res = (A @ u - Bt.T @ p - F)[free]
        for (gG, ph) in rank1:
            res = res - (ph * (gG @ p))[free]
        return np.linalg.norm(res) / max(np.linalg.norm(b0[free]), np.linalg.norm((A @ u)[free]), 1e-300)
    for it in range(maxit):
        rhs = b0 + Bt.T @ p
        for (gG, ph) in rank1:
            rhs = rhs + ph * (gG @ p)
        u[free] = solveK(rhs[free])
        Bu = B @ u
        p = p - r * (Minv @ Bu)
        dres = np.sqrt(abs(Bu @ (Minv @ Bu)))
        relres = momres(u, p)
        hist.append(dres); rhist.append(relres)
        if verbose:
            print(f"   ipm it {it} ||div u||={dres:.3e} relres={relres:.3e}", flush=True)
        if relres < rtol and dres < tol:
            break
        if it >= 8 and min(rhist[-5:]) > 0.7 * min(rhist[:-5]) and min(hist[-5:]) > 0.7 * min(hist[:-5]):
            break                                              # stagnated at round-off
    slv.free_memory(everything=True)
    return u, p, hist, relres


# ------------------------------------------------------------------------------ errors
def errors(S, u, ex, mu, nq=None):
    R = S.R; nt = len(S.T)
    Xq, Wq = (R.QX, R.QW) if nq is None else tet_quad(nq)
    PHI = R.basis(Xq); DR = R.dbasis(Xq)
    H1 = L2 = 0.0
    for s in range(0, nt, CH):
        t = np.arange(s, min(nt, s + CH))
        X = S.phys(t, Xq)
        U = u[S.dofs[t]]                                       # (nt,nl,3)
        gr = S.pgrads(t, DR)
        guh = np.einsum('tqam,tac->tqcm', gr, U)
        G = ex["gu"](X[..., 0], X[..., 1], X[..., 2])
        w = Wq[None] * S.det[t][:, None]
        H1 += (w * ((G - guh) ** 2).sum((-1, -2))).sum()
        uh = np.einsum('qa,tac->tqc', PHI, U)
        L2 += (w * ((ex["u"](X[..., 0], X[..., 1], X[..., 2]) - uh) ** 2).sum(-1)).sum()
    out = dict(H1=np.sqrt(H1), L2=np.sqrt(L2))
    if not S.info.get("hole"):
        return out
    bl2 = bdn = 0.0; area = 0.0; pint = 0.0; recs = []
    for fd in S.face_data(S.fGamma):
        t, w, n, X = fd["t"], fd["w"], fd["n"], fd["X"]
        ph = R.basis(fd["ref"]); gr = S.pgrads(t, R.dbasis(fd["ref"]))
        U = u[S.dofs[t]]
        uh = np.einsum('qa,tac->tqc', ph, U)
        dnh = np.einsum('tqam,tm,tac->tqc', gr, n, U)
        ue = ex["u"](X[..., 0], X[..., 1], X[..., 2]); Ge = ex["gu"](X[..., 0], X[..., 1], X[..., 2])
        dne = np.einsum('tqcm,tm->tqc', Ge, n)
        pe = ex["p"](X[..., 0], X[..., 1], X[..., 2])
        bl2 += (w * ((ue - uh) ** 2).sum(-1)).sum()
        bdn += (fd["diam"][:, None] * w * ((dne - dnh) ** 2).sum(-1)).sum()
        area += w.sum(); pint += (w * pe).sum()
        recs.append((w, pe, np.einsum('tqc,tc->tq', uh, n), np.einsum('tqc,tc->tq', ue - uh, n)))
    pbar = pint / area
    G2 = sum((w * (pe - pbar) ** 2).sum() for w, pe, _, _ in recs)
    un2 = sum((w * un ** 2).sum() for w, _, un, _ in recs)
    leak = sum((w * un * (pe - pbar)).sum() for w, pe, un, _ in recs)
    en = np.sqrt(H1 + (mu / S.h) * bl2 + bdn)
    lam = S.h / mu; G = np.sqrt(G2)
    out.update(bL2=np.sqrt(bl2), bdn=np.sqrt(bdn), energy=en, G=G, pbar=pbar, area=area,
               unL2=np.sqrt(un2))
    if G > 1e-12:
        out.update(slip_ratio=np.sqrt(bl2) / (lam * G), un_ratio=np.sqrt(un2) / (lam * G),
                   leak_corr=leak / (lam * G2), H1_ratio=np.sqrt(H1) / (lam * G),
                   energy_ratio=en * np.sqrt(mu / S.h) / G)
    return out


def div_norm(S, sysm, u):
    Bu = sysm["B"] @ u
    return np.sqrt(abs(Bu @ (sysm["Minv"] @ Bu)))


# ------------------------------------------------------------------------------ penalty threshold
def local_threshold(S, faces, flux="grad", tol=1e-3):
    """Sufficient coercivity threshold: smallest mu such that on EVERY boundary micro-tet T with face F
    a_T(v,v) - 2<flux(v), v>_F + (mu/h)||v||_F^2 >= 0 for all v in P_k(T)^3 (batched bisection).
    Since a_h = sum_T a_T with a_T >= 0 and each Alfeld micro-tet has at most one boundary face,
    mu >= max_T mu_T implies N_h(v,v) >= 0 on V_h.  Returns the per-face thresholds."""
    R = S.R; nl = R.nl; n3 = 3 * nl
    res = []
    for fd in S.face_data(faces):
        t, w, n = fd["t"], fd["w"], fd["n"]
        G = np.einsum('tim,tjn,ijab->tmnab', S.Jinv[t], S.Jinv[t], R.K) * S.det[t][:, None, None, None, None]
        Aloc = np.zeros((len(t), nl, 3, nl, 3)); tr = G[:, 0, 0] + G[:, 1, 1] + G[:, 2, 2]
        for c in range(3):
            Aloc[:, :, c, :, c] += tr
            for d in range(3):
                Aloc[:, :, c, :, d] += G[:, d, c]
        ph = R.basis(fd["ref"]); gr = S.pgrads(t, R.dbasis(fd["ref"]))
        dn = np.einsum('tqam,tm->tqa', gr, n)
        M = np.einsum('tq,qa,qb->tab', w, ph, ph); Dn = np.einsum('tq,qa,tqb->tab', w, ph, dn)
        T4 = np.zeros_like(Aloc)
        for c in range(3):
            T4[:, :, c, :, c] += Dn
        if flux == "sym":
            for c in range(3):
                for d in range(3):
                    T4[:, :, c, :, d] += np.einsum('tq,qa,tqb->tab', w, ph, gr[..., c]) * n[:, d][:, None, None]
        Mf = np.zeros_like(Aloc)
        for c in range(3):
            Mf[:, :, c, :, c] = M
        K0 = (Aloc - T4 - np.transpose(T4, (0, 3, 4, 1, 2))).reshape(-1, n3, n3)
        Mf = Mf.reshape(-1, n3, n3) / S.h
        lo = np.zeros(len(t)); hi = np.full(len(t), 1e7)
        scale = np.abs(K0).max((1, 2))[:, None, None]
        while np.max(hi / np.maximum(lo, 1e-9)) > 1 + tol:
            mid = np.sqrt(np.maximum(lo, 1e-3) * hi)
            ev = np.linalg.eigvalsh((K0 + mid[:, None, None] * Mf) / scale)[:, 0]
            ok = ev > -1e-11
            hi = np.where(ok, mid, hi); lo = np.where(ok, lo, mid)
        res.append(hi)
    return np.concatenate(res)
