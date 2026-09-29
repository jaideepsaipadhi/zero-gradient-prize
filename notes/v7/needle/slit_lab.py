"""Interior slit lab: P4/P3disc SV on [-1,1]^2 (v=0 on the square), uniform diagonal mesh, one interior
vertex w split along the horizontal edge line z-w-f into x_up, x_dn = w +- (delta/2) e_y  (two needles
N=[z,x_dn,x_up] (apex angle eps at z) and N'=[f,x_up,x_dn] sharing the short edge).
Functions: build(n, eps), schur(S) -> (Sd, M, one), beta/mode, functional tests."""
import os, sys
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "v6", "strong2"))
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, scipy.linalg as sla
import svn
from infsup import laplace

def square_mesh(n, L=1.0):
    xs = np.linspace(-L, L, n + 1)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    verts = np.stack([X.ravel(), Y.ravel()], 1)
    vid = np.arange((n + 1) ** 2).reshape(n + 1, n + 1)
    tris = []
    for i in range(n):
        for j in range(n):
            a, b, c, d = vid[i, j], vid[i + 1, j], vid[i + 1, j + 1], vid[i, j + 1]
            if True:  # uniform diagonal: every interior vertex has 6 triangles on 3 lines
                tris += [[a, b, c], [a, c, d]]
            else:
                tris += [[a, b, d], [b, c, d]]
    return verts, np.array(tris), vid

def orient(verts, tris):
    tris = np.array(tris); P = verts[tris]
    ar = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]); tris[ar < 0] = tris[ar < 0][:, [0, 2, 1]]
    return tris

def split(verts, tris, w, z, f, eps, shift=0.0):
    """split w along line z-w-f; delta chosen so the apex angle at z is eps. shift moves the split point
    along the line (asymmetric kite)."""
    d = verts[f] - verts[z]; d = d / np.linalg.norm(d); nrm = np.array([-d[1], d[0]])
    r = np.linalg.norm(verts[w] - verts[z])
    delta = 2 * r * np.tan(eps / 2)
    up = len(verts); dn = w
    V = np.vstack([verts, verts[w]])
    newt = []
    for tri in tris:
        tri = list(tri)
        if w in tri:
            c = verts[tri].mean(0) - verts[w]
            if np.dot(c, nrm) > 0:
                tri = [up if v == w else v for v in tri]
        newt.append(tri)
    V[up] = verts[w] + shift * d + 0.5 * delta * nrm
    V[dn] = verts[w] + shift * d - 0.5 * delta * nrm
    newt += [[z, dn, up], [f, up, dn]]
    return V, orient(V, newt), up, dn

def build(n, eps, shift=0.0, which=None):
    verts, tris, vid = square_mesh(n)
    i0 = n // 2 if which is None else which
    w = vid[i0, n // 2]; z = vid[i0 - 1, n // 2]; f = vid[i0 + 1, n // 2]
    V, T, up, dn = split(verts, tris, w, z, f, eps, shift)
    S = svn.Space(V, T, {int(vid[0, 0]), int(vid[1, 0])}, None)   # dummy Gamma edge on the square (unused)
    return S, dict(z=z, f=f, up=up, dn=dn)

def schur(S):
    sysm = svn.assemble(S, 1.0, lambda X, Y: np.zeros(X.shape + (2,)), None, want=("A", "B", "M"))
    A = laplace(S); B = sysm["B"]
    free = np.setdiff1d(np.arange(A.shape[0]), svn.dofs(S.dnodes).ravel())
    Af = A[free][:, free].tocsc(); Bf = B[:, free].tocsr()
    lu = spl.splu(Af)
    X = lu.solve(Bf.T.toarray()); Sd = Bf @ X; Sd = .5 * (Sd + Sd.T)
    M = sp.block_diag(list(sysm["Ml"]), format="csr").toarray()
    one = np.zeros(M.shape[0]); one[0::10] = 1
    return Sd, M, one, sysm, (lu, Af, Bf, free)

def constrained_beta(Sd, M, one, cons=(), k=2):
    Y = np.array([M @ one * 0 + one] + list(cons)).T
    # constraints are given as representers r (q -> r^T M q); build basis of {q: (r,q)_M=0}
    Q, _ = np.linalg.qr(M @ Y, mode="complete"); Z = Q[:, Y.shape[1]:]
    ev, Vv = sla.eigh(Z.T @ Sd @ Z, Z.T @ M @ Z, subset_by_index=[0, k - 1])
    return np.sqrt(np.maximum(ev, 0)), Z @ Vv

def dual_norm_ratio(Sd, M, r):
    """for a pressure rho (coefficient vector r): sup_v (rho, div v)/(|rho| |grad v|) = sqrt(r^T Sd r / r^T M r)"""
    return float(np.sqrt(max(r @ Sd @ r, 0) / (r @ M @ r)))

def tri_of(S, a, b, c):
    s = {a, b, c}
    for t, tri in enumerate(S.tris):
        if set(tri.tolist()) == s:
            return t
    raise KeyError

def point_rep(S, sysm, t, x):
    """L2(T) Riesz representer of q -> q|_T(x) (x physical)"""
    ref = np.linalg.solve(S.J[t], x - S.P0[t])
    e = svn._mono(svn.MON3, np.array([ref[0]]), np.array([ref[1]]))[0]
    r = np.zeros(10 * len(S.tris)); r[10 * t:10 * t + 10] = np.linalg.solve(sysm["Ml"][t], e)
    return r

def eval_p(S, q, t, x):
    ref = np.linalg.solve(S.J[t], x - S.P0[t])
    return float(svn._mono(svn.MON3, np.array([ref[0]]), np.array([ref[1]]))[0] @ q[10 * t:10 * t + 10])

def interior_basis(S, layers=1):
    """pressure dofs on triangles not touching the square boundary (removes the exact corner/boundary
    spurious modes of the uniform mesh; no mean constraint needed then)"""
    L = np.abs(S.verts).max()
    onb = (np.abs(np.abs(S.verts[:, 0]) - L) < 1e-12) | (np.abs(np.abs(S.verts[:, 1]) - L) < 1e-12)
    bad = onb.copy()
    for _ in range(layers - 1):
        bad = bad | np.isin(np.arange(len(S.verts)), S.tris[bad[S.tris].any(1)].ravel())
    keep = np.where(~bad[S.tris].any(1))[0]
    cols = (10 * keep[:, None] + np.arange(10)[None, :]).ravel()
    E = np.zeros((10 * len(S.tris), len(cols))); E[cols, np.arange(len(cols))] = 1
    return E

def beta_sub(Sd, M, E, cons=(), k=2):
    """inf-sup over q = E y subject to (r_i, q)_M = 0"""
    if len(cons):
        C = np.array(cons) @ M @ E          # (nc, ny)
        _, s, Vt = np.linalg.svd(C); Nn = Vt[len(cons):].T
        E = E @ Nn
    ev, Vv = sla.eigh(E.T @ Sd @ E, E.T @ M @ E, subset_by_index=[0, k - 1])
    return np.sqrt(np.maximum(ev, 0)), E @ Vv
