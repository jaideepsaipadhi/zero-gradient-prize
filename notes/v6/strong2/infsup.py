"""Discrete inf-sup constants for (V_O, Pi_h cap L^2_0), P4-P3disc Scott-Vogelius, strong no-slip space
V_O = {v in V_h : v = 0 on dR and on Gamma_h}.

  beta^2 = min_{q in Q, q != 0}  (B A^{-1} B^T q, q) / (M q, q),    A = vector Laplacian (grad:grad) on V_O

computed on three pressure subspaces Q (all mean-free):
  'full'  : Pi_h cap L^2_0
  'lockw' : Q_lock  = {q : q|T1(z) = q|T2(z)  at every locked z}     (Theorem L', one condition per lock)
  'locks' : Q_lock0 = {q : q|T1(z) = q|T2(z) = 0 at every locked z}  (hypothesis (L) of notes/v6/strong)
A locked vertex = Gamma_h vertex lying in exactly 2 triangles.
Meshes: std (every wall vertex in 3 triangles), one (one flipped diagonal: 1 locked vertex),
        alt (mesh (a): every other wall vertex locked).
Usage: python3 infsup.py MESH N [N ...]        (dense eigensolver; N <= 32 is a few seconds..minutes)
       python3 infsup.py MESH N --lobpcg [--nev=8]   (iterative, two random starts, for larger N)
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, scipy.linalg as sla
import svn, strong_bc

make_mesh_alt = strong_bc.make_mesh_alt


def make_mesh_one(N, L=2.5):        # identical to notes/v6/strong/one_lock.py
    verts, tris, circ, outer = svn.make_mesh(N, L=L)
    tris = tris.copy(); m = N // 4
    vid = np.arange((m + 1) * N).reshape(m + 1, N)
    i = 1; i1 = 2
    a, b, c, d = vid[0, i], vid[0, i1], vid[1, i1], vid[1, i]
    tris[2 * i] = [a, b, d]; tris[2 * i + 1] = [b, c, d]
    P = verts[tris]
    ar = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    tris[ar < 0] = tris[ar < 0][:, [0, 2, 1]]
    return verts, tris, circ, outer


def make_mesh_fan(N, L=2.5, i=2):
    """two CONSECUTIVE locked wall vertices z_i, z_{i+1} sharing the apex w = d_{i+1} (layer-1 vertex):
    quads i-1, i, i+1 of layer 0 are retriangulated as
      (z_{i-1}, z_i, w), (z_i, z_{i+1}, w), (z_{i+1}, z_{i+2}, w)            [fan: z_i, z_{i+1} in 2 triangles]
      (z_{i-1}, w, d_i), (z_{i-1}, d_i, d_{i-1}), (z_{i+2}, d_{i+2}, w)."""
    verts, tris, circ, outer = svn.make_mesh(N, L=L)
    tris = tris.copy(); m = N // 4
    vid = np.arange((m + 1) * N).reshape(m + 1, N)
    z = lambda j: vid[0, j % N]; d = lambda j: vid[1, j % N]
    w = d(i + 1)
    new = [[z(i - 1), z(i), w], [z(i), z(i + 1), w], [z(i + 1), z(i + 2), w],
           [z(i - 1), w, d(i)], [z(i - 1), d(i), d(i - 1)], [z(i + 2), d(i + 2), w]]
    for k, q in enumerate([i - 1, i, i + 1]):
        tris[2 * q] = new[2 * k]; tris[2 * q + 1] = new[2 * k + 1]
    P = verts[tris]
    ar = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    tris[ar < 0] = tris[ar < 0][:, [0, 2, 1]]
    return verts, tris, circ, outer


def min_angle(verts, tris):
    P = verts[tris]; out = []
    for a in range(3):
        u = P[:, (a + 1) % 3] - P[:, a]; v = P[:, (a + 2) % 3] - P[:, a]
        out.append(np.degrees(np.arccos(np.sum(u * v, 1) / np.linalg.norm(u, axis=1) / np.linalg.norm(v, axis=1))))
    return float(np.min(out))


MESHES = dict(std=lambda N: svn.make_mesh(N), one=make_mesh_one, alt=make_mesh_alt, fan=make_mesh_fan)
REFV = np.array([[0., 0.], [1., 0.], [0., 1.]])


def laplace(S):
    nt = len(S.tris)
    gx, gy = S.grads(svn.QX, svn.QY)
    wq = svn.QW[None, :] * np.abs(S.detJ)[:, None]
    K = np.einsum('tq,tqa,tqb->tab', wq, gx, gx) + np.einsum('tq,tqa,tqb->tab', wq, gy, gy)
    loc = np.zeros((nt, 15, 2, 15, 2)); loc[:, :, 0, :, 0] = K; loc[:, :, 1, :, 1] = K
    D = svn.dofs(S.ids).reshape(nt, 30)
    rows = np.repeat(D, 30, axis=1).ravel(); cols = np.tile(D, (1, 30)).ravel()
    return sp.coo_matrix((loc.reshape(nt, 30, 30).ravel(), (rows, cols)), shape=(2 * S.nn,) * 2).tocsr()


def locked(S, circ):
    """list of (z, [(t, a), (t', a')]) for Gamma_h vertices in exactly two triangles"""
    inc = {}
    for t, tri in enumerate(S.tris):
        for a, v in enumerate(tri):
            if v in circ:
                inc.setdefault(int(v), []).append((t, a))
    return [(z, L) for z, L in sorted(inc.items()) if len(L) == 2], inc


def setup(mesh, N):
    verts, tris, circ, outer = MESHES[mesh](N)
    S = svn.Space(verts, tris, circ, outer)
    sysm = svn.assemble(S, 1.0, lambda X, Y: np.zeros(X.shape + (2,)), None)
    return S, sysm, circ


def constraint_vectors(S, sysm, circ):
    """Riesz vectors (w.r.t. M) of the vertex functionals: returns (Y_const, list of (rT1, rT2)) """
    nt = len(S.tris); npr = 10 * nt
    lk, inc = locked(S, circ)
    Ml = sysm["Ml"]
    one = np.zeros(npr); one[0::10] = 1.0
    reps = []
    for z, L in lk:
        pair = []
        for (t, a) in L:
            e = svn._mono(svn.MON3, REFV[a:a + 1, 0], REFV[a:a + 1, 1])[0]
            r = np.zeros(npr); r[10 * t:10 * t + 10] = np.linalg.solve(Ml[t], e)
            pair.append(r)
        reps.append(pair)
    return one, reps, lk


def infsup(mesh, N, method="dense", nev=4):
    t0 = time.time()
    S, sysm, circ = setup(mesh, N)
    A = laplace(S); B = sysm["B"]
    ndof = A.shape[0]
    gn = strong_bc.gamma_nodes(S)
    dD = np.union1d(svn.dofs(S.dnodes).ravel(), svn.dofs(gn).ravel())
    free = np.setdiff1d(np.arange(ndof), dD)
    Af = A[free][:, free].tocsc(); Bf = B[:, free].tocsr()
    M = sp.block_diag(list(sysm["Ml"]), format="csr")
    one, reps, lk = constraint_vectors(S, sysm, circ)
    Ys = dict(full=[one], lockw=[one] + [r1 - r2 for r1, r2 in reps],
              locks=[one] + [r for pr in reps for r in pr])
    lu = spl.splu(Af, permc_spec="COLAMD")
    out = dict(mesh=mesh, N=N, hG=S.hGamma, h=S.h, nlocked=len(lk), npr=B.shape[0], ndof=len(free), minang=min_angle(S.verts, S.tris))
    # turning angle at locked vertices
    V = S.verts
    if lk:
        phis = []
        for z, _ in lk:
            nb = [v for v in range(N) if abs(v - z) % N in (1, N - 1)]
            a, b = V[nb[0]] - V[z], V[nb[1]] - V[z]
            phis.append(np.pi - np.arccos(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b)))
        out["phi_min"] = float(min(phis))
    if method == "dense":
        X = lu.solve(Bf.T.toarray())
        Sd = Bf @ X; Sd = 0.5 * (Sd + Sd.T)
        Md = M.toarray()
        for name, Y in Ys.items():
            Y = np.array(Y).T                      # npr x c
            Q, _ = np.linalg.qr(Md @ Y, mode="complete")
            Z = Q[:, Y.shape[1]:]                  # M-orthogonal complement of span Y
            ev = sla.eigh(Z.T @ Sd @ Z, Z.T @ Md @ Z, eigvals_only=True, subset_by_index=[0, nev - 1])
            out["beta_" + name] = float(np.sqrt(max(ev[0], 0)))
            out["ev_" + name] = [float(e) for e in ev]
    else:
        Sop = spl.LinearOperator(M.shape, matvec=lambda q: Bf @ lu.solve(Bf.T @ q), dtype=float)
        Minv = sysm["Minv"]
        for name, Y in Ys.items():
            Y = np.array(Y).T
            evs = []
            for seed in (0, 1):                    # two random starts; keep the smaller spectrum
                X0 = np.random.default_rng(seed).standard_normal((M.shape[0], nev))
                e, _ = spl.lobpcg(Sop, X0, B=M, M=Minv, Y=Y, largest=False, tol=1e-9, maxiter=3000)
                evs.append(np.sort(e))
            ev = min(evs, key=lambda e: e[0])
            out["beta_" + name] = float(np.sqrt(max(ev[0], 0))); out["ev_" + name] = [float(e) for e in ev]
    out["secs"] = time.time() - t0
    return out


if __name__ == "__main__":
    mesh = sys.argv[1]
    meth = "lobpcg" if "--lobpcg" in sys.argv else "dense"
    nev = int(([a.split("=")[1] for a in sys.argv if a.startswith("--nev=")] or ["4"])[0])
    for N in [int(a) for a in sys.argv[2:] if not a.startswith("--")]:
        r = infsup(mesh, N, meth, nev)
        r["method"] = meth; r["nev"] = nev
        print(json.dumps(r), flush=True)
