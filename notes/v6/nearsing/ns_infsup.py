"""Discrete inf-sup constants of (V_O, Pi_h cap L^2_0) (P4/P3disc, strong no-slip) on nearly singular meshes,
on constrained pressure spaces (dense; small N).  At every modified wall vertex z (star: two regular triangles
Tf, Tl and the needle T*, identified as the triangle with the smallest angle at z):
   full   : Pi_h cap L^2_0
   wired  : q|Tf(z) = q|Tl(z)                       (one condition per z; the lock limit of Q_lock)
   wirefn : q|Tf(z) = q|T*(z)                       (a WRONG wiring, control)
   zero3  : q|Tf(z) = q|T*(z) = q|Tl(z) = 0         (the space the velocity proof needs, cf. (L))
For each space two constants:
   beta  = min_{q in Q} sup_v (q, div v)/(|q| |grad v|)          (dual, = sqrt(min eig (S q,q)/(M q,q)))
   gamma = min_{q in Q} |q| / min{|grad v| : div v = q}           (primal / right-inverse constant)
Usage: python3 ns_infsup.py MESH N eps [eps ...]
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
sys.path.insert(0, os.path.join(HERE, "..", "strong2"))
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "2")
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, scipy.linalg as sla
import svn, strong_bc
import ns_mesh
from infsup import laplace

REFV = np.array([[0., 0.], [1., 0.], [0., 1.]])


def vertex_rep(S, sysm, t, v):
    a = int(np.where(S.tris[t] == v)[0][0])
    e = svn._mono(svn.MON3, REFV[a:a + 1, 0], REFV[a:a + 1, 1])[0]
    r = np.zeros(10 * len(S.tris)); r[10 * t:10 * t + 10] = np.linalg.solve(sysm["Ml"][t], e)
    return r


def run(mesh, N, eps):
    t0 = time.time()
    if mesh == "std":
        verts, tris, circ, outer = svn.make_mesh(N); lk = []
    else:
        verts, tris, circ, outer, lk = ns_mesh.make(mesh, N, eps)
    S = svn.Space(verts, tris, circ, outer)
    sysm = svn.assemble(S, 1.0, lambda X, Y: np.zeros(X.shape + (2,)), None)
    A = laplace(S); B = sysm["B"]
    gn = strong_bc.gamma_nodes(S)
    dD = np.union1d(svn.dofs(S.dnodes).ravel(), svn.dofs(gn).ravel())
    free = np.setdiff1d(np.arange(A.shape[0]), dD)
    Af = A[free][:, free].tocsc(); Bf = B[:, free].tocsr()
    M = sp.block_diag(list(sysm["Ml"]), format="csr").toarray()
    npr = M.shape[0]
    one = np.zeros(npr); one[0::10] = 1.0
    ang = ns_mesh.angles(S.verts, S.tris)
    cons = dict(full=[], wired=[], wirefn=[], zero3=[], wiredslit=[], zero3slit=[])
    # slit wiring: at both endpoints of every short edge (length < 0.3 x the longest edge of its two triangles)
    # wire the two triangles sharing it
    edges = {}
    for t, tri in enumerate(S.tris):
        for a in range(3):
            e = tuple(sorted((int(tri[a]), int(tri[(a + 1) % 3]))))
            edges.setdefault(e, []).append(t)
    slit = []
    for e, T in edges.items():
        if len(T) == 2:
            L = np.linalg.norm(S.verts[e[0]] - S.verts[e[1]])
            Lmax = max(np.linalg.norm(S.verts[S.tris[t][a]] - S.verts[S.tris[t][(a + 1) % 3]]) for t in T for a in range(3))
            opp = [ang[t, [a for a in range(3) if int(S.tris[t][a]) not in e][0]] for t in T]
            if L < 0.3 * Lmax and max(opp) < np.radians(10):
                for v in e:
                    slit.append(vertex_rep(S, sysm, T[0], v) - vertex_rep(S, sysm, T[1], v))
    info = []; nslit = 0
    for z in lk:
        T = [t for t in range(len(S.tris)) if z in S.tris[t]]
        az = [ang[t, int(np.where(S.tris[t] == z)[0][0])] for t in T]
        if len(T) == 2:                          # plain lock (mesh 'one'/'alt')
            tf, tl = T; tn = None
        else:
            assert len(T) == 3, (z, T)
            k = int(np.argmin(az)); tn = T[k]; tf, tl = [T[j] for j in range(3) if j != k]
        rf, rl = vertex_rep(S, sysm, tf, z), vertex_rep(S, sysm, tl, z)
        cons["wired"].append(rf - rl)
        if tn is not None:
            rn = vertex_rep(S, sysm, tn, z)
            cons["wirefn"].append(rf - rn); cons["zero3"] += [rf, rn, rl]
        else:
            cons["wirefn"].append(rf - rl); cons["zero3"] += [rf, rl]
        info.append(float(np.degrees(min(az))))
    cons["wiredslit"] = cons["wired"] + slit
    cons["zero3slit"] = cons["zero3"] + slit
    nslit = len(slit)
    lu = spl.splu(Af, permc_spec="COLAMD")
    X = lu.solve(Bf.T.toarray())
    Sd = Bf @ X; Sd = 0.5 * (Sd + Sd.T)
    out = dict(mesh=mesh, N=N, eps=eps, hG=S.hGamma, nmod=len(lk), nslitcons=nslit, minang_z=min(info) if info else None,
               minang=float(np.degrees(ang.min())))
    for name, Y in cons.items():
        if mesh in ("std",) and name != "full":
            continue
        Y = np.array([one] + Y).T
        Q, _ = np.linalg.qr(M @ Y, mode="complete")
        Z = Q[:, Y.shape[1]:]
        MZ = Z.T @ M @ Z
        ev = sla.eigh(Z.T @ Sd @ Z, MZ, eigvals_only=True, subset_by_index=[0, 2])
        out["beta_" + name] = float(np.sqrt(max(ev[0], 0)))
        # primal: max_q (Mq)^T S^+ (Mq) / (q^T M q) over q = Z y ; S^+ on one^perp via the full mean-free basis
        if name in ("full", "zero3", "wired", "zero3slit", "wiredslit"):
            Qf, _ = np.linalg.qr(one[:, None], mode="complete"); U = Qf[:, 1:]       # Euclidean complement of 'one'
            SU = U.T @ Sd @ U
            G = M @ Z                                                               # columns M q
            W = U @ np.linalg.solve(SU, U.T @ G)                                    # S^+ M q
            K = G.T @ W; K = 0.5 * (K + K.T)
            evp = sla.eigh(K, MZ, eigvals_only=True, subset_by_index=[K.shape[0] - 1, K.shape[0] - 1])
            out["gamma_" + name] = float(1 / np.sqrt(evp[-1]))
    out["secs"] = time.time() - t0
    return out


if __name__ == "__main__":
    mesh, N = sys.argv[1], int(sys.argv[2])
    for e in sys.argv[3:]:
        r = run(mesh, N, float(e))
        print(json.dumps({k: (float(f"{v:.5g}") if isinstance(v, float) else v) for k, v in r.items()}), flush=True)
