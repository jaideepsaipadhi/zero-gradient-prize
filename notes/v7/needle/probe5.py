"""Theorem S on the actual wall meshes (notes/v6/nearsing/ns_mesh.py): inf-sup of the SLIT-DOMAIN pair
  V = {v in V_O : v = 0 on every needle (all P4 nodes of needle triangles)},
  Q = {q in Pi_h cap L^2_0, q = 0 on needles, + corner/lock constraints of Omega^K}:
    family II ('split'): q|T1(z) = q|T3(z) = 0 (T1, T3 are corner triangles of Omega^K at z)
    family I  ('wall') : q|T1(z) = q|T2(z) (lock of Omega^K at z)  and  q|T2(y) = 0 (T2 is the only triangle at y)
Compared with the lock-free constant of mesh 'std' (0.10229 at N=16).  Usage: python3 probe5.py MESH N eps..."""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "v6", "nearsing"))
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, scipy.linalg as sla
import ns_infsup as I
from ns_infsup import svn, strong_bc, ns_mesh, laplace, vertex_rep

def run(mesh, N, eps):
    t0 = time.time()
    verts, tris, circ, outer, lk = ns_mesh.make(mesh, N, eps)
    S = svn.Space(verts, tris, circ, outer)
    sysm = svn.assemble(S, 1.0, lambda X, Y: np.zeros(X.shape + (2,)), None)
    A = laplace(S); B = sysm["B"]
    ang = ns_mesh.angles(S.verts, S.tris)
    needles = [t for t in range(len(S.tris)) if ang[t].min() < np.radians(min(5.0, np.degrees(eps) * 3 + 1e-9) if eps < 0.05 else 25)]
    gn = strong_bc.gamma_nodes(S)
    fixed = np.union1d(np.union1d(svn.dofs(S.dnodes).ravel(), svn.dofs(gn).ravel()), svn.dofs(S.ids[needles].ravel()).ravel())
    free = np.setdiff1d(np.arange(A.shape[0]), fixed)
    Af = A[free][:, free].tocsc(); Bf = B[:, free].tocsr()
    X = spl.splu(Af).solve(Bf.T.toarray()); Sd = Bf @ X; Sd = .5 * (Sd + Sd.T)
    M = sp.block_diag(list(sysm["Ml"]), format="csr").toarray()
    npr = M.shape[0]; one = np.zeros(npr); one[0::10] = 1
    cons = []
    for t in needles:
        for j in range(10):
            e = np.zeros(npr); e[10 * t + j] = 1; cons.append(np.linalg.solve(M, e) if False else e)   # coefficient constraints
    # corner / lock constraints (as M-representers)
    for z in lk:
        T = [t for t in range(len(S.tris)) if z in S.tris[t] and t not in needles]
        if mesh.endswith("split"):
            for t in T:
                cons.append(M @ vertex_rep(S, sysm, t, z))
        else:
            assert len(T) == 2, T
            cons.append(M @ (vertex_rep(S, sysm, T[0], z) - vertex_rep(S, sysm, T[1], z)))
    if mesh.endswith("wall"):     # y: interior vertex whose only non-needle triangle is T2
        for x in range(len(S.verts)):
            T = [t for t in range(len(S.tris)) if x in S.tris[t]]
            nn = [t for t in T if t not in needles]
            if len(T) == 3 and len(nn) == 1 and x not in circ:
                cons.append(M @ vertex_rep(S, sysm, nn[0], x))
    C = np.array([M @ one] + cons)        # rows: linear functionals on coefficient vectors
    _, s, Vt = np.linalg.svd(C, full_matrices=True)
    rank = int((s > 1e-10 * s.max()).sum()); Z = Vt[rank:].T
    ev = sla.eigh(Z.T @ Sd @ Z, Z.T @ M @ Z, eigvals_only=True, subset_by_index=[0, 1])
    print(json.dumps(dict(mesh=mesh, N=N, eps=eps, nneedles=len(needles), ncons=len(cons), beta_slit=float(f"{np.sqrt(max(ev[0],0)):.5g}"),
                          beta2=float(f"{np.sqrt(max(ev[1],0)):.5g}"), secs=round(time.time() - t0, 1))), flush=True)

if __name__ == "__main__":
    mesh, N = sys.argv[1], int(sys.argv[2])
    for e in sys.argv[3:]:
        run(mesh, N, float(e))
