"""Where does the smallest inf-sup mode live?  Dense, small N.  Usage: python3 locate_mode.py MESH N [space]"""
import sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, scipy.linalg as sla
from infsup import *
mesh, N = sys.argv[1], int(sys.argv[2]); space = sys.argv[3] if len(sys.argv) > 3 else "full"
S, sysm, circ = setup(mesh, N)
A = laplace(S); B = sysm["B"]
gn = strong_bc.gamma_nodes(S)
dD = np.union1d(svn.dofs(S.dnodes).ravel(), svn.dofs(gn).ravel()); free = np.setdiff1d(np.arange(A.shape[0]), dD)
lu = spl.splu(A[free][:, free].tocsc()); Bf = B[:, free]
Sd = Bf @ lu.solve(Bf.T.toarray()); Sd = .5 * (Sd + Sd.T); Md = sp.block_diag(list(sysm["Ml"])).toarray()
one, reps, lk = constraint_vectors(S, sysm, circ)
Y = np.array([one] + ([r1 - r2 for r1, r2 in reps] if space == "lockw" else [])).T
Q, _ = np.linalg.qr(Md @ Y, mode="complete"); Z = Q[:, Y.shape[1]:]
ev, V = sla.eigh(Z.T @ Sd @ Z, Z.T @ Md @ Z, subset_by_index=[0, 5])
print("ev", ev)
for j in range(2):
    q = Z @ V[:, j]; nt = len(S.tris)
    e = np.array([q[10*t:10*t+10] @ sysm["Ml"][t] @ q[10*t:10*t+10] for t in range(nt)]); e /= e.sum()
    idx = np.argsort(-e)[:6]
    cen = S.verts[S.tris[idx]].mean(1)
    print("mode", j, "top-6 triangle energy fractions", np.round(e[idx], 3), "centroids", np.round(cen, 2).tolist(),
          "radius", np.round(np.linalg.norm(cen, axis=1), 2).tolist())
