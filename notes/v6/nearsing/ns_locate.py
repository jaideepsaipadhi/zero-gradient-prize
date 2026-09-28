"""Localise the smallest inf-sup mode (full mean-free pressure space) on a nearly singular mesh.
Usage: python3 ns_locate.py MESH N eps      Prints the triangles carrying most of the M-mass of the mode."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, scipy.linalg as sla, scipy.sparse as sp, scipy.sparse.linalg as spl
import ns_infsup as I
from ns_infsup import svn, strong_bc, ns_mesh, laplace

def mode(mesh, N, eps, extra=None):
    if mesh == "interior":
        verts, tris, circ, outer = svn.make_mesh(N); lk = []
        verts, tris = interior_split(verts, tris, N, eps)
        verts, tris = ns_mesh._orient(verts, np.array(tris))
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
    one = np.zeros(M.shape[0]); one[0::10] = 1
    X = spl.splu(Af).solve(Bf.T.toarray()); Sd = Bf @ X; Sd = .5 * (Sd + Sd.T)
    Q, _ = np.linalg.qr((M @ one)[:, None], mode="complete"); Z = Q[:, 1:]
    ev, V = sla.eigh(Z.T @ Sd @ Z, Z.T @ M @ Z, subset_by_index=[0, 1])
    q = Z @ V[:, 0]
    mass = np.array([q[10*t:10*t+10] @ sysm["Ml"][t] @ q[10*t:10*t+10] for t in range(len(S.tris))])
    mass /= mass.sum()
    ang = np.degrees(ns_mesh.angles(S.verts, S.tris))
    order = np.argsort(-mass)[:8]
    print(json.dumps(dict(mesh=mesh, N=N, eps=eps, beta=float(np.sqrt(ev[0])), beta2=float(np.sqrt(ev[1])))))
    for t in order:
        print(f"  tri {t:5d} mass {mass[t]:.4f}  angles {np.round(ang[t],2)}  verts {S.tris[t]}  wall? {[int(v) in circ for v in S.tris[t]]}")
    return S, q

def interior_split(verts, tris, N, eps, layer=None):
    """split an interior vertex d (layer N//8) along its radial edge d->f, with apex-side vertex c (layer-1):
    exactly the slit of ns_mesh.split_apex but away from the wall"""
    m = N // 4; vid = np.arange((m + 1) * N).reshape(m + 1, N)
    j = layer or max(2, m // 2); i = 1
    z = vid[j - 1, i]; w = vid[j, i]           # radial edge z->w, w split (slit needles at z and at vid[j+1,i])
    return ns_mesh.split_apex(verts, np.array(tris), z, eps, w=w)

if __name__ == "__main__":
    mode(sys.argv[1], int(sys.argv[2]), float(sys.argv[3]))
