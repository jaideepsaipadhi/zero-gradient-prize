"""Discrete inf-sup constant of P_k / P_{k-1}^disc on the Alfeld mesh (IS3):
   beta_h^2 = min eig of  (B A^{-1} B^T, M_p)  on L^2_0,  V_0 = velocities vanishing on the whole boundary,
   A = vector Laplacian (grad u, grad v).  Also counts spurious pressure modes (eigenvalues < 1e-10).
  python3 infsup3d.py box K N [N ...]   |   python3 infsup3d.py sph K N [N ...]"""
import os, sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, scipy.linalg as sl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svn3d

zero = lambda a, b, c: np.zeros(np.shape(a) + (3,))


def infsup(S):
    R = S.R; nl = R.nl; nt = len(S.T); n3 = 3 * nl
    G = np.einsum('tim,tjn,ijab->tmnab', S.Jinv, S.Jinv, R.K) * S.det[:, None, None, None, None]
    tr = G[:, 0, 0] + G[:, 1, 1] + G[:, 2, 2]
    loc = np.zeros((nt, nl, 3, nl, 3))
    for c in range(3):
        loc[:, :, c, :, c] = tr
    D = S.dofs.reshape(nt, n3)
    A = sp.coo_matrix((loc.ravel(), (np.repeat(D, n3, 1).ravel(), np.tile(D, (1, n3)).ravel())), (3 * S.nn,) * 2).tocsc()
    s = svn3d.assemble(S, 0.0, zero, [], [])
    B = s["B"]
    # boundary nodes: nodes on any boundary face
    bn = set()
    FL = [(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)]
    lam = np.array(R.LBARY)
    bt = np.concatenate([S.fBox[0], S.fGamma[0]]); bi = np.concatenate([S.fBox[1], S.fGamma[1]])
    for i in range(4):
        loc_nodes = np.where(lam[:, i] == 0)[0]
        tt = bt[bi == i]
        bn.update(S.ids[tt][:, loc_nodes].ravel().tolist())
    bn = np.array(sorted(bn))
    dD = (3 * bn[:, None] + np.arange(3)).ravel(); free = np.setdiff1d(np.arange(3 * S.nn), dD)
    Af = A[free][:, free]; Bf = B[:, free].tocsc()
    lu = spl.splu(Af)
    X = lu.solve(Bf.T.toarray())
    Sm = Bf @ X; Sm = (Sm + Sm.T) / 2
    Mp = sp.block_diag([R.Mp * d for d in S.det]).toarray()
    ev = sl.eigh(Sm, Mp, eigvals_only=True)
    return ev, len(free), B.shape[0]


if __name__ == "__main__":
    kind, k = sys.argv[1], int(sys.argv[2])
    for N in map(int, sys.argv[3:]):
        S = svn3d.build(kind, N, k, **({"L": 1.0} if kind == "box" else {"L": 2.0}))
        ev, nv, npr = infsup(S)
        print(f"{kind} k={k} N={N} vel_dofs(V0)={nv} p_dofs={npr} smallest eigs: {np.array2string(ev[:4], precision=4)}  "
              f"#(<1e-10)={int((ev < 1e-10).sum())}  beta_h(on L2_0)={np.sqrt(ev[1]):.4f}  max={ev[-1]:.4f}", flush=True)
