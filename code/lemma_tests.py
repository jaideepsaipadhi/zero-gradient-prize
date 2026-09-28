"""Numerical stress tests of individual lemmas via exact discrete dual norms.
  dual_Z(l)  = sup_{v in Z_h}  l(v)/||grad v||     (saddle solve with pointwise div-free constraint)
  dual_R(l)  = sup_{v in V_h^R} l(v)/||grad v||
Tests
  L4.1  normal load   S = q(x) x/|x|      : <S, d_n v>   -> should stay BOUNDED     (Lemma 4.1)
        tangential    S = q(x) x^perp/|x| : <S, d_n v>   -> should grow ~ h^{-1/2}  (remark after Thm C)
        normal load on V_h^R (no div-free) :             -> should grow (incompressibility is what saves it)
  L1.5  <(grad u)^T n_e, v> for Scott's shear flow and test B: dual_R ~ h^{3/2} (corrected Lemma 1.5)
"""
import os, sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, warnings
warnings.filterwarnings("ignore")
import svn

q = lambda X, Y: X**2 * Y - Y + X / 2


def gram_grad(S):
    nt = len(S.tris)
    gx, gy = S.grads(svn.QX, svn.QY)
    wq = svn.QW[None, :] * np.abs(S.detJ)[:, None]
    K = np.einsum('tq,tqa,tqb->tab', wq, gx, gx) + np.einsum('tq,tqa,tqb->tab', wq, gy, gy)
    loc = np.zeros((nt, 15, 2, 15, 2)); loc[:, :, 0, :, 0] = K; loc[:, :, 1, :, 1] = K
    D = svn.dofs(S.ids).reshape(nt, 30)
    rows = np.repeat(D, 30, axis=1).ravel(); cols = np.tile(D, (1, 30)).ravel()
    return sp.coo_matrix((loc.reshape(nt, 30, 30).ravel(), (rows, cols)), shape=(2 * S.nn, 2 * S.nn)).tocsr()


def edge_functional(S, kind, fn):
    """kind 'dn': l(v)=sum_e int S . d_n v ;  kind 'val': l(v)=sum_e int W . v"""
    l = np.zeros(2 * S.nn)
    for (t, ref, Xe, we, n) in S.edge_data():
        ph = svn.basis(ref[:, 0], ref[:, 1])
        dX, dY = svn.dbasis(ref[:, 0], ref[:, 1]); Ji = S.Jinv[t]
        gxe = Ji[0, 0] * dX + Ji[1, 0] * dY; gye = Ji[0, 1] * dX + Ji[1, 1] * dY
        dn = n[0] * gxe + n[1] * gye
        W = fn(Xe[:, 0], Xe[:, 1], n)                 # (ng,2)
        B = dn if kind == "dn" else ph
        d = svn.dofs(S.ids[t])
        for c in range(2):
            np.add.at(l, d[:, c], (we[:, None] * W[:, c:c + 1] * B).sum(0))
    return l


def duals(S, sysm, ls):
    H = gram_grad(S)
    dD = svn.dofs(S.dnodes).ravel(); free = np.setdiff1d(np.arange(2 * S.nn), dD)
    Hf = H[free][:, free].tocsc(); Bf = sysm["B"][:, free]
    luR = spl.splu(Hf)
    K = sp.bmat([[Hf, Bf.T], [Bf, None]], format="csc")
    luZ = spl.splu(K, permc_spec="COLAMD", diag_pivot_thresh=1.0)
    out = []
    for l in ls:
        lf = l[free]
        xR = luR.solve(lf)
        xZ = luZ.solve(np.concatenate([lf, np.zeros(Bf.shape[0])]))[:len(free)]
        out.append((np.sqrt(abs(lf @ xR)), np.sqrt(abs(lf @ xZ))))
    return out


if __name__ == "__main__":
    exA = svn.make_exact("A"); exB = svn.make_exact("B1")
    Snorm = lambda X, Y, n: np.stack([q(X, Y) * X, q(X, Y) * Y], -1) / np.hypot(X, Y)[:, None]
    Stang = lambda X, Y, n: np.stack([-q(X, Y) * Y, q(X, Y) * X], -1) / np.hypot(X, Y)[:, None]

    def transposed(ex):
        def f(X, Y, n):
            G = ex["gu"](X, Y)                                  # G[...,c,d] = d_d u_c
            return np.einsum('gcd,c->gd', G, n)                 # ((grad u)^T n)_d = sum_c d_d u_c n_c
        return f

    rows = []
    for N in [int(a) for a in (sys.argv[1] if len(sys.argv) > 1 else "16,24,32,48,64").split(",")]:
        v, t, c, o = svn.make_mesh(N); S = svn.Space(v, t, c, o)
        sysm = svn.assemble(S, 100.0, exA["f"], None)
        ls = [edge_functional(S, "dn", Snorm), edge_functional(S, "dn", Stang),
              edge_functional(S, "val", transposed(exA)), edge_functional(S, "val", transposed(exB))]
        d = duals(S, sysm, ls)
        rows.append((N, S.hGamma, d))
        print(f"N={N:3d} hG={S.hGamma:.4f} | normal S: Z {d[0][1]:.4e}  R {d[0][0]:.4e} | tangential S: Z {d[1][1]:.4e}"
              f" | (grad u)^T n  A: R {d[2][0]:.4e}  B: R {d[3][0]:.4e}", flush=True)
    print("\nobserved rates vs h_Gamma (positive = decays, negative = grows):")
    names = ["normal S on Z_h", "normal S on V_h^R", "tangential S on Z_h", "(grad u)^T n, test A", "(grad u)^T n, test B"]
    pick = [lambda d: d[0][1], lambda d: d[0][0], lambda d: d[1][1], lambda d: d[2][0], lambda d: d[3][0]]
    for nm, f in zip(names, pick):
        r = [np.log(f(rows[i][2]) / f(rows[i - 1][2])) / np.log(rows[i][1] / rows[i - 1][1]) for i in range(1, len(rows))]
        print(f"  {nm:24s}: " + "  ".join(f"{x:5.2f}" for x in r))
