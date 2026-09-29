"""
notes/v7/nsconst: the ACTUAL discrete stability constant of (L) on the production meshes (REPORT.tex, Section 2.5).

  beta_L(h) := inf_{delta in Z_h} sup_{v in Z_h} A(delta, v) / (|||delta||| |||v|||),   A = nu N_h + K,  K(w,v) = c1(w;u,v) + c1(u;w,v),

GS / Scott--Vogelius P4-P3disc + Nitsche (code/svn_k.py, code/svn_gen.py, unchanged), test B velocity u (nodal interpolant),
mu = 100 (the paper's production penalty), |||.||| the paper's energy norm (02_setting.tex l.73).
Also the same quantity for the Stokes operator (K = 0), whose ratio to nu is the discrete Nitsche coercivity scale.

python discrete_L.py N nu
"""
import os, sys, json, time
for v_ in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v_, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
import svn_k, svn_gen as sg
from consts import exact_B


def energy_gram(S, mu, sysm):
    R = S.R; nt = len(S.tris); nl = R.nl
    gx, gy = S.grads(R.QX, R.QY)
    wq = R.QW[None, :] * np.abs(S.detJ)[:, None]
    loc = np.einsum('tq,tqa,tqb->tab', wq, gx, gx) + np.einsum('tq,tqa,tqb->tab', wq, gy, gy)
    rows = np.repeat(S.ids, nl, axis=1).ravel(); cols = np.tile(S.ids, (1, nl)).ravel()
    Ks = sp.coo_matrix((loc.ravel(), (rows, cols)), shape=(S.nn, S.nn)).tocsr()
    rE, cE, vE = [], [], []
    for (t, ref, Xe, we, n) in S.edge_data():
        dX, dY = R.dbasis(ref[:, 0], ref[:, 1]); Ji = S.Jinv[t]
        dn = n[0] * (Ji[0, 0] * dX + Ji[1, 0] * dY) + n[1] * (Ji[0, 1] * dX + Ji[1, 1] * dY)
        le = we.sum()
        blk = le * np.einsum('g,ga,gb->ab', we, dn, dn)
        rE.append(np.repeat(S.ids[t], nl)); cE.append(np.tile(S.ids[t], nl)); vE.append(blk.ravel())
    Ke = sp.coo_matrix((np.concatenate(vE), (np.concatenate(rE), np.concatenate(cE))), shape=(S.nn, S.nn)).tocsr()
    Tsc = (Ks + Ke).tocoo()
    Tv = sp.coo_matrix((np.concatenate([Tsc.data, Tsc.data]),
                        (np.concatenate([2 * Tsc.row, 2 * Tsc.row + 1]), np.concatenate([2 * Tsc.col, 2 * Tsc.col + 1]))),
                       shape=(2 * S.nn, 2 * S.nn)).tocsr()
    return (Tv + mu * sysm["Pen"]).tocsr()          # Pen = boundary mass / h  ->  (mu/h) ||v||^2_Gamma


def infsup(Kmat, Tmat, B, free):
    Kf = Kmat[free][:, free].tocsc(); Bf = B[:, free].tocsc(); Tf = Tmat[free][:, free].tocsr()
    nf, npr = Kf.shape[0], B.shape[0]
    M = sp.bmat([[Kf, -Bf.T], [Bf, None]], format="csc")
    lu = spl.splu(M, permc_spec="COLAMD")

    def sol(r, trans):
        rhs = np.concatenate([r, np.zeros(npr)])
        x = lu.solve(rhs, trans="T" if trans else "N")
        x = x + lu.solve(rhs - (M.T @ x if trans else M @ x), trans="T" if trans else "N")
        return x[:nf]
    Op = spl.LinearOperator((nf, nf), matvec=lambda y: sol(Tf @ sol(Tf @ y, False), True), dtype=float)
    w = spl.eigs(Op, k=1, which="LM", return_eigenvectors=False, tol=1e-9, maxiter=5000)
    return float(1 / np.sqrt(np.abs(w).max()))


def main(N, nu, mu=100.0):
    t0 = time.time()
    S = sg.circle_space(N, 4)
    sysm = svn_k.assemble(S, mu, lambda X, Y: np.zeros(np.shape(X) + (2,)))
    U, G = exact_B()
    u = U(S.xy[:, 0], S.xy[:, 1]).ravel()
    _, Jc = sg.conv_system(S, u, "c1")
    T = energy_gram(S, mu, sysm)
    free = np.setdiff1d(np.arange(2 * S.nn), svn_k.dofs(S.dnodes).ravel())
    bS = infsup(nu * sysm["A"], T, sysm["B"], free)
    bL = infsup(nu * sysm["A"] + Jc, T, sysm["B"], free)
    print(json.dumps(dict(N=N, nu=nu, mu=mu, h=S.h, hG=S.hGamma, gamma=mu * S.hGamma / S.h,
                          beta_Stokes_over_nu=bS / nu, beta_L=bL, beta_L_over_nu=bL / nu,
                          ndof=int(len(free)), secs=round(time.time() - t0, 1))), flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), float(sys.argv[2]))
