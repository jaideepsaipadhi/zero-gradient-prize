"""
notes/v7/nsconst: the cell prediction for lim d_h / hG^{1/2} on the production (LP) meshes (REPORT.tex, Section 4).

A(gamma; a) := cell value of notes/v6/unifmu2/floquet_cell.py (one-edge sector of the regular N-gon annulus, first-layer
aspect a, penalty gamma in edge units, mode m), i.e. ||grad W|| / (G l^{1/2}).  Two-scale conjecture (REPORT, Conj. 4.3):

    d_h^2 / hG  ->  D^2 := ||grad U||^{-2} int_Gamma A(gamma l(x); a(x))^2 l(x) (p - pbar)^2 ds ,

l(x) = |e(x)|/hG (production: l = max(cos^2, sin^2)), a(x) the local first-layer aspect, gamma = mu hG / h.
Test P:  p - pbar = -3/4 sin t + 1/4 sin 3t + 1/2 cos t on Gamma, ||grad U|| = 4.57024450 (paper, Section hc).

Usage: python leak_cell.py table N          A(gamma; a) on an (a, gamma) grid
       python leak_cell.py predict N        D for mu = 1e8 and mu = 100 (uses lim h/hG of the production family)
"""
import os, sys, json, time
for v_ in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v_, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "..")
sys.path.insert(0, os.path.join(ROOT, "code"))
sys.path.insert(0, os.path.join(ROOT, "notes", "v6", "unifmu2"))
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
import floquet_cell as fc
from cell_limit import a_of_theta

GRADU = 4.57024450


def A_values(N, a, gams, m=2):
    """copy of the solve loop of floquet_cell.cmd_phi, returning A for each gamma."""
    cell = fc.Cell(N, a); S = cell.S
    T = cell.T(m * cell.al); TH = T.conj().T.tocsr(); BT = (cell.B @ T).tocsr()
    u, gu = fc.leak_flow(m)
    Wd = cell.wall
    th = np.arctan2(Wd["X"][:, 1], Wd["X"][:, 0])
    g = -np.exp(1j * m * th)[:, None] * Wd["n"][None, :]
    F = np.zeros(2 * S.nn, complex); d = fc.svn.dofs(S.ids[Wd["t"]])
    for c in range(2):
        np.add.at(F, d[:, c], (Wd["w"][:, None] * Wd["ph"] * g[:, c:c + 1]).sum(0))
    Fr = TH @ F; nB = BT.shape[0]; G = np.sqrt(2 * np.pi)
    out = []
    for gam in gams:
        lam = gam / cell.ell
        Kr = (TH @ cell.K(lam) @ T).tocsc()
        M = sp.bmat([[Kr, -BT.conj().T], [BT, None]], format="csc")
        rhs = np.concatenate([lam * Fr, np.zeros(nB)])
        lu = spl.splu(M, permc_spec="COLAMD", diag_pivot_thresh=1.0)
        x = lu.solve(rhs)
        for _ in range(2):
            x = x + lu.solve(rhs - M @ x)
        X = T @ x[:Kr.shape[0]]
        out.append(float(np.sqrt(cell.N) * fc.grad_err(cell, X, gu) / (G * np.sqrt(cell.ell))))
    return out


AGRID = [0.75, 0.8, 0.9, 1.0, 1.15, 1.35, 1.6, 1.9, 2.25, 2.6]


def cmd_table(N, gams=(20, 30, 100, 1e3, 1e10)):
    rows = {}
    for a in AGRID:
        t0 = time.time()
        rows[a] = A_values(N, a, list(gams))
        print(json.dumps(dict(N=N, a=a, gammas=list(gams), A=[round(v, 5) for v in rows[a]],
                              secs=round(time.time() - t0, 1))), flush=True)
    return rows


def ratio_h_hG(N):
    import svn
    S = svn.Space(*svn.make_mesh(N))
    return S.h / S.hGamma


def cmd_predict(N):
    rh = {n: ratio_h_hG(n) for n in (64, 128, 256, 512)}
    print(json.dumps(dict(h_over_hG={k: round(v, 5) for k, v in rh.items()})), flush=True)
    rinf = rh[512]
    th = np.linspace(0, 2 * np.pi, 2001)[:-1]
    thr = np.abs(((th + np.pi / 4) % (np.pi / 2)) - np.pi / 4)        # reduced angle in [0, pi/4]
    ell = np.cos(thr) ** 2
    a = a_of_theta(thr)
    pp = -0.75 * np.sin(th) + 0.25 * np.sin(3 * th) + 0.5 * np.cos(th)
    for mu in (1e8, 100.0):
        gam = mu / rinf
        # A on the needed (a, gamma*ell) pairs: tabulate on AGRID at gamma*ell values, interpolate
        gl = sorted(set(np.round(gam * np.cos(np.linspace(0, np.pi / 4, 6)) ** 2, 3))) if mu < 1e6 else [1e10]
        tab = np.array([A_values(N, aa, gl) for aa in AGRID])           # (na, ng)
        Aint = np.empty_like(th)
        for i in range(len(th)):
            col = np.array([np.interp(a[i], AGRID, tab[:, j]) for j in range(len(gl))])
            Aint[i] = col[0] if len(gl) == 1 else np.interp(gam * ell[i], gl, col)
        D2 = np.sum(Aint ** 2 * ell * pp ** 2) * (2 * np.pi / len(th)) / GRADU ** 2
        print(json.dumps(dict(N_cell=N, mu=mu, gamma=round(gam, 3), gammas_tab=gl, D=round(float(np.sqrt(D2)), 4))),
              flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "table":
        cmd_table(int(sys.argv[2]))
    elif sys.argv[1] == "predict":
        cmd_predict(int(sys.argv[2]))
