"""Slit lab numerics (interior needle pair, P4/P3disc, v=0 on [-1,1]^2).
 (a) beta and the 4 smallest inf-sup eigenvalues (interior pressures) / eps
 (b) gamma_off: right-inverse constant on {q = 0 on the two needles}   [Theorem S analogue]
 (c) dual ratio of the single-vertex functional  q -> q_N(x_up) - q_N'(x_up)   [Prop. sqrt(eps)]
 (d) best approximation in Z_h of a smooth div-free u = curl psi, vs eps and vs the unsplit mesh
Usage: python3 probe2.py n eps1 eps2 ..."""
import sys, json, time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, scipy.linalg as sla
from slit_lab import *

def psi_derivs(X, Y):
    # psi = (1-x^2)^2 (1-y^2)^2 * cos(1.3 x + 0.7 y);  u = (psi_y, -psi_x)
    import sympy as sy
    return None

import sympy as sy
xs, ys = sy.symbols("x y")
PSI = (1 - xs**2)**2 * (1 - ys**2)**2 * sy.cos(1.3 * xs + 0.7 * ys + 0.4)
U = [sy.diff(PSI, ys), -sy.diff(PSI, xs)]
GU = [[sy.lambdify((xs, ys), sy.diff(U[i], v), "numpy") for v in (xs, ys)] for i in range(2)]
UF = [sy.lambdify((xs, ys), U[i], "numpy") for i in range(2)]

def best_approx(S, sysm, parts):
    lu, Af, Bf, free = parts
    gx, gy = S.grads(svn.QX, svn.QY)
    X = S.phys(svn.QX, svn.QY)
    wq = svn.QW[None, :] * np.abs(S.detJ)[:, None]
    nt = len(S.tris); F = np.zeros(2 * S.nn)
    D = svn.dofs(S.ids)
    E0 = 0.0
    for c in range(2):
        dux = GU[c][0](X[..., 0], X[..., 1]) * np.ones_like(X[..., 0]); duy = GU[c][1](X[..., 0], X[..., 1]) * np.ones_like(X[..., 0])
        Fl = np.einsum('tq,tqa,tq->ta', wq, gx, dux) + np.einsum('tq,tqa,tq->ta', wq, gy, duy)
        np.add.at(F, D[:, :, c].ravel(), Fl.ravel())
        E0 += np.sum(wq * (dux**2 + duy**2))
    Ff = F[free]
    nq = Bf.shape[0]
    M = sp.block_diag(list(sysm["Ml"]), format="csr")
    delta = 1e-11
    A = sp.bmat([[Af, Bf.T], [Bf, -delta * M]], format="csc")
    sol = spl.splu(A).solve(np.concatenate([Ff, np.zeros(nq)]))
    z = sol[:len(Ff)]
    Az = Af @ z
    err2 = E0 - 2 * Ff @ z + z @ Az
    divres = np.linalg.norm(Bf @ z) / max(np.sqrt(z @ Az), 1e-300)
    # also plain best approximation in V_h (no div constraint)
    zv = lu.solve(Ff); errv2 = E0 - 2 * Ff @ zv + zv @ (Af @ zv)
    return float(np.sqrt(max(err2, 0))), float(np.sqrt(max(errv2, 0))), float(divres)

def run(n, eps):
    t0 = time.time()
    S, g = build(n, eps) if eps > 0 else (None, None)
    if eps == 0:          # unsplit reference mesh
        verts, tris, vid = square_mesh(n)
        S = svn.Space(verts, orient(verts, tris), {int(vid[0, 0]), int(vid[1, 0])}, None)
    Sd, M, one, sysm, parts = schur(S)
    E = interior_basis(S)
    out = dict(n=n, eps=eps)
    b, V = beta_sub(Sd, M, E, k=4)
    out["beta"] = [float(f"{x:.4g}") for x in b]
    if eps > 0:
        out["beta_over_eps"] = [float(f"{x/eps:.4g}") for x in b]
        tN = tri_of(S, g["z"], g["dn"], g["up"]); tNp = tri_of(S, g["f"], g["up"], g["dn"])
        # (b) right inverse constant on interior pressures vanishing on both needles
        keep = [j for j in range(E.shape[1]) if not (E[:, j].nonzero()[0][0] // 10 in (tN, tNp))]
        Eo = E[:, keep]
        MZ = Eo.T @ M @ Eo
        G = M @ Eo
        W = np.linalg.lstsq(Sd, G, rcond=1e-13)[0]
        K = G.T @ W; K = .5 * (K + K.T)
        evp = sla.eigh(K, MZ, eigvals_only=True)
        out["gamma_off"] = float(f"{1/np.sqrt(evp[-1]):.4g}")
        b_off, _ = beta_sub(Sd, M, Eo, k=1); out["beta_off"] = float(f"{b_off[0]:.4g}")
        # (c) single-vertex functional at x_up
        r = point_rep(S, sysm, tN, S.verts[g["up"]]) - point_rep(S, sysm, tNp, S.verts[g["up"]])
        out["vertex_fn_ratio"] = float(f"{dual_norm_ratio(Sd, M, r):.4g}")
        out["vertex_fn_ratio_over_sqrt_eps"] = float(f"{dual_norm_ratio(Sd, M, r)/np.sqrt(eps):.4g}")
        # mode mass on needles
        q = V[:, 0]
        mass = np.array([q[10*t:10*t+10] @ sysm["Ml"][t] @ q[10*t:10*t+10] for t in range(len(S.tris))]); mass /= mass.sum()
        out["mode1_needle_mass"] = float(f"{mass[tN] + mass[tNp]:.5f}")
        q = V[:, 1]
        mass = np.array([q[10*t:10*t+10] @ sysm["Ml"][t] @ q[10*t:10*t+10] for t in range(len(S.tris))]); mass /= mass.sum()
        out["mode2_needle_mass"] = float(f"{mass[tN] + mass[tNp]:.5f}")
    e_z, e_v, dr = best_approx(S, sysm, parts)
    out.update(best_Zh=float(f"{e_z:.6g}"), best_Vh=float(f"{e_v:.6g}"), divres=float(f"{dr:.2g}"), secs=round(time.time() - t0, 1))
    print(json.dumps(out), flush=True)

if __name__ == "__main__":
    n = int(sys.argv[1])
    for e in sys.argv[2:]:
        run(n, float(e))
