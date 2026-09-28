"""Pressure at locked vertices: raw vs filtered.  Solves method (S) (strong no-slip, strong_bc.solve_strong)
on mesh one/alt, then removes the Q_lock^perp component of p_h:
    p_f = p_h - P_perp p_h,   P_perp = M-orthogonal projection onto span{rho_z},
    rho_z = R_{T1}(z) - R_{T2}(z)  (Riesz representers of the vertex values; local on the locked star).
Theorem P(ii) of REPORT.tex: ||p_f - Pi(p~ - pbar)|| <= C (||grad(u~-u_h)|| + h_G^2) (velocity rate),
while the raw p_h error is only bounded by C (...) / phi_min, phi_min ~ h_G.
Usage: python3 filter.py MESH N [N ...]   (MESH in one, alt; lu solver unless pypardiso importable)
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code")); sys.path.insert(0, HERE)
import numpy as np, scipy.sparse as sp
import svn, strong_bc, pressure_err
from infsup import MESHES, constraint_vectors


def run(mesh, N, solver):
    t0 = time.time()
    ex = svn.make_exact("A")
    verts, tris, circ, outer = MESHES[mesh](N)
    S = svn.Space(verts, tris, circ, outer)
    sysm = svn.assemble(S, 100.0, ex["f"], None)
    u, p, info = strong_bc.solve_strong(S, sysm, ex["u"], "zero", solver)
    E = svn.errors(S, u, ex, 100.0)
    one, reps, lk = constraint_vectors(S, sysm, circ)
    M = sp.block_diag(list(sysm["Ml"]), format="csr")
    Y = np.array([r1 - r2 for r1, r2 in reps]).T
    G = Y.T @ (M @ Y)
    coef = np.linalg.solve(G, Y.T @ (M @ p))
    pf = p - Y @ coef
    raw = pressure_err.pressure_errors(S, p, ex["p"]); fil = pressure_err.pressure_errors(S, pf, ex["p"])
    return dict(mesh=mesh, N=N, hG=S.hGamma, nlocked=len(lk), H1=E["H1"], pL2=raw["pL2opt"], pL2_f=fil["pL2opt"],
                pLinf_layer=raw["pLinf_layer"], pLinf_layer_f=fil["pLinf_layer"], pL2_off=raw["pL2_off"],
                pL2_off_f=fil["pL2_off"], perp_norm=float(np.sqrt(abs((Y @ coef) @ (M @ (Y @ coef))))),
                divres=info["divres"], secs=time.time() - t0)


if __name__ == "__main__":
    mesh = sys.argv[1]
    try:
        import pypardiso; solver = "pardiso"
    except Exception:
        solver = "lu"
    for N in [int(a) for a in sys.argv[2:]]:
        print(json.dumps(run(mesh, N, solver)), flush=True)
