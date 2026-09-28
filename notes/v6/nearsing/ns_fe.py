"""Method (S) (strong no-slip SV, P4/P3disc, code/strong_bc.py) on the nearly singular meshes of ns_mesh.py.
Test case A (Scott's shear flow), mu = 100 (only enters the reported energy norm).
Usage: python3 ns_fe.py MESH EPSMODE N [N ...]
   MESH    in {onesplit, onewall, altsplit, altwall, one, alt, std}
   EPSMODE 'c0.3' (eps = 0.3), 'h1.0' (eps = 1.0*hG), 'q1.0' (eps = 1.0*hG^2)   [hG = 2 sin(pi/N) (chords)]
Prints one JSON line per N with H1 error, pressure errors, divres, the local star data (eps, phi, w_pred).
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import numpy as np
import svn, strong_bc, pressure_err
import ns_mesh


def eps_of(mode, N):
    hG = 2 * np.sin(np.pi / N)
    c = float(mode[1:])
    return {"c": c, "h": c * hG, "q": c * hG ** 2}[mode[0]]


def run(mesh, mode, N):
    t0 = time.time()
    eps = eps_of(mode, N) if mesh[3:] else 0.0
    if mesh == "std":
        verts, tris, circ, outer = svn.make_mesh(N); lk = []
    else:
        verts, tris, circ, outer, lk = ns_mesh.make(mesh, N, eps)
    S = svn.Space(verts, tris, circ, outer)
    ex = svn.make_exact("A")
    sysm = svn.assemble(S, 100.0, ex["f"], None)
    u, p, info = strong_bc.solve_strong(S, sysm, ex["u"], "zero", "lu")
    E = svn.errors(S, u, ex, 100.0)
    E.update(pressure_err.pressure_errors(S, p, ex["p"]))
    E.update(info)
    # |d_nu u(z)| at the modified vertices (exact gradient Frobenius norm)
    dn = [float(np.linalg.norm(ex["gu"](verts[z:z + 1, 0], verts[z:z + 1, 1])[0])) for z in lk]
    chk = ns_mesh.check(verts, tris, circ)
    rec = dict(mesh=mesh, mode=mode, N=N, eps=eps, hG=S.hGamma, nmod=len(lk), dnu_mean=float(np.mean(dn)) if dn else 0.0,
               H1=E["H1"], pL2=E["pL2"], divres=info["divres"], relres=info["relres"], vgrad=strong_bc.vertex_grad(S, u, ex),
               minang=chk["minang"], maxang=chk["maxang"], secs=time.time() - t0)
    return rec


if __name__ == "__main__":
    mesh, mode = sys.argv[1], sys.argv[2]
    for N in [int(a) for a in sys.argv[3:]]:
        r = run(mesh, mode, N)
        print(json.dumps({k: (float(f"{v:.5g}") if isinstance(v, float) else v) for k, v in r.items()}), flush=True)
