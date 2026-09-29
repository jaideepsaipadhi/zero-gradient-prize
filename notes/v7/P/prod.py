"""notes/v7/P: slim CNS* run on the production meshes (code/svn.make_mesh), boundary moments P_h(phi) only.
usage: python3 prod.py KIND MU N1,N2    (KIND in S1, B1)"""
import os, sys, json, time
for v_ in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v_, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
import svn
from annulus import exact_S

kind, MU = sys.argv[1], float(sys.argv[2])
ex = exact_S("S1") if kind == "S1" else svn.make_exact(kind)
for N in [int(s) for s in sys.argv[3].split(",")]:
    t0 = time.time()
    verts, tris, circ, outer = svn.make_mesh(N); S = svn.Space(verts, tris, circ, outer); nt = len(tris)
    sysm = svn.assemble(S, MU, ex["f"], None); A, B, F = sysm["A"], sysm["B"], sysm["F"]
    Bt = svn.coupling(sysm, 1); ndof = A.shape[0]
    dD = svn.dofs(S.dnodes).ravel(); free = np.setdiff1d(np.arange(ndof), dD)
    ug = np.zeros(ndof); ug[dD] = ex["u"](S.xy[S.dnodes, 0], S.xy[S.dnodes, 1]).ravel()
    K = sp.bmat([[A[free][:, free], -Bt[:, free].T], [B[:, free], None]], format="csc")
    rhs = np.concatenate([F[free] - A[free][:, dD] @ ug[dD], -B[:, dD] @ ug[dD]])
    lu = spl.splu(K, permc_spec="COLAMD", diag_pivot_thresh=1.0); x = lu.solve(rhs)
    for _ in range(3):
        x = x + lu.solve(rhs - K @ x)
    p = x[len(free):]
    ED = S.edge_data(); num = per = 0.0
    for (t, ref, Xe, we, n) in ED:
        num += (we * ex["p"](Xe[:, 0], Xe[:, 1])).sum(); per += we.sum()
    pbar = num / per; eds = []; epG = 0.0
    for (t, ref, Xe, we, n) in ED:
        epe = ex["p"](Xe[:, 0], Xe[:, 1]) - pbar - svn._mono(svn.MON3, ref[:, 0], ref[:, 1]) @ p[10 * t:10 * t + 10]
        eds.append((Xe, we, epe)); epG += (we * epe).sum()
    epG /= per
    rec = dict(kind=kind, mu=MU, N=N, h=S.h, hG=S.hGamma)
    for name, f in {"cos1": np.cos, "cos2": lambda t: np.cos(2 * t), "cos4": lambda t: np.cos(4 * t)}.items():
        rec["P_" + name + "/hG"] = float(sum((we * (epe - epG) * f(np.arctan2(Xe[:, 1], Xe[:, 0]))).sum()
                                             for (Xe, we, epe) in eds) / S.hGamma)
    rec["secs"] = round(time.time() - t0, 1)
    print(json.dumps(rec), flush=True)
