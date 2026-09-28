"""One job = one unit of work for the degree-general solver svn_k (own process; an OOM kill loses one job).

  python3 job_k.py study  K MESH N KIND MU   -> results/v3/study/{MESH}_k{K}_N{N}_{KIND}_mu{MU}.jsonl
                                                (GS theta=0, CNS* theta=1, and D = GS - CNS*)
  python3 job_k.py thresh K MESH N           -> results/v3/thresh/{MESH}_k{K}_N{N}.json  (coercivity threshold mu*)

MESH is "std" (svn.make_mesh) or "ct" (barycentric / Clough-Tocher refinement of it).
Env SVN_OUT overrides the output root (default results/v3).  Env SVN_SOLVER = "lu" (SuperLU, default,
identical to job.py) or "pardiso" (MKL PARDISO, pressure block regularised by 1e-12*M, + refinement; much less memory
at large N; agrees with "lu" to ~1e-10).
"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time, warnings, resource
warnings.filterwarnings("ignore")
import numpy as np
import svn_k

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ.get("SVN_OUT", os.path.join(ROOT, "..", "results", "v3"))


def _maxrss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2


def _write(path, recs):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + ".tmp", "w") as f:
        for r in recs:
            f.write(json.dumps(r) + "\n")
    os.replace(path + ".tmp", path)


def study(k, mesh, N, kind, mu):
    t0 = time.time()
    S = svn_k.build(N, k, mesh)
    ex = svn_k.make_exact(kind)
    sysm = svn_k.assemble(S, mu, ex["f"])
    zero = dict(u=lambda X, Y: np.zeros(np.shape(X) + (2,)), gu=lambda X, Y: np.zeros(np.shape(X) + (2, 2)))
    base = dict(k=k, mesh=mesh, N=N, kind=kind, mu=mu, h=S.h, hG=S.hGamma, ntri=len(S.tris),
                ndof=2 * S.nn, npres=int(sysm["B"].shape[0]))
    solver = os.environ.get("SVN_SOLVER", "lu")          # "lu" (SuperLU, reference) or "pardiso"
    base["solver"] = solver
    solve = dict(lu=svn_k.solve_lu, pardiso=svn_k.solve_pardiso)[solver]
    recs, sol = [], {}
    for th in (0, 1):
        t1 = time.time()
        u, p, hist, relres, nnzlu = solve(S, sysm, ex["u"], th)
        E = svn_k.errors(S, u, ex, mu)
        E.update(base); E.update(theta=th, divres=float(hist[-1]), relres=float(relres), nnzLU=int(nnzlu),
                                 secs=time.time() - t1)
        recs.append({kk: (float(v) if isinstance(v, (np.floating, float)) else v) for kk, v in E.items()})
        sol[th] = u
    Dn = svn_k.errors(S, sol[0] - sol[1], zero, mu)
    d = dict(base); d.update(theta="D", H1=float(Dn["H1"]), bL2=float(Dn["bL2"]), energy=float(Dn["energy"]),
                             bdn=float(Dn["bdn"]), secs=time.time() - t0, maxrss_gb=_maxrss_gb())
    recs.append(d)
    _write(os.path.join(OUT, "study", f"{mesh}_k{k}_N{N}_{kind}_mu{mu:g}.jsonl"), recs)


def thresh(k, mesh, N):
    import penalty_k
    t0 = time.time()
    S = svn_k.build(N, k, mesh)
    mu = penalty_k.threshold(S)
    rec = dict(k=k, mesh=mesh, N=N, ntri=len(S.tris), hG=S.hGamma, hOmega=S.h, rho=S.h / S.hGamma,
               mu_star=float(mu), gamma_star=float(mu * S.hGamma / S.h), secs=time.time() - t0,
               maxrss_gb=_maxrss_gb())
    _write(os.path.join(OUT, "thresh", f"{mesh}_k{k}_N{N}.json"), [rec])


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "study":
        study(int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), sys.argv[5], float(sys.argv[6]))
    elif mode == "thresh":
        thresh(int(sys.argv[2]), sys.argv[3], int(sys.argv[4]))
