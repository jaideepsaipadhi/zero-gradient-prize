"""One job = one unit of work, run as its own process (so an OOM kill loses only this job).

  python3 job.py study   N KIND MU        -> results/study/N{N}_{KIND}_mu{MU}.jsonl  (GS, CNS, and D = GS - CNS)
  python3 job.py thresh  N L LAYERS       -> results/thresh/N{N}_L{L}_m{LAYERS}.json (coercivity threshold mu*)
"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
import svn

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ.get("SVN_OUT", os.path.join(ROOT, "..", "results", "rerun"))


def study(N, kind, mu):
    t0 = time.time()
    verts, tris, circ, outer = svn.make_mesh(N)
    S = svn.Space(verts, tris, circ, outer)
    ex = svn.make_exact(kind)
    sysm = svn.assemble(S, mu, ex["f"], None)
    zero = dict(u=lambda X, Y: np.zeros(np.shape(X) + (2,)), gu=lambda X, Y: np.zeros(np.shape(X) + (2, 2)))
    recs, sol = [], {}
    for th in (0, 1):
        t1 = time.time()
        u, p, hist, relres = svn.solve_lu(S, sysm, ex["u"], th)
        E = svn.errors(S, u, ex, mu)
        E.update(N=N, kind=kind, mu=mu, theta=th, h=S.h, hG=S.hGamma, ntri=len(tris), ndof=2 * S.nn,
                 divres=float(hist[-1]), relres=float(relres), secs=time.time() - t1)
        recs.append({k: (float(v) if isinstance(v, (np.floating, float)) else v) for k, v in E.items()})
        sol[th] = u
    Dn = svn.errors(S, sol[0] - sol[1], zero, mu)
    recs.append(dict(N=N, kind=kind, mu=mu, theta="D", h=S.h, hG=S.hGamma,
                     H1=float(Dn["H1"]), bL2=float(Dn["bL2"]), energy=float(Dn["energy"]), bdn=float(Dn["bdn"]),
                     secs=time.time() - t0))
    os.makedirs(os.path.join(OUT, "study"), exist_ok=True)
    fn = os.path.join(OUT, "study", f"N{N}_{kind}_mu{mu:g}.jsonl")
    with open(fn + ".tmp", "w") as f:
        for r in recs:
            f.write(json.dumps(r) + "\n")
    os.replace(fn + ".tmp", fn)


def thresh(N, L, m):
    import penalty
    t0 = time.time()
    S = penalty.build(N, L=L, layers=m)
    mu = penalty.threshold(S)
    rec = dict(N=N, L=L, layers=m, ntri=len(S.tris), hG=S.hGamma, hOmega=S.h, rho=S.h / S.hGamma,
               mu_star=float(mu), gamma_star=float(mu * S.hGamma / S.h), secs=time.time() - t0)
    os.makedirs(os.path.join(OUT, "thresh"), exist_ok=True)
    fn = os.path.join(OUT, "thresh", f"N{N}_L{L:g}_m{m}.json")
    with open(fn + ".tmp", "w") as f:
        json.dump(rec, f)
    os.replace(fn + ".tmp", fn)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "study":
        study(int(sys.argv[2]), sys.argv[3], float(sys.argv[4]))
    elif mode == "thresh":
        thresh(int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]))
