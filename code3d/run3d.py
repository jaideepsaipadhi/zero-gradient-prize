"""Drivers for the 3D tests.  Usage:
  python3 run3d.py box   K N [N ...]            -> manufactured solution on (-1,1)^3, consistent Nitsche on all faces
  python3 run3d.py sph   K KIND MU N [N ...]    -> Omega_h = box minus inscribed polyhedron; KIND in rot, sph, P; GS and CNS*
  python3 run3d.py interp K KIND MU N [N ...]   -> only the nodal-interpolation error of the exact solution (no solve)
Results appended to ../notes/v6/num3d/results.jsonl (one JSON record per solve)."""
import sys, os, time, json, resource
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svn3d, exact3d

OUT = os.environ.get("NUM3D_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "notes", "v6", "num3d", "results.jsonl"))


def rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2


def dump(rec):
    rec = {k: (float(v) if isinstance(v, (np.floating, float)) else (int(v) if isinstance(v, (np.integer,)) else v))
           for k, v in rec.items()}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(json.dumps(rec), flush=True)


def run_box(k, N, mu=None, r=1e4):
    t0 = time.time()
    ex = exact3d.make("mms")
    S = svn3d.build("box", N, k, L=1.0)
    mu = mu or float(os.environ.get("NUM3D_MUBOX", "400"))
    allf = S.fBox
    sysm = svn3d.assemble(S, mu, ex["f"], [(allf, ex["u"], "sym")], [(allf, True)])
    u, p, hist, relres = svn3d.solve_ipm(S, sysm, None, r=r)
    E = svn3d.errors(S, u, ex, mu, nq=k + 3)
    dump(dict(test="box_mms", k=k, N=N, mu=mu, h=S.h, ntet=len(S.T), ndof=3 * S.nn, npdof=sysm["B"].shape[0],
              H1=E["H1"], L2=E["L2"], div=svn3d.div_norm(S, sysm, u), ipm_its=len(hist), relres=relres,
              secs=time.time() - t0, rss_gb=rss()))


def interp_record(S, ex, k, kind, mu, N, L):
    """nodal interpolant of the (extended) exact solution: approximation part of the error, no solve."""
    uI = ex["u"](S.xyz[:, 0], S.xyz[:, 1], S.xyz[:, 2]).ravel()
    E = svn3d.errors(S, uI, ex, mu, nq=k + 3)
    rec = dict(test=f"sphere_{kind}", method="interp", k=k, N=N, m=S.info["m"], L=L, mu=mu, h=S.h, hG=S.hGamma,
               ntet=len(S.T), ndof=3 * S.nn)
    rec.update({kk: E[kk] for kk in ("H1", "L2", "bL2", "bdn", "energy")})
    dump(rec)


def run_sph(k, kind, mu, N, L=2.0, r=1e4, methods=("GS", "CNS"), symK=True):
    ex = exact3d.make(kind)
    t0 = time.time()
    S = svn3d.build("sph", N, k, L=L)
    sols = {}
    interp_record(S, ex, k, kind, mu, N, L)
    for meth in methods:
        t1 = time.time()
        pc = [(S.fGamma, True)] if meth == "CNS" else []
        sysm = svn3d.assemble(S, mu, ex["f"], [(S.fGamma, None, "grad")], pc)
        u, p, hist, relres = svn3d.solve_ipm(S, sysm, ex["u"], r=r, symK=symK)
        sols[meth] = u
        E = svn3d.errors(S, u, ex, mu, nq=k + 3)
        rec = dict(test=f"sphere_{kind}", method=meth, k=k, N=N, m=S.info["m"], L=L, mu=mu, h=S.h, hG=S.hGamma,
                   ntet=len(S.T), ndof=3 * S.nn, npdof=sysm["B"].shape[0], div=svn3d.div_norm(S, sysm, u),
                   ipm_its=len(hist), relres=relres, secs=time.time() - t1, rss_gb=rss())
        rec.update(E)
        dump(rec)
        del sysm
    if len(sols) == 2:
        # D = u_GS - u_CNS*: theory D ~ (h/mu) v_phi, v_phi|Gamma ~ -(p - pbar) n  ->  D.n ~ (h/mu)(p - pbar)
        zero = dict(u=lambda a, b, c: np.zeros(np.shape(a) + (3,)), gu=lambda a, b, c: np.zeros(np.shape(a) + (3, 3)),
                    p=ex["p"])
        E = svn3d.errors(S, sols["GS"] - sols["CNS"], zero, mu, nq=k + 3)   # norms of D; leak_corr uses D.n
        rec = dict(test=f"sphere_{kind}", method="D=GS-CNS", k=k, N=N, m=S.info["m"], L=L, mu=mu, h=S.h, hG=S.hGamma,
                   ntet=len(S.T), ndof=3 * S.nn)
        rec.update(E)
        dump(rec)


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "box":
        for N in a[2:]:
            run_box(int(a[1]), int(N))
    elif a[0] == "interp":
        for N in a[4:]:
            S = svn3d.build("sph", int(N), int(a[1]), L=2.0)
            interp_record(S, exact3d.make(a[2]), int(a[1]), a[2], float(a[3]), int(N), 2.0)
    elif a[0] == "sph":
        L = float(os.environ.get("NUM3D_L", "2.0"))
        for N in a[4:]:
            run_sph(int(a[1]), a[2], float(a[3]), int(N), L=L)
