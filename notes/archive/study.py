"""Refinement study: GS (theta=0, as printed) vs consistent (theta=1) on
  A      : Scott's shear flow (p = 0)
  B1     : non-constant pressure, amplitude 1
  B100   : same velocity, pressure amplitude 100 (pressure-dominated, drag-like)
Also measures D = u_GS - u_CNS, which isolates the response to the missing -p n term.
"""
import json, sys, time, numpy as np, warnings
warnings.filterwarnings("ignore")
import svn

NS = [int(a) for a in sys.argv[1].split(",")] if len(sys.argv) > 1 else [16, 24, 32, 48, 64, 96, 128]
MUS = [float(a) for a in sys.argv[2].split(",")] if len(sys.argv) > 2 else [100.0, 1000.0, 10000.0]
KINDS = sys.argv[3].split(",") if len(sys.argv) > 3 else ["A", "B1", "B100"]
out = open(sys.argv[4] if len(sys.argv) > 4 else "study.jsonl", "a")

zero = dict(u=lambda X, Y: np.zeros(np.shape(X) + (2,)), gu=lambda X, Y: np.zeros(np.shape(X) + (2, 2)))
for N in NS:
    verts, tris, circ, outer = svn.make_mesh(N)
    S = svn.Space(verts, tris, circ, outer)
    for kind in KINDS:
        ex = svn.make_exact(kind)
        for mu in MUS:
            t0 = time.time()
            sysm = svn.assemble(S, mu, ex["f"], None)
            sol = {}
            for th in (0, 1):
                u, w, hist, relres = svn.solve_lu(S, sysm, ex["u"], th)
                E = svn.errors(S, u, ex, mu)
                E.update(N=N, kind=kind, mu=mu, theta=th, h=S.h, hG=S.hGamma, ntri=len(tris),
                         divres=float(hist[-1]), relres=float(relres))
                sol[th] = u
                out.write(json.dumps({k: (float(v) if isinstance(v, (np.floating, float)) else v)
                                      for k, v in E.items()}) + "\n"); out.flush()
            Dn = svn.errors(S, sol[0] - sol[1], zero, mu)
            rec = dict(N=N, kind=kind, mu=mu, theta="D", h=S.h, hG=S.hGamma,
                       H1=float(Dn["H1"]), bL2=float(Dn["bL2"]), energy=float(Dn["energy"]), bdn=float(Dn["bdn"]))
            out.write(json.dumps(rec) + "\n"); out.flush()
            print(f"N={N} {kind} mu={mu:g}  done in {time.time()-t0:.1f}s", flush=True)
