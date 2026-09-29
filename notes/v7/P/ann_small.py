"""notes/v7/P: annulus.py with a smaller outer radius (ROUT, default 1.4) so that N = 256 fits in ~1 min.
The exact solution is imposed on the outer ring, so the limit of P_h/h_G does not depend on ROUT.
usage: python3 ann_small.py DELTA A LAMHAT N1,N2 [S1|S2] [ROUT]"""
import sys, json, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import annulus
annulus.ROUT = float(sys.argv[6]) if len(sys.argv) > 6 else 1.4
_orig = annulus.annulus_mesh
annulus.annulus_mesh = lambda N, a, delta, rout=None: _orig(N, a, delta, annulus.ROUT)
delta, a, lamhat = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
for N in [int(s) for s in sys.argv[4].split(",")]:
    r = annulus.run(N, a, delta, lamhat, sys.argv[5] if len(sys.argv) > 5 else "S1"); r["rout"] = annulus.ROUT
    print(json.dumps({k: (float(f"{v:.6g}") if isinstance(v, float) else v) for k, v in r.items()}), flush=True)
