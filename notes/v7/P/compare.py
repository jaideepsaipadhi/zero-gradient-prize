"""notes/v7/P: measured (notes/v6/misc run_*.jsonl, CNS*) vs cell-predicted boundary moments P_h(phi)/h_Gamma on the
production meshes.  'pred' = full formula, 'no-k0' = without the global flux-multiplier term k0*Ac.
usage: python3 compare.py      (output: compare.log)"""
import os, sys, json, glob
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from predict import predict

MISC = os.path.join(HERE, "..", "..", "v6", "misc")
meas = {}
for f in glob.glob(os.path.join(MISC, "run_*.jsonl")):
    for l in open(f):
        try:
            r = json.loads(l)
        except Exception:
            continue
        if r.get("method") == "CNS*":
            meas[(r["kind"], r["mu"], r["N"])] = r
cases = [("S1", 100.0, (16, 32, 64)), ("B1", 100.0, (16, 32, 64)), ("S1", 1000.0, (32, 64)), ("B1", 1000.0, (64,))]
print("kind mu N | moment: measured  predicted (no-k0)  ratio meas/pred")
for kind, mu, Ns in cases:
    for N in Ns:
        m = meas.get((kind, mu, N))
        if m is None:
            continue
        p = predict(kind, mu, N, "limit")
        hG = m["hG"]
        parts = []
        for ph in ("cos1", "cos2", "sin3", "cos4"):
            if "P_" + ph not in m or m["P_" + ph] is None or not np.isfinite(m["P_" + ph]):
                continue
            mv = m["P_" + ph] / hG; pv = p["P_" + ph + "/hG"]; pnk = pv - p["Pc_" + ph + "/hG"]
            if abs(pv) < 1e-8 and abs(mv) < 1e-8:
                continue
            parts.append(f"{ph}: {mv:+.3f} {pv:+.3f} ({pnk:+.3f}) {mv / pv if abs(pv) > 1e-8 else float('nan'):.2f}")
        print(f"{kind} {mu:6.0f} {N:3d} | " + " | ".join(parts), flush=True)
for kind, mu in (("S1", 100.0), ("B1", 100.0), ("S1", 1000.0), ("B1", 1000.0)):
    p = predict(kind, mu, 512, "limit")
    print(f"limit (N=512 geometry) {kind} mu={mu:.0f}: " + " ".join(f"{ph}={p['P_' + ph + '/hG']:+.3f}" for ph in ("cos1", "cos2", "sin3", "cos4")),
          f"lam/lam_res in [{p['lam_over_res_min']:.3f},{p['lam_over_res_max']:.3f}]", flush=True)
