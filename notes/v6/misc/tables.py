"""Collect the misc_numerics runs into the tables of REPORT.tex (CNS* rows; GS where noted)."""
import json, glob, math
recs = []
for f in sorted(glob.glob("run_*.jsonl")):
    for l in open(f):
        recs.append(json.loads(l))
seen = {}
for r in recs:
    seen[(r["kind"], r["mu"], r["N"], r["method"])] = r      # later files override
def rows(kind, mu, method="CNS*"):
    return sorted([r for (k, m, N, me), r in seen.items() if k == kind and m == mu and me == method], key=lambda r: r["N"])
def rate(a, b, ka, key="hG"):
    return math.log(abs(a[ka]) / abs(b[ka])) / math.log(a[key] / b[key])
for kind in ("S1", "B1", "A"):
    for mu in (100.0, 1000.0):
        R = rows(kind, mu)
        if not R: continue
        print(f"\n== {kind} mu={mu:g} CNS*:  N hG | ||e_p-mean|| (rate) | ||div u^L|| | |||u_h-u^L||| | ratio ep/div  ep/E | vol_y (rate) | P_cos1/hG P_cos4/hG | (mu/h)S+P cos1 (rate)")
        prev = None
        for r in R:
            s = f"{r['N']:4d} {r['hG']:.4f} | {r['ep_L2opt']:.3e}"
            s += f" ({rate(prev, r, 'ep_L2opt'):.2f})" if prev else "       "
            s += f" | {r.get('divL', float('nan')):.3e} | {r.get('dEnergy', float('nan')):.3e} | {r['ep_L2opt']/r.get('divL', float('nan')):.2f} {r['ep_L2opt']/r.get('dEnergy', float('nan')):.2f}"
            s += f" | {r['vol_y']:+.3e}" + (f" ({rate(prev, r, 'vol_y'):.2f})" if prev else "       ")
            s += f" | {r['P_cos1']/r['hG']:+.3f} {r.get('P_cos4', float('nan'))/r['hG']:+.3f}"
            s += f" | {r['muS+P_cos1']:+.2e}" + (f" ({rate(prev, r, 'muS+P_cos1'):.2f})" if prev else "")
            print(s); prev = r
print("\n== torque: -omega(zeta^1) / predicted (h/mu) int sigma n . d_n psi   and  (J^chi - J^N)(e1)")
for kind in ("S1", "B1", "A", "Q1"):
    for mu in (100.0, 1000.0, 10000.0, 100000.0):
        R = rows(kind, mu)
        for r in R:
            if "-omega(zeta1)_rot" in r:
                print(f"{kind:3s} mu={mu:8g} N={r['N']:3d} gamma={mu*r['hG']/r['h']:7.1f} gamma*hG={mu*r['hG']**2/r['h']:6.2f}  "
                      f"-omega(z1)={r['-omega(zeta1)_rot']:+.4e} pred={r['pred_rot']:+.4e} ratio={r['-omega(zeta1)_rot']/r['pred_rot']:.3f}"
                      f"  Jchi-JN(rot)={r['Jchi-JN_rot']:+.3e}  Jchi-JN(e1)={r['Jchi-JN_e1']:+.3e}")
