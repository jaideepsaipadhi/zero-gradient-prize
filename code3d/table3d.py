"""Print convergence tables from results.jsonl (last record wins per key).  python3 table3d.py [results.jsonl ...]"""
import sys, os, json, math
paths = sys.argv[1:] or [os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "notes", "v6", "num3d", "results.jsonl")]
recs = {}
for p in paths:
    for l in open(p):
        r = json.loads(l)
        key = (r["test"], r.get("method", "-"), r["k"], r.get("mu"), r.get("L", 1.0), r["N"])
        recs[key] = r
groups = {}
for key, r in sorted(recs.items(), key=lambda kv: (kv[0][:5], kv[0][5])):
    groups.setdefault(key[:5], []).append(r)


def rate(a, b, ha, hb):
    try:
        return math.log(a / b) / math.log(ha / hb)
    except Exception:
        return float("nan")


for (test, meth, k, mu, L), rs in groups.items():
    print(f"\n== {test}  method={meth}  k={k}  mu={mu:g}  L={L:g}   [NUMERICAL]")
    if test == "box_mms":
        print(f"{'N':>3} {'vel dofs':>9} {'p dofs':>8} {'h':>7} {'H1 err':>10} {'rate':>5} {'L2 err':>10} {'rate':>5} {'|div uh|':>9} {'ipm':>4} {'relres':>8}")
        for i, r in enumerate(rs):
            rh = rate(rs[i - 1]["H1"], r["H1"], rs[i - 1]["h"], r["h"]) if i else float("nan")
            rl = rate(rs[i - 1]["L2"], r["L2"], rs[i - 1]["h"], r["h"]) if i else float("nan")
            print(f"{r['N']:3d} {r['ndof']:9d} {r['npdof']:8d} {r['h']:7.4f} {r['H1']:10.3e} {rh:5.2f} {r['L2']:10.3e} {rl:5.2f} "
                  f"{r['div']:9.1e} {r['ipm_its']:4d} {r['relres']:8.1e}")
        continue
    hasG = "slip_ratio" in rs[0]
    print(f"{'N':>3} {'m':>2} {'vel dofs':>9} {'p dofs':>8} {'hG':>7} {'h':>6} {'H1 err':>10} {'rate':>5} {'r(h)':>5} {'|e|_Gh':>9} {'rate':>5} {'energy':>9} {'rate':>5} "
          f"{'|div uh|':>9}" + (f" {'slip':>6} {'un':>6} {'leak':>6} {'H1rat':>6}" if hasG else ""))
    for i, r in enumerate(rs):
        pr = rs[i - 1] if i else None
        f = lambda key: rate(pr[key], r[key], pr["hG"], r["hG"]) if pr else float("nan")
        line = (f"{r['N']:3d} {r['m']:2d} {r['ndof']:9d} {r.get('npdof', 0):8d} {r['hG']:7.4f} {r['h']:6.3f} {r['H1']:10.3e} {f('H1'):5.2f} {(rate(pr['H1'], r['H1'], pr['h'], r['h']) if pr else float('nan')):5.2f} "
                f"{r['bL2']:9.2e} {f('bL2'):5.2f} {r['energy']:9.2e} {f('energy'):5.2f} {r.get('div', float('nan')):9.1e}")
        if hasG:
            line += f" {r['slip_ratio']:6.3f} {r['un_ratio']:6.3f} {r['leak_corr']:6.3f} {r['H1_ratio']:6.2f}"
        print(line)
print("\nrates are w.r.t. h (box) or hG (sphere) between consecutive levels; r(h) = H1 rate w.r.t. h.  slip = ||e||_{Gamma_h}/((h/mu)G),"
      "\nun = ||u_h.n||_{Gamma_h}/((h/mu)G), leak = <u_h.n, p-pbar>_{Gamma_h}/((h/mu)G^2), H1rat = ||grad e||/((h/mu)G);"
      "\nG = ||p-pbar||_{Gamma_h}.  For method D=GS-CNS the 'error' is D itself (u_h := D).")
