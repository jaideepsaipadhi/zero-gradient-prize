"""Collect existing L2(Omega_h) velocity errors (||u~ - u_h||) from results/pod/study and results/local/study
and print tables with rates in h_Gamma.  Theta 0 = GS, 1 = CNS*.  Usage: python3 collect_l2.py"""
import json, glob, os, math, collections
ROOT = os.path.join(os.path.dirname(__file__), "../../../results")
rows = collections.defaultdict(dict)
for d in ["pod/study", "local/study"]:
    for f in glob.glob(os.path.join(ROOT, d, "*.jsonl")):
        for line in open(f):
            r = json.loads(line)
            if r.get("theta") not in (0, 1) or "L2" not in r:
                continue
            key = (r["kind"], r["mu"], r["theta"])
            rows[key][r["N"]] = r
for key in sorted(rows):
    kind, mu, th = key
    print(f"\n== {kind} mu={mu:g} {'GS' if th == 0 else 'CNS*'}")
    print(f"{'N':>4} {'hG':>8} {'L2':>11} {'rate':>6} {'H1':>11} {'rate':>6} {'L2*mu/(h)':>10} {'L2/hG^2':>9}")
    prev = None
    for N in sorted(rows[key]):
        r = rows[key][N]
        s = f"{N:4d} {r['hG']:8.4f} {r['L2']:11.4e} "
        if prev:
            lr = math.log(prev['hG'] / r['hG'])
            s += f"{math.log(prev['L2'] / r['L2']) / lr:6.2f} {r['H1']:11.4e} {math.log(prev['H1'] / r['H1']) / lr:6.2f}"
        else:
            s += f"{'':6} {r['H1']:11.4e} {'':6}"
        s += f" {r['L2'] * mu / r['h']:10.4f} {r['L2'] / r['hG']**2:9.4f}"
        print(s); prev = r
