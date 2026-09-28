import json, sys, numpy as np
from collections import defaultdict
rows = [json.loads(l) for l in open(sys.argv[1] if len(sys.argv) > 1 else "study.jsonl")]
G1 = np.sqrt(7 * np.pi / 8)
tab = defaultdict(list)
for r in rows:
    tab[(r["kind"], r["mu"], str(r["theta"]))].append(r)


def amp(kind):
    return 0.0 if kind == "A" else float(kind[1:])


def rates(h, e):
    return [np.nan] + [np.log(e[i] / e[i - 1]) / np.log(h[i] / h[i - 1]) for i in range(1, len(e))]


for kind in ("A", "B1", "B100"):
    for mu in (100.0, 1000.0, 10000.0):
        print(f"\n=== {kind}  mu={mu:g}  (G = {amp(kind)*G1:.3f})")
        for th, lab in (("0", "GS (printed)"), ("1", "consistent"), ("D", "D = GS - CNS")):
            R = sorted(tab[(kind, mu, th)], key=lambda r: r["N"])
            if not R:
                continue
            h = np.array([r["h"] for r in R])
            print(f"  {lab}")
            head = f"    {'N':>4} {'h':>7} {'H1':>10} {'rt':>5} {'energy':>10} {'rt':>5} {'slip':>10} {'rt':>5} {'bdn':>10} {'rt':>5}"
            if th == "D" and kind != "A":
                head += f" {'H1*mu/(hG)':>11} {'slip*mu/(hG)':>12} {'en*sqrt(mu/h)/G':>15}"
            print(head)
            cols = {k: np.array([r[k] for r in R]) for k in ("H1", "energy", "bL2", "bdn")}
            rt = {k: rates(h, v) for k, v in cols.items()}
            for i, r in enumerate(R):
                s = f"    {r['N']:>4} {h[i]:7.4f}"
                for k in ("H1", "energy", "bL2", "bdn"):
                    s += f" {cols[k][i]:10.3e} {rt[k][i]:5.2f}"
                if th == "D" and kind != "A":
                    G = amp(kind) * G1
                    s += f" {cols['H1'][i]*mu/(h[i]*G):11.4f} {cols['bL2'][i]*mu/(h[i]*G):12.4f} {cols['energy'][i]*np.sqrt(mu/h[i])/G:15.4f}"
                print(s)
