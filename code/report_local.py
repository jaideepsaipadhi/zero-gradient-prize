"""Tables for the local runs (tests P and C): python3 report_local.py [results dir, default ../results/local]"""
import json, glob, os, sys, math
ROOT = os.path.dirname(os.path.abspath(__file__))
RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "..", "results", "local")
POD = os.path.join(ROOT, "..", "results", "pod")
G = math.sqrt(7 * math.pi / 8)
def load(kind, mu, base=RES):
    out = []
    for fn in glob.glob(os.path.join(base, "study", f"N*_{kind}_mu{mu:g}.jsonl")):
        out.append([json.loads(l) for l in open(fn)])
    return sorted(out, key=lambda r: r[0]["N"])
print("Test P (u = 0, f = grad p, p as in B1): GS output = exact traction response; CNS* must return u_h = 0.")
for mu in (100.0, 1000.0, 10000.0):
    print(f"\n  mu = {mu:g}")
    print("     N       h   GS H1*mu/(hG)  slip*mu/(hG)  en*(mu/h)^.5/G |  CNS* H1    CNS* slip  divres")
    for R in load("P1", mu):
        g, c = R[0], R[1]; h = g["h"]
        print(f"  {g['N']:4d}  {h:.4f}   {g['H1']*mu/(h*G):10.4f}  {g['bL2']*mu/(h*G):12.4f}  {g['energy']*math.sqrt(mu/h)/G:14.4f} |  {c['H1']:.1e}   {c['bL2']:.1e}  {g['divres']:.1e}")
print("\nTest C (psi as in B, p = exp(x/2) cos y), mu = 100: CNS* vs B1 (polynomial pressure) and GS rate")
print("     N       h    CNS* H1       rate   rel.diff vs B1   GS H1        rate")
prev = None
B = {r[0]["N"]: r for r in load("B1", 100.0, POD)}
for R in load("C", 100.0):
    g, c = R[0], R[1]; h = g["h"]; N = g["N"]
    rc = rg = float("nan")
    if prev:
        rc = math.log(c["H1"] / prev[0]) / math.log(h / prev[2]); rg = math.log(g["H1"] / prev[1]) / math.log(h / prev[2])
    rd = abs(c["H1"] / B[N][1]["H1"] - 1) if N in B else float("nan")
    print(f"  {N:4d}  {h:.4f}  {c['H1']:.6e}  {rc:5.2f}   {rd:.1e}        {g['H1']:.4e}  {rg:5.2f}")
    prev = (c["H1"], g["H1"], h)
