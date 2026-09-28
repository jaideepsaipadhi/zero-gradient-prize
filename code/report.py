"""Summarize results/: refinement tables with observed rates and theory-normalized ratios, and threshold studies."""
import json, glob, os, numpy as np
from collections import defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
rows = []
import sys
RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "..", "results", "pod")
for fn in glob.glob(os.path.join(RES, "study", "*.jsonl")):
    rows += [json.loads(l) for l in open(fn)]
G1 = np.sqrt(7 * np.pi / 8)          # ||p - pbar||_{L2(Gamma)} for test B with amplitude 1
tab = defaultdict(list)
for r in rows:
    tab[(r["kind"], r["mu"], str(r["theta"]))].append(r)


def amp(kind):
    return 0.0 if kind == "A" else float(kind[1:])


def rates(h, e):
    return [np.nan] + [np.log(e[i] / e[i - 1]) / np.log(h[i] / h[i - 1]) for i in range(1, len(e))]


print("PPL 115 numerics — Scott–Vogelius P4, Nitsche on inscribed polygon, strong Dirichlet on square")
print("GS = Gjerde–Scott (27)-(28) as printed; CNS = pressure-consistent <p - pbar_Gamma(p), v.n>;")
print("D = u_GS - u_CNS isolates the response to the missing -p n traction.  h = h_Omega (penalty h).")
print("Theory: GS/D:  ||grad D|| ~ c1 h/mu G,  slip ~ c2 h/mu G,  energy ~ c3 (h/mu)^(1/2) G   (G = ||p - pbar||_Gamma)")
print("        CNS:   rate 3/2 in H1 (gamma^(1/2) h_Gamma^(3/2) + h^4), independent of pressure amplitude")
bad = [r for r in rows if r.get("theta") in (0, 1) and (r.get("relres", 0) > 1e-8 or r.get("divres", 0) > 1e-6)]
print(f"\nsolver health: {len(rows)} records; {len(bad)} with relres>1e-8 or divres>1e-6")
for r in bad:
    print("   CHECK", r["N"], r["kind"], r["mu"], r["theta"], r.get("relres"), r.get("divres"))

for kind in ("A", "B1", "B100"):
    for mu in (100.0, 1000.0, 10000.0):
        if not tab[(kind, mu, "0")]:
            continue
        print(f"\n=== test {kind}  mu={mu:g}   G = {amp(kind) * G1:.3f}")
        for th, lab in (("0", "GS (printed)"), ("1", "CNS (consistent)"), ("D", "D = GS - CNS")):
            R = sorted(tab[(kind, mu, th)], key=lambda r: r["N"])
            if not R:
                continue
            h = np.array([r["h"] for r in R])
            print(f"  {lab}")
            head = f"    {'N':>4} {'h':>7} {'H1':>10} {'rt':>5} {'energy':>10} {'rt':>5} {'slip':>10} {'rt':>5} {'bdn':>10} {'rt':>5}"
            if th == "D" and kind != "A":
                head += f" {'H1*mu/hG':>9} {'slip*mu/hG':>10} {'en*(mu/h)^.5/G':>14}"
            print(head)
            cols = {k: np.array([r[k] for r in R]) for k in ("H1", "energy", "bL2", "bdn")}
            rt = {k: rates(h, v) for k, v in cols.items()}
            for i, r in enumerate(R):
                s = f"    {r['N']:>4} {h[i]:7.4f}"
                for k in ("H1", "energy", "bL2", "bdn"):
                    s += f" {cols[k][i]:10.3e} {rt[k][i]:5.2f}"
                if th == "D" and kind != "A":
                    G = amp(kind) * G1
                    s += f" {cols['H1'][i] * mu / (h[i] * G):9.4f} {cols['bL2'][i] * mu / (h[i] * G):10.4f} {cols['energy'][i] * np.sqrt(mu / h[i]) / G:14.4f}"
                print(s)

T = [json.load(open(fn)) for fn in glob.glob(os.path.join(RES, "thresh", "*.json"))]
if T:
    print("\n=== Coercivity threshold mu* on Z_h (sign of pivots of A_mu + r B^T M^-1 B)")
    print("  fixed rho, refinement (L = 2.5, layers = N/4):")
    for t in sorted([t for t in T if t["L"] == 2.5 and t["layers"] == t["N"] // 4], key=lambda t: t["N"]):
        print(f"    N={t['N']:3d} ntri={t['ntri']:6d} rho={t['rho']:7.3f} mu*={t['mu_star']:10.4f} gamma*=mu*/rho={t['gamma_star']:8.4f}")
    for N in sorted({t["N"] for t in T if t["L"] != 2.5}):
        print(f"  rho sweep at N={N} (enlarged outer square, inner layers kept shape-similar):")
        for t in sorted([t for t in T if t["N"] == N and (t["L"] != 2.5 or t["layers"] == N // 4)], key=lambda t: t["L"]):
            print(f"    L={t['L']:5.1f} layers={t['layers']:4d} rho={t['rho']:8.3f} mu*={t['mu_star']:11.4f} gamma*={t['gamma_star']:8.4f}")
