"""Tables for the v3 runs (degree-general solver svn_k):  python3 report_v3.py [results dir, default ../results/v3]
Writes <dir>/report_v3.txt (and prints it)."""
import json, glob, os, sys, io
import numpy as np
from collections import defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "..", "results", "v3")
G1 = np.sqrt(7 * np.pi / 8)          # ||p - pbar||_{L2(Gamma)} for amplitude 1 (tests B1, P1)
out = io.StringIO()


def P(*a):
    print(*a, file=out)


rows = []
for fn in sorted(glob.glob(os.path.join(RES, "study", "*.jsonl"))):
    rows += [json.loads(l) for l in open(fn)]
T = [json.loads(open(fn).read().splitlines()[0]) for fn in glob.glob(os.path.join(RES, "thresh", "*.json"))]
thr = defaultdict(dict)
for t in T:
    thr[(t["mesh"], t["k"])][t["N"]] = t["mu_star"]

tab = defaultdict(list)
for r in rows:
    tab[(r["mesh"], r["k"], r["kind"], r["mu"], str(r["theta"]))].append(r)


def rates(h, e):
    with np.errstate(all="ignore"):
        return [np.nan] + [np.log(e[i] / e[i - 1]) / np.log(h[i] / h[i - 1]) for i in range(1, len(e))]


def mustar_note(mesh, k, mu):
    d = thr.get((mesh, k))
    if not d:
        return ""
    s = ", ".join(f"mu*(N={N})={v:.1f}" for N, v in sorted(d.items()))
    worst = max(d.values())
    flag = "   <<< mu BELOW coercivity threshold: results not covered by the theory" if mu <= worst else ""
    return f"   [{s}]{flag}"


P("PPL 115 numerics v3 — degree-general Scott–Vogelius (P_k velocity / P_{k-1}-disc pressure) + Nitsche on the")
P("inscribed polygon.  mesh std = svn.make_mesh; ct = barycentric (Clough–Tocher/Alfeld) refinement of it.")
P("h = h_Omega (max element diameter; for ct equal to the macro h), G = ||p - pbar||_Gamma = sqrt(7 pi/8) (amp. 1).")
P("GS = theta 0 (as printed); CNS* = theta 1 (mean-corrected pressure coupling); D = u_GS - u_CNS*.")
P("Theory: CNS* H1 rate -> 3/2 for k >= 2 (geometric term), independent of pressure amplitude;")
P("        GS H1 rate -> 1 when G > 0;  test P: GS slip*mu/(hG) -> 1, energy*(mu/h)^(1/2)/G -> 1;  CNS* on P: u_h = 0.")

solvers = sorted({r.get("solver", "?") for r in rows})
bad = [r for r in rows if r.get("theta") in (0, 1) and
       (not np.isfinite(r.get("relres", 0)) or r.get("relres", 0) > 1e-8 or r.get("divres", 0) > 1e-6)]
P(f"\nsolver health: {len(rows)} records (solvers: {', '.join(solvers)}); "
  f"{len(bad)} with relres > 1e-8 or divres > 1e-6")
for r in bad:
    P(f"   CHECK {r['mesh']} k={r['k']} N={r['N']} {r['kind']} mu={r['mu']:g} theta={r['theta']} "
      f"relres={r.get('relres'):.1e} divres={r.get('divres'):.1e}")
dmax = max((r["divres"] for r in rows if r.get("theta") in (0, 1)), default=np.nan)
rmax = max((r["relres"] for r in rows if r.get("theta") in (0, 1)), default=np.nan)
P(f"   max divres = {dmax:.1e}, max relres = {rmax:.1e}")

keys = sorted({(r["mesh"], r["k"], r["kind"], r["mu"]) for r in rows},
              key=lambda k: (k[0] != "std", k[1], k[2], k[3]))
for (mesh, k, kind, mu) in keys:
    P(f"\n=== mesh={mesh} k={k} test {kind} mu={mu:g}{mustar_note(mesh, k, mu)}")
    for th, lab in (("0", "GS"), ("1", "CNS*"), ("D", "D = GS - CNS*")):
        R = sorted(tab[(mesh, k, kind, mu, th)], key=lambda r: r["N"])
        if not R:
            continue
        h = np.array([r["h"] for r in R])
        cols = {c: np.array([r[c] for r in R]) for c in ("H1", "energy", "bL2", "bdn")}
        rt = {c: rates(h, v) for c, v in cols.items()}
        P(f"  {lab}")
        head = f"    {'N':>4} {'h':>7} {'H1':>10} {'rt':>5} {'energy':>10} {'rt':>5} {'slip':>10} {'rt':>5} {'bdn':>10} {'rt':>5}"
        ratio = (kind.startswith("P") and th == "0") or (kind.startswith("B") and th == "D")
        if ratio:
            head += f" {'H1*mu/hG':>9} {'slip*mu/hG':>10} {'en*(mu/h)^.5/G':>14}"
        if th in ("0", "1"):
            head += f" {'divres':>8} {'relres':>8}"
        P(head)
        for i, r in enumerate(R):
            s = f"    {r['N']:>4} {h[i]:7.4f}"
            for c in ("H1", "energy", "bL2", "bdn"):
                s += f" {cols[c][i]:10.3e} {rt[c][i]:5.2f}"
            if ratio:
                s += (f" {cols['H1'][i] * mu / (h[i] * G1):9.4f} {cols['bL2'][i] * mu / (h[i] * G1):10.4f}"
                      f" {cols['energy'][i] * np.sqrt(mu / h[i]) / G1:14.4f}")
            if th in ("0", "1"):
                s += f" {r['divres']:8.1e} {r['relres']:8.1e}"
            P(s)

if T:
    P("\n=== Coercivity threshold mu* on Z_h (sign of pivots of A_mu + r B^T M^-1 B, r = 1e6), L = 2.5, layers = N/4")
    P("    (k = 4 std reference from results/pod: N=16 59.27, N=32 66.02, N=64 70.72; gamma* ~ 18-21)")
    for t in sorted(T, key=lambda t: (t["mesh"], t["k"], t["N"])):
        P(f"    {t['mesh']:>3} k={t['k']} N={t['N']:3d} ntri={t['ntri']:6d} rho={t['rho']:6.3f} "
          f"mu*={t['mu_star']:9.3f} gamma*=mu*/rho={t['gamma_star']:8.3f}")

txt = out.getvalue()
open(os.path.join(RES, "report_v3.txt"), "w").write(txt)
print(txt)
