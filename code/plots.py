"""Figures for the paper: python3 plots.py  (reads ../results/{pod,local}, writes ../paper/figures/*.pdf)"""
import json, glob, os, math
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.abspath(__file__))
POD = os.path.join(ROOT, "..", "results", "pod"); LOC = os.path.join(ROOT, "..", "results", "local")
FIG = os.path.join(ROOT, "..", "paper", "figures"); os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 9, "figure.figsize": (3.3, 2.6), "savefig.bbox": "tight"})

def series(base, kind, mu, th, key="H1"):
    out = []
    for fn in glob.glob(os.path.join(base, "study", f"N*_{kind}_mu{mu:g}.jsonl")):
        for l in open(fn):
            r = json.loads(l)
            if str(r["theta"]) == str(th):
                out.append((r["h"], r[key]))
    return sorted(out)

# Figure 1: H1 error, GS vs CNS*
fig, ax = plt.subplots()
for kind, mu, th, lab, mk in [("B100", 100.0, 0, r"GS, $B_{100}$", "o"), ("B1", 100.0, 0, r"GS, $B_1$", "s"),
                               ("B1", 100.0, 1, r"CNS$^*$, $B_1$ (= $B_{100}$)", "^"), ("A", 100.0, 1, r"CNS$^*$, $A$", "v")]:
    s = series(POD, kind, mu, th); ax.loglog([a for a, _ in s], [b for _, b in s], mk + "-", ms=3, lw=1, label=lab)
h = [0.11, 0.6]
ax.loglog(h, [6 * x for x in h], "k--", lw=0.7); ax.text(0.3, 2.4, "slope 1", fontsize=7)
ax.loglog(h, [0.02 * x ** 1.5 for x in h], "k:", lw=0.7); ax.text(0.3, 0.0022, "slope 3/2", fontsize=7)
ax.set_xlabel(r"$h_\Omega$"); ax.set_ylabel(r"$\|\nabla(\tilde u-u_h)\|$"); ax.legend(fontsize=6, loc="lower right")
ax.set_title(r"$\mu=100$", fontsize=8); fig.savefig(os.path.join(FIG, "h1_rates.pdf"))

# Figure 2: Corollary B' ratios on test P
fig, ax = plt.subplots(); G = math.sqrt(7 * math.pi / 8)
for mu, mk in [(100.0, "o"), (1000.0, "s"), (10000.0, "^")]:
    sl = series(LOC, "P1", mu, 0, "bL2"); en = series(LOC, "P1", mu, 0, "energy")
    ax.plot([a for a, _ in sl], [b * mu / (a * G) for a, b in sl], mk + "-", ms=3, lw=1, label=rf"slip, $\mu={mu:g}$")
    ax.plot([a for a, _ in en], [b * math.sqrt(mu / a) / G for a, b in en], mk + "--", ms=3, lw=1, mfc="none", label=rf"energy, $\mu={mu:g}$")
ax.axhline(1, color="k", lw=0.6); ax.set_xlim(0, 1.6); ax.set_xlabel(r"$h_\Omega$"); ax.set_ylabel("normalised GS error"); ax.legend(fontsize=6)
ax.set_title("test P: slip$\\cdot\\mu/(hG)$, energy$\\cdot(\\mu/h)^{1/2}/G$", fontsize=7); fig.savefig(os.path.join(FIG, "leading_term.pdf"))

# Figure 3: mu* vs rho
fig, ax = plt.subplots(); T = [json.load(open(f)) for f in glob.glob(os.path.join(POD, "thresh", "*.json"))]
for N, mk in [(16, "o"), (32, "s")]:
    t = sorted([r for r in T if r["N"] == N and (r["L"] != 2.5 or r["layers"] == N // 4)], key=lambda r: r["rho"])
    ax.loglog([r["rho"] for r in t], [r["mu_star"] for r in t], mk + "-", ms=3, lw=1, label=f"N = {N}")
rr = [3, 45]; ax.loglog(rr, [20 * x for x in rr], "k--", lw=0.7, label=r"$\mu=20\rho$")
ax.set_xlabel(r"$\rho$"); ax.set_ylabel(r"$\mu^*$"); ax.legend(fontsize=7); fig.savefig(os.path.join(FIG, "threshold.pdf"))
print("wrote", sorted(os.listdir(FIG)))
