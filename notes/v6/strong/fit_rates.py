"""T3: asymptotic rate of strong imposition on mesh (s).  Reads results/v5/strong (N<=128) and
notes/v6/strong/runs/strong (N=160,192), plus GS/CNS* from results/pod/study (N<=256)."""
import os, json, glob, numpy as np
ROOT = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.join(ROOT, "..", "..", "..")
def load(pattern):
    d = {}
    for f in glob.glob(pattern):
        for l in open(f):
            r = json.loads(l); d[r["N"]] = r
    return d
S = load(os.path.join(REPO, "results/v5/strong/N*_A_zero.jsonl")); S.update(load(os.path.join(ROOT, "runs/strong/N*_A_zero.jsonl")))
ONE = load(os.path.join(ROOT, "runs/strong/N*_A_zero_one.jsonl"))
ALT = load(os.path.join(REPO, "results/v5/strong/N*_A_zero_alt.jsonl"))
GS = {}; CN = {}
for f in glob.glob(os.path.join(REPO, "results/pod/study/N*_A_mu100.jsonl")):
    for l in open(f):
        r = json.loads(l)
        if r.get("theta") == 0: GS[r["N"]] = r
        if r.get("theta") == 1: CN[r["N"]] = r
out = []
def P(s=""): print(s); out.append(s)
def table(name, D, key="H1"):
    Ns = sorted(D); P(f"== {name}: N, hG, {key}, rate(hG), rate(h)")
    for i, N in enumerate(Ns):
        r = D[N]; s = f"{N:4d} {r['hG']:.5f} {r[key]:.5e}"
        if i:
            q = D[Ns[i - 1]]
            s += f"  {np.log(q[key] / r[key]) / np.log(q['hG'] / r['hG']):.3f}  {np.log(q[key] / r[key]) / np.log(q['h'] / r['h']):.3f}"
        P(s)
table("strong, mesh (s)", S); table("strong, mesh (s), vertex gradient error", S, "vgrad_err")
table("strong, mesh (s), ||p_h|| (p = 0)", S, "ph_L2")
table("GS, mesh (s)", GS); table("CNS*, mesh (s)", CN)
if ONE: table("strong, one flipped quad (one locked vertex)", ONE); table("  its vgrad", ONE, "vgrad_err")
if ALT: table("strong, mesh (a)", ALT)
P("== ratio strong/GS on mesh (s)")
for N in sorted(S):
    if N in GS: P(f"{N:4d} {S[N]['H1'] / GS[N]['H1']:.4f}")
# least-squares fits on N >= 32
Ns = [N for N in sorted(S) if N >= 32]
h = np.array([S[N]["hG"] for N in Ns]); E = np.array([S[N]["H1"] for N in Ns])
P("== fits of strong, mesh (s), N>=32, weighted by 1/E (relative residuals)")
for name, pw in [("a h^1.5 + b h^2", (1.5, 2)), ("a h^1.5 + b h^2.5", (1.5, 2.5)), ("a h^1.5 + b h^2 + c h^2.5", (1.5, 2, 2.5)),
                 ("a h^1 + b h^1.5", (1, 1.5)), ("a h^1 + b h^2", (1, 2)), ("a h^2 + b h^2.5", (2, 2.5))]:
    A = np.stack([h ** p for p in pw], 1) / E[:, None]
    c, *_ = np.linalg.lstsq(A, np.ones_like(E), rcond=None)
    rel = A @ c - 1
    P(f"  {name:28s} coef={np.array2string(c, precision=4)}  max rel resid={np.abs(rel).max():.2e}")
P("== E/hG^1.5 and first differences against hG (model E = a hG^1.5 (1 + beta hG))")
prev = None
for N in sorted(S):
    q = S[N]["H1"] / S[N]["hG"] ** 1.5
    s_ = f"{N:4d} E/hG^1.5 = {q:.5f}"
    if prev: s_ += f"   d(E/hG^1.5)/d(hG) = {(prev[1] - q) / (prev[0] - S[N]['hG']):.4f}"
    P(s_); prev = (S[N]["hG"], q)
P("  predicted local rate 1.5 + beta*hG/(1+beta*hG) with beta = b/a of the 2-term fit a h^1.5 + b h^2.5:")
A_ = np.stack([h ** 1.5, h ** 2.5], 1) / E[:, None]; c_, *_ = np.linalg.lstsq(A_, np.ones_like(E), rcond=None); beta = c_[1] / c_[0]
Ns_all = sorted(S)
for i in range(1, len(Ns_all)):
    h1, h2 = S[Ns_all[i - 1]]["hG"], S[Ns_all[i]]["hG"]
    pred = 1.5 + np.log((1 + beta * h1) / (1 + beta * h2)) / np.log(h1 / h2)
    obs = np.log(S[Ns_all[i - 1]]["H1"] / S[Ns_all[i]]["H1"]) / np.log(h1 / h2)
    P(f"  {Ns_all[i - 1]:4d}->{Ns_all[i]:4d}: predicted {pred:.3f}  observed {obs:.3f}")
for Nn in (256, 384, 512, 1024):
    hh_ = S[128]["hG"] * 128 / Nn
    P(f"  extrapolated local rate near N={Nn}: {1.5 + beta * hh_ / (1 + beta * hh_):.3f}")
if ONE:
    P("== one locked vertex: locked part L = sqrt(E_one^2 - E_std^2) and its rate")
    prev = None
    for N in sorted(ONE):
        if N in S:
            L_ = np.sqrt(max(ONE[N]["H1"] ** 2 - S[N]["H1"] ** 2, 0)); s_ = f"{N:4d} L = {L_:.5e}  L/hG = {L_ / S[N]['hG']:.4f}"
            if prev: s_ += f"  rate {np.log(prev[1] / L_) / np.log(prev[0] / S[N]['hG']):.3f}"
            P(s_); prev = (S[N]["hG"], L_)
P("== solver health of the v6 runs")
for D, nm in ((S, "strong"), (ONE, "one")):
    for N in sorted(D):
        r = D[N]
        if "divres" in r: P(f"  {nm} N={N}: relres={r['relres']:.1e} divres={r['divres']:.1e} flux={r.get('flux', float('nan')):.1e} secs={r['secs']:.0f} rss={r['maxrss_gb']:.2f}GB")
P("== free-exponent fit E = a h^p (1 + b h) over 4 consecutive levels (grid on p)")
for lo in range(0, len(Ns) - 3):
    hh, EE = h[lo:lo + 4], E[lo:lo + 4]; best = None
    for p in np.linspace(0.8, 2.5, 1701):
        A = np.stack([hh ** p, hh ** (p + 1)], 1) / EE[:, None]
        c, *_ = np.linalg.lstsq(A, np.ones_like(EE), rcond=None)
        r = np.linalg.norm(A @ c - 1)
        if best is None or r < best[0]: best = (r, p, c)
    P(f"  N in {Ns[lo:lo + 4]}: p = {best[1]:.3f}, a = {best[2][0]:.4f}, b = {best[2][1]:.4f}, resid = {best[0]:.1e}")
open(os.path.join(ROOT, "fit_rates.log"), "w").write("\n".join(out) + "\n")
