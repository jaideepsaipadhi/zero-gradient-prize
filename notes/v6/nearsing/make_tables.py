"""Build tab_fe.tex (errors, rates, lower-bound quantities) and tab_infsup.tex from runs/ and ns_bound.log."""
import os, json, glob, math
HERE = os.path.dirname(os.path.abspath(__file__))
B = {}
for l in open(os.path.join(HERE, "ns_bound.log")):
    r = json.loads(l); B[(r["mesh"], r["mode"], r["N"])] = r
PRED = {("alt", "c0"): "1/2", ("altsplit", "c0.3"): "3/2", ("altsplit", "h1.0"): "1", ("altsplit", "q1.0"): "1/2",
        ("altwall", "c0.3"): "3/2", ("altwall", "h1.0"): "1", ("altwall", "q1.0"): "1/2",
        ("onesplit", "q1.0"): "1", ("onesplit", "h1.0"): "3/2", ("onewall", "q1.0"): "1", ("std", "c0"): "3/2"}
EPS = {"c0.3": r"$0.3$", "h1.0": r"$\hG$", "q1.0": r"$\hG^2$", "c0": "--"}
rows = []
for (mesh, mode), pred in PRED.items():
    f = os.path.join(HERE, "runs", f"{mesh}_{mode}.jsonl")
    if not os.path.exists(f):
        continue
    R = sorted([json.loads(l) for l in open(f) if l.startswith("{")], key=lambda r: r["N"])
    Es, rates = [], []
    for i, r in enumerate(R):
        rate = "" if i == 0 else "%.2f" % (math.log(R[i - 1]["H1"] / r["H1"]) / math.log(R[i - 1]["hG"] / r["hG"]))
        b = B.get((mesh, mode, r["N"]))
        extra = ("%.3f & %.2f & %.2f" % (b["gam_sqrtD"], b["D_over_ind"], b["E_over_gamD"])) if b else "-- & -- & --"
        rows.append(f"{mesh} & {EPS[mode]} & {r['N']} & {r['H1']:.4f} & {rate} & {extra} & {r['divres']:.0e} & {pred if i == len(R) - 1 else ''}\\\\")
    rows.append(r"\midrule")
with open(os.path.join(HERE, "tab_fe.tex"), "w") as fh:
    fh.write(r"""\begin{table}[ht]\centering\small
\begin{tabular}{llrrrrrrrl}\toprule
mesh & $\varepsilon$ & $N$ & $\norm{\nabla(\tilde u-u_h)}$ & rate & $\gamma_4(\sum d_z^2)^{1/2}$ & $\frac{(\sum d_z^2)^{1/2}}{\mathcal I}$ & $\frac{E}{\gamma_4(\sum d_z^2)^{1/2}}$ & divres & pred.\\\midrule
""" + "\n".join(rows[:-1]) + r"""
\bottomrule\end{tabular}
\caption{Method (S), $k=4$, test A. $\mathcal I=(\sum_{z}h_z^2|\nabla u(z)|^2\min(1,\phi_z/\sqrt{\varepsilon_z})^2)^{1/2}$ (Theorem~\ref{thm:N2} indicator);
$d_z$ computed exactly on the actual stars (\texttt{ns\_bound.py}); ``pred.'' = limiting rate predicted by Theorem~\ref{thm:N3}.}\label{tab:fe}
\end{table}
""")
# inf-sup
rows = []
for f in sorted(glob.glob(os.path.join(HERE, "runs", "infsup_*.jsonl"))):
    for l in open(f):
        if not l.startswith("{"):
            continue
        r = json.loads(l)
        g = lambda k: ("%.4f" % r[k]) if k in r else "--"
        rows.append(f"{r['mesh']} & {r['N']} & {r['eps']:g} & {r['minang_z']:.2f} & {g('beta_full')} & {g('beta_wired')} & {g('gamma_wired')} & "
                    f"{g('gamma_zero3')} & {g('beta_wiredslit')} & {g('gamma_wiredslit')}\\\\")
with open(os.path.join(HERE, "tab_infsup.tex"), "w") as fh:
    fh.write(r"""\begin{table}[ht]\centering\small
\begin{tabular}{lrrrrrrrrr}\toprule
mesh & $N$ & $\varepsilon$ & $\min\angle_z$ ($^\circ$) & $\beta_{\rm full}$ & $\beta_{\rm wired}$ & $\gamma_{\rm wired}$ & $\gamma_{\rm zero3}$ & $\beta_{\rm wired+slit}$ & $\gamma_{\rm wired+slit}$\\\midrule
""" + "\n".join(rows) + r"""
\bottomrule\end{tabular}
\caption{Dense inf-sup constants (\texttt{ns\_infsup.py}); lock-free reference $\beta(\texttt{std})=0.1023$ ($N=16$), $0.0913$ ($N=24$);
plain lock ($\texttt{one}$, $N=16$): $\beta_{\rm full}=0.0583$, $\beta_{\rm wired}=\gamma_{\rm wired}=0.1023$.}\label{tab:infsup}
\end{table}
""")
print("ok")
