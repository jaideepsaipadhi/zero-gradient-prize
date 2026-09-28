"""Build tab_infsup.tex and tab_filter.tex from the logs.  Usage: python3 make_tables.py"""
import json, glob, math
from infsup import MESHES, min_angle
def rd(pat):
    out = []
    for f in sorted(glob.glob(pat)):
        for l in open(f):
            try: out.append(json.loads(l))
            except Exception: pass
    return out
inf = rd("infsup_lobpcg_*.log") + rd("infsup_dense*.log")   # dense entries win
seen = {}
for r in inf:
    seen.setdefault((r["mesh"], r["N"]), []).append(r)
order = {"std": 0, "one": 1, "alt": 2, "fan": 3}
rows = []
for (m, N), rs in sorted(seen.items(), key=lambda kv: (order[kv[0][0]], kv[0][1])):
    r = rs[-1]
    alt = "" if len(rs) == 1 else " (%s)" % ",".join("%.4f" % x["beta_full"] for x in rs)
    phi = r.get("phi_min")
    rows.append("%s & %d & %.4f & %d & %s & %.4f%s & %s & %.4f & %.4f & %.1f\\\\" % (
        m, N, r["hG"], r["nlocked"], "--" if phi is None else "%.4f" % phi, r["beta_full"], alt,
        "--" if phi is None else "%.3f" % (r["beta_full"] / phi), r["beta_lockw"], r["beta_locks"], min_angle(*MESHES[m](N)[:2])))
T = ["\\begin{table}[h]\\centering\\small", "\\caption{Inf-sup constants ($k=4$). $\\beta_{\\rm full}$ on $\\Pi_h\\cap L^2_0$, "
     "$\\beta_{\\Ql}$ on $\\Ql$ (Theorem~\\ref{thm:L}), $\\beta_{(L)}$ on the (L)-space. Values in parentheses: repeated "
     "computations (dense / LOBPCG).}\\label{tab:infsup}",
     "\\begin{tabular}{lrrrrlrrrr}\\toprule",
     "mesh & $N$ & $h_\\Gamma$ & $\\#\\mathcal L$ & $\\phi_{\\min}$ & $\\beta_{\\rm full}$ & $\\beta_{\\rm full}/\\phi_{\\min}$ & $\\beta_{\\Ql}$ & $\\beta_{(L)}$ & min angle\\\\\\midrule"]
T += rows + ["\\bottomrule\\end{tabular}\\end{table}"]
open("tab_infsup.tex", "w").write("\n".join(T) + "\n")
fl = rd("filter.log") + rd("filter_fan.log")
rows = []; prev = {}
for r in fl:
    m = r["mesh"]
    def rate(key):
        if m not in prev: return "--"
        p = prev[m]; return "%.2f" % (math.log(p[key] / r[key]) / math.log(p["hG"] / r["hG"]))
    rows.append("%s & %d & %.4f & %d & %.3e & %s & %.3f & %s & %.4f & %s & %.1f & %.2f & %.3f\\\\" % (
        m, r["N"], r["hG"], r["nlocked"], r["H1"], rate("H1"), r["pL2"], rate("pL2"), r["pL2_f"], rate("pL2_f"),
        r["pLinf_layer"], r["pLinf_layer_f"], r["perp_norm"]))
    prev[m] = r
T = ["\\begin{table}[h]\\centering\\small", "\\caption{Method (S), test case A, $\\mu=100$, $k=4$. $\\|\\pi_h\\|$: raw pressure error; "
     "$\\|\\pi_h^f\\|$: filtered pressure error (Theorem~\\ref{thm:P}(iii)); $L^\\infty$ errors on the boundary layer; last column "
     "$\\|(I-P_{\\rm lock})p_h\\|$.}\\label{tab:filter}",
     "\\begin{tabular}{lrrrrrrrrrrrr}\\toprule",
     "mesh & $N$ & $h_\\Gamma$ & $\\#\\mathcal L$ & $H^1$ err & rate & $\\|\\pi_h\\|$ & rate & $\\|\\pi_h^f\\|$ & rate & $L^\\infty$ & $L^\\infty_f$ & perp\\\\\\midrule"]
T += rows + ["\\bottomrule\\end{tabular}\\end{table}"]
open("tab_filter.tex", "w").write("\n".join(T) + "\n")
print(open("tab_infsup.tex").read()); print(open("tab_filter.tex").read())
