"""v5: tables for the Applied Mathematics Letters version.

  python3 report_v5.py thresh    # (once) recompute h/h_Gamma and true rho = h/min|e| for every
                                 #   results/pod/thresh run  ->  results/v5/thresh_rho.jsonl
  python3 report_v5.py           # regenerate results/v5/letter_tables.txt from the JSONL only

Inputs (JSONL):  results/pod/study/N*_{A,B1,B100}_mu100.jsonl   (H1 errors, GS/CNS*; no rerun)
                 results/local/study/N*_C_mu100.jsonl            (velocity cross-check for test C)
                 results/v5/pressure/*.jsonl                      (code/pressure_err.py)
                 results/v5/strong/*.jsonl                        (code/strong_bc.py)
                 results/v5/thresh_rho.jsonl                      (this script, mode 'thresh')
Rates are log(e_{i-1}/e_i)/log(h_{i-1}/h_i) with h = h_Omega (the h of mu/h).
"""
import os, sys, json, glob, math
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(ROOT, "..", "results")
V5 = os.path.join(RES, "v5")


def load(pattern):
    out = []
    for fn in sorted(glob.glob(pattern)):
        with open(fn) as f:
            out += [json.loads(l) for l in f if l.strip()]
    return out


def build_thresh():
    import svn
    recs = []
    for fn in sorted(glob.glob(os.path.join(RES, "pod", "thresh", "*.json"))):
        r = json.load(open(fn))
        v, t, c, o = svn.make_mesh(r["N"], L=r["L"], layers=r["layers"])
        S = svn.Space(v, t, c, o)
        el = np.array([np.linalg.norm(v[t[k, a]] - v[t[k, b]]) for (k, a, b) in S.bedges])
        assert abs(S.h - r["hOmega"]) < 1e-12 and abs(S.hGamma - r["hG"]) < 1e-12 and len(t) == r["ntri"]
        recs.append(dict(N=r["N"], L=r["L"], layers=r["layers"], ntri=r["ntri"], h=S.h, hG=S.hGamma,
                         emin=float(el.min()), emax=float(el.max()), h_over_hG=S.h / S.hGamma,
                         rho_true=S.h / float(el.min()), mu_star=r["mu_star"],
                         gamma_star=r["mu_star"] * S.hGamma / S.h,
                         mu_star_over_rho=r["mu_star"] / (S.h / float(el.min())), src=os.path.relpath(fn, RES)))
    recs.sort(key=lambda r: (r["L"] != 2.5, r["N"], r["L"]))
    with open(os.path.join(V5, "thresh_rho.jsonl"), "w") as f:
        for r in recs:
            f.write(json.dumps(r) + "\n")


def rates(rows, key):
    out = [None]
    for a, b in zip(rows[:-1], rows[1:]):
        ea, eb = a.get(key), b.get(key)
        out.append(math.log(ea / eb) / math.log(a["h"] / b["h"]) if ea and eb and ea > 0 and eb > 0 else None)
    return out


def by_N(recs, **sel):
    rs = [r for r in recs if all(r.get(k) == v for k, v in sel.items())]
    d = {}
    for r in rs:
        d[r["N"]] = r
    return d


def e(x):
    if x is None:
        return "--"
    m, ex = f"{x:.2e}".split("e")
    return f"{m}e{int(ex)}"


def rt(x):
    return "--" if x is None else f"{x:.2f}"


def ltx_e(x):
    if x is None:
        return "--"
    m, ex = f"{x:.2e}".split("e")
    return f"${m}\\mathrm{{e}}{{{int(ex)}}}$"


def ltx_r(x):
    return "--" if x is None else f"{x:.2f}"


def series(d, Ns, key):
    rows = [d[N] if N in d else None for N in Ns]
    vals, rts, prev = [], [], None
    for r in rows:
        if r is None:
            vals.append(None); rts.append(None); continue
        v = r.get(key)
        vals.append(v)
        rts.append(math.log(prev[1] / v) / math.log(prev[0] / r["h"]) if prev and v and prev[1] else None)
        prev = (r["h"], v)
    return vals, rts


def plain_name(x):
    for a, b in (("$", ""), ("^*", "*"), ("H^1", "H1"), ("L^2", "L2"), ("\\Gamma", "Gamma"), ("|\\cdot|_1", "H1")):
        x = x.replace(a, b)
    return x


def table(title, Ns, hs, cols, latex_caption, label):
    """cols: list of (LaTeX header, vals, rates) ; returns (plain, latex)"""
    head = f"{'N':>5} {'h':>7} " + " ".join(f"{plain_name(c[0]):>18}" for c in cols)
    lines = [title, head]
    for i, N in enumerate(Ns):
        lines.append(f"{N:>5} {hs[i]:>7.4f} " + " ".join(f"{e(c[1][i]):>11} {rt(c[2][i]):>6}" for c in cols))
    L = ["\\begin{table}[t]", "\\centering", f"\\caption{{{latex_caption}}}\\label{{{label}}}", "\\small",
         "\\begin{tabular}{rr" + "rr" * len(cols) + "}", "\\toprule",
         "$N$ & $h$ & " + " & ".join(f"\\multicolumn{{2}}{{c}}{{{c[0]}}}" for c in cols) + "\\\\",
         " ".join(f"\\cmidrule(lr){{{3 + 2 * j}-{4 + 2 * j}}}" for j in range(len(cols))),
         " & & " + " & ".join("error & rate" for _ in cols) + "\\\\", "\\midrule"]
    for i, N in enumerate(Ns):
        L.append(f"{N} & {hs[i]:.4f} & " + " & ".join(f"{ltx_e(c[1][i])} & {ltx_r(c[2][i])}" for c in cols) + "\\\\")
    L += ["\\bottomrule", "\\end{tabular}", "\\end{table}"]
    return "\n".join(lines), "\n".join(L)


def main():
    pod = load(os.path.join(RES, "pod", "study", "N*_mu100.jsonl"))
    loc = load(os.path.join(RES, "local", "study", "N*_C_mu100.jsonl"))
    prs = load(os.path.join(V5, "pressure", "*.jsonl"))
    stg = load(os.path.join(V5, "strong", "*.jsonl"))
    for r in stg:
        r.setdefault("mesh", "std"); r.setdefault("solver", "lu")
    thr = load(os.path.join(V5, "thresh_rho.jsonl"))
    out_plain, out_ltx, notes = [], [], []

    # ------------------------------------------------ health
    allrec = [r for r in prs + stg]
    bad = [r for r in allrec if r.get("relres", 0) > 1e-8 or r.get("divres", 0) > 1e-6]
    podm = [r for r in pod if r["theta"] in (0, 1) and r["kind"] in ("A", "B100")]
    badp = [r for r in podm if r.get("relres", 0) > 1e-8 or r.get("divres", 0) > 1e-6]
    # velocity cross-check v5 pressure runs vs stored pod/local runs
    ref = {(r["N"], r["kind"], r["theta"]): r for r in pod + loc if r["theta"] in (0, 1)}
    diffs = []
    for r in prs:
        k = (r["N"], r["kind"], r["theta"])
        if k in ref:
            diffs.append((abs(r["H1"] - ref[k]["H1"]) / ref[k]["H1"], k))
    md = max(diffs) if diffs else (float("nan"), None)
    health = [f"solver health: v5 records {len(allrec)} (pressure {len(prs)}, strong {len(stg)}); "
              f"flagged (relres>1e-8 or divres>1e-6): {len(bad)}",
              f"  v5 max relres {max((r['relres'] for r in allrec), default=float('nan')):.2e}, "
              f"max divres {max((r['divres'] for r in allrec), default=float('nan')):.2e}, "
              f"max RSS {max((r.get('maxrss_gb', 0) for r in allrec), default=0):.2f} GB, "
              f"solvers used: {sorted({r.get('solver', 'lu') for r in allrec})}",
              f"  pod A/B100 mu=100 records used: {len(podm)}, flagged: {len(badp)}; "
              f"max relres {max(r['relres'] for r in podm):.2e}, max divres {max(r['divres'] for r in podm):.2e}",
              f"  velocity H1 of the v5 reruns vs stored pod/local JSONL: max rel. diff {md[0]:.1e} at {md[1]} "
              f"({len(diffs)} pairs)"]

    # ------------------------------------------------ T1: test A, strong vs GS vs CNS*
    NsA = [16, 24, 32, 48, 64, 96, 128, 192, 256]
    gsA = by_N(pod, kind="A", theta=0); cnA = by_N(pod, kind="A", theta=1)
    stA = by_N(stg, kind="A", method="strong-zero", mesh="std"); siA = by_N(stg, kind="A", method="strong-interp", mesh="std")
    saA = by_N(stg, kind="A", method="strong-zero", mesh="alt")
    NsA = [N for N in NsA if N in gsA and N in cnA]
    hs = [gsA[N]["h"] for N in NsA]
    cols = [("strong, mesh (a)", *series(saA, NsA, "H1")), ("strong, mesh (s)", *series(stA, NsA, "H1")),
            ("GS", *series(gsA, NsA, "H1")), ("CNS$^*$", *series(cnA, NsA, "H1"))]
    p, l = table("T1  test A (Scott's shear flow, p = 0), mu = 100: |u - u_h|_{1,Omega_h}",
                 NsA, hs, cols,
                 "Test A ($u=(1-r^{-2})(-y,x)$, $p\\equiv0$), $k=4$, $\\mu=100$: $H^1$ seminorm error "
                 "$|\\tilde u-u_h|_{1,\\Omega_h}$. Strong: $u_h=0$ at all nodes of $\\Gamma_h$, no Nitsche terms, "
                 "on (s) the mesh of GS/CNS$^*$ (every polygon vertex in 3 triangles) and (a) the same mesh with "
                 "alternating first-layer diagonals (polygon vertices in 2 or 4 triangles).",
                 "tab:A")
    out_plain.append(p); out_ltx.append(l)
    # supplement: strong with interpolated data, energy norms, vertex gradient defect
    sup = ["T1'  strong-imposition diagnostics (test A): H1 of the strong runs; 'interp' = u_h = u~ at the Gamma_h nodes;",
           "     vgrad = max_{z vertex of Gamma_h, T ∋ z} |grad u_h|_T(z) - grad u(z)|;  pressure ||p_h||_{L2} (p = 0)",
           f"{'N':>5} {'h':>7} {'zero H1':>10} {'rt':>5} {'vgrad':>8} {'||p_h||':>9} {'rt':>5}   {'interp H1':>10} {'rt':>5} {'vgrad':>8}"]
    Ns_s = sorted(set(stA) | set(siA))
    zH, zR = series(stA, Ns_s, "H1"); iH, iR = series(siA, Ns_s, "H1"); zP, zPR = series(stA, Ns_s, "pL2opt")
    for i, N in enumerate(Ns_s):
        r = stA.get(N) or siA.get(N)
        sup.append(f"{N:>5} {r['h']:>7.4f} {e(zH[i]):>10} {rt(zR[i]):>5} {stA[N]['vgrad_err'] if N in stA else float('nan'):>8.3f} "
                   f"{e(zP[i]):>9} {rt(zPR[i]):>5}   {e(iH[i]):>10} {rt(iR[i]):>5} "
                   f"{siA[N]['vgrad_err'] if N in siA else float('nan'):>8.1e}")
    sup.append("  mesh (a) = alternating first-layer diagonals: strong u_h = 0, and GS / CNS* on the same mesh")
    sup.append(f"{'N':>5} {'h':>7} {'strong H1':>10} {'rt':>5} {'vgrad':>8} {'||p_h||':>9} {'rt':>5}   "
               f"{'GS H1':>10} {'rt':>5} {'CNS* H1':>10} {'rt':>5} {'GS ||p||':>9} {'CNS* ||p||':>10}")
    gaA = by_N(stg, kind="A", method="GS", mesh="alt"); caA = by_N(stg, kind="A", method="CNS*", mesh="alt")
    Ns_a = sorted(set(saA) | set(gaA))
    aH, aR = series(saA, Ns_a, "H1"); aP, aPR = series(saA, Ns_a, "pL2opt")
    gH, gR = series(gaA, Ns_a, "H1"); cH, cR = series(caA, Ns_a, "H1")
    for i, N in enumerate(Ns_a):
        r = saA.get(N) or gaA.get(N)
        sup.append(f"{N:>5} {r['h']:>7.4f} {e(aH[i]):>10} {rt(aR[i]):>5} "
                   f"{saA[N]['vgrad_err'] if N in saA else float('nan'):>8.3f} {e(aP[i]):>9} {rt(aPR[i]):>5}   "
                   f"{e(gH[i]):>10} {rt(gR[i]):>5} {e(cH[i]):>10} {rt(cR[i]):>5} "
                   f"{e(gaA[N]['pL2'] if N in gaA else None):>9} {e(caA[N]['pL2'] if N in caA else None):>10}")
    out_plain.append("\n".join(sup))

    # ------------------------------------------------ T2: B100, H1 + pressure
    NsB = [16, 24, 32, 48, 64, 96, 128, 192]
    gsB = by_N(pod, kind="B100", theta=0); cnB = by_N(pod, kind="B100", theta=1)
    NsB = [N for N in NsB if N in gsB and N in cnB]
    hs = [gsB[N]["h"] for N in NsB]
    pgB = by_N(prs, kind="B100", theta=0); pcB = by_N(prs, kind="B100", theta=1)
    cols = [("GS, $|\\cdot|_1$", *series(gsB, NsB, "H1")), ("CNS$^*$, $|\\cdot|_1$", *series(cnB, NsB, "H1")),
            ("GS, $L^2$ pressure", *series(pgB, NsB, "pL2")), ("CNS$^*$, $L^2$ pressure", *series(pcB, NsB, "pL2"))]
    p, l = table("T2  test B100 (p = 100(x^2y - y + x/2)), mu = 100: |u - u_h|_1 and ||p* - p_h||_{L2(Omega_h)}, "
                 "p* = p~ - mean_{Gamma_h} p~",
                 NsB, hs, cols,
                 "Test B$_{100}$ (large pressure variation, $p=100(x^2y-y+x/2)$), $k=4$, $\\mu=100$: velocity error "
                 "$|\\tilde u-u_h|_{1,\\Omega_h}$ and pressure error $\\|p^\\ast-p_h\\|_{L^2(\\Omega_h)}$, "
                 "$p^\\ast=\\tilde p-\\bar p_{\\Gamma_h}(\\tilde p)$.", "tab:B100")
    out_plain.append(p); out_ltx.append(l)

    # ------------------------------------------------ T3: penalty thresholds
    ref_ = [r for r in thr if r["L"] == 2.5]
    swp = [r for r in thr if r["N"] in (16, 32)]
    P = ["T3  coercivity threshold mu* (results/pod/thresh), gamma* = mu* h_Gamma/h; h/h_Gamma = h/max|e|, "
         "rho = h/min|e| (e in E_Gamma)",
         "(a) refinement at L = 2.5 (fixed h/h_Gamma ~ 3.3)",
         f"{'N':>5} {'h/hG':>7} {'rho':>7} {'max/min|e|':>10} {'mu*':>8} {'gamma*':>7} {'mu*/rho':>8}"]
    for r in ref_:
        P.append(f"{r['N']:>5} {r['h_over_hG']:>7.3f} {r['rho_true']:>7.3f} {r['emax']/r['emin']:>10.3f} "
                 f"{r['mu_star']:>8.2f} {r['gamma_star']:>7.2f} {r['mu_star_over_rho']:>8.2f}")
    P.append("(b) rho-sweep through the outer box (-L,L)^2")
    P.append(f"{'N':>5} {'L':>5} {'layers':>6} {'h/hG':>7} {'rho':>7} {'mu*':>8} {'gamma*':>7} {'mu*/rho':>8}")
    for N in (16, 32):
        for r in sorted([r for r in swp if r["N"] == N], key=lambda r: r["L"]):
            P.append(f"{N:>5} {r['L']:>5g} {r['layers']:>6} {r['h_over_hG']:>7.3f} {r['rho_true']:>7.3f} "
                     f"{r['mu_star']:>8.2f} {r['gamma_star']:>7.2f} {r['mu_star_over_rho']:>8.2f}")
    out_plain.append("\n".join(P))
    L = ["\\begin{table}[t]", "\\centering",
         "\\caption{Coercivity threshold $\\mu^\\ast$ of the Nitsche form on $Z_h$ ($k=4$) and "
         "$\\gamma^\\ast=\\mu^\\ast\\hG/h$. $h/\\hG=h/\\max_e|e|$; $\\rho=h/\\min_e|e|$, $e\\in\\mathcal E_\\Gamma$. "
         "(a) refinement in the box $L=2.5$; (b) $h$ varied through the box size $L$.}\\label{tab:thresh}", "\\small",
         "\\textit{(a)}\\\\[2pt]",
         "\\begin{tabular}{l" + "c" * len(ref_) + "}", "\\toprule",
         "$N$ & " + " & ".join(str(r["N"]) for r in ref_) + "\\\\", "\\midrule",
         "$h/\\hG$ & " + " & ".join(f"{r['h_over_hG']:.2f}" for r in ref_) + "\\\\",
         "$\\rho$ & " + " & ".join(f"{r['rho_true']:.2f}" for r in ref_) + "\\\\",
         "$\\mu^\\ast$ & " + " & ".join(f"{r['mu_star']:.1f}" for r in ref_) + "\\\\",
         "$\\gamma^\\ast$ & " + " & ".join(f"{r['gamma_star']:.1f}" for r in ref_) + "\\\\",
         "\\bottomrule", "\\end{tabular}\\\\[6pt]", "\\textit{(b)}\\\\[2pt]"]
    Ls = sorted({r["L"] for r in swp})
    L += ["\\begin{tabular}{ll" + "c" * len(Ls) + "}", "\\toprule",
          " & $L$ & " + " & ".join(f"{x:g}" for x in Ls) + "\\\\"]
    for N in (16, 32):
        rr = {r["L"]: r for r in swp if r["N"] == N}
        L += ["\\midrule",
              f"$N={N}$ & $h/\\hG$ & " + " & ".join(f"{rr[x]['h_over_hG']:.2f}" for x in Ls) + "\\\\",
              " & $\\rho$ & " + " & ".join(f"{rr[x]['rho_true']:.2f}" for x in Ls) + "\\\\",
              " & $\\mu^\\ast$ & " + " & ".join(f"{rr[x]['mu_star']:.1f}" for x in Ls) + "\\\\",
              " & $\\gamma^\\ast$ & " + " & ".join(f"{rr[x]['gamma_star']:.1f}" for x in Ls) + "\\\\"]
    L += ["\\bottomrule", "\\end{tabular}", "\\end{table}"]
    out_ltx.append("\n".join(L))

    # ------------------------------------------------ T4: pressure errors, all tests (supplement)
    NsP = sorted({r["N"] for r in prs})
    for kind in ("A", "B1", "B100", "C"):
        for th, name in ((0, "GS"), (1, "CNS*")):
            d = by_N(prs, kind=kind, theta=th)
            if not d:
                continue
            Ns = [N for N in NsP if N in d]
            keys = [("pL2", "L2 (pbar_G)"), ("pL2opt", "L2 opt.const"), ("pL2G", "L2 both G-norm."),
                    ("pL2_layer", "first layer"), ("pL2_off", "off layer"), ("pLinf_layer", "Linf layer")]
            lines = [f"T4  pressure, test {kind}, {name}, mu = 100 "
                     f"(lam_h = Gamma_h-mean of p_h; err_mean = Omega_h-mean of p* - p_h)",
                     f"{'N':>5} {'h':>7} " + " ".join(f"{k[1]:>17}" for k in keys) + f" {'lam_h':>10} {'err_mean':>10}"]
            S = {k[0]: series(d, Ns, k[0]) for k in keys}
            for i, N in enumerate(Ns):
                lines.append(f"{N:>5} {d[N]['h']:>7.4f} " + " ".join(f"{e(S[k[0]][0][i]):>10} {rt(S[k[0]][1][i]):>6}" for k in keys)
                             + f" {d[N]['lam_h']:>10.2e} {d[N]['err_mean']:>10.2e}")
            out_plain.append("\n".join(lines))
    # LaTeX for the pressure supplement (B1 and C, GS vs CNS*, pL2)
    for kind, cap in (("B1", "Test B$_1$"), ("C", "Test C (non-polynomial pressure $e^{x/2}\\cos y$)"), ("A", "Test A ($p\\equiv0$)")):
        dg = by_N(prs, kind=kind, theta=0); dc = by_N(prs, kind=kind, theta=1)
        Ns = [N for N in NsP if N in dg and N in dc]
        if not Ns:
            continue
        cols = [("GS", *series(dg, Ns, "pL2")), ("GS, off first layer", *series(dg, Ns, "pL2_off")),
                ("CNS$^*$", *series(dc, Ns, "pL2"))]
        _, l = table("", Ns, [dg[N]["h"] for N in Ns], cols,
                     f"{cap}, $\\mu=100$: pressure error $\\|p^\\ast-p_h\\|_{{L^2(\\Omega_h)}}$ "
                     "($p^\\ast=\\tilde p-\\bar p_{\\Gamma_h}(\\tilde p)$); the off-layer column removes the optimal constant.",
                     f"tab:p{kind}")
        out_ltx.append(l)

    txt = ["PPL 115 / letter numerics (v5).  k = 4, Scott-Vogelius P4/P3disc, meshes of svn.make_mesh (L = 2.5).",
           "Generated by code/report_v5.py from JSONL; do not edit by hand.", ""] + health + [""]
    txt += ["=" * 100, "PLAIN TEXT", "=" * 100]
    for p in out_plain:
        txt += [p, ""]
    txt += ["=" * 100, "LaTeX (booktabs; needs \\usepackage{booktabs,amsmath}; \\hG = h_\\Gamma macro of the paper)", "=" * 100]
    for l in out_ltx:
        txt += [l, ""]
    with open(os.path.join(V5, "letter_tables.txt"), "w") as f:
        f.write("\n".join(txt))
    print("\n".join(txt[:6 + len(health)]))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "thresh":
        sys.path.insert(0, ROOT)
        build_thresh()
    else:
        main()
