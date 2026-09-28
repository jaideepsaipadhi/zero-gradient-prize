"""Tables for the v4 runs:   python3 report_v4.py [RESULTS_DIR]  -> RESULTS_DIR/report_v4.txt (default ../results/v4)
Rates are computed in hG between consecutive refinement levels that were run."""
import os, sys, json, glob
import numpy as np

D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "v4")
out = []


def P(s=""):
    out.append(s)


def load(sub):
    recs = []
    for f in sorted(glob.glob(os.path.join(D, sub, "*.jsonl"))):
        for l in open(f):
            recs.append(json.loads(l))
    return recs


def rate(e, h):
    r = [float("nan")]
    for i in range(1, len(e)):
        r.append(np.log(e[i - 1] / e[i]) / np.log(h[i - 1] / h[i]) if e[i] > 0 and e[i - 1] > 0 else float("nan"))
    return r


def f(x, w=9, p=3):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return " " * (w - 3) + "  -"
    if isinstance(x, str):
        return x.rjust(w)
    if abs(x) != 0 and (abs(x) < 1e-3 or abs(x) >= 1e4):
        return f"{x:{w}.{p - 1}e}"
    return f"{x:{w}.{p + 1}f}" if abs(x) < 10 else f"{x:{w}.{p}f}"


def fr(x):
    return "   -" if x is None or np.isnan(x) else f"{x:5.2f}"


health = {}


def note_health(fam, r):
    if "relres" not in r:
        return
    h = health.setdefault(fam, dict(n=0, relres=0.0, divres=0.0, notconv=0))
    h["n"] += 1; h["relres"] = max(h["relres"], r["relres"]); h["divres"] = max(h["divres"], r["divres"])
    if r.get("converged") is False:
        h["notconv"] += 1


P("v4 numerics: general domains (Section 10), steady Navier-Stokes (Section 11), cylinder drag (Schafer-Turek 2D-1)")
P("k = 4 Scott-Vogelius + Nitsche throughout.  GS = printed method (theta = 0); CNS = CNS* (theta = 1, GLOBAL boundary")
P("mean of p over all Gamma_h components); PER = per-component means; NAIVE = theta = 1 without a mean.")
P("Rates are in hG between consecutive levels.  Ratios: slip = ||e||_{Gamma_h} nu mu/(h G), energy = |||e||| nu (mu/h)^{1/2}/G,")
P("H1 = ||grad e|| nu mu/(h G) (nu = 1 for Stokes); e = u~ - u_h (= -u_h when u = 0), D = u_GS - u_CNS.")
P()

# ------------------------------------------------------------------ thresholds / meshes
P("=" * 110)
P("0. Meshes (mesh_general.py) and coercivity thresholds mu* (penalty_k.threshold on the general mesh, all Sigma Dirichlet)")
P("=" * 110)
th = {}
for fn in glob.glob(os.path.join(D, "thresh", "*.json")):
    r = json.loads(open(fn).read())
    th[(r["geom"], r["N"])] = r
P(f"{'geom':8s} {'N':>4s} {'ntri':>6s} {'h':>8s} {'hG':>8s} {'h/hG':>6s} {'minang':>6s} {'Theta':>6s} {'#fix':>4s} {'mu*':>8s} {'gamma*':>7s}")
for (g, N) in sorted(th):
    r = th[(g, N)]; m = r["mesh"]
    P(f"{g:8s} {N:4d} {r['ntri']:6d} {r['hOmega']:8.4f} {r['hG']:8.4f} {r['h_over_hG']:6.2f} {m['minangle']:6.1f} "
      f"{m['theta_min']:6.3f} {sum(m['fixes'].values()):4d} {r['mu_star']:8.1f} {r['gamma_star']:7.2f}")
P("(#fix = local repairs to enforce (M1),(M1'),(M2); every mesh has >= 3 triangles at each Gamma_h vertex and at each")
P(" straight-side vertex of Sigma, >= 2 at each corner, Theta >= 0.15.)")
P()

# ------------------------------------------------------------------ (A)(i) two disks
gd = load("gd")
for r in gd:
    note_health("general domains (Stokes)", r)
P("=" * 110)
P("A(i). Two disks (Prop. gd:prop:Gglobal): u = 0, p = +1 near disk 1 (r = 0.5), -1 near disk 2 (r = 0.8), f = grad p")
P("      G_per = 0 but G (global mean) = 2.7809.  Prediction for GS: flux_i = int_{Gamma_h,i} u_h.n ~ (h/mu)|Gamma_i|(pbar_i - pbar)")
P("=" * 110)
for mu in sorted(set(r["mu"] for r in gd if r["geom"] == "twodisk")):
    rs = sorted([r for r in gd if r["geom"] == "twodisk" and r["mu"] == mu and r["method"] != "meta"], key=lambda r: r["N"])
    P(f"mu = {mu:g}")
    P(f"{'meth':5s} {'N':>4s} {'hG':>7s} {'h':>7s} {'flux1/pr':>9s} {'flux2/pr':>9s} {'slip':>7s} {'energy':>7s} {'H1rat':>7s} "
      f"{'misfit':>7s} {'|u_h|_1':>9s} {'rH1(h)':>6s} {'rE(h)':>6s} {'relres':>8s} {'divres':>8s}")
    for m in ("GS", "PER", "CNS"):
        rr = [r for r in rs if r["method"] == m]
        if not rr:
            continue
        hG = [r["h"] for r in rr]                         # leak quantities scale with h/mu: rates in h
        rh = rate([r["H1"] for r in rr], hG); re_ = rate([r["energy"] for r in rr], hG)
        for i, r in enumerate(rr):
            fl = r["flux"]; pr = r["flux_pred_GS"]
            if m == "CNS":
                P(f"{m:5s} {r['N']:4d} {r['hG']:7.4f} {r['h']:7.4f} {'':>9s} {'':>9s} {'':>7s} {'':>7s} {'':>7s} {'':>7s} "
                  f"{r['H1']:9.2e} {'':>6s} {'':>6s} {r['relres']:8.1e} {r['divres']:8.1e}   |flux| {max(abs(x) for x in fl):.1e}")
            else:
                P(f"{m:5s} {r['N']:4d} {r['hG']:7.4f} {r['h']:7.4f} {fl[0]/pr[0]:9.4f} {fl[1]/pr[1]:9.4f} {r['slip_ratio']:7.4f} "
                  f"{r['energy_ratio']:7.4f} {r['H1_ratio']:7.3f} {r['leak_misfit']:7.4f} {r['H1']:9.2e} {fr(rh[i]):>6s} {fr(re_[i]):>6s} "
                  f"{r['relres']:8.1e} {r['divres']:8.1e}")
    P()

# ------------------------------------------------------------------ (A)(ii) ellipse / nonconvex
for geom, title in (("ellipse", "ellipse x^2/1.2^2 + y^2/0.7^2 = 1"), ("star", "nonconvex curve r = 1 + 0.3 cos 3t")):
    rs = [r for r in gd if r["geom"] == geom and r["method"] != "meta"]
    if not rs:
        continue
    P("=" * 110)
    P(f"A(ii). {title}; u = curl(phi^2 (x + y^2/2)/10), p = a (x^2 y - y + x/2); B = a 1, B100 = a 100; P: u = 0")
    P("=" * 110)
    for mu in sorted(set(r["mu"] for r in rs)):
        P(f"--- mu = {mu:g}")
        # velocity: GS on B100 and B, CNS (B), D ratios
        def grab(kind, m):
            return sorted([r for r in rs if r["mu"] == mu and r["kind"] == kind and r["method"] == m], key=lambda r: r["N"])
        cns = grab("B", "CNS"); cns100 = grab("B100", "CNS"); gs = grab("B", "GS"); gs100 = grab("B100", "GS")
        dd = grab("B100", "D")
        if cns:
            hG = [r["hG"] for r in cns]
            P("  (GS B100 and D rates are in h = h_Omega, since the leak is (h/mu)(p - pbar); all other rates in hG;")
            P("   CNS100-CNS = relative H1 difference of the CNS velocities for a = 100 and a = 1)")
            P(f"{'N':>4s} {'hG':>7s} {'GS B100 H1':>11s} {'r(h)':>5s} {'GS B1 H1':>9s} {'rate':>5s} {'CNS H1':>9s} {'rate':>5s} "
              f"{'CNS en':>9s} {'rate':>5s} {'CNS100-CNS':>10s} {'D100 slip':>9s} {'energy':>7s} {'Dr(h)':>6s}")
            r1 = rate([r["H1"] for r in gs100], [r["h"] for r in gs100]) if gs100 else []
            r2 = rate([r["H1"] for r in gs], [r["hG"] for r in gs]) if gs else []
            r3 = rate([r["H1"] for r in cns], hG); r4 = rate([r["energy"] for r in cns], hG)
            r5 = rate([r["H1"] for r in dd], [r["h"] for r in dd]) if dd else []
            for i, c in enumerate(cns):
                def at(lst, key, j=i):
                    x = [r for r in lst if r["N"] == c["N"]]
                    return x[0][key] if x else None
                ii = lambda lst: [k for k, r in enumerate(lst) if r["N"] == c["N"]]
                g100 = at(gs100, "H1"); g1 = at(gs, "H1"); d100 = [r for r in dd if r["N"] == c["N"]]
                c100 = at(cns100, "H1")
                rel = abs(c100 - c["H1"]) / c["H1"] if c100 is not None else None
                P(f"{c['N']:4d} {c['hG']:7.4f} {f(g100, 11)} {fr(r1[ii(gs100)[0]]) if ii(gs100) else '   -':>5s} "
                  f"{f(g1)} {fr(r2[ii(gs)[0]]) if ii(gs) else '   -':>5s} {c['H1']:9.3e} {fr(r3[i]):>5s} {c['energy']:9.3e} {fr(r4[i]):>5s} "
                  f"{f(rel, 10)} {f(d100[0]['slip_ratio'] if d100 else None)} {f(d100[0]['energy_ratio'] if d100 else None, 7)} "
                  f"{fr(r5[ii(dd)[0]]) if ii(dd) else '   -':>6s}")
            # pressure
            P("  pressure L2 error (modulo constants):")
            P(f"  {'N':>4s} {'GS B1':>9s} {'rate':>5s} {'CNS B1':>9s} {'rate':>5s} {'GS B100':>9s} {'rate':>5s} {'CNS B100':>9s} {'rate':>5s}")
            cols = [grab("B", "GS"), grab("B", "CNS"), grab("B100", "GS"), grab("B100", "CNS")]
            rts = [rate([r["pL2"] for r in c_], [r["hG"] for r in c_]) for c_ in cols]
            for i, c in enumerate(cns):
                row = f"  {c['N']:4d}"
                for c_, rt in zip(cols, rts):
                    j = [k for k, r in enumerate(c_) if r["N"] == c["N"]]
                    row += f" {c_[j[0]]['pL2']:9.3e} {fr(rt[j[0]]):>5s}" if j else f" {'-':>9s} {'-':>5s}"
                P(row)
        pp = grab("P", "GS"); pc = grab("P", "CNS")
        if pp:
            P("  test P (u = 0): GS output = -e is exactly the leak response; CNS returns u_h = 0 (p in Pi_h)")
            P(f"  {'N':>4s} {'hG':>7s} {'slip':>7s} {'energy':>7s} {'H1rat':>7s} {'misfit':>7s} {'tang/nrm':>8s} {'H1 rate':>7s} "
              f"{'CNS |u|_1':>9s} {'relres':>8s}")
            rt = rate([r["H1"] for r in pp], [r["h"] for r in pp])
            for i, r in enumerate(pp):
                c = [x for x in pc if x["N"] == r["N"]]
                P(f"  {r['N']:4d} {r['hG']:7.4f} {r['slip_ratio']:7.4f} {r['energy_ratio']:7.4f} {r['H1_ratio']:7.3f} "
                  f"{r['leak_misfit']:7.4f} {r['tang_over_norm']:8.4f} {fr(rt[i]):>7s} {c[0]['H1'] if c else float('nan'):9.1e} "
                  f"{r['relres']:8.1e}")
        P()

# ------------------------------------------------------------------ (B) Navier-Stokes
ns = load("ns")
for r in ns:
    note_health("Navier-Stokes (circle)", r)
if ns:
    P("=" * 110)
    P("B. Steady Navier-Stokes, circle geometry (svn_k std mesh), k = 4, mu = 100; test B (psi, p of test B1), test P (u = 0)")
    P("   nu N_h + c(u;u,v): c1 = standard, cs = skew.  nu < 1 reached by continuation in nu (exact solution nu-independent)")
    P("=" * 110)
    for nu in sorted(set(r["nu"] for r in ns), reverse=True):
        for conv in ("c1", "cs"):
            rs = [r for r in ns if r["nu"] == nu and r["conv"] == conv and r["method"] != "meta"]
            if not rs:
                continue
            gs = sorted([r for r in rs if r["kind"] == "B" and r["method"] == "GS"], key=lambda r: r["N"])
            cn = sorted([r for r in rs if r["kind"] == "B" and r["method"] == "CNS"], key=lambda r: r["N"])
            dd = sorted([r for r in rs if r["kind"] == "B" and r["method"] == "D"], key=lambda r: r["N"])
            if gs:
                P(f"--- nu = {nu:g}, {conv}, test B")
                P(f"{'N':>4s} {'hG':>7s} {'GS H1':>9s} {'rate':>5s} {'CNS H1':>9s} {'rate':>5s} {'CNS en':>9s} {'rate':>5s} "
                  f"{'D slip':>7s} {'D en':>7s} {'D H1':>7s} {'D rate':>6s} {'GS pL2':>9s} {'CNS pL2':>9s} {'rate':>5s} {'Newton its':>14s}")
                hG = [r["hG"] for r in gs]
                a = rate([r["H1"] for r in gs], hG); b = rate([r["H1"] for r in cn], hG); c = rate([r["energy"] for r in cn], hG)
                d = rate([r["H1"] for r in dd], hG); e = rate([r["pL2"] for r in cn], hG)
                for i in range(len(gs)):
                    its = f"{sum(gs[i]['newton_its_path'])}/{sum(cn[i]['newton_its_path'])}"
                    P(f"{gs[i]['N']:4d} {gs[i]['hG']:7.4f} {gs[i]['H1']:9.3e} {fr(a[i]):>5s} {cn[i]['H1']:9.3e} {fr(b[i]):>5s} "
                      f"{cn[i]['energy']:9.3e} {fr(c[i]):>5s} {dd[i]['slip_ratio']:7.4f} {dd[i]['energy_ratio']:7.4f} "
                      f"{dd[i]['H1_ratio']:7.3f} {fr(d[i]):>6s} {gs[i]['pL2']:9.3e} {cn[i]['pL2']:9.3e} {fr(e[i]):>5s} {its:>14s}")
            pp = sorted([r for r in rs if r["kind"] == "P" and r["method"] == "GS"], key=lambda r: r["N"])
            pc = sorted([r for r in rs if r["kind"] == "P" and r["method"] == "CNS"], key=lambda r: r["N"])
            if pp:
                P(f"--- nu = {nu:g}, {conv}, test P (u = 0; leak (h/(nu mu))(p - pbar))")
                P(f"{'N':>4s} {'hG':>7s} {'slip':>7s} {'energy':>7s} {'H1rat':>7s} {'misfit':>7s} {'CNS |u|_1':>9s}")
                for r, c in zip(pp, pc):
                    P(f"{r['N']:4d} {r['hG']:7.4f} {r['slip_ratio']:7.4f} {r['energy_ratio']:7.4f} {r['H1_ratio']:7.3f} "
                      f"{r['leak_misfit']:7.4f} {c['H1']:9.1e}")
            P()

# ------------------------------------------------------------------ (C) Schafer-Turek
st = load("st")
for r in st:
    note_health("Schafer-Turek", r)
if st:
    REF = dict(cD=5.57953523384, cL=0.010618948146, dp=0.11752016697)
    P("=" * 110)
    P("C. Schafer-Turek 2D-1 (Re = 20, nu = 1e-3, U_max = 0.3), k = 4, Newton, do-nothing outflow")
    P(f"   reference cD = {REF['cD']}, cL = {REF['cL']}, dp = {REF['dp']}.  err = computed - reference.")
    P("   bdry = -int_{Gamma_h} (-p_h n + nu D(u_h) n);  vol = -omega(chi_h) (Gjerde-Scott (7)), chi_h = e_j at the Gamma_h nodes")
    P("   flux = int_{Gamma_h} u_h.n; heuristic leak prediction (h/(nu mu)) int_{Gamma_h} p_h (for GS and CNS*, see README)")
    P("=" * 110)
    for mu in sorted(set(r["mu"] for r in st)):
        for m in ("GS", "CNS", "NAIVE"):
            rr = sorted([r for r in st if r["mu"] == mu and r["method"] == m], key=lambda r: r["N"])
            if not rr:
                continue
            P(f"--- mu = {mu:g}, {m}")
            P(f"{'N':>4s} {'ntri':>6s} {'hG':>8s} {'cD vol':>12s} {'err':>10s} {'rate':>5s} {'cD bdry':>12s} {'err':>10s} {'rate':>5s} "
              f"{'cL vol':>10s} {'err':>9s} {'cL bdry':>10s} {'err':>9s} {'dp err':>9s} {'flux':>9s} {'pred':>9s} {'its':>3s}")
            hG = [r["hG"] for r in rr]
            a = rate([abs(r["cD_vol"] - REF["cD"]) for r in rr], hG); b = rate([abs(r["cD_bdry"] - REF["cD"]) for r in rr], hG)
            for i, r in enumerate(rr):
                pred = r["h"] / (r["nu"] * r["mu"]) * r["int_p_Gamma"]
                P(f"{r['N']:4d} {r['ntri']:6d} {r['hG']:8.5f} {r['cD_vol']:12.7f} {r['cD_vol'] - REF['cD']:10.2e} {fr(a[i]):>5s} "
                  f"{r['cD_bdry']:12.7f} {r['cD_bdry'] - REF['cD']:10.2e} {fr(b[i]):>5s} {r['cL_vol']:10.6f} "
                  f"{r['cL_vol'] - REF['cL']:9.1e} {r['cL_bdry']:10.6f} {r['cL_bdry'] - REF['cL']:9.1e} {r['dp'] - REF['dp']:9.1e} "
                  f"{r['flux']:9.2e} {pred:9.2e} {len(r['newton_hist']):3d}")
        P()

# ------------------------------------------------------------------ health / runtime
P("=" * 110)
P("Solver health (relres = relative residual of the UNregularised system; divres = ||div u_h||_{M^-1})")
P("=" * 110)
for k, h in health.items():
    P(f"{k:28s} solves {h['n']:4d}   max relres {h['relres']:.1e}   max divres {h['divres']:.1e}   not converged {h['notconv']}")
P()
P("Runtimes (wall seconds per job, 1 thread) and peak RSS")
for sub in ("gd", "ns", "st"):
    meta = [r for r in load(sub) if r.get("method") == "meta"]
    if meta:
        s = [r["secs"] for r in meta]; m = [r["maxrss_gb"] for r in meta]
        P(f"  {sub:3s} jobs {len(meta):4d}  total {sum(s)/3600:6.2f} h  max {max(s):7.0f} s  peak RSS {max(m):.2f} GB")
ths = list(th.values())
if ths:
    P(f"  thresh jobs {len(ths):4d}  total {sum(r['secs'] for r in ths)/3600:6.2f} h  peak RSS {max(r['maxrss_gb'] for r in ths):.2f} GB")

txt = "\n".join(out) + "\n"
open(os.path.join(D, "report_v4.txt"), "w").write(txt)
print(txt)
