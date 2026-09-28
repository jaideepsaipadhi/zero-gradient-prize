"""One job = one unit of v4 work (own process).  Output: JSONL under results/v4 (env SVN_OUT overrides).

  python3 job_v4.py gd GEOM N KIND MU      general domains (Section 10), k = 4
        GEOM = twodisk | ellipse | star ;  KIND = X (two-disk counterexample, u = 0, p = c_i near disk i)
                                                   | B (manufactured u = curl(phi^2 g), p cubic) | P (u = 0, same p)
        methods: GS, CNS (global mean), and PER (per-component means) for twodisk; D = GS - CNS for B
  python3 job_v4.py thresh GEOM N          coercivity threshold mu* on the general mesh (penalty_k)
  python3 job_v4.py ns N KIND MU NU CONV   steady Navier-Stokes on the paper's circle mesh, k = 4,
        KIND = B | P (circle test B / P), CONV = c1 | cs; NU < 1 is reached by continuation in nu
        (the exact solution does not depend on nu).  methods GS, CNS; D = GS - CNS for B
  python3 job_v4.py st N MU                Schafer-Turek 2D-1 (Re = 20), k = 4, do-nothing outflow;
        methods GS, CNS (global mean, as printed in the paper), NAIVE (theta = 1 without a mean: the
        consistent method when Sigma has an outflow part); drag/lift by boundary traction and volume formula
"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time, warnings, resource
warnings.filterwarnings("ignore")
import numpy as np
import svn_k, svn_gen as sg, mesh_general as mg, exact_v4 as ev

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ.get("SVN_OUT", os.path.join(ROOT, "..", "results", "v4"))
K = 4

GEOMS = {
    "twodisk": dict(box=(-3.0, 3.0, -2.0, 2.0), curves=lambda: [mg.Circle(-1.4, 0.0, 0.5), mg.Circle(1.4, 0.0, 0.8)],
                    extra=[((-0.3, -2.0), (-0.3, 2.0)), ((0.3, -2.0), (0.3, 2.0))]),
    "ellipse": dict(box=(-2.5, 2.5, -2.5, 2.5), curves=lambda: [mg.Ellipse(0.0, 0.0, 1.2, 0.7)], extra=[]),
    "star": dict(box=(-2.5, 2.5, -2.5, 2.5), curves=lambda: [mg.PolarCurve(0.0, 0.0, 1.0, 0.3, 3)], extra=[]),
}
MESHOPT = dict(R_cap=4.0, grad=3.0, area_fac=3.0)
ST = dict(box=(0.0, 2.2, 0.0, 0.41), R_cap=8.0, grad=40.0, area_fac=3.0, Um=0.3, H=0.41, nu=1e-3,
          cD=5.57953523384, cL=0.010618948146, dp=0.11752016697)


def _maxrss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2


def _write(path, recs):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + ".tmp", "w") as f:
        for r in recs:
            f.write(json.dumps(r) + "\n")
    os.replace(path + ".tmp", path)


def _clean(d):
    out = {}
    for k, v in d.items():
        if isinstance(v, (np.floating,)):
            v = float(v)
        elif isinstance(v, np.ndarray):
            v = v.tolist()
        elif isinstance(v, (np.integer,)):
            v = int(v)
        out[k] = v
    return out


def gen_mesh(geom, N):
    G = GEOMS[geom]
    return mg.make_mesh(G["box"], G["curves"](), N, extra_lines=G["extra"], **MESHOPT)


def leak_block(S, ext, u, pex, gs, nu, mu, e_errs):
    """normalised leak diagnostics for a velocity error e (u_h itself when u = 0)."""
    G = gs["G"]; lam = S.h / (nu * mu)
    d = dict(G=G, pbar=gs["pbar"], pbar_i=gs["pbar_i"], len_i=gs["len_i"])
    if G > 0:
        d.update(slip_ratio=e_errs["bL2"] / (lam * G),
                 energy_ratio=e_errs["energy"] * nu * np.sqrt(mu / S.h) / G,
                 H1_ratio=e_errs["H1"] / (lam * G))
    return d


def gd(geom, N, kind, mu):
    t0 = time.time()
    M = gen_mesh(geom, N)
    S = sg.GSpace(M, svn_k.Ref(K))
    curves = M.curves
    if kind == "X":
        ex = ev.twodisk(1.0, -1.0, 0.3); methods = ("GS", "CNS", "PER")
    else:
        ex = ev.levelset(curves[0], kind); methods = ("GS", "CNS")
    sysm = svn_k.assemble(S, mu, ex["f"]); ext = sg.gamma_extras(S)
    gs = sg.gamma_stats(S, ex["p"], curves)
    base = dict(geom=geom, N=N, kind=kind, mu=mu, k=K, h=S.h, hG=S.hGamma, ntri=len(S.tris), ndof=2 * S.nn,
                npres=int(sysm["B"].shape[0]), mesh=M.stats, t_setup=time.time() - t0)
    zero = dict(u=lambda X, Y: np.zeros(np.shape(X) + (2,)), gu=lambda X, Y: np.zeros(np.shape(X) + (2, 2)))
    recs, sol = [], {}
    lam = S.h / mu
    for m in methods:
        t1 = time.time()
        u, p, info = sg.solve(S, sysm, ext, m, ex["u"])
        E = svn_k.errors(S, u, ex, mu)
        E.update(sg.pressure_errors(S, p, ex["p"]))
        flux = sg.fluxes(ext, u)
        pred = [lam * L * (pb - gs["pbar"]) for L, pb in zip(gs["len_i"], gs["pbar_i"])]
        r = dict(base); r.update(E); r.update(method=m, flux=flux.tolist(), flux_pred_GS=pred,
                                    relres=info["relres"], divres=info["divres"], secs=time.time() - t1)
        if kind in ("X", "P"):                           # u = 0: e = -u_h is exactly the leak response
            r.update(leak_block(S, ext, u, ex["p"], gs, 1.0, mu, E))
            r.update(sg.leak_profile(S, u, ex["p"], gs["pbar"], lam))
        recs.append(_clean(r)); sol[m] = u
    if kind.startswith("B"):
        Dn = svn_k.errors(S, sol["GS"] - sol["CNS"], zero, mu)
        d = dict(base); d.update(method="D", **Dn); d.update(leak_block(S, ext, None, ex["p"], gs, 1.0, mu, Dn))
        d.update(sg.leak_profile(S, sol["GS"] - sol["CNS"], ex["p"], gs["pbar"], lam))
        d["flux"] = sg.fluxes(ext, sol["GS"] - sol["CNS"]).tolist()
        recs.append(_clean(d))
    recs.append(dict(method="meta", secs=time.time() - t0, maxrss_gb=_maxrss_gb(), **{k: base[k] for k in ("geom", "N", "kind", "mu")}))
    _write(os.path.join(OUT, "gd", f"{geom}_N{N}_{kind}_mu{mu:g}.jsonl"), recs)


def thresh(geom, N):
    import penalty_k
    t0 = time.time()
    if geom == "st":
        M = mg.make_mesh(ST["box"], [mg.Circle(0.2, 0.2, 0.05)], N, R_cap=ST["R_cap"], grad=ST["grad"],
                         area_fac=ST["area_fac"])
    else:
        M = gen_mesh(geom, N)
    S = sg.GSpace(M, svn_k.Ref(K))                        # all of Sigma Dirichlet (symmetric form)
    mu = penalty_k.threshold(S)
    rec = dict(geom=geom, N=N, k=K, ntri=len(S.tris), hG=S.hGamma, hOmega=S.h, h_over_hG=S.h / S.hGamma,
               mu_star=float(mu), gamma_star=float(mu * S.hGamma / S.h), mesh=M.stats, secs=time.time() - t0,
               maxrss_gb=_maxrss_gb())
    _write(os.path.join(OUT, "thresh", f"{geom}_N{N}.json"), [_clean(rec)])


def ns(N, kind, mu, nu, conv):
    t0 = time.time()
    S = sg.circle_space(N, K)
    ext = sg.gamma_extras(S)
    path = [1.0] if nu >= 1 else list(np.geomspace(1.0, nu, int(np.ceil(np.log(1 / nu) / np.log(2.2))) + 1))
    fine_path = [1.0] if nu >= 1 else list(np.geomspace(1.0, nu, int(np.ceil(np.log(1 / nu) / np.log(1.3))) + 1))
    gsx = sg.gamma_stats(S, ev.circle(kind, nu, 1)["p"], [mg.Circle(0, 0, 1.0)])
    base = dict(N=N, kind=kind, mu=mu, nu=nu, conv=conv, k=K, h=S.h, hG=S.hGamma, ntri=len(S.tris),
                ndof=2 * S.nn, nu_path=[float(x) for x in path])
    zero = dict(u=lambda X, Y: np.zeros(np.shape(X) + (2,)), gu=lambda X, Y: np.zeros(np.shape(X) + (2, 2)))
    recs, sol = [], {}
    for m in ("GS", "CNS"):
        t1 = time.time()
        for pth in (path, fine_path):                     # retry with a finer nu-path if Newton fails
            x0 = None; its = []; conv_ok = True
            for nuj in pth:
                ex = ev.circle(kind, nuj, 1)
                sysm = svn_k.assemble(S, mu, ex["f"])
                u, p, info = sg.solve(S, sysm, ext, m, ex["u"], nu=nuj, conv=conv, x0=x0)
                x0 = (u, p); its.append(info["newton_its"]); conv_ok &= info["converged"]
                if not info["converged"]:
                    break
            if conv_ok:
                break
        E = svn_k.errors(S, u, ex, mu)
        E.update(sg.pressure_errors(S, p, ex["p"]))
        r = dict(base); r.update(E); r.update(method=m, relres=info["relres"], divres=info["divres"],
                                    newton_hist=info["newton"], newton_its_path=its, nu_path_used=[float(v) for v in pth],
                                    converged=bool(conv_ok),
                                    flux=sg.fluxes(ext, u).tolist(), secs=time.time() - t1)
        if kind == "P":
            r.update(leak_block(S, ext, u, ex["p"], gsx, nu, mu, E))
            r.update(sg.leak_profile(S, u, ex["p"], gsx["pbar"], S.h / (nu * mu)))
        recs.append(_clean(r)); sol[m] = u
    if kind == "B":
        Dn = svn_k.errors(S, sol["GS"] - sol["CNS"], zero, mu)
        d = dict(base); d.update(method="D", **Dn); d.update(leak_block(S, ext, None, ex["p"], gsx, nu, mu, Dn))
        d.update(sg.leak_profile(S, sol["GS"] - sol["CNS"], ex["p"], gsx["pbar"], S.h / (nu * mu)))
        recs.append(_clean(d))
    recs.append(dict(method="meta", secs=time.time() - t0, maxrss_gb=_maxrss_gb(), N=N, kind=kind, mu=mu, nu=nu, conv=conv))
    _write(os.path.join(OUT, "ns", f"circle_N{N}_{kind}_mu{mu:g}_nu{nu:g}_{conv}.jsonl"), recs)


def st(N, mu, methods=("GS", "CNS", "NAIVE")):
    t0 = time.time()
    M = mg.make_mesh(ST["box"], [mg.Circle(0.2, 0.2, 0.05)], N, R_cap=ST["R_cap"], grad=ST["grad"],
                     area_fac=ST["area_fac"])
    S = sg.GSpace(M, svn_k.Ref(K), dirichlet_sides=(1, 3, 4), outflow_sides=(2,))
    Um, H, nu = ST["Um"], ST["H"], ST["nu"]
    g = lambda X, Y: np.stack([4 * Um * Y * (H - Y) / H**2 * (X < 1e-12), 0 * X], -1)
    f = lambda X, Y: np.zeros(np.shape(X) + (2,))
    sysm = svn_k.assemble(S, mu, f); ext = sg.gamma_extras(S); Tout = nu * sg.outflow_matrix(S)
    base = dict(N=N, mu=mu, nu=nu, k=K, h=S.h, hG=S.hGamma, ntri=len(S.tris), ndof=2 * S.nn,
                npres=int(sysm["B"].shape[0]), mesh=M.stats)
    recs = []
    for m in methods:
        t1 = time.time()
        u, p, info = sg.solve(S, sysm, ext, m, g, nu=nu, conv="c1", Tout=Tout)
        F = sg.forces(S, sysm, u, p, nu)
        pf, pb = sg.point_pressure(S, p, [(0.15, 0.2), (0.25, 0.2)])
        # discrete pressure on Gamma_h: mean and int p ds (for the heuristic leak predictions)
        pint = 0.0; plen = 0.0
        for (t, ref, Xe, we, n) in S.edge_data():
            pv = S.R.pbasis(ref[:, 0], ref[:, 1]) @ p[t * S.R.npl:(t + 1) * S.R.npl]
            pint += we @ pv; plen += we.sum()
        # outflow pressure mean (do-nothing normalises p there)
        r = dict(base)
        r.update(method=m, cD_bdry=500 * F["Fb"][0], cD_vol=500 * F["Fv"][0], cL_bdry=500 * F["Fb"][1],
                 cL_vol=500 * F["Fv"][1], dp=pf - pb, flux=float(sg.fluxes(ext, u)[0]),
                 int_p_Gamma=float(pint), len_Gamma_h=float(plen),
                 relres=info["relres"], divres=info["divres"], newton_hist=info["newton"],
                 converged=info["converged"], secs=time.time() - t1)
        recs.append(_clean(r))
    recs.append(dict(method="meta", secs=time.time() - t0, maxrss_gb=_maxrss_gb(), N=N, mu=mu))
    _write(os.path.join(OUT, "st", f"st_N{N}_mu{mu:g}.jsonl"), recs)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "gd":
        gd(sys.argv[2], int(sys.argv[3]), sys.argv[4], float(sys.argv[5]))
    elif mode == "thresh":
        thresh(sys.argv[2], int(sys.argv[3]))
    elif mode == "ns":
        ns(int(sys.argv[2]), sys.argv[3], float(sys.argv[4]), float(sys.argv[5]), sys.argv[6])
    elif mode == "st":
        st(int(sys.argv[2]), float(sys.argv[3]))
    else:
        raise SystemExit(__doc__)
