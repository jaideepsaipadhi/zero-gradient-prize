"""v5: pressure errors for GS (theta = 0) and CNS* (theta = 1), k = 4, meshes of svn.make_mesh.

One job = one (N, KIND, MU): both methods are solved on the same assembled system and the pressure
error is measured against

    p* = p~ - pbar_{Gamma_h}(p~)          (p~ = the analytic pressure, extended by its formula)

in L^2(Omega_h), for both methods (for CNS* this is the p* of Theorem D; for GS it is the target of
Prop. 4.1 of notes/v4/pressure_drag.tex).  Because either method might carry a different natural
constant, three normalisations are reported:
    pL2      = ||p* - p_h||                        (p* normalised by the Gamma_h mean, p_h as computed)
    pL2G     = ||p* - (p_h - lam_h)||              (both normalised to zero Gamma_h mean; lam_h = pbar_Gamma(p_h))
    pL2opt   = min_c ||p* - p_h - c||              (optimal constant removed = Omega_h mean of the error)
plus the first-layer / off-layer split (after optimal constant), L^inf on the first layer, lam_h, and
the Omega_h mean of the error.  Velocity norms (svn.errors) are recorded too, for a consistency check
against results/pod and results/local.

Usage:  python3 pressure_err.py N KIND MU [lu|pardiso]      -> $SVN_OUT/pressure/N{N}_{KIND}_mu{MU}.jsonl
Solver: 'lu' = svn.solve_lu (SuperLU, the reference); 'pardiso' = svn_k.solve_pardiso's algorithm on svn's assembly.
"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time, resource, warnings
warnings.filterwarnings("ignore")
import numpy as np
import svn

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ.get("SVN_OUT", os.path.join(ROOT, "..", "results", "v5"))


def gamma_mean(S, fn):
    """mean over Gamma_h of a function given as fn(t, ref, Xe) -> values at the Gauss points"""
    num = per = 0.0
    for (t, ref, Xe, we, n) in S.edge_data():
        num += (we * fn(t, ref, Xe)).sum(); per += we.sum()
    return num / per


def ph_edge(p):
    def f(t, ref, Xe):
        return svn._mono(svn.MON3, ref[:, 0], ref[:, 1]) @ p[10 * t:10 * t + 10]
    return f


def pressure_errors(S, p, pex):
    nt = len(S.tris)
    X = S.phys(svn.QX, svn.QY)
    wq = svn.QW[None, :] * np.abs(S.detJ)[:, None]
    ph = np.einsum('qj,tj->tq', svn.Q3_Q, p.reshape(nt, 10))
    pbar_ex = gamma_mean(S, lambda t, ref, Xe: pex(Xe[:, 0], Xe[:, 1]))
    lam_h = gamma_mean(S, ph_edge(p))
    pstar = pex(X[..., 0], X[..., 1]) - pbar_ex
    err = pstar - ph
    area = wq.sum()
    emean = (wq * err).sum() / area
    L2 = lambda e, m=slice(None): float(np.sqrt((wq[m] * e[m] ** 2).sum()))
    layer = np.zeros(nt, bool); layer[[t for (t, _, _) in S.bedges]] = True
    eo = err - emean
    return dict(pL2=L2(err), pL2G=L2(err + lam_h), pL2opt=L2(eo),
                pL2_layer=L2(eo, layer), pL2_off=L2(eo, ~layer),
                pLinf_layer=float(np.abs(eo[layer]).max()),
                pLinf_layer_G=float(np.abs(err[layer] + lam_h).max()),
                lam_h=float(lam_h), err_mean=float(emean), pbar_ex=float(pbar_ex),
                pstar_L2=L2(pstar), ph_L2=L2(ph))


def solve_pardiso(S, sysm, g, theta, eps=1e-12, maxref=10):
    """svn_k.solve_pardiso (MKL PARDISO, pressure block regularised by eps*M, iterative refinement;
    relres of the UNregularised system) transcribed onto svn's k = 4 assembly and coupling."""
    import scipy.sparse as sp, pypardiso
    A, B, F = sysm["A"], sysm["B"], sysm["F"]
    Bt = svn.coupling(sysm, theta)
    M = sp.block_diag(list(sysm["Ml"]), format="csr")
    ndof = A.shape[0]
    dD = svn.dofs(S.dnodes).ravel(); free = np.setdiff1d(np.arange(ndof), dD)
    ug = np.zeros(ndof); ug[dD] = g(S.xy[S.dnodes, 0], S.xy[S.dnodes, 1]).ravel()
    Aff = A[free][:, free]; Btf = Bt[:, free]; Bf = B[:, free]
    K = sp.bmat([[Aff, -Btf.T], [Bf, eps * M]], format="csr"); K.sort_indices()
    rhs = np.concatenate([F[free] - A[free][:, dD] @ ug[dD], -B[:, dD] @ ug[dD]])
    slv = pypardiso.PyPardisoSolver(mtype=11)
    slv.set_iparm(11, 1); slv.set_iparm(13, 1)
    slv.factorize(K)
    x = slv.solve(K, rhs); nb = np.linalg.norm(rhs)
    for _ in range(maxref):
        res = rhs - K @ x
        if np.linalg.norm(res) < 1e-15 * nb:
            break
        x = x + slv.solve(K, res)
    slv.free_memory(everything=True)
    nf = len(free)
    u = ug.copy(); u[free] = x[:nf]; p = x[nf:]
    Bu = B @ u
    dres = np.sqrt(abs(Bu @ (sysm["Minv"] @ Bu)))
    r0 = np.concatenate([Aff @ x[:nf] - Btf.T @ p, Bf @ x[:nf]]) - rhs
    return u, p, [dres], np.linalg.norm(r0) / nb


def run(N, kind, mu, solver="lu"):
    t0 = time.time()
    ex = svn.make_exact(kind)
    verts, tris, circ, outer = svn.make_mesh(N)
    S = svn.Space(verts, tris, circ, outer)
    sysm = svn.assemble(S, mu, ex["f"], None)
    recs = []
    for th in (0, 1):
        t1 = time.time()
        if solver == "lu":
            u, p, hist, relres = svn.solve_lu(S, sysm, ex["u"], th)
        else:
            u, p, hist, relres = solve_pardiso(S, sysm, ex["u"], th)
        E = svn.errors(S, u, ex, mu)
        E.update(pressure_errors(S, p, ex["p"]))
        E.update(N=N, kind=kind, mu=mu, theta=th, method=("GS", "CNS*")[th], solver=solver,
                 h=S.h, hG=S.hGamma, ntri=len(tris), ndof=2 * S.nn, divres=float(hist[-1]),
                 relres=float(relres), secs=time.time() - t1)
        recs.append({k: (float(v) if isinstance(v, (np.floating, float)) else v) for k, v in E.items()})
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2
    for r in recs:
        r["maxrss_gb"] = rss
    os.makedirs(os.path.join(OUT, "pressure"), exist_ok=True)
    fn = os.path.join(OUT, "pressure", f"N{N}_{kind}_mu{mu:g}.jsonl")
    with open(fn + ".tmp", "w") as f:
        for r in recs:
            f.write(json.dumps(r) + "\n")
    os.replace(fn + ".tmp", fn)
    return recs


if __name__ == "__main__":
    N, kind, mu = int(sys.argv[1]), sys.argv[2], float(sys.argv[3])
    solver = sys.argv[4] if len(sys.argv) > 4 else "lu"
    for r in run(N, kind, mu, solver):
        print(json.dumps({k: (float(f"{v:.5g}") if isinstance(v, float) else v) for k, v in r.items()}))
