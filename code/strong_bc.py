"""v5: strongly imposed no-slip on the polygon Gamma_h (no Nitsche), Scott-Vogelius P4/P3-disc,
reusing svn.py's mesh, space, volume assembly and error routines.

    find u_h in V_h, u_h = I_h g on the square, u_h = bc on the Gamma_h nodes,  p_h in Pi_h:
        1/2 (D u_h, D v) - (p_h, div v) = (f, v)     v = 0 on the square and on Gamma_h
        (q, div u_h) = 0

bc = "zero":   u_h = 0 at every P4 node of Gamma_h (Scott's setting: the no-slip value of the curved
               wall transferred to the polygon);
bc = "interp": u_h = u~ (the analytic u extended by its formula) at the Gamma_h nodes (a control run:
               the boundary data are then exact at the nodes, but the zero-gradient constraints remain
               because u~ != 0 there ... only at the vertices is u~ = 0).
With all velocity boundary values prescribed, the pressure is determined up to a constant; one pressure
unknown is pinned and p_h is then shifted to zero Omega_h-mean.  Divergence-free fields that
vanish on two non-collinear polygon edges have zero gradient at the common vertex; this is the
mechanism behind the h^{1/2} rate (paper, Section 1).

Usage: python3 strong_bc.py N KIND MU [zero|interp|nitsche] [lu|pardiso] [std|alt]
       -> $SVN_OUT/strong/N{N}_{KIND}_{bc}[_alt].jsonl
bc = "nitsche" runs GS and CNS* (svn.py) on the chosen mesh, for a same-mesh comparison.
mesh = "alt": first-layer diagonals alternated (see make_mesh_alt), polygon vertices in 2 or 4 triangles.
(MU only enters the reported energy norm, which uses the Nitsche weights mu/h and h for comparability.)
"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time, resource, warnings
warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
import svn
import pressure_err

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ.get("SVN_OUT", os.path.join(ROOT, "..", "results", "v5"))


def make_mesh_alt(N, L=2.5):
    """svn.make_mesh with the diagonals of the FIRST layer alternated: quad i (between polygon vertices
    i, i+1) is split along (z_i, z'_{i+1}) for even i and along (z_{i+1}, z'_i) for odd i.  Polygon
    vertices then alternately lie in 2 and in 4 triangles (svn.make_mesh: always 3).  At a vertex in 2
    triangles, a divergence-free P_k field vanishing on both Gamma_h edges has zero gradient
    (8 linear conditions on the 8 gradient entries), which is Scott's locking mechanism."""
    verts, tris, circ, outer = svn.make_mesh(N, L=L)
    tris = tris.copy()
    m = N // 4
    vid = np.arange((m + 1) * N).reshape(m + 1, N)
    for i in range(1, N, 2):                            # odd quads of layer 0: flip the diagonal
        i1 = (i + 1) % N
        a, b, c, d = vid[0, i], vid[0, i1], vid[1, i1], vid[1, i]
        tris[2 * i] = [a, b, d]; tris[2 * i + 1] = [b, c, d]
    P = verts[tris]
    ar = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    tris[ar < 0] = tris[ar < 0][:, [0, 2, 1]]
    return verts, tris, circ, outer


def build(N, mesh="std"):
    return svn.Space(*(svn.make_mesh(N) if mesh == "std" else make_mesh_alt(N)))


def gamma_nodes(S):
    """all P4 nodes lying on Gamma_h edges (vertices + 3 edge nodes per edge)"""
    loc = {(0, 1): [], (1, 2): [], (2, 0): []}
    for k, (l1, l2, l0) in enumerate(svn.LBARY):     # local vertex 0 <-> l0, 1 <-> l1, 2 <-> l2
        lam = (l0, l1, l2)
        for (a, b) in loc:
            c = 3 - a - b
            if lam[c] == 0:
                loc[(a, b)].append(k)
    nodes = set()
    for (t, u, w) in S.bedges:
        key = (u, w) if (u, w) in loc else (w, u)
        nodes.update(S.ids[t, loc[key]].tolist())
    return np.array(sorted(nodes))


def solve_strong(S, sysm, g, bc="zero", solver="lu"):
    A, B, F = sysm["Avol"], sysm["B"], sysm["F"]            # no Nitsche terms
    ndof = A.shape[0]; npr = B.shape[0]; nt = npr // 10
    gn = gamma_nodes(S)
    dS = svn.dofs(S.dnodes).ravel(); dG = svn.dofs(gn).ravel()
    dD = np.union1d(dS, dG)
    ug = np.zeros(ndof)
    ug[dS] = g(S.xy[S.dnodes, 0], S.xy[S.dnodes, 1]).ravel()
    if bc == "interp":
        ug[dG] = g(S.xy[gn, 0], S.xy[gn, 1]).ravel()
    else:
        ug[dG] = 0.0
    free = np.setdiff1d(np.arange(ndof), dD)
    # The pressure is determined up to a constant: pin the constant mode of element 0 (drop that pressure
    # unknown and its continuity row, which is implied by the others since the data have zero net flux),
    # then shift p_h to zero Omega_h-mean.  (A dense mean-zero bordering row caused heavy SuperLU fill.)
    keep = np.arange(1, npr)
    Bk = B[keep]
    K = sp.bmat([[A[free][:, free], -Bk[:, free].T], [Bk[:, free], None]], format="csc")
    rhs = np.concatenate([F[free] - A[free][:, dD] @ ug[dD], -Bk[:, dD] @ ug[dD]])
    if solver == "lu":
        lu = spl.splu(K, permc_spec="COLAMD", diag_pivot_thresh=1.0)
        x = lu.solve(rhs)
        for _ in range(3):
            x = x + lu.solve(rhs - K @ x)
        nnzLU = int(lu.L.nnz + lu.U.nnz); del lu
    else:     # MKL PARDISO (as svn_k.solve_pardiso): pressure block regularised by 1e-12*M, refinement
        import pypardiso
        Mk = sp.block_diag(list(sysm["Ml"]), format="csr")[keep][:, keep]
        K0 = K.tocsr()
        K = sp.bmat([[A[free][:, free], -Bk[:, free].T], [Bk[:, free], 1e-12 * Mk]], format="csr"); K.sort_indices()
        slv = pypardiso.PyPardisoSolver(mtype=11)
        slv.set_iparm(11, 1); slv.set_iparm(13, 1)
        slv.factorize(K)
        x = slv.solve(K, rhs); nb = np.linalg.norm(rhs)
        for _ in range(10):
            res = rhs - K @ x
            if np.linalg.norm(res) < 1e-15 * nb:
                break
            x = x + slv.solve(K, res)
        slv.free_memory(everything=True); nnzLU = -1
        K = K0                                        # relres of the UNregularised system
    relres = np.linalg.norm(K @ x - rhs) / np.linalg.norm(rhs)
    nf = len(free)
    u = ug.copy(); u[free] = x[:nf]
    p = np.zeros(npr); p[keep] = x[nf:]
    wq = svn.QW[None, :] * np.abs(S.detJ)[:, None]
    cvec = np.einsum('tq,qj->tj', wq, svn.Q3_Q).ravel()
    p[0::10] -= (cvec @ p) / wq.sum()                   # zero Omega_h-mean
    flux = float(np.sum(B[0::10] @ u))                  # int div u_h = net boundary flux of the data
    Bu = B @ u
    dres = np.sqrt(abs(Bu @ (sysm["Minv"] @ Bu)))
    return u, p, dict(divres=float(dres), relres=float(relres), flux=flux, nnzLU=nnzLU,
                      n_gamma_nodes=int(len(gn)))


def vertex_grad(S, u, ex):
    """max over Gamma_h vertices z and triangles T containing z of |grad u_h|_T(z) - grad u(z)|_max"""
    gx, gy = svn.dbasis(np.array([0., 1., 0.]), np.array([0., 0., 1.]))   # rows: local vertices 0,1,2
    oncirc = np.abs(np.linalg.norm(S.verts, axis=1) - 1) < 1e-12
    worst = 0.0
    for t in np.where(oncirc[S.tris].any(1))[0]:
        Ji = S.Jinv[t]
        Ut = u[svn.dofs(S.ids[t])]
        for a in range(3):
            v = S.tris[t, a]
            if not oncirc[v]:
                continue
            dx = Ji[0, 0] * gx[a] + Ji[1, 0] * gy[a]; dy = Ji[0, 1] * gx[a] + Ji[1, 1] * gy[a]
            Gh = np.stack([dx @ Ut, dy @ Ut], -1)            # Gh[c, d] = d_d u_c
            Ge = ex["gu"](S.verts[v:v + 1, 0], S.verts[v:v + 1, 1])[0]
            worst = max(worst, float(np.abs(Gh - Ge).max()))
    return worst


def run(N, kind, mu, bc="zero", solver="lu", mesh="std"):
    t0 = time.time()
    ex = svn.make_exact(kind)
    S = build(N, mesh); tris = S.tris
    sysm = svn.assemble(S, mu, ex["f"], None)
    if bc == "nitsche":                                  # GS and CNS* on the same mesh, for comparison
        recs = []
        for th in (0, 1):
            if solver == "lu":
                u, p, hist, relres = svn.solve_lu(S, sysm, ex["u"], th)
            else:
                u, p, hist, relres = pressure_err.solve_pardiso(S, sysm, ex["u"], th)
            E = svn.errors(S, u, ex, mu); E.update(pressure_err.pressure_errors(S, p, ex["p"]))
            E.update(N=N, kind=kind, mu=mu, method=("GS", "CNS*")[th], theta=th, mesh=mesh, solver=solver,
                     h=S.h, hG=S.hGamma, ntri=len(tris), ndof=2 * S.nn, divres=float(hist[-1]),
                     relres=float(relres), vgrad_err=vertex_grad(S, u, ex), secs=time.time() - t0,
                     maxrss_gb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2)
            recs.append({k: (float(v) if isinstance(v, (np.floating, float)) else v) for k, v in E.items()})
        _write(recs, N, kind, bc, mesh)
        return recs
    u, p, info = solve_strong(S, sysm, ex["u"], bc, solver)
    E = svn.errors(S, u, ex, mu)
    E.update(pressure_err.pressure_errors(S, p, ex["p"]))
    E.update(info)
    E.update(N=N, kind=kind, mu=mu, method="strong-" + bc, mesh=mesh, solver=solver, h=S.h, hG=S.hGamma, ntri=len(tris),
             ndof=2 * S.nn, vgrad_err=vertex_grad(S, u, ex), secs=time.time() - t0,
             maxrss_gb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2)
    rec = {k: (float(v) if isinstance(v, (np.floating, float)) else v) for k, v in E.items()}
    _write([rec], N, kind, bc, mesh)
    return [rec]


def _write(recs, N, kind, bc, mesh):
    os.makedirs(os.path.join(OUT, "strong"), exist_ok=True)
    fn = os.path.join(OUT, "strong", f"N{N}_{kind}_{bc}" + ("" if mesh == "std" else "_" + mesh) + ".jsonl")
    with open(fn + ".tmp", "w") as f:
        for r in recs:
            f.write(json.dumps(r) + "\n")
    os.replace(fn + ".tmp", fn)


if __name__ == "__main__":
    N, kind, mu = int(sys.argv[1]), sys.argv[2], float(sys.argv[3])
    bc = sys.argv[4] if len(sys.argv) > 4 else "zero"
    solver = sys.argv[5] if len(sys.argv) > 5 else "lu"
    mesh = sys.argv[6] if len(sys.argv) > 6 else "std"
    for r in run(N, kind, mu, bc, solver, mesh):
        print(json.dumps({k: (float(f"{v:.5g}") if isinstance(v, float) else v) for k, v in r.items()}))
