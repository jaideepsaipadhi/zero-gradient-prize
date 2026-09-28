"""Rescue studies for the printed method (notes/v4/rescue.tex, labels rs:).  New file; reuses svn_k unchanged.

Task A ("rescue B"): GS with a mesh-dependent penalty mu = c (h_ref/h)^alpha, h_ref = h(N=16).
Task B: flux-corrected outer data  g~_I = g_I - (m_R / (8 L^2)) x,  m_R = int_{dR} g_I . nu  (exact for the
        P_k interpolant), so that int_{dR} g~_I . nu = 0 up to round-off.  x (the position field) is linear, hence
        reproduced exactly by the P_k interpolant, and int_{dR} x . nu = int_R div x = 8 L^2.

  python3 rescue.py study  KIND N ALPHA C          -> results/v4/rescue/study/{KIND}_N{N}_a{ALPHA}_c{C}.jsonl
  python3 rescue.py musweep KIND N MU[,MU...]      -> results/v4/rescue/musweep/{KIND}_N{N}.jsonl  (GS and CNS*)
  python3 rescue.py flux   KIND N MU[,MU...]       -> results/v4/rescue/flux/{KIND}_N{N}.jsonl
                                                     (CNS*, uncorrected vs flux-corrected data)
  python3 rescue.py cond   N MU[,MU...]            -> results/v4/rescue/cond/N{N}.jsonl
                                                     (extreme eigenvalues of the velocity block and of the
                                                      symmetric GS saddle matrix; pressure Schur complement)
  python3 rescue.py fluxsize KIND N[,N...]         -> prints m_R, |Gamma_h|, ||u~||_{Gamma_h}

KIND: A, B1, B100, P1 (svn.make_exact) or F, F2 (non-polynomial boundary flux, defined below).
All runs: k = 4, standard meshes, L = 2.5, 1 thread.  Solver: SuperLU (svn_k.solve_lu), or MKL PARDISO
(svn_k.solve_pardiso, pressure block regularised by 1e-12 M) if env RESCUE_SOLVER=pardiso (used for N = 96,
where SuperLU needs 5.6 GB).  Every record stores the solver and the relative residual of the unregularised system.
"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time, warnings, resource
warnings.filterwarnings("ignore")
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import sympy as sy
import svn_k

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ.get("RESCUE_OUT", os.path.join(ROOT, "..", "results", "v4", "rescue"))
K_DEG = 4
L_BOX = 2.5
H_REF = None          # h at N = 16, filled lazily


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
        if isinstance(v, (np.floating, float)):
            out[k] = float(v)
        elif isinstance(v, (np.integer,)):
            out[k] = int(v)
        else:
            out[k] = v
    return out


# ---------------------------------------------------------------- exact solutions
def make_exact(kind):
    """Tests A, B*, P*, C from svn.make_exact; F and F2 are new.

    F : psi = (r^2-1)^3 sin(2x+1) cos(3y/2+7/10) / 400,  p = x^2 y - y + x/2.
        The cube makes grad u = 0 on Gamma (wall shear W = 0), so u~ = O(h_Gamma^4) on Gamma_h and the
        geometric error is small; g.nu on dR is not a polynomial, so m_R != 0.
    F2: psi = (r^2-1)^2 sin(2x+1) cos(3y/2+7/10) / 40, same p (W != 0; realistic)."""
    if kind not in ("F", "F2"):
        return svn_k.make_exact(kind)
    x, y = sy.symbols('x y', real=True)
    r2 = x**2 + y**2
    if kind == "F":
        psi = (r2 - 1)**3 * sy.sin(2 * x + 1) * sy.cos(sy.Rational(3, 2) * y + sy.Rational(7, 10)) / 400
    else:
        psi = (r2 - 1)**2 * sy.sin(2 * x + 1) * sy.cos(sy.Rational(3, 2) * y + sy.Rational(7, 10)) / 40
    p = x**2 * y - y + x / 2
    u = [sy.diff(psi, y), -sy.diff(psi, x)]
    f = [-(sy.diff(u[i], x, 2) + sy.diff(u[i], y, 2)) + sy.diff(p, [x, y][i]) for i in range(2)]
    gu = [[sy.diff(u[i], v) for v in (x, y)] for i in range(2)]
    U = sy.lambdify((x, y), u, 'numpy'); Fl = sy.lambdify((x, y), f, 'numpy')
    GU = sy.lambdify((x, y), gu, 'numpy'); P = sy.lambdify((x, y), p, 'numpy')

    def vec(fn):
        def g(X, Y):
            vals = fn(X, Y)
            return np.stack([np.broadcast_to(np.asarray(v, float), np.shape(X)) for v in vals], -1)
        return g

    def gmat(X, Y):
        vals = GU(X, Y)
        return np.stack([np.stack([np.broadcast_to(np.asarray(vals[i][j], float), np.shape(X))
                                   for j in range(2)], -1) for i in range(2)], -2)
    return dict(u=vec(U), f=vec(Fl), gu=gmat,
                p=lambda X, Y: np.broadcast_to(np.asarray(P(X, Y), float), np.shape(X)))


# ---------------------------------------------------------------- boundary fluxes
def outer_edges(S):
    """(t, iu, iw) for triangle edges on the square."""
    L = np.abs(S.verts).max()
    on = lambda P: abs(abs(P[0]) - L) < 1e-10 or abs(abs(P[1]) - L) < 1e-10
    out = []
    for t in range(len(S.tris)):
        for (a, b) in [(0, 1), (1, 2), (2, 0)]:
            A = S.verts[S.tris[t, a]]; B = S.verts[S.tris[t, b]]
            if on(A) and on(B) and (abs(A[0] - B[0]) < 1e-12 or abs(A[1] - B[1]) < 1e-12):
                out.append((t, a, b))
    return out


def flux_outer(S, g):
    """m_R = int_{dR} g_I . nu for the P_k Lagrange interpolant g_I of g (Gauss rule exact for P_k)."""
    R = S.R
    Rr = np.array([[0, 0], [1, 0], [0, 1]], float)
    m = 0.0
    for (t, a, b) in outer_edges(S):
        ref = Rr[a][None, :] + R.GL[:, None] * (Rr[b] - Rr[a])[None, :]
        ph = R.basis(ref[:, 0], ref[:, 1])                      # (q, nl)
        nodes = S.xy[S.ids[t]]
        gv = g(nodes[:, 0], nodes[:, 1])                        # (nl, 2)
        gi = ph @ gv
        A = S.verts[S.tris[t, a]]; B = S.verts[S.tris[t, b]]
        Lg = np.linalg.norm(B - A); M = (A + B) / 2
        nu = np.sign(M) * (np.abs(M) >= np.abs(M).max() - 1e-12)   # outward normal of the square side
        m += (R.GLW * Lg) @ (gi @ nu)
    return m


def flux_gamma(S, u):
    """int_{Gamma_h} u_h . n (n outward from Omega_h, i.e. into the polygon) and |Gamma_h|."""
    R = S.R
    fl, per = 0.0, 0.0
    for (t, ref, Xe, we, n) in S.edge_data():
        ph = R.basis(ref[:, 0], ref[:, 1])
        ue = ph @ u[svn_k.dofs(S.ids[t])]
        fl += we @ (ue @ n); per += we.sum()
    return fl, per


def trace_norm_exact(S, ex):
    """||u~||_{L2(Gamma_h)} of the exact velocity (extension) on the chords."""
    s = 0.0
    for (t, ref, Xe, we, n) in S.edge_data():
        s += (we * (ex["u"](Xe[:, 0], Xe[:, 1]) ** 2).sum(-1)).sum()
    return np.sqrt(s)


def corrected(S, g):
    """Flux-corrected Dirichlet function: g~ = g - c x, c = m_R / (8 L^2)."""
    mR = flux_outer(S, g)
    L = np.abs(S.verts).max()
    c = mR / (8 * L * L)
    gt = lambda X, Y: g(X, Y) - c * np.stack([X, Y], -1)
    return gt, mR, c


# ---------------------------------------------------------------- one solve
def _href():
    global H_REF
    if H_REF is None:
        H_REF = svn_k.build(16, K_DEG).h
    return H_REF


def solve_one(S, ex, mu, theta, g=None, sysm=None):
    if sysm is None:
        sysm = svn_k.assemble(S, mu, ex["f"])
    solve = dict(lu=svn_k.solve_lu, pardiso=svn_k.solve_pardiso)[os.environ.get("RESCUE_SOLVER", "lu")]
    u, p, hist, relres, _ = solve(S, sysm, g if g is not None else ex["u"], theta)
    E = svn_k.errors(S, u, ex, mu)
    fl, per = flux_gamma(S, u)
    E.update(solver=os.environ.get("RESCUE_SOLVER", "lu"), theta=theta, mu=mu, divres=float(hist[-1]), relres=float(relres), fluxG=fl, perG=per,
             gamma=mu * S.hGamma / S.h)
    return E, u


def study(kind, N, alpha, c):
    t0 = time.time()
    S = svn_k.build(N, K_DEG)
    mu = c * (_href() / S.h) ** alpha
    ex = make_exact(kind)
    sysm = svn_k.assemble(S, mu, ex["f"])
    recs = []
    for th in (0, 1):
        E, _ = solve_one(S, ex, mu, th, sysm=sysm)
        E.update(kind=kind, N=N, alpha=alpha, c=c, h=S.h, hG=S.hGamma, secs=time.time() - t0,
                 maxrss_gb=_maxrss_gb())
        recs.append(_clean(E))
    _write(os.path.join(OUT, "study", f"{kind}_N{N}_a{alpha:g}_c{c:g}.jsonl"), recs)


def musweep(kind, N, mus):
    S = svn_k.build(N, K_DEG); ex = make_exact(kind)
    recs = []
    for mu in mus:
        t0 = time.time()
        sysm = svn_k.assemble(S, mu, ex["f"])
        for th in (0, 1):
            E, _ = solve_one(S, ex, mu, th, sysm=sysm)
            E.update(kind=kind, N=N, h=S.h, hG=S.hGamma, secs=time.time() - t0, maxrss_gb=_maxrss_gb())
            recs.append(_clean(E))
    _write(os.path.join(OUT, "musweep", f"{kind}_N{N}.jsonl"), recs)


def flux(kind, N, mus):
    S = svn_k.build(N, K_DEG); ex = make_exact(kind)
    gt, mR, c = corrected(S, ex["u"])
    mR_after = flux_outer(S, gt)
    utr = trace_norm_exact(S, ex)
    recs = []
    for mu in mus:
        t0 = time.time()
        sysm = svn_k.assemble(S, mu, ex["f"])
        for corr in (0, 1):
            E, u = solve_one(S, ex, mu, 1, g=(gt if corr else None), sysm=sysm)
            floor = np.sqrt(mu / S.h) * abs(mR) / np.sqrt(E["perG"])
            E.update(kind=kind, N=N, h=S.h, hG=S.hGamma, corrected=corr, mR=mR, mR_after=mR_after, c_corr=c,
                     floor=floor, slip_floor=abs(mR) / np.sqrt(E["perG"]), utrace=utr,
                     secs=time.time() - t0, maxrss_gb=_maxrss_gb())
            recs.append(_clean(E))
    _write(os.path.join(OUT, "flux", f"{kind}_N{N}.jsonl"), recs)


def fluxsize(kind, Ns):
    ex = make_exact(kind)
    for N in Ns:
        S = svn_k.build(N, K_DEG)
        gt, mR, c = corrected(S, ex["u"])
        per = sum(we.sum() for (_, _, _, we, _) in S.edge_data())
        print(f"{kind} N={N:4d} h={S.h:.4f} hG={S.hGamma:.4f} m_R={mR:.3e} after={flux_outer(S, gt):.1e} "
              f"|G_h|={per:.5f} slipfloor={abs(mR)/np.sqrt(per):.3e} ||u~||_Gh={trace_norm_exact(S, ex):.3e}",
              flush=True)


# ---------------------------------------------------------------- conditioning
def cond(N, mus):
    """Extreme eigenvalues (Euclidean nodal basis, free dofs) of
         A_mu  (velocity block, SPD for mu >= mu*),
         K_mu = [[A_mu, B^T], [B, 0]]  (symmetric form of the GS saddle matrix; theta = 0),
       with the pressure in an L2-orthonormal element basis, and of the pressure Schur complement
       S_mu = B A_mu^{-1} B^T (L2 metric): its smallest eigenvalues are squared inf-sup constants in the A_mu-norm."""
    S = svn_k.build(N, K_DEG)
    zero_f = lambda X, Y: np.zeros(np.shape(X) + (2,))
    s0 = svn_k.assemble(S, 0.0, zero_f)
    dD = svn_k.dofs(S.dnodes).ravel(); free = np.setdiff1d(np.arange(2 * S.nn), dD)
    A0 = s0["A"][free][:, free]; P = s0["Pen"][free][:, free]
    B0 = s0["B"][:, free].tocsr()
    # L2-orthonormal pressure basis on each element (local Cholesky of the pressure mass matrix): B = Linv B0,
    # so the pressure metric is Euclidean and S = B A^{-1} B^T is the Schur complement in the L2 metric.
    Linv = sp.block_diag([np.linalg.inv(np.linalg.cholesky(m)) for m in s0["Ml"]], format="csr")
    B = (Linv @ B0).tocsr()
    npres = B.shape[0]
    recs = []
    for mu in mus:
        t0 = time.time()
        A = (A0 + mu * P).tocsc(); A = ((A + A.T) / 2).tocsc()
        lmaxA = spl.eigsh(A, k=1, which="LA", return_eigenvectors=False, tol=1e-8)[0]
        luA = spl.splu(A)
        opA = spl.LinearOperator(A.shape, matvec=luA.solve, dtype=float)
        lminA = 1.0 / spl.eigsh(opA, k=1, which="LA", return_eigenvectors=False, tol=1e-8)[0]
        K = sp.bmat([[A, B.T], [B, None]], format="csc")
        luK = spl.splu(K)
        nv = A.shape[0]
        opK = spl.LinearOperator(K.shape, matvec=luK.solve, dtype=float)
        smallK = spl.eigsh(opK, k=2, which="LM", return_eigenvectors=False, tol=1e-8)
        minabsK = 1.0 / np.max(np.abs(smallK))
        bigK = spl.eigsh(K, k=2, which="LM", return_eigenvectors=False, tol=1e-8)
        maxabsK = np.max(np.abs(bigK))
        Sop = spl.LinearOperator((npres,) * 2, matvec=lambda q: B @ luA.solve(B.T @ q), dtype=float)
        smax = spl.eigsh(Sop, k=1, which="LA", return_eigenvectors=False, tol=1e-6)[0]

        def Sinv(q):                      # S^{-1} q = -[K^{-1}(0, q)]_p
            return -luK.solve(np.concatenate([np.zeros(nv), q]))[nv:]
        inv = spl.eigsh(spl.LinearOperator((npres,) * 2, matvec=Sinv, dtype=float), k=2, which="LA",
                        return_eigenvectors=False, tol=1e-6)
        smin = 1.0 / np.sort(inv)[::-1]
        rec = dict(N=N, mu=mu, h=S.h, hG=S.hGamma, gamma=mu * S.hGamma / S.h, nfree=nv, npres=npres,
                   lmaxA=lmaxA, lminA=lminA, condA=lmaxA / lminA, maxabsK=maxabsK, minabsK=minabsK,
                   condK=maxabsK / minabsK, schur_max=smax, schur_min=float(smin[0]), schur_min2=float(smin[1]),
                   secs=time.time() - t0, maxrss_gb=_maxrss_gb())
        print(_clean(rec), flush=True)
        recs.append(_clean(rec))
        del luA, luK
    _write(os.path.join(OUT, "cond", f"N{N}.jsonl"), recs)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "study":
        study(sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))
    elif mode == "musweep":
        musweep(sys.argv[2], int(sys.argv[3]), [float(m) for m in sys.argv[4].split(",")])
    elif mode == "flux":
        flux(sys.argv[2], int(sys.argv[3]), [float(m) for m in sys.argv[4].split(",")])
    elif mode == "fluxsize":
        fluxsize(sys.argv[2], [int(n) for n in sys.argv[3].split(",")])
    elif mode == "cond":
        cond(int(sys.argv[2]), [float(m) for m in sys.argv[3].split(",")])
