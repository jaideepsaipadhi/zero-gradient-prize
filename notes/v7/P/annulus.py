"""notes/v7/P: CNS* on locally periodic ANNULUS meshes with a prescribed first-layer cell (REPORT.tex, Section 5).

Omega_h = (polygon of radius R_out) minus (regular N-gon inscribed in the unit circle).
Rings r_0 = 1, r_{j+1} = r_j (1 + a*l), l = 2 sin(pi/N) (all cells similar, aspect a), last ring = R_out.
Vertex (j,i) at angle (i + j*delta) * 2pi/N; quads split forward: [v(j,i), v(j,i+1), v(j+1,i+1)], [v(j,i), v(j+1,i+1), v(j+1,i)].
  delta = 0   : rectangular lattice, apex of T_e over the counterclockwise endpoint (production type);
  delta = -1  : its mirror image;
  delta = -1/2: isosceles lattice, mirror symmetric about every radial line through a vertex or an edge midpoint.
Penalty: mu/h := lamhat / l exactly (S.h = 1), so the cell penalty is lamhat for every edge.
Test S: u = curl((r^2-1)^2 (1+x)/4) in P_4, p = x^2 y - y + x/2 in P_3 (interior discretisation exact).
Outputs: P_h(phi) = <e_p - pbar_G(e_p), phi(theta)> for phi in {cos, sin, cos2, cos3}, ||e_p - mean||_{L2},
  ||(T_h - sigma n) . n||_{L2(Gamma_h)} with T_h = D(u_h) n - p_h n (Lemma C1 form, D = grad + grad^T).
usage: python3 annulus.py DELTA A LAMHAT N1,N2,... [S1|S2]   (S2: psi = (r^2-1)^2 (1+x+y)/4)
"""
import os, sys, json, time
for v_ in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v_, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, sympy as sy
import svn

ROUT = 2.5


def annulus_mesh(N, a, delta, rout=ROUT):
    al = 2 * np.pi / N; ell = 2 * np.sin(np.pi / N)
    r = [1.0]
    while r[-1] * (1 + a * ell) < rout:
        r.append(r[-1] * (1 + a * ell))
    if rout - r[-1] < 0.5 * a * ell * r[-1]:
        r[-1] = rout
    else:
        r.append(rout)
    J = len(r) - 1
    verts = np.array([[rj * np.cos((i + j * delta) * al), rj * np.sin((i + j * delta) * al)]
                      for j, rj in enumerate(r) for i in range(N)])
    vid = lambda j, i: j * N + (i % N)
    tris = []
    for j in range(J):
        for i in range(N):
            tris.append([vid(j, i), vid(j, i + 1), vid(j + 1, i + 1)])
            tris.append([vid(j, i), vid(j + 1, i + 1), vid(j + 1, i)])
    tris = np.array(tris)
    P = verts[tris]; ar = (P[:, 1, 0] - P[:, 0, 0]) * (P[:, 2, 1] - P[:, 0, 1]) - (P[:, 1, 1] - P[:, 0, 1]) * (P[:, 2, 0] - P[:, 0, 0])
    tris[ar < 0] = tris[ar < 0][:, [0, 2, 1]]
    circ = set(range(N)); outer = set(range(J * N, (J + 1) * N))
    return verts, tris, circ, outer, ell, J


def exact_S(kind="S1"):
    x, y = sy.symbols('x y', real=True); r2 = x ** 2 + y ** 2
    psi = (r2 - 1) ** 2 * ((1 + x) if kind == "S1" else (1 + x + y)) / 4
    u = [sy.expand(sy.diff(psi, y)), sy.expand(-sy.diff(psi, x))]
    p = x ** 2 * y - y + x / 2
    f = [sy.expand(-(sy.diff(u[i], x, 2) + sy.diff(u[i], y, 2)) + sy.diff(p, [x, y][i])) for i in range(2)]
    gu = [[sy.diff(u[i], v) for v in (x, y)] for i in range(2)]
    U = sy.lambdify((x, y), u, 'numpy'); Fl = sy.lambdify((x, y), f, 'numpy')
    GU = sy.lambdify((x, y), gu, 'numpy'); P = sy.lambdify((x, y), p, 'numpy')
    vec = lambda fn: (lambda X, Y: np.stack([np.broadcast_to(np.asarray(v, float), np.shape(X)) for v in fn(X, Y)], -1))
    gmat = lambda X, Y: np.stack([np.stack([np.broadcast_to(np.asarray(GU(X, Y)[i][j], float), np.shape(X))
                                            for j in range(2)], -1) for i in range(2)], -2)
    return dict(u=vec(U), f=vec(Fl), gu=gmat, p=lambda X, Y: np.broadcast_to(np.asarray(P(X, Y), float), np.shape(X)))


def run(N, a, delta, lamhat, kind="S1"):
    t0 = time.time()
    ex = exact_S(kind)
    verts, tris, circ, outer, ell, J = annulus_mesh(N, a, delta)
    S = svn.Space(verts, tris, circ, outer)
    # Dirichlet nodes: all P4 nodes on the outer ring edges
    dn = set()
    for t in range(len(tris)):
        g = tris[t]
        for (u_, w_) in [(0, 1), (1, 2), (2, 0)]:
            if g[u_] in outer and g[w_] in outer:
                for k, (l1, l2, l0) in enumerate(svn.LBARY):
                    if {0: l0, 1: l1, 2: l2}[3 - u_ - w_] == 0:
                        dn.add(S.ids[t, k])
    S.dnodes = np.array(sorted(dn))
    hG = S.hGamma; S.h = 1.0
    mu = lamhat / ell
    sysm = svn.assemble(S, mu, ex["f"], None)
    A, B, F = sysm["A"], sysm["B"], sysm["F"]
    Bt = svn.coupling(sysm, 1)
    ndof = A.shape[0]; nt = len(tris)
    dD = svn.dofs(S.dnodes).ravel(); free = np.setdiff1d(np.arange(ndof), dD)
    ug = np.zeros(ndof); ug[dD] = ex["u"](S.xy[S.dnodes, 0], S.xy[S.dnodes, 1]).ravel()
    K = sp.bmat([[A[free][:, free], -Bt[:, free].T], [B[:, free], None]], format="csc")
    rhs = np.concatenate([F[free] - A[free][:, dD] @ ug[dD], -B[:, dD] @ ug[dD]])
    lu = spl.splu(K, permc_spec="COLAMD", diag_pivot_thresh=1.0)
    x = lu.solve(rhs)
    for _ in range(3):
        x = x + lu.solve(rhs - K @ x)
    relres = float(np.linalg.norm(K @ x - rhs) / np.linalg.norm(rhs))
    u = ug.copy(); u[free] = x[:len(free)]; p = x[len(free):]
    X = S.phys(svn.QX, svn.QY); wq = svn.QW[None, :] * np.abs(S.detJ)[:, None]
    ED = S.edge_data()
    num = per = 0.0
    for (t, ref, Xe, we, n) in ED:
        num += (we * ex["p"](Xe[:, 0], Xe[:, 1])).sum(); per += we.sum()
    pbarex = num / per
    pst = ex["p"](X[..., 0], X[..., 1]) - pbarex
    ep = pst - np.einsum('qj,tj->tq', svn.Q3_Q, p.reshape(nt, 10))
    area = wq.sum(); epm = (wq * ep).sum() / area
    rec = dict(kind=kind, N=N, a=a, delta=delta, lamhat=lamhat, layers=J, hG=hG, ell=ell, mu=mu, relres=relres,
               ep_L2=float(np.sqrt((wq * (ep - epm) ** 2).sum())))
    # boundary quantities
    eds = []; epG = 0.0; tr2 = 0.0; trn2 = 0.0
    for (t, ref, Xe, we, n) in ED:
        q3 = svn._mono(svn.MON3, ref[:, 0], ref[:, 1])
        ph_ = q3 @ p[10 * t:10 * t + 10]
        epe = ex["p"](Xe[:, 0], Xe[:, 1]) - pbarex - ph_
        dX, dY = svn.dbasis(ref[:, 0], ref[:, 1]); Ji = S.Jinv[t]
        Ut = u[svn.dofs(S.ids[t])]
        gxe = (Ji[0, 0] * dX + Ji[1, 0] * dY) @ Ut; gye = (Ji[0, 1] * dX + Ji[1, 1] * dY) @ Ut
        Gh = np.stack([gxe, gye], -1)                     # [g, c, d] = d_d u_c
        Ge = ex["gu"](Xe[:, 0], Xe[:, 1])
        eD = (Ge - Gh) + np.transpose(Ge - Gh, (0, 2, 1))
        trerr = np.einsum('gcd,d->gc', eD, n) - epe[:, None] * n[None, :]   # (sigma(u~,p*) - T_h) n  (up to pbar shift)
        tr2 += (we * (trerr ** 2).sum(-1)).sum(); trn2 += (we * (trerr @ n) ** 2).sum()
        eds.append((Xe, we, epe)); epG += (we * epe).sum()
    epG /= per
    rec["trac_L2"] = float(np.sqrt(tr2)); rec["tracn_L2"] = float(np.sqrt(trn2))
    for name, f in {"cos1": np.cos, "sin1": np.sin, "cos2": lambda t: np.cos(2 * t), "cos3": lambda t: np.cos(3 * t)}.items():
        Pm = sum((we * (epe - epG) * f(np.arctan2(Xe[:, 1], Xe[:, 0]))).sum() for (Xe, we, epe) in eds)
        rec["P_" + name] = float(Pm); rec["P_" + name + "/hG"] = float(Pm / hG)
    rec["ndof"] = int(K.shape[0]); rec["secs"] = round(time.time() - t0, 1)
    return rec


if __name__ == "__main__":
    delta, a, lamhat = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    for N in [int(s) for s in sys.argv[4].split(",")]:
        r = run(N, a, delta, lamhat, sys.argv[5] if len(sys.argv) > 5 else "S1")
        print(json.dumps({k: (float(f"{v:.6g}") if isinstance(v, float) else v) for k, v in r.items()}), flush=True)
