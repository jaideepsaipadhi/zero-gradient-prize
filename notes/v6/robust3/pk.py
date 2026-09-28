"""
Local condition (P_k) for the pressure-blind star certificate (REPORT.tex, Sec. 4).

For a boundary star om_z (the triangles at a vertex z of Gamma_h; e, e' its two Gamma_h edges) build
  V#(om_z) = { v in V_h(om_z): v = 0 on d(om_z) \ Gamma_h, B*(q, v) = 0 for all q in P_{k-1}(T), T in om_z,
               int_{Gamma} v.n = 0, int_{e u e'} v = 0, int_{e u e'} d_n v = 0 }
and the scale-invariant quadratic forms (s = |e_z|, w ranges over V_h|_om modulo constants):
  |grad v|^2,  Phi0(v)^2 = sup_w N0(w,v)^2/|grad w|^2,  Phi1(v)^2 = sup_w (s^-1 <w,v>_Gamma)^2/|grad w|^2,
  N0(w,v) = a_h(w,v) - <d_n w, v> - <w, d_n v>.     PhiQ^2 = sum of the three.
sigma_z(H)^2 = sup_v <(I-pi)H, v.n>^2 / (s^{2k+2} PhiQ(v)^2) for H in Hom_k (Frobenius norm of its tensor).
Reports the eigenvalues of that 5x5 (k+1 x k+1) form and
  kappa_z = min { sigma_z(H) : dist_F(H, R (n_e.x)^k) = 1 }   (Schur complement)  -- (P_k) asks kappa_z >= kappa_P > 0.

usage: python pk.py mesh N1,N2,...      production stars (all boundary vertices)
       python pk.py flat                exactly collinear e, e' (the h -> 0 limit shapes) + perturbations
       python pk.py random n            random admissible 3- and 4-triangle stars
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../code"))
import numpy as np
import scipy.linalg as sl
from math import comb
import svn_k as K
import svn


def star_forms(P, fan, k=4):
    """P: (nv,2) vertices, P[0] = z; fan: list of triangles (index triples) all containing 0.
    Gamma_h edges: (0, g1) and (0, g2) given as the first/last outer vertex of the fan (fan is ordered).
    Returns dict with M (k+1 x k+1, Frobenius coords), kappa, eigen info."""
    R = K.Ref(k)
    verts = np.asarray(P, float); tris = np.array(fan)
    P_ = verts[tris]; d1 = P_[:, 1] - P_[:, 0]; d2 = P_[:, 2] - P_[:, 0]; ar = d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]
    tris[ar < 0] = tris[ar < 0][:, [0, 2, 1]]
    g1, g2 = GAMMA
    circ = {0, g1, g2}
    S = K.Space(verts, tris, circ, set(), R)
    # Gamma_h edges are exactly (0,g1), (0,g2): drop any other edge with both ends in circ (none if g1,g2 not adjacent)
    S.bedges = [(t, u, w) for (t, u, w) in S.bedges if 0 in (S.tris[t, u], S.tris[t, w])]
    s = np.linalg.norm(verts[g1] - verts[0])
    nn = S.nn; ndof = 2 * nn
    # boundary edges of the star mesh that are not Gamma_h edges -> fixed nodes
    from collections import Counter
    ecount = Counter()
    for t in range(len(tris)):
        for (a, b) in [(0, 1), (1, 2), (2, 0)]:
            ecount[tuple(sorted((tris[t, a], tris[t, b])))] += 1
    fixed = set()
    for t in range(len(tris)):
        for kk, (l1, l2, l0) in enumerate(R.LBARY):
            lam = {0: l0, 1: l1, 2: l2}
            nz = [tris[t, v] for v in range(3) if lam[v] > 0]
            # node lies on the closed edge between the vertices with lam=0 complement...
            for (a, b) in [(0, 1), (1, 2), (2, 0)]:
                ga, gb = tris[t, a], tris[t, b]
                c = 3 - a - b
                if lam[c] == 0 and ecount[tuple(sorted((ga, gb)))] == 1 and not (0 in (ga, gb) and ({ga, gb} & {g1, g2})):
                    fixed.add(int(S.ids[t, kk]))
    free = np.array(sorted(set(range(nn)) - fixed))
    fdof = np.array([[2 * g, 2 * g + 1] for g in free]).ravel()
    sysm = K.assemble(S, 0.0, lambda X, Y: np.zeros(np.shape(X) + (2,)))
    N0 = sysm["A"].toarray()                       # mu = 0: a_h - <dn.,.> - <.,dn .>
    # boundary mass (unscaled) and edge data
    Mb = np.zeros((ndof, ndof)); rows_pb = []
    Mv = np.zeros((2, ndof)); Dv = np.zeros((2, ndof)); flux = np.zeros(ndof)
    edges = []
    for (t, ref, Xe, we, n) in S.edge_data():
        ph = R.basis(ref[:, 0], ref[:, 1]); dX, dY = R.dbasis(ref[:, 0], ref[:, 1]); Ji = S.Jinv[t]
        dn = n[0] * (Ji[0, 0] * dX + Ji[1, 0] * dY) + n[1] * (Ji[0, 1] * dX + Ji[1, 1] * dY)
        d = K.dofs(S.ids[t])
        Mphi = np.einsum('g,ga,gb->ab', we, ph, ph)
        for c in range(2):
            Mb[np.ix_(d[:, c], d[:, c])] += Mphi
            Mv[c, d[:, c]] += we @ ph; Dv[c, d[:, c]] += we @ dn
            flux[d[:, c]] += (we @ ph) * n[c]
        edges.append((t, ref, Xe, we, n, ph, d))
    Bm = (sysm["B"] - sysm["C"]).toarray()
    Cm = np.vstack([Bm[:, fdof], flux[None, fdof], Mv[:, fdof], Dv[:, fdof]])
    ns = sl.null_space(Cm, rcond=1e-10)
    # plain-grad Gram on all dofs
    gx, gy = S.grads(R.QX, R.QY); wq = R.QW[None, :] * np.abs(S.detJ)[:, None]
    Gt = np.einsum('tq,tqa,tqb->tab', wq, gx, gx) + np.einsum('tq,tqa,tqb->tab', wq, gy, gy)
    G = np.zeros((ndof, ndof))
    for t in range(len(tris)):
        d = K.dofs(S.ids[t])
        for c in range(2):
            G[np.ix_(d[:, c], d[:, c])] += Gt[t]
    Gp = np.linalg.pinv(G, rcond=1e-11, hermitian=True)
    A0 = N0[:, fdof] @ ns; A1 = (Mb[:, fdof] / s) @ ns
    Gv = ns.T @ G[np.ix_(fdof, fdof)] @ ns
    Q0 = A0.T @ Gp @ A0; Q1 = A1.T @ Gp @ A1
    Q = Gv + Q0 + Q1
    # consistency: constants annihilated
    one = np.zeros((ndof, 2)); one[0::2, 0] = 1; one[1::2, 1] = 1
    cres = max(np.abs(one.T @ N0[:, fdof] @ ns).max(), np.abs(one.T @ Mb[:, fdof] @ ns).max() / s)
    # pairings with H_j = ((x-z)/s)^j ((y-z)/s)^(k-j), divided by s
    z = verts[0]
    Pi = np.zeros((k + 1, ns.shape[1]))
    for j in range(k + 1):
        f = lambda X, Y, j=j: ((X - z[0]) / s)**j * ((Y - z[1]) / s)**(k - j)
        for (t, ref, Xe, we, n, ph, d) in edges:
            Xq = S.phys(R.QX, R.QY)[t]; wqt = R.QW * abs(S.detJ[t])
            Mp = np.einsum('q,qi,qj->ij', wqt, R.QP_Q, R.QP_Q)
            cpi = np.linalg.solve(Mp, R.QP_Q.T @ (wqt * f(Xq[:, 0], Xq[:, 1])))
            eta = f(Xe[:, 0], Xe[:, 1]) - R.pbasis(ref[:, 0], ref[:, 1]) @ cpi
            wv = ph.T @ (we * eta)
            vec = np.zeros(ndof)
            for c in range(2):
                vec[d[:, c]] += wv * n[c]
            Pi[j] += (vec[fdof] @ ns) / s
    Sm = Pi @ np.linalg.solve(Q, Pi.T)
    D = np.diag([np.sqrt(comb(k, j)) for j in range(k + 1)])
    M = D @ Sm @ D
    lam, U = np.linalg.eigh(M)
    # normal power of e_z = (0,g1)
    tv = (verts[g1] - z) / s; n = np.array([tv[1], -tv[0]])
    ell = np.array([comb(k, j) * n[0]**j * n[1]**(k - j) for j in range(k + 1)])
    ell = np.linalg.solve(D, ell); ell /= np.linalg.norm(ell)
    # Schur complement kappa^2 = min over c perp ell of min_a (a ell + c)^T M (a ell + c)
    B = sl.null_space(ell[None, :])
    Mll = ell @ M @ ell; Mlb = ell @ M @ B; Mbb = B.T @ M @ B
    Sch = Mbb - np.outer(Mlb, Mlb) / Mll if Mll > 1e-14 * np.trace(M) else Mbb
    kappa = float(np.sqrt(max(np.linalg.eigvalsh(Sch).min(), 0)))
    align = float(abs(U[:, 0] @ ell))
    return dict(dimV=int(ns.shape[1]), sig=[float(np.sqrt(max(x, 0))) for x in lam], kappa=kappa,
                align_min_eigvec_normal=align, const_res=float(cres))


GAMMA = (None, None)


def production_star(N, i, L=2.5):
    v, t, circ, outer = svn.make_mesh(N, L=L)
    return v, t


def mesh_stars(N, L=2.5):
    global GAMMA
    v, t, circ, outer = svn.make_mesh(N, L=L)
    m = N // 4
    out = []
    for i in range(N):
        z = i                      # vid[0, i] = i
        fan = [tr for tr in t if z in tr]
        loc = sorted({int(a) for tr in fan for a in tr})
        mp = {z: 0}
        for a in loc:
            if a != z:
                mp[a] = len(mp)
        P = np.zeros((len(mp), 2))
        for a, b in mp.items():
            P[b] = v[a]
        g1 = mp[(i + 1) % N]; g2 = mp[(i - 1) % N]
        GAMMA = (g1, g2)
        fl = [[mp[int(a)] for a in tr] for tr in fan]
        r = star_forms(P, fl)
        ang = np.degrees(np.arccos(np.clip(np.dot(P[g1] - P[0], P[g2] - P[0]) /
                                           (np.linalg.norm(P[g1] - P[0]) * np.linalg.norm(P[g2] - P[0])), -1, 1)))
        r.update(N=N, i=i, angle_deg=float(ang))
        out.append(r)
    return out


def flat_star(a=1.0, H=1.2, x3=0.0, delta=0.0, xs=(1.0, -1.0)):
    """z = 0, e = [0,(1,0)], e' = [0,(-1,0)] rotated by +-delta/2 (delta = 0: exactly collinear);
    outer vertices at height H: p1 above e's far end (x = xs[0]), p2 above z shifted by x3, p3 above e' far end."""
    global GAMMA
    c, s_ = np.cos(delta / 2), np.sin(delta / 2)
    g1 = np.array([c, -s_]) * a; g2 = np.array([-c, -s_]) * a
    # interior on the +y side
    P = np.array([[0, 0], g1, [xs[0], H], [x3, H], [xs[1], H], g2], float)
    fan = [[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5]]
    GAMMA = (1, 5)
    return star_forms(P, fan)


if __name__ == "__main__":
    if sys.argv[1] == "mesh":
        for N in [int(x) for x in sys.argv[2].split(",")]:
            rs = mesh_stars(N)
            ks = [r["kappa"] for r in rs]; s5 = [r["sig"][0] for r in rs]; s4 = [r["sig"][1] for r in rs]
            print(json.dumps(dict(N=N, nstars=len(rs), dimV=sorted({r["dimV"] for r in rs}),
                                  kappa_min=min(ks), kappa_max=max(ks), sig_min_max=max(s5), sig_min_min=min(s5),
                                  sig2_min=min(s4), align_min=min(r["align_min_eigvec_normal"] for r in rs),
                                  angle=[min(r["angle_deg"] for r in rs), max(r["angle_deg"] for r in rs)],
                                  const_res=max(r["const_res"] for r in rs), sample=rs[0])), flush=True)
    elif sys.argv[1] == "flat":
        for kw in [dict(), dict(x3=0.3), dict(H=0.6), dict(H=2.0, x3=-0.4), dict(xs=(1.0, -0.8)),
                   dict(delta=0.05), dict(delta=0.2), dict(delta=-0.2)]:
            r = flat_star(**kw); r.update(kw)
            print(json.dumps(r), flush=True)


def flat3(H=1.0, delta=0.0, xw=(1.0, 0.0), k=4):
    """production-type 3-triangle star, flat limit: z=0, e=[0,(1,0)], e'=[0,(-1,0)] (rotated by delta/2 each),
    layer-1 vertices w_{i+1} = (xw[0], H), w_i = (xw[1], H); triangles (z,g1,w1), (z,w1,w0), (g2,z,w0)."""
    global GAMMA
    c, s_ = np.cos(delta / 2), np.sin(delta / 2)
    P = np.array([[0, 0], [c, -s_], [xw[0], H], [xw[1], H], [-c, -s_]], float)
    GAMMA = (1, 4)
    return star_forms(P, [[0, 1, 2], [0, 2, 3], [4, 0, 3]], k)


def kernel_dirs(r_M=None):
    pass


def random_flat(n, seed=1):
    """random 3- and 4-triangle stars with exactly collinear Gamma_h edges (delta = 0) and with delta small."""
    rng = np.random.default_rng(seed)
    out = []
    for i in range(n):
        H = rng.uniform(0.6, 1.6)
        if i % 2 == 0:
            xw = (rng.uniform(0.5, 1.3), rng.uniform(-0.4, 0.4))
            r = flat3(H=H, xw=xw); r.update(kind="3tri", H=H, xw=xw)
        else:
            x3 = rng.uniform(-0.4, 0.4); xs = (rng.uniform(0.7, 1.3), -rng.uniform(0.7, 1.3))
            r = flat_star(H=H, x3=x3, xs=xs); r.update(kind="4tri", H=H, x3=x3, xs=xs)
        out.append(r)
    return out


if __name__ == "__main__" and sys.argv[1] == "random":
    rs = random_flat(int(sys.argv[2]))
    for kind in ("3tri", "4tri"):
        q = [r for r in rs if r["kind"] == kind]
        print(json.dumps(dict(kind=kind, n=len(q), dimV=sorted({r["dimV"] for r in q}),
                              sig1_max=max(r["sig"][0] for r in q), sig2_min=min(r["sig"][1] for r in q),
                              kappa_min=min(r["kappa"] for r in q), align_min=min(r["align_min_eigvec_normal"] for r in q))))
    worst = min(rs, key=lambda r: r["kappa"]); print(json.dumps(worst))
