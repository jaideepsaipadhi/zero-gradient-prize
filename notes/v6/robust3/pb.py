"""
Pressure-blind lower-bound certificate for ||grad W_*(phi)||  (notes/v6/robust3/REPORT.tex).

Pressure-blind test fields: v in V_R with B*(q,v) = 0 for all q in Pi_h.  For such v and the CNS* velocity W
of the hydrostatic data (grad phi, 0):   N_h(W, v) = <eta, v.n>,   eta = phi - pi_h phi on Gamma_h.
Patch space V#(om): v supported in the patch om (v = 0 on the part of d(om) off Gamma_h), pressure-blind,
  sum_{e in om} int_e v = 0,  sum_{e in om} int_e d_n v = 0   (mode "sum")   or per edge (mode "edge").
Then N_h(., v) kills constants on om and |<eta, v.n>| <= ||l_v||_* ||grad W||_om  (l_v(w) = N_h(w, v)),
with ||l_v||_* the exact dual norm over V_h|_om modulo constants.

usage: python pb.py dims N [k]                  -> dim V#(star) per boundary vertex, pairing rank on Hom_k
       python pb.py cert N phi mu1,mu2 [k] [mode] [solve]  -> certificate LB vs measured ||grad W||
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../code"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../robust2"))
import numpy as np
import scipy.linalg as sl
import svn_k as K
import leak as LK


def setup(N, k, phi, gphi, mu):
    S = K.build(N, k)
    sysm = K.assemble(S, mu, gphi)
    E = LK.edge_info(S, phi)
    R = S.R
    gx, gy = S.grads(R.QX, R.QY)
    wq = R.QW[None, :] * np.abs(S.detJ)[:, None]
    Gt = np.einsum('tq,tqa,tqb->tab', wq, gx, gx) + np.einsum('tq,tqa,tqb->tab', wq, gy, gy)   # plain grad stiffness
    # per-edge extra data: dn basis, edge index of triangle
    ed = {}
    for (e, (t, ref, Xe, we, n)) in zip(E, S.edge_data()):
        dX, dY = R.dbasis(ref[:, 0], ref[:, 1]); Ji = S.Jinv[t]
        gxe = Ji[0, 0] * dX + Ji[1, 0] * dY; gye = Ji[0, 1] * dX + Ji[1, 1] * dY
        e["dn"] = n[0] * gxe + n[1] * gye
        ed[t] = e
    # node -> triangles
    n2t = {}
    for t in range(len(S.tris)):
        for g in S.ids[t]:
            n2t.setdefault(int(g), set()).add(t)
    return S, sysm, E, ed, Gt, n2t


def patch_space(S, sysm, ed, n2t, tris, mode="sum"):
    """basis (global dof indices, matrix Vb) of V#(patch)."""
    tris = set(tris)
    nodes = sorted({int(g) for t in tris for g in S.ids[t] if n2t[int(g)] <= tris})
    dof = np.array([[2 * g, 2 * g + 1] for g in nodes]).ravel()
    alld = np.array(sorted({int(d) for t in tris for d in K.dofs(S.ids[t]).ravel()}))
    npl = S.R.npl
    rows = []
    Bm = (sysm["B"] - sysm["C"]).tocsr()
    for t in sorted(tris):
        rows.append(Bm[t * npl:(t + 1) * npl][:, dof].toarray())
    flux = np.zeros(len(dof)); Mv = np.zeros((2, len(dof))); Dv = np.zeros((2, len(dof)))
    perrows = []
    pos = {d: i for i, d in enumerate(dof)}
    for t in sorted(tris):
        if t not in ed:
            continue
        e = ed[t]
        me = e["w"] @ e["ph"]; md = e["w"] @ e["dn"]
        Me = np.zeros((2, len(dof))); De = np.zeros((2, len(dof)))
        for a, g in enumerate(e["d"]):
            for c in range(2):
                if int(g[c]) in pos:
                    j = pos[int(g[c])]
                    flux[j] += me[a] * e["n"][c]
                    Me[c, j] += me[a]; De[c, j] += md[a]
        Mv += Me; Dv += De; perrows += [Me, De]
    rows.append(flux[None, :])
    if mode == "sum":
        rows += [Mv, Dv]
    else:
        rows += perrows
    Cm = np.vstack(rows)
    if Cm.shape[1] == 0:
        return dof, alld, np.zeros((0, 0))
    ns = sl.null_space(Cm, rcond=1e-10)
    return dof, alld, ns


def dual_matrix(S, sysm, Gt, tris, dof, alld, mu):
    """Q with v^T Q v = ||l_v||_*^2 = sup_w N_h(w,v)^2 / ||grad w||_om^2 (w in V_h|_om mod constants)."""
    A = sysm["A"]                                       # full N_h at the assembly mu
    Nm = A[alld][:, dof].toarray()
    # local plain-gradient Gram over patch triangles only; but a_h(w,v) involves only tris in supp v = patch
    pos = {d: i for i, d in enumerate(alld)}
    G = np.zeros((len(alld), len(alld)))
    for t in tris:
        d = K.dofs(S.ids[t])
        for c in range(2):
            ix = [pos[int(x)] for x in d[:, c]]
            G[np.ix_(ix, ix)] += Gt[t]
    # a_h in sysm["A"] also contains contributions of triangles outside patch between alld dofs and dof:
    # supp(v) is inside the patch, so those vanish.  Pseudo-inverse modulo constants:
    Gp = np.linalg.pinv(G, rcond=1e-11, hermitian=True)
    return Nm.T @ Gp @ Nm, Nm, G


MU0 = [None]


def pairing_vec(E, ed, tris, dof):
    pos = {d: i for i, d in enumerate(dof)}
    p = np.zeros(len(dof))
    for t in tris:
        if t not in ed:
            continue
        e = ed[t]
        wv = e["ph"].T @ (e["w"] * e["eta"])
        for a, g in enumerate(e["d"]):
            for c in range(2):
                if int(g[c]) in pos:
                    p[pos[int(g[c])]] += wv[a] * e["n"][c]
    return p


def stars(S):
    circ_tris = {}
    bverts = set()
    for (t, u, w) in S.bedges:
        bverts |= {int(S.tris[t, u]), int(S.tris[t, w])}
    out = {}
    for t in range(len(S.tris)):
        for v in S.tris[t]:
            if int(v) in bverts:
                out.setdefault(int(v), []).append(t)
    return out


def eta_on_edge(S, e, f):
    """eta = f - pi_T f on the Gamma_h edge record e (f vectorised in X,Y)."""
    R = S.R; t = e["t"]
    Xq = S.phys(R.QX, R.QY)[t]; wq = R.QW * abs(S.detJ[t])
    M = np.einsum('q,qi,qj->ij', wq, R.QP_Q, R.QP_Q)
    cpi = np.linalg.solve(M, R.QP_Q.T @ (wq * f(Xq[:, 0], Xq[:, 1])))
    Rr = np.array([[0, 0], [1, 0], [0, 1]], float)
    # recover ref points of the edge from X via inverse map
    ref = np.linalg.solve(S.J[t], (e["X"] - S.P0[t]).T).T
    return f(e["X"][:, 0], e["X"][:, 1]) - R.pbasis(ref[:, 0], ref[:, 1]) @ cpi


def hom_pairings(S, ed, tris, dof, z, k):
    """matrix P (k+1, ndof): row j = <(I-pi)H_j, v.n> with H_j = (x-z)^j (y-z)^(k-j) / h^k."""
    pos = {d: i for i, d in enumerate(dof)}
    zc = S.verts[z]; hs = S.hGamma
    P = np.zeros((k + 1, len(dof)))
    for j in range(k + 1):
        f = lambda X, Y, j=j: ((X - zc[0]) / hs)**j * ((Y - zc[1]) / hs)**(k - j)
        for t in tris:
            if t not in ed:
                continue
            e = ed[t]
            et = eta_on_edge(S, e, f)
            wv = e["ph"].T @ (e["w"] * et)
            for a, g in enumerate(e["d"]):
                for c in range(2):
                    if int(g[c]) in pos:
                        P[j, pos[int(g[c])]] += wv[a] * e["n"][c]
    return P


def dims(N, k=4, mode="sum", nshow=4):
    phi, gphi = LK.PHIS["sin"]
    MU0[0] = 100.0
    S, sysm, E, ed, Gt, n2t = setup(N, k, phi, gphi, 100.0)
    st = stars(S)
    res = []
    for i, (z, tris) in enumerate(sorted(st.items())):
        dof, alld, ns = patch_space(S, sysm, ed, n2t, tris, mode)
        dimV = ns.shape[1]
        P = hom_pairings(S, ed, tris, dof, z, k) @ ns if dimV else np.zeros((k + 1, 0))
        sv = np.linalg.svd(P, compute_uv=False) if dimV else np.array([])
        rank = int((sv > 1e-9 * max(sv.max() if len(sv) else 1, 1e-300)).sum()) if dimV else 0
        # kernel of H -> pairing, compare with edge normal powers of the patch's boundary triangles
        info = dict(z=z, m=len(tris), ndof=len(dof), dimV=dimV, rank=rank, sv=[float(s) for s in sv[:k + 1]])
        if dimV:
            U, s_, Vt = np.linalg.svd(P.T)       # rows of P correspond to H_j
            ker = Vt[rank:]                      # coefficient vectors of H in kernel
            # test: is (n.x)^k in the kernel for the Gamma_h edges' normals?
            tests = []
            for t in tris:
                if t in ed:
                    n = ed[t]["n"]
                    from math import comb
                    cvec = np.array([comb(k, j) * n[0]**j * n[1]**(k - j) for j in range(k + 1)])
                    tests.append(float(np.linalg.norm(P.T @ cvec) / np.linalg.norm(cvec) / max(sv.max(), 1e-300)))
            info["normal_power_residual"] = tests
        res.append(info)
        if i < nshow:
            print(json.dumps(info), flush=True)
    dv = [r["dimV"] for r in res]; rk = [r["rank"] for r in res]
    print(json.dumps(dict(N=N, k=k, mode=mode, dimV_min=min(dv), dimV_max=max(dv), rank_min=min(rk), rank_max=max(rk))), flush=True)


def dims_tri(N, k):
    phi, gphi = LK.PHIS["sin"]
    MU0[0] = 100.0
    S, sysm, E, ed, Gt, n2t = setup(N, k, phi, gphi, 100.0)
    out = []
    for t in sorted(ed)[:3]:
        dof, alld, ns = patch_space(S, sysm, ed, n2t, [t], "edge")
        z = int(S.tris[t][0])
        dimV = ns.shape[1]
        if dimV:
            P = hom_pairings(S, ed, [t], dof, z, k) @ ns
            sv = np.linalg.svd(P, compute_uv=False)
            n = ed[t]["n"]
            from math import comb
            cvec = np.array([comb(k, j) * n[0]**j * n[1]**(k - j) for j in range(k + 1)])
            res_n = float(np.linalg.norm(P.T @ cvec) / np.linalg.norm(cvec) / sv.max())
        else:
            sv = []; res_n = None
        out.append(dict(k=k, t=t, ndof=len(dof), dimV=dimV, sv=[float(s) for s in sv], normal_res=res_n))
        print(json.dumps(out[-1]), flush=True)


def cert(N, name, mus, k=4, mode="sum", do_solve=True, check=True):
    phi, gphi = LK.PHIS[name]
    S = K.build(N, k)
    E = LK.edge_info(S, phi)
    eta_n = np.sqrt(sum(e["w"] @ e["eta"]**2 for e in E))
    out = []
    for mu in mus:
        t0 = time.time()
        S2, sysm, E2, ed, Gt, n2t = setup(N, k, phi, gphi, mu)
        st = stars(S2)
        zs = sorted(st)
        LB2 = {}; vbest = {}
        for z in zs:
            tris = st[z]
            dof, alld, ns = patch_space(S2, sysm, ed, n2t, tris, mode)
            Q, Nm, G = dual_matrix(S2, sysm, Gt, tris, dof, alld, mu)
            p = pairing_vec(E2, ed, tris, dof)
            Qv = ns.T @ Q @ ns; pv = ns.T @ p
            c = np.linalg.lstsq(Qv, pv, rcond=1e-12)[0]
            LB2[z] = float(pv @ c)
            vbest[z] = (dof, ns @ c)
        allLB = np.sqrt(sum(LB2.values()) / 2)
        # disjoint parity sub-families (stars at every other boundary vertex, in angular order)
        ang = sorted(zs, key=lambda z: np.arctan2(S2.verts[z][1], S2.verts[z][0]))
        par = [np.sqrt(sum(LB2[z] for z in ang[i::2])) for i in (0, 1)]
        gam = mu * S2.hGamma / S2.h
        norm = S2.hGamma**0.5 / gam * eta_n
        rec = dict(N=N, k=k, phi=name, mu=mu, gamma=gam, mode=mode, LB=float(max(allLB, *par)), LB_all=float(allLB),
                   LB_par=[float(x) for x in par], LB_n=float(max(allLB, *par) / norm), eta=eta_n, hG=S2.hGamma)
        if do_solve:
            Bt = K.coupling(S2, sysm, 1)
            u, pr, rel = LK.solve(S2, sysm["A"], Bt, sysm["B"], sysm["F"])
            Er = K.errors(S2, u, LK.ZERO, mu)
            rec.update(H1=Er["H1"], H1_n=Er["H1"] / norm, ratio=Er["H1"] / rec["LB"], relres=rel)
            if check:
                # identity N_h(W, v) = <eta, v.n> for the assembled best v (sum over stars of one parity)
                v = np.zeros(2 * S2.nn)
                for z in ang[0::2]:
                    dof, vv = vbest[z]; v[dof] += vv
                lhs = u @ (sysm["A"] @ v)
                rhs = sum(e["w"] @ (e["eta"] * ((e["ph"] @ v[e["d"]]) @ e["n"])) for e in E2)
                pbres = np.abs(Bt @ v).max() / max(np.abs(sysm["B"]).max() * np.abs(v).max(), 1e-300)
                rec.update(identity_lhs=float(lhs), identity_rhs=float(rhs), pb_residual=float(pbres))
        rec["time"] = time.time() - t0
        print(json.dumps(rec), flush=True)
        out.append(rec)
    return out


if __name__ == "__main__":
    if sys.argv[1] == "cert":
        N = int(sys.argv[2]); name = sys.argv[3]; mus = [float(x) for x in sys.argv[4].split(",")]
        k = int(sys.argv[5]) if len(sys.argv) > 5 else 4
        mode = sys.argv[6] if len(sys.argv) > 6 else "sum"
        cert(N, name, mus, k, mode, do_solve=(len(sys.argv) <= 7 or sys.argv[7] != "nosolve"))
    elif sys.argv[1] == "dims":
        N = int(sys.argv[2]); k = int(sys.argv[3]) if len(sys.argv) > 3 else 4
        mode = sys.argv[4] if len(sys.argv) > 4 else "sum"
        dims(N, k, mode)


def patch_forms(S, sys0, ed, n2t, Gt, tris, z, k, mode="sum"):
    """scale-invariant 5x5 form M (Frobenius coords) of sigma_om(H)^2 for a patch (as pk.star_forms), in the
    global mesh; sys0 assembled with mu = 0 so sys0['A'] = N0."""
    from math import comb
    dof, alld, ns = patch_space(S, sys0, ed, n2t, tris, mode)
    s = min(np.linalg.norm(S.verts[S.tris[t, u]] - S.verts[S.tris[t, w]]) for (t, u, w) in S.bedges if t in tris)
    N0 = sys0["A"][alld][:, dof].toarray() @ ns
    M1 = (sys0["Pen"] * S.h)[alld][:, dof].toarray() @ ns / s
    pos = {d: i for i, d in enumerate(alld)}
    G = np.zeros((len(alld), len(alld)))
    for t in tris:
        d = K.dofs(S.ids[t])
        for c in range(2):
            ix = [pos[int(x)] for x in d[:, c]]
            G[np.ix_(ix, ix)] += Gt[t]
    Gp = np.linalg.pinv(G, rcond=1e-11, hermitian=True)
    posd = [pos[int(d)] for d in dof]
    Q = ns.T @ G[np.ix_(posd, posd)] @ ns + N0.T @ Gp @ N0 + M1.T @ Gp @ M1
    zc = S.verts[z]
    Pi = np.zeros((k + 1, ns.shape[1]))
    for j in range(k + 1):
        f = lambda X, Y, j=j: ((X - zc[0]) / s)**j * ((Y - zc[1]) / s)**(k - j)
        posv = {d: i for i, d in enumerate(dof)}
        vec = np.zeros(len(dof))
        for t in tris:
            if t not in ed:
                continue
            e = ed[t]; et = eta_on_edge(S, e, f); wv = e["ph"].T @ (e["w"] * et)
            for a, g in enumerate(e["d"]):
                for c in range(2):
                    if int(g[c]) in posv:
                        vec[posv[int(g[c])]] += wv[a] * e["n"][c]
        Pi[j] = vec @ ns / s
    Sm = Pi @ np.linalg.lstsq(Q, Pi.T, rcond=1e-13)[0]
    D = np.diag([np.sqrt(comb(k, j)) for j in range(k + 1)])
    return D @ Sm @ D, ns.shape[1]


def patch_scan(N, k=4, kind="star2"):
    phi, gphi = LK.PHIS["sin"]
    S, sys0, E, ed, Gt, n2t = setup(N, k, phi, gphi, 0.0)
    st = stars(S)
    ang = sorted(st, key=lambda z: np.arctan2(S.verts[z][1], S.verts[z][0]))
    ev = []
    for i, z in enumerate(ang):
        if kind == "star":
            tris = st[z]
        elif kind == "star2":
            tris = sorted(set(st[z]) | set(st[ang[(i + 1) % len(ang)]]))
        elif kind == "star3":
            tris = sorted(set(st[z]) | set(st[ang[(i + 1) % len(ang)]]) | set(st[ang[(i - 1) % len(ang)]]))
        M, dimV = patch_forms(S, sys0, ed, n2t, Gt, tris, z, k)
        lam = np.linalg.eigvalsh(M)
        ev.append(np.sqrt(np.clip(lam, 0, None)))
    ev = np.array(ev)
    print(json.dumps(dict(N=N, kind=kind, dimV=dimV, sig_min=[float(x) for x in ev.min(0)], sig_max=[float(x) for x in ev.max(0)])), flush=True)
