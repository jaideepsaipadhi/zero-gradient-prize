"""Numerical checks for the H^1 lower bound (notes/v3/lower_bound.tex).

    python3 lower_bound_tests.py stars   16,32,64,128,256,512    # (M3): star constants kappa_z for every boundary vertex
    python3 lower_bound_tests.py field   16,32,64,128,256        # assembled test field v = sum_z W(z) v_z, tests A and B1
    python3 lower_bound_tests.py global  16,24,32,48,64          # exact dual norms on Z_h, Y_h, Y_h^1 and the CNS*/GS errors

Notation (see the note):
  Y_h   = Z_h cap V_h^0                           (discrete div-free, zero on Gamma_h and on dR)
  Y_h^1 = { v in Y_h : int_e d_n v . tau_e = 0 on every chord e }
  star constant  kappa_z = sup_{v in Y^1(omega_z)} |sum_e int_e s_e^2 d_n v . tau_e| / (|e_z|^2 ||grad v||),
  omega_z = the (three) triangles at the boundary vertex z, v = 0 on the boundary of omega_z.
All solves are small (<= 3 GB, one core).
"""
import sys, os, time
from collections import Counter
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, scipy.linalg as sla
import svn
from lemma_tests import gram_grad, edge_functional

MU = 100.0


# ------------------------------------------------------------------ helpers
def sub_space(N, nlayers):
    """first nlayers layers of the production mesh (same geometry as svn.make_mesh(N))"""
    v, t, c, o = svn.make_mesh(N)
    return svn.Space(v, t[:2 * N * nlayers], c, o)


def edge_tables(S):
    out = []
    for (t, ref, Xe, we, n) in S.edge_data():
        dX, dY = svn.dbasis(ref[:, 0], ref[:, 1]); Ji = S.Jinv[t]
        gxe = Ji[0, 0] * dX + Ji[1, 0] * dY; gye = Ji[0, 1] * dX + Ji[1, 1] * dY
        dn = n[0] * gxe + n[1] * gye
        L = we.sum(); s = (svn.GL - 0.5) * L
        out.append(dict(t=t, dn=dn, we=we, s=s, tau=np.array([-n[1], n[0]]), n=n, L=L, X=Xe))
    return out


def gamma_nodes(S):
    gam = set()
    for (t, u_, w_) in S.bedges:
        for k, (l1, l2, l0) in enumerate(svn.LBARY):
            lam = {0: l0, 1: l1, 2: l2}
            if lam[3 - u_ - w_] == 0:
                gam.add(S.ids[t, k])
    return gam


def dn_tau_row(S, E, weight=None):
    """row r with r.v = int_e weight(s) d_n v . tau_e"""
    r = np.zeros(2 * S.nn); d = svn.dofs(S.ids[E["t"]])
    w = E["we"] if weight is None else E["we"] * weight
    for c in range(2):
        np.add.at(r, d[:, c], E["tau"][c] * (w @ E["dn"]))
    return r


def star_space(S, B, H, ET, tri_set, gam):
    """orthonormal (in grad) basis of Y^1(omega) for the patch tri_set; returns (dofs, basis) """
    tri_set = np.array(sorted(tri_set))
    ecount = Counter()
    for t in tri_set:
        g = S.tris[t]
        for (a, b) in [(0, 1), (1, 2), (2, 0)]:
            ecount[tuple(sorted((g[a], g[b])))] += 1
    bnodes = set(gam)
    for t in tri_set:
        g = S.tris[t]
        for (a, b) in [(0, 1), (1, 2), (2, 0)]:
            if ecount[tuple(sorted((g[a], g[b])))] == 1:
                for k, (l1, l2, l0) in enumerate(svn.LBARY):
                    if {0: l0, 1: l1, 2: l2}[3 - a - b] == 0:
                        bnodes.add(S.ids[t, k])
    nodes = sorted(set(S.ids[tri_set].ravel()) - bnodes)
    dof = svn.dofs(np.array(nodes)).ravel()
    prow = (tri_set[:, None] * 10 + np.arange(10)[None, :]).ravel()
    rows = [B[prow][:, dof].toarray()]
    tset = set(tri_set.tolist())
    for E in ET:
        if E["t"] in tset:
            d = svn.dofs(S.ids[E["t"]]); r = Counter()
            for c in range(2):
                for a, val in zip(d[:, c], E["tau"][c] * (E["we"] @ E["dn"])):
                    r[a] += val
            pos = {g: j for j, g in enumerate(dof)}
            row = np.zeros(len(dof))
            for a, val in r.items():
                if a in pos:
                    row[pos[a]] += val
            rows.append(row[None, :])
    Z = sla.null_space(np.vstack(rows), rcond=1e-10)
    if Z.shape[1] == 0:
        return dof, Z
    G = Z.T @ H[dof][:, dof].toarray() @ Z
    w, U = np.linalg.eigh(G)
    return dof, Z @ (U / np.sqrt(w))           # columns have ||grad||=1 and are grad-orthogonal


def star_tris(N, i):
    """boundary star of circle vertex c_i in the structured mesh: T_a(i), T_b(i), T_a(i-1)"""
    return [2 * i, 2 * i + 1, 2 * ((i - 1) % N)]


def star_fields(N):
    """for every boundary vertex z_i: unit-gradient v_z in Y^1(omega_z) maximising |M_z|, M_z := -1/2 sum_e int s^2 d_n v.tau,
    oriented so that M_z > 0; returns S, ET, list of (v_z full vector, M_z/|e|^2 = kappa_z)"""
    S = sub_space(N, 1)
    sysm = svn.assemble(S, MU, lambda X, Y: np.zeros(X.shape + (2,)), 1)
    B = sysm["B"]; H = gram_grad(S).tocsr(); ET = edge_tables(S); gam = gamma_nodes(S)
    mrow = np.zeros(2 * S.nn)
    for E in ET:
        mrow += -0.5 * dn_tau_row(S, E, E["s"] ** 2)
    edge_of_tri = {E["t"]: E for E in ET}
    out = []
    for i in range(N):
        tri = star_tris(N, i)
        dof, Q = star_space(S, B, H, ET, tri, gam)
        m = Q.T @ mrow[dof]
        c = m / np.linalg.norm(m)
        vals = Q @ c
        M = mrow[dof] @ vals
        e0 = edge_of_tri[2 * i]["L"]
        out.append(((dof, vals), M / e0 ** 2, Q.shape[1]))     # v_z stored sparsely
    return S, ET, H, out


# ------------------------------------------------------------------ (M3): star constants
def cmd_stars(Ns):
    print("N    h_Gamma   dim Y1(star)   min kappa_z   mean kappa_z   max kappa_z")
    for N in Ns:
        t0 = time.time()
        S, ET, H, st = star_fields(N)
        k = np.array([x[1] for x in st]); dims = sorted(set(x[2] for x in st))
        print(f"{N:4d}  {S.hGamma:.5f}   {dims}   {k.min():.5f}      {k.mean():.5f}      {k.max():.5f}   ({time.time()-t0:.0f}s)", flush=True)


# ------------------------------------------------------------------ assembled test field
def wall_shear(ex):
    """W(x) = d_n u . tau on Gamma, n = -e_r (into the disk), tau = n rotated (+90 deg) -- sign convention irrelevant here"""
    def W(X, Y):
        G = ex["gu"](X, Y); r = np.hypot(X, Y)
        n = -np.stack([X / r, Y / r], -1); tau = np.stack([-n[..., 1], n[..., 0]], -1)
        dnu = np.einsum('...cd,...d->...c', G, n)
        return (dnu * tau).sum(-1)
    return W


def cmd_field(Ns):
    print("test  N    h_Gamma   l(v)/||grad v||    /h^1.5     leading-term prediction /h^1.5")
    rows = {}
    for N in Ns:
        S, ET, H, st = star_fields(N)
        V = S.verts
        circ = [i for i in range(N)]              # circle vertices are the first N vertices (layer 0)
        for kind in ("A", "B1"):
            ex = svn.make_exact(kind); W = wall_shear(ex)
            Wz = np.array([W(V[i, 0], V[i, 1]) for i in circ])
            v = np.zeros(2 * S.nn)
            for i in range(N):
                np.add.at(v, st[i][0][0], Wz[i] * st[i][0][1])
            l = edge_functional(S, "dn", lambda X, Y, n: ex["u"](X, Y))     # <u~, d_n v>
            g = np.sqrt(v @ (H @ v))
            # leading term  -1/2 sum_e W_e int_e s^2 d_n v.tau
            lead = 0.0
            for E in ET:
                xm = E["X"].mean(0); We = W(xm[0], xm[1])
                lead += -0.5 * We * (dn_tau_row(S, E, E["s"] ** 2) @ v)
            r = abs(l @ v) / g
            rows.setdefault(kind, []).append((N, S.hGamma, r))
            print(f"{kind:4s} {N:4d}  {S.hGamma:.5f}   {r:.5e}   {r/S.hGamma**1.5:.5f}    {abs(lead)/g/S.hGamma**1.5:.5f}", flush=True)
    for kind, rr in rows.items():
        print(kind, "rates:", " ".join(f"{np.log(rr[i][2]/rr[i-1][2])/np.log(rr[i][1]/rr[i-1][1]):.3f}" for i in range(1, len(rr))))


# ------------------------------------------------------------------ global dual norms + errors
def dual(H, C, l, free):
    Hf = H[free][:, free].tocsc(); Cf = C[:, free]
    K = sp.bmat([[Hf, Cf.T], [Cf, None]], format="csc")
    lu = spl.splu(K, permc_spec="COLAMD", diag_pivot_thresh=1.0)
    rhs = np.concatenate([l[free], np.zeros(Cf.shape[0])])
    x = lu.solve(rhs)
    for _ in range(2):
        x = x + lu.solve(rhs - K @ x)
    v = np.zeros(H.shape[0]); v[free] = x[:len(free)]
    return np.sqrt(abs(l @ v)), v


def energy_gram(S, H):
    rM, cM, vM, vD = [], [], [], []
    for (t, ref, Xe, we, n) in S.edge_data():
        ph = svn.basis(ref[:, 0], ref[:, 1]); dX, dY = svn.dbasis(ref[:, 0], ref[:, 1]); Ji = S.Jinv[t]
        dn = n[0] * (Ji[0, 0] * dX + Ji[1, 0] * dY) + n[1] * (Ji[0, 1] * dX + Ji[1, 1] * dY)
        Mm = np.einsum('g,ga,gb->ab', we, ph, ph); Dm = we.sum() * np.einsum('g,ga,gb->ab', we, dn, dn)
        d = svn.dofs(S.ids[t])
        for c in range(2):
            rM.append(np.repeat(d[:, c], 15)); cM.append(np.tile(d[:, c], 15)); vM.append(Mm.ravel()); vD.append(Dm.ravel())
    r = np.concatenate(rM); c = np.concatenate(cM); sh = H.shape
    Mb = sp.coo_matrix((np.concatenate(vM), (r, c)), shape=sh).tocsr()
    Dn = sp.coo_matrix((np.concatenate(vD), (r, c)), shape=sh).tocsr()
    return H + (MU / S.h) * Mb + Dn


def cmd_global(Ns, kind="A"):
    ex = svn.make_exact(kind)
    keys = ["ZgradG", "ZenG", "ZgradT", "ZenT", "ZgradD", "ZenD", "ZenP", "Y", "Y1", "CNS_H1", "CNS_slip", "GS_H1"]
    print("N  hG  " + "  ".join(keys))
    rows = []
    for N in Ns:
        vv, tt, cc, oo = svn.make_mesh(N); S = svn.Space(vv, tt, cc, oo)
        H = gram_grad(S).tocsr(); Hen = energy_gram(S, H)
        uex = lambda X, Y, n: ex["u"](X, Y)
        lT = edge_functional(S, "val", lambda X, Y, n: np.einsum('gcd,c->gd', ex["gu"](X, Y), n))
        lD = -edge_functional(S, "dn", uex)
        lP = (MU / S.h) * edge_functional(S, "val", uex)
        sysm = svn.assemble(S, MU, ex["f"], 1); B = sysm["B"]
        ndof = 2 * S.nn; freeR = np.setdiff1d(np.arange(ndof), svn.dofs(S.dnodes).ravel())
        free0 = np.setdiff1d(freeR, svn.dofs(np.array(sorted(gamma_nodes(S)))).ravel())
        o = dict(N=N, hG=S.hGamma)
        o["ZgradG"] = dual(H, B, lT + lD + lP, freeR)[0]; o["ZenG"] = dual(Hen, B, lT + lD + lP, freeR)[0]
        o["ZgradT"] = dual(H, B, lT, freeR)[0]; o["ZenT"] = dual(Hen, B, lT, freeR)[0]
        o["ZgradD"] = dual(H, B, lD, freeR)[0]; o["ZenD"] = dual(Hen, B, lD, freeR)[0]
        o["ZenP"] = dual(Hen, B, lP, freeR)[0]
        B0 = B[1:, :]                                          # div V_h^0 has zero mean: drop one dependent row
        o["Y"] = dual(H, B0, lD, free0)[0]
        Em = sp.csr_matrix(np.array([dn_tau_row(S, E) for E in edge_tables(S)]))
        o["Y1"] = dual(H, sp.vstack([B0, Em]).tocsr(), lD, free0)[0]
        for th, nm in [(1, "CNS"), (0, "GS")]:
            st = svn.assemble(S, MU, ex["f"], th)
            u, p, dres, rr = svn.solve_lu(S, st, ex["u"], th)
            E = svn.errors(S, u, ex, MU)
            o[nm + "_H1"] = E["H1"]; o[nm + "_slip"] = E["bL2"]
        rows.append(o)
        print(f"{N:3d} {S.hGamma:.4f} " + " ".join(f"{o[k]:.4e}" for k in keys), flush=True)
    print("rates in h_Gamma:")
    for k in keys:
        print(f"  {k:9s} " + " ".join(f"{np.log(rows[i][k]/rows[i-1][k])/np.log(rows[i]['hG']/rows[i-1]['hG']):6.3f}" for i in range(1, len(rows))))
    print("Y1 / CNS_H1: " + " ".join(f"{o['Y1']/o['CNS_H1']:.4f}" for o in rows))


if __name__ == "__main__":
    cmd = sys.argv[1]; Ns = [int(a) for a in sys.argv[2].split(",")]
    {"stars": cmd_stars, "field": cmd_field, "global": cmd_global}[cmd](Ns)
