"""NUMERICAL test of the smooth-data stability hypothesis (S) (robust5 REPORT, Sec. 5).
CNS* (k = 4) velocity W^a driven by the pure boundary load  v -> <a, v.n>_{Gamma_h},  a = cos(2 theta) (smooth, mean 0).
(S) predicts ||grad W^a|| ~ gamma^{-1} h_Gamma; the Nitsche pressure layer argument predicts an extra
gamma^{-3/2} h_Gamma^{1/2} term whose leading part is the load (20/gamma) (h_Gamma/H_e) a: smooth when the first-layer
height ratio H_e/h_Gamma is constant ('u' meshes), oscillating when it alternates ('alt' meshes).
Meshes: robust4 mesh_u, first-layer height rho_i*|e| at vertex i (rho_i = rho for 'u'; alternating r1, r2 for 'alt').
usage: python run_S.py N1,N2 mu1,mu2 mesh(u|alt) [r1 r2]
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../code"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../robust2"))
import numpy as np
import svn_k as K
import leak as LK


def mesh(N, rhos, L=2.5):
    m = N // 4
    th = 2 * np.pi * np.arange(N) / N
    a = np.stack([np.cos(th), np.sin(th)], 1)
    b = a * (L / np.max(np.abs(a), 1))[:, None]
    X = np.zeros((m + 1, N, 2))
    for j in range(m + 1):
        X[j] = (1 - j / m) * a + (j / m) * b
    X[1] = a * (1 + np.asarray(rhos)[:, None] * 2 * np.sin(np.pi / N))
    assert np.all(np.linalg.norm(X[1], axis=1) < np.linalg.norm(X[2], axis=1))
    vid = np.arange((m + 1) * N).reshape(m + 1, N)
    verts = X.reshape(-1, 2); tris = []
    for j in range(m):
        for i in range(N):
            i1 = (i + 1) % N
            tris.append([vid[j, i], vid[j, i1], vid[j + 1, i1]])
            tris.append([vid[j, i], vid[j + 1, i1], vid[j + 1, i]])
    tris = np.array(tris)
    P = verts[tris]; ar = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    tris[ar < 0] = tris[ar < 0][:, [0, 2, 1]]
    return verts, tris, set(vid[0].tolist()), set(vid[m].tolist())


def minangle(verts, tris):
    P = verts[tris]; out = 180.0
    for i in range(3):
        u = P[:, (i + 1) % 3] - P[:, i]; w = P[:, (i + 2) % 3] - P[:, i]
        c = np.einsum('ij,ij->i', u, w) / np.linalg.norm(u, axis=1) / np.linalg.norm(w, axis=1)
        out = min(out, np.degrees(np.arccos(np.clip(c, -1, 1))).min())
    return out


def run(N, mus, kind, r1=0.7, r2=1.4, k=4, afun=lambda th: np.cos(2 * th)):
    rhos = np.ones(N) if kind == "u" else np.where(np.arange(N) % 2 == 0, r1, r2)
    v, t, c, o = mesh(N, rhos)
    S = K.Space(v, t, c, o, K.Ref(k))
    R = S.R
    for mu in mus:
        t0 = time.time()
        sysm = K.assemble(S, mu, lambda X, Y: np.zeros(np.shape(X) + (2,)))
        F = np.zeros(2 * S.nn)
        for (tt, ref, Xe, we, n) in S.edge_data():
            ph = R.basis(ref[:, 0], ref[:, 1]); d = K.dofs(S.ids[tt])
            a = afun(np.arctan2(Xe[:, 1], Xe[:, 0]))
            for cc in range(2):
                F[d[:, cc]] += ph.T @ (we * a) * n[cc]
        Bt = K.coupling(S, sysm, 1)
        u, p, rel = LK.solve(S, sysm["A"], Bt, sysm["B"], F)
        Er = K.errors(S, u, LK.ZERO, mu)
        gam = mu * S.hGamma / S.h; hG = S.hGamma
        rec = dict(N=N, mesh=kind, mu=mu, gamma=gam, hG=hG, minangle=minangle(v, t), H1=Er["H1"],
                   S_ratio=Er["H1"] * gam / hG, layer_ratio=Er["H1"] * gam**1.5 / hG**0.5, relres=rel,
                   time=round(time.time() - t0, 1))
        print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    Ns = [int(x) for x in sys.argv[1].split(",")]; mus = [float(x) for x in sys.argv[2].split(",")]
    kind = sys.argv[3]; r = [float(x) for x in sys.argv[4:6]] or [0.7, 1.4]
    for N in Ns:
        run(N, mus, kind, *r)
