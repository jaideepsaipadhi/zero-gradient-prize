"""Rigorous lower bound (S3 of notes/v6/strong + Theorem N2) evaluated on the nearly singular meshes, test A.

For every modified wall vertex z (star T_1..T_m) compute, with the ACTUAL star geometry,
   d_z^2 = min_{(G_T) in A_z} sum_T |T| |grad u~(z) - G_T|_F^2,
   A_z  = {tr G_T = 0, G_T t = 0 on the Gamma_h edges, (G_T - G_T')s = 0 across interior edges s}.
S3: for every div-free v in V_h vanishing on Gamma_h (in particular u_h),
   ||grad(u~ - v)|| >= LB := gamma_k (sum_z d_z^2)^{1/2} - (sum_z R_z^2)^{1/2},   gamma_4 = 1/10,
   R_z^2 = sum_{T in star z} ||grad u~ - grad u~(z)||^2_T      (stars of distinct modified z are disjoint).
Also prints the Theorem-N2 indicator  I := (sum_z h_z^2 |grad u(z)|^2 min(1, phi_z/sqrt(eps_z))^2)^{1/2}
and the ratio sqrt(sum d_z^2)/I (Theorem N2 predicts it to stay in a fixed interval).
Usage: python3 ns_bound.py      (reads runs/*.jsonl for the measured errors)
"""
import os, sys, json, glob
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
import numpy as np
import svn, ns_mesh
from ns_fe import eps_of

k = 4; gam = 2.0 / (k * (k + 1))
ex = svn.make_exact("A")
g, w = np.polynomial.legendre.leggauss(8)
g = (g + 1) / 2; w = w / 2
A_, B_ = np.meshgrid(g, g, indexing="ij"); WA, WB = np.meshgrid(w, w, indexing="ij")
QX = (A_ * (1 - B_)).ravel(); QY = B_.ravel(); QW = (WA * WB * (1 - B_)).ravel()


def star_data(verts, tris, z):
    T = [t for t in range(len(tris)) if z in tris[t]]
    oncirc = lambda v: abs(np.linalg.norm(verts[v]) - 1) < 1e-12
    # order the star: rays from z
    rays = {}
    for t in T:
        for v in tris[t]:
            if v != z:
                rays.setdefault(int(v), []).append(t)
    wallnb = [v for v in rays if oncirc(v) and len(rays[v]) == 1]
    assert len(wallnb) == 2
    start = wallnb[0]; order = [start]; tlist = []
    cur, prev_t = start, None
    while True:
        t = [tt for tt in rays[cur] if tt != prev_t][0]
        tlist.append(t)
        nxt = [int(v) for v in tris[t] if v != z and v != cur][0]
        order.append(nxt); prev_t = t; cur = nxt
        if nxt == wallnb[1]:
            break
    return order, tlist


def dz(verts, tris, z):
    order, tl = star_data(verts, tris, z)
    m = len(tl); Z = verts[z]
    rows = []
    def blk(i):
        e = np.zeros((4, 4 * m)); e[:, 4 * i:4 * i + 4] = np.eye(4); return e   # G = [[a,b],[c,d]] -> (a,b,c,d)
    def Gt(i, v):              # rows of (G_i v)
        M = np.zeros((2, 4 * m)); M[0, 4 * i:4 * i + 2] = v; M[1, 4 * i + 2:4 * i + 4] = v; return M
    for i in range(m):
        r = np.zeros((1, 4 * m)); r[0, 4 * i] = 1; r[0, 4 * i + 3] = 1; rows.append(r)          # trace
    rows.append(Gt(0, verts[order[0]] - Z)); rows.append(Gt(m - 1, verts[order[-1]] - Z))       # wall edges
    for i in range(m - 1):
        rows.append(Gt(i, verts[order[i + 1]] - Z) - Gt(i + 1, verts[order[i + 1]] - Z))      # continuity
    C = np.vstack(rows)
    u, s, vt = np.linalg.svd(C); Nsp = vt[int((s > 1e-10 * s[0]).sum()):].T
    Gz = ex["gu"](Z[None, 0], Z[None, 1])[0]
    area = [0.5 * abs(np.cross(verts[tris[t][1]] - verts[tris[t][0]], verts[tris[t][2]] - verts[tris[t][0]])) for t in tl]
    Amat = np.vstack([np.sqrt(area[i]) * blk(i) @ Nsp for i in range(m)])
    b = np.concatenate([np.sqrt(area[i]) * Gz.ravel() for i in range(m)])
    y, *_ = np.linalg.lstsq(Amat, b, rcond=None)
    d2 = float(np.sum((Amat @ y - b) ** 2))
    R2 = 0.0
    for t in tl:
        P = verts[tris[t]]
        J = np.array([P[1] - P[0], P[2] - P[0]]).T; det = abs(np.linalg.det(J))
        X = P[0][0] + J[0, 0] * QX + J[0, 1] * QY; Y = P[0][1] + J[1, 0] * QX + J[1, 1] * QY
        R2 += det * np.sum(QW[:, None, None] * (ex["gu"](X, Y) - Gz) ** 2)
    # indicator data
    angs = []
    for t in tl:
        o = [verts[v] - Z for v in tris[t] if v != z]
        angs.append(np.arccos(np.dot(o[0], o[1]) / np.linalg.norm(o[0]) / np.linalg.norm(o[1])))
    eps = min(angs); phi = sum(angs) - np.pi
    hz = max(np.linalg.norm(verts[v] - Z) for v in order)
    wz = 1.0 if m == 2 else min(1.0, phi / np.sqrt(eps))      # plain lock: w_z = 1
    ind2 = hz ** 2 * np.sum(Gz ** 2) * wz ** 2
    return d2, R2, ind2, m


if __name__ == "__main__":
    out = []
    for f in sorted(glob.glob(os.path.join(HERE, "runs", "*.jsonl"))):
        if "infsup" in f:
            continue
        for line in open(f):
            if not line.startswith("{"):
                continue
            r = json.loads(line)
            if r["mesh"] in ("std",):
                continue
            eps = r["eps"]
            verts, tris, circ, outer, lk = ns_mesh.make(r["mesh"], r["N"], eps)
            tris = np.array(tris)
            D2 = R2 = I2 = 0.0
            for z in lk:
                d2, rr, i2, m = dz(verts, tris, z)
                D2 += d2; R2 += rr; I2 += i2
            LB = gam * np.sqrt(D2) - np.sqrt(R2)
            rec = dict(mesh=r["mesh"], mode=r["mode"], N=r["N"], hG=r["hG"], eps=eps, E=r["H1"], sqrtD=np.sqrt(D2),
                       gam_sqrtD=gam * np.sqrt(D2), R=np.sqrt(R2), LB=LB, ind=np.sqrt(I2),
                       D_over_ind=np.sqrt(D2 / I2) if I2 > 0 else float("nan"), E_over_gamD=r["H1"] / (gam * np.sqrt(D2)))
            out.append(rec)
            print(json.dumps({a: (float(f"{b:.4g}") if isinstance(b, float) else b) for a, b in rec.items()}))
    with open(os.path.join(HERE, "ns_bound.log"), "w") as fh:
        for rec in out:
            fh.write(json.dumps(rec) + "\n")
