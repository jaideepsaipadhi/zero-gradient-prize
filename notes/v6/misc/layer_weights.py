"""(d): first-layer weights w_e = |e|/H_e on the production meshes (code/svn.py make_mesh), their
symmetry, and the Riemann-sum limit of h_G * sum_e w_e q(m_e*) n_e  (Prop. D1).  usage: python3 layer_weights.py N1,N2,..."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../code"))
import numpy as np, svn
for N in [int(s) for s in sys.argv[1].split(",")]:
    v, t, c, o = svn.make_mesh(N); S = svn.Space(v, t, c, o)
    th, w, le = [], [], []
    for (tt, a, b) in S.bedges:
        A, B = v[t[tt, a]], v[t[tt, b]]; C = v[t[tt, 3 - a - b]]
        L = np.linalg.norm(B - A); area = abs(np.cross(B - A, C - A)) / 2; H = 2 * area / L
        m = (A + B) / 2; th.append(np.arctan2(m[1], m[0])); w.append(L / H); le.append(L)
    th, w, le = map(np.array, (th, w, le)); o_ = np.argsort(th); th, w, le = th[o_], w[o_], le[o_]
    # D4 symmetry test: compare w at theta and at reflected / rotated angles
    def at(x):
        x = np.mod(x + np.pi, 2 * np.pi) - np.pi
        return w[np.argmin(np.abs(np.angle(np.exp(1j * (th - x)))))]
    rot = max(abs(at(x + np.pi / 2) - wi) for x, wi in zip(th, w))
    ref = max(abs(at(-x) - wi) for x, wi in zip(th, w))
    hG = le.max()
    # q = cos(theta) (dipole) and q = cos(2 theta): V_h = hG * sum w q n  (n = -x*, into D)
    n = -np.stack([np.cos(th), np.sin(th)], 1)
    out = []
    for name, q in [("cos1", np.cos(th)), ("sin1", np.sin(th)), ("cos2", np.cos(2 * th)), ("cos3", np.cos(3 * th))]:
        V = hG * (w[:, None] * q[:, None] * n).sum(0)
        out.append(f"{name}: hG*V=({V[0]:+.4f},{V[1]:+.4f})")
    print(f"N={N:4d} #edges={len(w)} w in [{w.min():.3f},{w.max():.3f}] mean {w.mean():.3f}  "
          f"max|w(th+pi/2)-w|={rot:.2e} max|w(-th)-w|={ref:.2e}  len ratio max/min={le.max()/le.min():.3f}  " + "  ".join(out))
