"""Star seminorm sigma_z (robust3 pk.star_forms) on flat 3-triangle stars with prescribed boundary triangles.
z = (0,0), e = [z, g1 = (1,0)], e' = [g2 = (-l,0), z]; T = (z, g1, w1), T' = (g2, z, w0), middle (z, w1, w0).
Checks (robust4 REPORT Sec. 3):
  kernel dimension vs the parallel-edge criterion; kernel = span{n^4, l_L^4 - l_R^4} on homothetic pairs;
  kappa_S := min{ sigma_z(H) : dist_F(H, S) = 1 } for S = Ker_T (3-dim), S = K_T (2-dim), S = R n^4.
usage: python stars.py
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../robust3"))
import numpy as np
import scipy.linalg as sl
from math import comb
import pk

k = 4
D = np.diag([np.sqrt(comb(k, j)) for j in range(k + 1)])


def star3(w1, w0, l=1.0):
    P = np.array([[0, 0], [1, 0], w1, w0, [-l, 0]], float)
    pk.GAMMA = (1, 4)
    return pk.star_forms(P, [[0, 1, 2], [0, 2, 3], [4, 0, 3]], k), P


import inspect
_s = inspect.getsource(pk.star_forms).replace("return dict(dimV", "return M, dict(dimV")
_ns = dict(pk.__dict__); exec(_s, _ns); _star_forms_M = _ns["star_forms"]


def Mform(w1, w0, l=1.0):
    """(M, info): the 5x5 form of sigma_z^2 in Frobenius coords, for the 3-triangle star."""
    P = np.array([[0, 0], [1, 0], w1, w0, [-l, 0]], float)
    _ns["GAMMA"] = (1, 4)
    return _star_forms_M(P, [[0, 1, 2], [0, 2, 3], [4, 0, 3]], k)


def report(name, w1, w0, l=1.0):
    M, info = Mform(w1, w0, l)
    T = np.array([[0, 0], [1, 0], w1], float); Tp = np.array([[-l, 0], [0, 0], w0], float)
    ev = np.sqrt(np.clip(np.linalg.eigvalsh(M), 0, None))
    sc = ev.max()
    out = dict(case=name, dimV=info["dimV"], sig_rel=[float(x / sc) for x in ev],
               kap_n=kappa(M, powvec(0, 1)[None, :]) / sc, kap_K=kappa(M, KT(T)) / sc,
               kap_Kp=kappa(M, KT(Tp)) / sc, kap_KerT=kappa(M, KerT(T)) / sc)
    # is the 2-dim near-kernel equal to K_T ?  principal angle
    lam, U = np.linalg.eigh(M)
    Q = sl.orth(KT(T).T); out["angle_K"] = float(np.linalg.norm(U[:, :2] - Q @ (Q.T @ U[:, :2])))
    print(json.dumps(out), flush=True)
    return out


def powvec(lx, ly):
    """Frobenius coords of (lx x + ly y)^k."""
    h = np.array([comb(k, j) * lx**j * ly**(k - j) for j in range(k + 1)], float)
    return np.linalg.solve(D, h)


def KerT(T):
    """3 edge-normal powers of triangle T (3x2) -> 3 x 5 matrix (Frobenius coords, unnormalised)."""
    out = []
    for i in range(3):
        a, b = T[(i + 1) % 3], T[(i + 2) % 3]
        t = b - a; out.append(powvec(t[1], -t[0]))
    return np.array(out)


def KT(T):
    """K_T = span{l_c^k, l_L^k - l_R^k} with barycentric gradients; T = (L, R, apex), e = [L, R]."""
    A = np.hstack([np.asarray(T, float), np.ones((3, 1))])
    Cf = np.linalg.solve(A, np.eye(3))          # column i: coefficients (a, b, c) of lam_i
    gL, gR, gC = Cf[:2, 0], Cf[:2, 1], Cf[:2, 2]
    return np.array([powvec(*gC), powvec(*gL) - powvec(*gR)])


def kappa(M, S):
    """min over unit H orthogonal to span(S) of min_s (H+s)^T M (H+s)  (sqrt)."""
    Q = sl.orth(S.T)                      # 5 x d
    B = sl.null_space(Q.T)               # 5 x (5-d)
    Mss = Q.T @ M @ Q; Msb = Q.T @ M @ B; Mbb = B.T @ M @ B
    w_, V_ = np.linalg.eigh(Mss); tol = 1e-12 * max(np.trace(M), 1e-300)
    Mp = (V_[:, w_ > tol] / w_[w_ > tol]) @ V_[:, w_ > tol].T
    Sch = Mbb - Msb.T @ Mp @ Msb
    return float(np.sqrt(max(np.linalg.eigvalsh(Sch).min(), 0)))


def prod_scan(N, L=2.5):
    """production stars (svn.make_mesh): relative singular values and kappa_S for S = n^4, K_{T_e}, K_{T_e'}, Ker_{T_e}."""
    import svn
    v, t, circ, outer = svn.make_mesh(N, L=L)
    res = []
    for i in range(N):
        z = i; fan = [tr for tr in t if z in tr]
        loc = sorted({int(a) for tr in fan for a in tr}); mp = {z: 0}
        for a in loc:
            if a != z: mp[a] = len(mp)
        P = np.zeros((len(mp), 2))
        for a, b in mp.items(): P[b] = v[a]
        g1 = mp[(i + 1) % N]; g2 = mp[(i - 1) % N]
        _ns["GAMMA"] = (g1, g2)
        fl = [[mp[int(a)] for a in tr] for tr in fan]
        M, info = _star_forms_M(P, fl, k)
        # rotate to the frame of star_forms (it uses global x,y centred at z scaled by s): K in global coords
        s = np.linalg.norm(P[g1] - P[0])
        # boundary triangles: T_e contains z, g1; T_e' contains z, g2
        Te = [tr for tr in fl if 0 in tr and g1 in tr][0]; Tp = [tr for tr in fl if 0 in tr and g2 in tr][0]
        apex = [a for a in Te if a not in (0, g1)][0]; apexp = [a for a in Tp if a not in (0, g2)][0]
        T = P[[0, g1, apex]] ; Tq = P[[g2, 0, apexp]]
        # orientation: e = [L,R] with interior on the left of L->R ... use z->g1 as given; KT is orientation-free
        # except the sign of l_L - l_R (irrelevant for a span)
        ev = np.sqrt(np.clip(np.linalg.eigvalsh(M), 0, None)); sc = ev.max()
        n = np.array([(P[g1] - P[0])[1], -(P[g1] - P[0])[0]]) / s
        res.append(dict(sig=ev / sc, kn=kappa(M, powvec(*n)[None, :]) / sc, kK=kappa(M, KT(T)) / sc,
                        kKp=kappa(M, KT(Tq)) / sc, kKer=kappa(M, KerT(T)) / sc, dimV=info["dimV"]))
    agg = dict(N=N, dimV=sorted({r["dimV"] for r in res}),
               sig2_min=float(min(r["sig"][1] for r in res)), sig3_min=float(min(r["sig"][2] for r in res)),
               kap_n_min=float(min(r["kn"] for r in res)), kap_K_min=float(min(r["kK"] for r in res)),
               kap_Kp_min=float(min(r["kKp"] for r in res)), kap_KerT_min=float(min(r["kKer"] for r in res)))
    print(json.dumps(agg), flush=True)
    return agg


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "prod":
    for N in [int(x) for x in sys.argv[2].split(",")]:
        prod_scan(N)
    sys.exit()

if __name__ == "__main__":
    report("periodic rho=1", (1, 1), (0, 1))
    report("periodic rho=0.6 apex 0.3", (1.3, 0.6), (0.3, 0.6))
    report("homothetic l=0.8", (1, 1), (0, 0.8), 0.8)
    report("homothetic l=1.2 apex", (0.7, 0.9), (-1.2 + 0.7 * 1.2, 0.9 * 1.2), 1.2)
    report("share radial only", (1, 1), (0, 0.6))
    report("share hyp only", (1, 1), (0.4, 1.4))
    report("mirror", (1, 1), (-1, 1))
    report("generic", (0.8, 1.1), (-0.3, 0.7))
    report("production-like N=64", (1.0 + 0.05, 1.0), (0.0, 1.0 * 1.02))


