"""H^1 order of the CNS* pressure defect W_*(phi) on a UNIFORM production-type mesh (robust4 REPORT, Sec. 5).
Mesh 'u': N boundary vertices equally spaced on the unit circle; first layer radial at height rho*|e| (all boundary
triangles congruent up to rotation: apex radially above the ccw-forward endpoint); further layers as svn.make_mesh.
Potentials (k = 4):
  sin : sin 2x cos y                                        generic control
  Q   : Re (8/35) r^3... (x^4 tangential power, cos 4theta)   visible control, x^4 in Ker_T but not in K
  B   : grad^4 phi/4! = cos(4 theta) H_*,  H_* = x^4-(x-y/rho)^4     (invisible direction, zero mean eta)
  A   : (r-1)^4 cos 4theta                                   pure normal, varying amplitude
  R   : (r-1)^4                                              pure normal, constant amplitude
usage: python run_u.py N1,N2 phi1,phi2 mu1,mu2 [rho]
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../code"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../robust2"))
import numpy as np
import sympy as sp
import svn_k as K
import leak as LK

Xs, Ys = sp.symbols('X Y', real=True)
_r = sp.sqrt(Xs**2 + Ys**2)
_e4 = sp.expand(((Xs + sp.I * Ys))**4) / _r**4
CB = [sp.Rational(2, 7) + sp.I * sp.Rational(12, 35), sp.Rational(44, 35) + sp.I * sp.Rational(10, 7),
      sp.Rational(86, 35) + sp.I * sp.Rational(57, 35), sp.Rational(52, 35) - sp.I * sp.Rational(16, 35), -1]
EXPR = {
    "sin": sp.sin(2 * Xs) * sp.cos(Ys),
    "Q": sp.re(sp.expand(sp.Rational(8, 35) * (Xs + sp.I * Ys)**4 / _r)),
    "B": sp.re(sp.expand(sum(CB[j] * (_r - 1)**j for j in range(5)) * _e4)),
    "A": (_r - 1)**4 * sp.re(sp.expand(_e4)),
    "R": (_r - 1)**4,
}


def phis(name):
    f = EXPR[name]
    F = sp.lambdify((Xs, Ys), f, "numpy")
    G = sp.lambdify((Xs, Ys), [sp.diff(f, Xs), sp.diff(f, Ys)], "numpy")
    phi = lambda X, Y: np.asarray(F(X, Y), float) + 0 * X
    gphi = lambda X, Y: np.stack([np.asarray(g, float) + 0 * X for g in G(X, Y)], -1)
    return phi, gphi


def mesh_u(N, rho=1.0, L=2.5):
    assert N % 8 == 0
    m = N // 4
    th = 2 * np.pi * np.arange(N) / N
    a = np.stack([np.cos(th), np.sin(th)], 1)
    b = a * (L / np.max(np.abs(a), 1))[:, None]
    X = np.zeros((m + 1, N, 2))
    for j in range(m + 1):
        X[j] = (1 - j / m) * a + (j / m) * b
    X[1] = a * (1 + rho * 2 * np.sin(np.pi / N))
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


def build_u(N, k=4, rho=1.0):
    return K.Space(*mesh_u(N, rho), K.Ref(k))


def run(N, name, mus, k=4, rho=1.0):
    phi, gphi = phis(name)
    S = build_u(N, k, rho)
    E = LK.edge_info(S, phi)
    eta_n = np.sqrt(sum(e["w"] @ e["eta"]**2 for e in E))
    # edge-mean / oscillation split of eta
    osc = np.sqrt(sum(e["w"] @ (e["eta"] - (e["w"] @ e["eta"]) / e["w"].sum())**2 for e in E))
    out = []
    for mu in mus:
        t0 = time.time()
        sysm = K.assemble(S, mu, gphi)
        Bt = K.coupling(S, sysm, 1)
        u, p, rel = LK.solve(S, sysm["A"], Bt, sysm["B"], sysm["F"])
        Er = K.errors(S, u, LK.ZERO, mu)
        gam = mu * S.hGamma / S.h; hG = S.hGamma
        rec = dict(N=N, phi=name, mu=mu, rho=rho, gamma=gam, hG=hG, eta=eta_n, eta_h4=eta_n / hG**4,
                   eta_osc_frac=osc / eta_n, H1=Er["H1"], H1n=Er["H1"] * gam / hG**4.5,
                   H1_eta=Er["H1"] / (hG**0.5 / gam * eta_n), relres=rel, time=time.time() - t0,
                   maxrss_MB=__import__('resource').getrusage(__import__('resource').RUSAGE_SELF).ru_maxrss / 1024)
        print(json.dumps(rec), flush=True); out.append(rec)
    return out


if __name__ == "__main__":
    Ns = [int(x) for x in sys.argv[1].split(",")]; names = sys.argv[2].split(",")
    mus = [float(x) for x in sys.argv[3].split(",")]; rho = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
    for N in Ns:
        for nm in names:
            run(N, nm, mus, rho=rho)
