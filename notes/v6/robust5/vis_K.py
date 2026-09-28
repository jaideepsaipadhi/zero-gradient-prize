"""NUMERICAL evidence for hypothesis (V) of the (P^K) reduction (robust5 REPORT, Sec. 4):
visibility of the vertex functional  v -> v(z).n  by the pressure-blind star space V#(omega_z) of robust3,
  vis(star) := sup_{v in V#} |v(z).n_bar| / Phi_Q(v)       (scale-invariant units, |e_z| = 1),
on exactly periodic flat m = 3 stars (T' = T - t) and on tilted / perturbed stars approaching them along several
directions (eps1, eps2) = rotation of e and of e' separately, plus random shape perturbations.
Also reports kappa_K (robust4 definition) for the same stars.
usage: python vis_K.py"""
import sys, os, json, inspect
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../robust3"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../robust4"))
import pk, stars as S

_s = inspect.getsource(pk.star_forms).replace("return dict(dimV", "return M, Q, ns, fdof, dict(dimV")
_ns = dict(pk.__dict__); exec(_s, _ns); SF = _ns["star_forms"]


def star(c, eps1=0.0, eps2=0.0, pert=None, l=1.0):
    """periodic m=3 star from T = (z, a, c): e = [z,(1,0)], e' = [(-l,0), z], middle (z, c, c - (l,0)).
    eps1 rotates a about z (e tilts), eps2 rotates g2 about z; pert: displacement of the interior vertices."""
    c = np.asarray(c, float); cp = c - np.array([l, 0.0])
    a = np.array([np.cos(eps1), -np.sin(eps1)]); g2 = -l * np.array([np.cos(eps2), np.sin(eps2)])
    P = np.array([[0, 0], a, c, cp, g2], float)
    if pert is not None:
        P[2] += pert[0]; P[3] += pert[1]
    return P, [[0, 1, 2], [0, 2, 3], [4, 0, 3]]


def vis(P, fan):
    _ns["GAMMA"] = (1, 4)
    M, Q, ns, fdof, info = SF(P, fan, 4)
    pos = {int(d): i for i, d in enumerate(fdof)}
    n = np.array([0.0, -1.0])                 # outward normal of the flat configuration
    f = n[0] * ns[pos[0]] + n[1] * ns[pos[1]]
    val = float(np.sqrt(f @ np.linalg.solve(Q, f)))
    T = P[fan[0]]
    return val, S.kappa(M, S.KT(T)), S.kappa(M, S.KerT(T)), info["dimV"]


out = []
for c in [(1.0, 1.0), (0.5, 0.9), (1.3, 0.6), (0.2, 1.2)]:
    v0, kK0, kT0, d0 = vis(*star(c))
    rec = dict(c=c, periodic=dict(vis=v0, kappa_K=kK0, kappa_T=kT0, dimV=d0), approach=[])
    for (u1, u2) in [(1, 1), (1, -1), (1, 0), (0, 1), (1, 3)]:
        for t in [1e-1, 1e-2, 1e-3]:
            v, kK, kT, d = vis(*star(c, eps1=t * u1, eps2=t * u2))
            rec["approach"].append(dict(dir=(u1, u2), t=t, vis=v, kappa_K=kK, dimV=d))
    rng = np.random.default_rng(0)
    for t in [1e-1, 1e-2, 1e-3]:
        for rep in range(3):
            pert = rng.normal(size=(2, 2)) * t
            v, kK, kT, d = vis(*star(c, eps1=t * rng.normal(), eps2=t * rng.normal(), pert=pert))
            rec["approach"].append(dict(dir="random", t=t, vis=v, kappa_K=kK, dimV=d))
    print(json.dumps(rec), flush=True)
