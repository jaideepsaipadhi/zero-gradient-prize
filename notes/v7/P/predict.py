"""notes/v7/P: cell-coefficient prediction of P_h(phi) on the production meshes (code/svn.make_mesh) and comparison
with the measured moments of notes/v6/misc (run_*.jsonl).

For each wall edge e (chord of the unit circle, fluid outside, kappa = 1):
  l_e = |e|, H_e = 2|T_e|/|e| (first-layer height), a_e = H_e/l_e, lamhat_e = mu*l_e/h (h = h_Omega),
  W_e = d_{e_in}(u.tau_c) at the chord midpoint's projection, tau_c = -e_theta (the cell x-direction).
Cell: 'limit' = rectangular lattice (delta = -1, apex over the left = counterclockwise endpoint), aspect a_e;
      'layer1' = same but the first row uses the actual apex offset of T_e (finite-N fan included in row 1 only).
Prediction (Theorem 3.1 of REPORT.tex):
  P(phi) = sum_e l_e [ W_e l_e A1_e + k0 Ac_e ] phi(m_e),   k0 = sum_e l_e W_e l_e A1_e / sum_e l_e (1 - Ac_e).
usage: python3 predict.py KIND MU N1,N2,...   (KIND in S1, B1, Q1, A)
"""
import os, sys, json, time
for v_ in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v_, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
import numpy as np
import svn
from cell import CellProblem, lattice_rows

K_ = 4; CK = K_ * (K_ + 1)


def exact_grad(kind):
    import sympy as sy
    x, y = sy.symbols('x y', real=True); r2 = x ** 2 + y ** 2
    if kind == "S1":
        psi = (r2 - 1) ** 2 * (1 + x) / 4; u = [sy.diff(psi, y), -sy.diff(psi, x)]
    elif kind == "Q1":
        u = [(r2 - 1) * y, -(r2 - 1) * x]
    elif kind.startswith("B"):
        psi = (r2 - 1) ** 2 * (x + y ** 2 / 2) / 10; u = [sy.diff(psi, y), -sy.diff(psi, x)]
    elif kind == "A":
        u = [(1 - 1 / r2) * (-y), (1 - 1 / r2) * x]
    G = sy.lambdify((x, y), [[sy.diff(u[i], v) for v in (x, y)] for i in range(2)], "numpy")
    return lambda X, Y: np.array(G(X, Y), float)


def mesh_geometry(N):
    verts, tris, circ, outer = svn.make_mesh(N)
    P = verts[tris]
    h = max(np.linalg.norm(P[:, i] - P[:, (i + 1) % 3], axis=1).max() for i in range(3))
    edges = []
    for t in range(len(tris)):
        g = tris[t]
        on = [v in circ for v in g]
        if sum(on) == 2:
            ia, ib = [k for k in range(3) if on[k]]; ic = 3 - ia - ib
            A, B, C = verts[g[ia]], verts[g[ib]], verts[g[ic]]
            if A[0] * B[1] - A[1] * B[0] < 0:        # A -> B counterclockwise
                A, B = B, A
            m = (A + B) / 2; l = np.linalg.norm(B - A)
            e_in = m / np.linalg.norm(m)              # into the fluid (outside the circle)
            tau_c = np.array([e_in[1], -e_in[0]])     # cell x-direction (clockwise), (tau_c, e_in) right-handed
            H = abs(np.dot(C - m, e_in * 0 + np.array([-(B - A)[1], (B - A)[0]]) / l))
            off = np.dot(C - m, tau_c) / l            # apex offset in the cell frame
            edges.append(dict(m=m, l=l, H=H, off=off, e_in=e_in, tau_c=tau_c, A=A, B=B))
    hG = max(e["l"] for e in edges)
    return edges, h, hG


_cache = {}


def cell_A(a, lam, off=None, J=None):
    key = (round(a, 10), round(lam, 8), None if off is None else round(off, 10))
    if key in _cache:
        return _cache[key]
    J = J or int(np.ceil(8.0 / a)) + 1
    rows = lattice_rows(a, -1.0, J)
    if off is not None:                               # row 1 apex offset = off  ->  R_1.x = 0.5 + off
        rows[1:, 0] += (0.5 + off) - 0.0              # shift all rows >= 1 rigidly (keeps rows>=1 a lattice)
    c = CellProblem(rows, 'f', lam)
    o, _, _ = c.run("bump"); oc, _, _ = c.run("trac")
    _cache[key] = (o["A"], oc["A"], o["PL2"], o["trac_L2"])
    return _cache[key]


def predict(kind, mu, N, variant="limit"):
    G = exact_grad(kind)
    edges, h, hG = mesh_geometry(N)
    rec = dict(kind=kind, mu=mu, N=N, h=h, hG=hG, variant=variant)
    S1 = S2 = 0.0; rows = []
    for e in edges:
        a = e["H"] / e["l"]; lam = mu * e["l"] / h
        A1, Ac, PL2, TL2 = cell_A(a, lam, e["off"] if variant == "layer1" else None)
        th = np.arctan2(e["m"][1], e["m"][0]); x0 = np.array([np.cos(th), np.sin(th)])
        Gm = G(x0[0], x0[1])                            # Gm[c, d] = d_d u_c
        W = e["tau_c"] @ Gm @ e["e_in"]                 # d_{e_in} (u . tau_c)
        rows.append((th, e["l"], W, A1, Ac, a, lam))
        S1 += e["l"] * W * e["l"] * A1; S2 += e["l"] * (1 - Ac)
    k0 = S1 / S2
    phis = {"cos1": np.cos, "sin1": np.sin, "cos2": lambda t: np.cos(2 * t), "sin3": lambda t: np.sin(3 * t),
            "cos4": lambda t: np.cos(4 * t)}
    for name, f in phis.items():
        P1 = sum(l * W * l * A1 * f(th) for (th, l, W, A1, Ac, a, lam) in rows)
        Pc = sum(l * k0 * Ac * f(th) for (th, l, W, A1, Ac, a, lam) in rows)
        rec["P_" + name + "/hG"] = (P1 + Pc) / hG
        rec["Pc_" + name + "/hG"] = Pc / hG
    arr = np.array([r[3:] for r in rows])
    rec.update(A1_min=float(arr[:, 0].min()), A1_max=float(arr[:, 0].max()), Ac_min=float(arr[:, 1].min()),
               Ac_max=float(arr[:, 1].max()), a_min=float(arr[:, 2].min()), a_max=float(arr[:, 2].max()),
               lam_min=float(arr[:, 3].min()), lam_max=float(arr[:, 3].max()),
               lam_over_res_min=float((arr[:, 3] * arr[:, 2] / CK).min()),
               lam_over_res_max=float((arr[:, 3] * arr[:, 2] / CK).max()), k0_over_hG=k0 / hG,
               off_mean=float(np.mean([e["off"] for e in edges])))
    return rec


if __name__ == "__main__":
    kind, mu = sys.argv[1], float(sys.argv[2]); Ns = [int(s) for s in sys.argv[3].split(",")]
    variants = sys.argv[4].split(",") if len(sys.argv) > 4 else ["limit", "layer1"]
    for N in Ns:
        for var in variants:
            t0 = time.time(); r = predict(kind, mu, N, var); r["secs"] = round(time.time() - t0, 1)
            print(json.dumps({k: (float(f"{v:.6g}") if isinstance(v, float) else v) for k, v in r.items()}), flush=True)
