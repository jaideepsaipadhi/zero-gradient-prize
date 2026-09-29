"""
notes/v7/nsconst: the limit gamma*_inf of the coercivity threshold (REPORT.tex, Section 3).

Uses the one-edge Floquet cell of notes/v6/unifmu2/floquet_cell.py (imported, not modified):
a sector of a regular N-gon annulus, J layers of aspect a (first-layer height / chord), quads split on the
forward diagonal, quasi-periodic identification with phase phi.  For N -> inf this is the flat half-plane
lattice of rectangles [0,1] x [0,a] (J rows, zero Dirichlet below row J).

  gamma_c(a, J; N) := sup_phi  l * (smallest lam with N_h(lam) >= 0 on Z_h(phi))      (edge units)

Modes:
  python cell_limit.py phase  N a J          gamma_c(phi) on a phase grid (where is the sup?)
  python cell_limit.py table                 gamma_c(0.75, J) for J = 1..8, N = 1024, 4096, 16384
  python cell_limit.py theta J               f(theta) = gamma_c(a(theta), J)/cos^2(theta), theta in [0, pi/4]
  python cell_limit.py mesh Ns J             J-layer truncated threshold of the production meshes (penalty bisection)
  python cell_limit.py aspect N              local chord / aspect of the production mesh (check of a(theta))
"""
import os, sys, json, time
for v_ in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v_, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "..")
sys.path.insert(0, os.path.join(ROOT, "code"))
sys.path.insert(0, os.path.join(ROOT, "notes", "v6", "unifmu2"))
import numpy as np, scipy.linalg as sla
from scipy.optimize import minimize_scalar
import floquet_cell as fc

L = 2.5


class Thr:
    """precomputed cell; gamma(phase) on Z_h."""
    def __init__(self, N, a, J):
        self.cell = fc.Cell(N, a, J=J)
        self.A0 = (self.cell.Avol + self.cell.Cn); self.Mg = self.cell.Mg

    def gamma(self, phase):
        c = self.cell
        T = c.T(phase).toarray()
        Z = sla.null_space(c.B @ T, rcond=1e-11)
        Q = T @ Z
        A0q = Q.conj().T @ (self.A0 @ Q); Mq = Q.conj().T @ (self.Mg @ Q)
        A0q = (A0q + A0q.conj().T) / 2; Mq = (Mq + Mq.conj().T) / 2
        w, U = np.linalg.eigh(Mq); keep = w > 1e-12 * w.max()
        Uk, U0 = U[:, keep], U[:, ~keep]
        A00 = U0.conj().T @ A0q @ U0
        if A00.shape[0]:
            e0 = np.linalg.eigvalsh((A00 + A00.conj().T) / 2).min()
            assert e0 > 0, e0                      # N_h > 0 on {v in Z: v = 0 on the wall edge}
            Ak0 = Uk.conj().T @ A0q @ U0; Akk = Uk.conj().T @ A0q @ Uk
            Sch = Akk - Ak0 @ np.linalg.solve(A00, Ak0.conj().T)
        else:
            Sch = Uk.conj().T @ A0q @ Uk
        ev = sla.eigh((Sch + Sch.conj().T) / 2, np.diag(w[keep]), eigvals_only=True)
        return float(-ev.min() * c.ell)

    def sup(self, ngrid=17):
        ph = np.linspace(0, np.pi, ngrid)
        g = [self.gamma(p) for p in ph]
        i = int(np.argmax(g))
        lo, hi = ph[max(i - 1, 0)], ph[min(i + 1, ngrid - 1)]
        if hi - lo < 1e-12:
            return g[i], ph[i]
        r = minimize_scalar(lambda p: -self.gamma(p), bounds=(lo, hi), method="bounded", options=dict(xatol=1e-5))
        if -r.fun > g[i]:
            return float(-r.fun), float(r.x)
        return g[i], float(ph[i])


def a_of_theta(th):
    """limit first-layer aspect of svn.make_mesh at polar angle th in [0, pi/4] (side x = L)."""
    return (L / np.cos(th) - 1) / (2 * np.cos(th) ** 2)


def cmd_phase(N, a, J):
    t = Thr(N, a, J)
    for p in np.linspace(0, np.pi, 25):
        print(json.dumps(dict(N=N, a=a, J=J, phase=round(p, 5), gamma=round(t.gamma(p), 6))), flush=True)
    s, p = t.sup()
    print(json.dumps(dict(N=N, a=a, J=J, sup=s, argmax=p)), flush=True)


def cmd_table():
    for J in (1, 2, 3, 4, 6, 8):
        for N in (1024, 4096, 16384):
            t0 = time.time(); s, p = Thr(N, 0.75, J).sup()
            print(json.dumps(dict(a=0.75, J=J, N=N, gamma_c=round(s, 6), argmax_phase=round(p, 5),
                                  secs=round(time.time() - t0, 1))), flush=True)


def cmd_theta(J, N=4096):
    best = (-1, None)
    for th in np.linspace(0, np.pi / 4, 13):
        a = a_of_theta(th)
        s, p = Thr(N, a, J).sup()
        f = s / np.cos(th) ** 2
        best = max(best, (f, th))
        print(json.dumps(dict(theta=round(th, 5), a=round(a, 5), J=J, gamma_c=round(s, 5), phase=round(p, 4),
                              f=round(f, 5))), flush=True)
    print(json.dumps(dict(J=J, max_f=best[0], at_theta=best[1])), flush=True)


def cmd_aspect(N):
    import svn
    verts, tris, circ, outer = svn.make_mesh(N)
    m = N // 4
    a0 = verts[:N]; a1 = verts[N:2 * N]
    ch = np.linalg.norm(np.roll(a0, -1, 0) - a0, axis=1)
    th = np.arctan2(a0[:, 1], a0[:, 0])
    ht = np.linalg.norm(a1 - a0, axis=1)
    hG = ch.max()
    for i in range(N // 8 + 1):
        j = (i + N // 8) % N  # vertex index: side x = L starts at corner (L,-L); vertex N/8 is theta = 0
        thm = np.arctan2((a0[j] + a0[(j + 1) % N])[1], (a0[j] + a0[(j + 1) % N])[0])
        print(json.dumps(dict(N=N, theta_mid=round(float(thm), 4), chord_over_hG=round(float(ch[j] / hG), 4),
                              cos2=round(float(np.cos(thm) ** 2), 4), aspect=round(float(ht[j] / ch[j]), 4),
                              a_limit=round(float(a_of_theta(abs(thm))), 4))), flush=True)


def cmd_mesh(Ns, J):
    import svn, penalty
    from lamstar_family import patch_dirichlet
    for N in Ns:
        t0 = time.time()
        verts, tris, circ, outer = svn.make_mesh(N)
        layer = np.repeat(np.arange(N // 4), 2 * N)
        sub = tris[layer < J]; used = np.unique(sub)
        remap = -np.ones(len(verts), int); remap[used] = np.arange(len(used))
        S = svn.Space(verts[used], remap[sub], {int(remap[i]) for i in circ if remap[i] >= 0},
                      {int(remap[i]) for i in outer if remap[i] >= 0})
        S.dnodes = patch_dirichlet(S)
        mu = penalty.threshold(S, lo=1e-3, hi=1e6, tol=1e-6)
        lam = mu / S.h * S.hGamma          # (mu/h)* hG with the patch's own h: gamma normalised by max|e|
        print(json.dumps(dict(N=N, J=J, gamma_J=round(lam, 5), hG=S.hGamma, secs=round(time.time() - t0, 1))),
              flush=True)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "phase":
        cmd_phase(int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]))
    elif cmd == "table":
        cmd_table()
    elif cmd == "theta":
        cmd_theta(int(sys.argv[2]))
    elif cmd == "aspect":
        cmd_aspect(int(sys.argv[2]))
    elif cmd == "mesh":
        cmd_mesh([int(x) for x in sys.argv[2].split(",")], int(sys.argv[3]))
