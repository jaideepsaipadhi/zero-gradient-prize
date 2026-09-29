"""Goal 2 (NUMERICAL, with checks of (M0)-(M2)): a global conforming mesh of Omega_h = R \\ P_h containing the
column-star counterexample S^(H) of Prop. pb:counter, and the coercivity threshold of N_h on Z_h on that mesh.

Mesh (polar column layer + dyadic rings + structured outer part):
  ring 0 : N vertices on the unit circle (P_h = regular N-gon, all chords l = 2 sin(pi/N));
  layer 0: radial columns of height H*l (outer radius r1 = 1 + H l), each split by the forward diagonal
           -> the star of EVERY boundary vertex is a small perturbation of S^(H) (rim (1,0),(1,H),(0,H),(-1,0));
  rings 1..J: cells over two cells of the previous ring, radial thickness ~ aspect * angular width,
           cut into 3 triangles from the hanging bottom midpoint (dyadic coarsening);
  outer  : N/2^J rays to the square R = (-L,L)^2, m_out structured layers, forward diagonal.
Computes: min angle, Theta_0 of (M2), singular vertices (M1), then patch thresholds lambda = (mu/h)^* h_G for
the boundary star, windows of consecutive stars, the first layer, and the GLOBAL threshold, for several L
(h grows with L at fixed h_G, so rho = h/h_G grows)."""
import sys, os, math, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "code"))
import svn, penalty
from lamstar_family import patch_dirichlet


def build_mesh(N, H, J, L, aspect=1.0, m_out=None):
    V = []; T = []; layer = []
    ang0 = 2 * np.pi * np.arange(N) / N
    d0 = np.stack([np.cos(ang0), np.sin(ang0)], 1)
    ell = 2 * math.sin(math.pi / N)
    ring0 = list(range(N)); V += list(d0)
    r1 = 1 + H * ell
    ring1 = list(range(len(V), len(V) + N)); V += list(r1 * d0)
    for i in range(N):
        i1 = (i + 1) % N
        T += [(ring0[i], ring0[i1], ring1[i1]), (ring0[i], ring1[i1], ring1[i])]; layer += [0, 0]
    prev, R, n = ring1, r1, N
    for j in range(1, J + 1):
        n2 = n // 2
        width = R * 2 * np.pi / n2
        R2 = R + aspect * width
        ang = 2 * np.pi * np.arange(n2) / n2
        new = list(range(len(V), len(V) + n2)); V += list(R2 * np.stack([np.cos(ang), np.sin(ang)], 1))
        for k in range(n2):
            bl, mid, br = prev[2 * k], prev[2 * k + 1], prev[(2 * k + 2) % n]
            tl, tr = new[k], new[(k + 1) % n2]
            T += [(bl, mid, tl), (mid, br, tr), (mid, tr, tl)]; layer += [j] * 3
        prev, R, n = new, R2, n2
    assert n % 8 == 0
    ang = 2 * np.pi * np.arange(n) / n
    dr = np.stack([np.cos(ang), np.sin(ang)], 1)
    b = L * dr / np.max(np.abs(dr), 1)[:, None]
    for k in range(n):                                  # snap to the square exactly
        c = np.argmax(np.abs(b[k])); b[k, c] = np.sign(b[k, c]) * L
    a = R * dr
    if m_out is None:
        m_out = max(2, int(math.ceil(np.max(np.linalg.norm(b - a, axis=1)) / (R * 2 * np.pi / n) / 1.5)))
    rings = [prev]
    for j in range(1, m_out + 1):
        t = j / m_out
        new = list(range(len(V), len(V) + n)); V += list((1 - t) * a + t * b); rings.append(new)
    outer = set(rings[-1])
    for j in range(m_out):
        lo, hi = rings[j], rings[j + 1]
        for k in range(n):
            k1 = (k + 1) % n
            T += [(lo[k], lo[k1], hi[k1]), (lo[k], hi[k1], hi[k])]; layer += [J + 1 + j] * 2
    V = np.array(V); T = np.array(T); layer = np.array(layer)
    Pp = V[T]; ar = np.cross(Pp[:, 1] - Pp[:, 0], Pp[:, 2] - Pp[:, 0])
    T[ar < 0] = T[ar < 0][:, [0, 2, 1]]
    return V, T, set(ring0), outer, layer


def mesh_checks(V, T, circ):
    ang = []
    for t in T:
        for i in range(3):
            u, w = V[t[(i + 1) % 3]] - V[t[i]], V[t[(i + 2) % 3]] - V[t[i]]
            ang.append(math.degrees(math.acos(np.clip(np.dot(u, w) / np.linalg.norm(u) / np.linalg.norm(w), -1, 1))))
    minang = min(ang)
    inc = [[] for _ in range(len(V))]
    for ti, t in enumerate(T):
        for i in range(3):
            inc[t[i]].append(ti)
    Theta_min, nsing = 9, 0
    for z in range(len(V)):
        items = []
        for ti in inc[z]:
            t = list(T[ti]); i = t.index(z)
            p, q = V[t[(i + 1) % 3]] - V[z], V[t[(i + 2) % 3]] - V[z]      # ccw: p then q
            a1 = math.atan2(p[1], p[0]); th = math.acos(np.clip(np.dot(p, q) / np.linalg.norm(p) / np.linalg.norm(q), -1, 1))
            items.append((a1, th, p, q))
        dirs = {round((math.degrees(math.atan2(e[1], e[0])) % 180), 6) for it in items for e in (it[2], it[3])}
        # cluster directions mod 180
        ds = sorted(dirs); lines = []
        for d in ds:
            if not lines or min(abs(d - lines[-1]), 180 - abs(d - lines[-1])) > 1e-6:
                lines.append(d)
        if lines and min(abs(lines[0] - lines[-1]), 180 - abs(lines[0] - lines[-1])) < 1e-6 and len(lines) > 1:
            lines = lines[:-1]
        onb = z in circ
        if len(inc[z]) >= 1 and len(lines) <= 2 and (onb or len(inc[z]) >= 3):
            nsing += 1
        if onb and len(inc[z]) == 1:
            nsing += 1
        # order triangles ccw around z: chain by shared edge (q of one == p of next)
        items.sort(key=lambda it: it[0])
        if onb:           # start the chain at the triangle whose p-edge is a boundary edge (no predecessor)
            qs = [tuple(np.round(it[3], 12)) for it in items]
            start = [k for k, it in enumerate(items) if tuple(np.round(it[2], 12)) not in qs]
            k0 = start[0]; order = []; cur = k0
            for _ in range(len(items)):
                order.append(cur)
                nxt = [k for k, it in enumerate(items) if np.allclose(it[2], items[cur][3])]
                if not nxt:
                    break
                cur = nxt[0]
            ths = [items[k][1] for k in order]
            pairs = [(ths[k], ths[k + 1]) for k in range(len(ths) - 1)]
        else:
            ths = [it[1] for it in items]
            pairs = [(ths[k], ths[(k + 1) % len(ths)]) for k in range(len(ths))]
        if pairs:
            Theta = max(abs(math.sin(a + b)) for a, b in pairs)
            Theta_min = min(Theta_min, Theta)
    return minang, Theta_min, nsing


def patch_lambda(V, T, circ, outer, mask, lo=1e-3):
    sub = T[mask]; used = np.unique(sub)
    remap = -np.ones(len(V), int); remap[used] = np.arange(len(used))
    S = svn.Space(V[used], remap[sub], {int(remap[i]) for i in circ if remap[i] >= 0},
                  {int(remap[i]) for i in outer if remap[i] >= 0})
    S.dnodes = patch_dirichlet(S)
    penalty.negpivots.__defaults__[-1].clear()   # the factor cache is keyed by id(S); ids of freed patches are reused
    mu = penalty.threshold(S, lo=lo, hi=1e6, tol=1e-5)
    return (mu / S.h * S.hGamma) if mu > lo * 1.001 else None


if __name__ == "__main__":
    H = float(sys.argv[1]) if len(sys.argv) > 1 else 5.0
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 128
    mode = sys.argv[3] if len(sys.argv) > 3 else "all"
    t0 = time.time()
    print(f"== polar column mesh, H = {H}, N = {N} boundary edges")
    Ls = (2.5, 5.0, 10.0) if mode == "all" else (2.5,)
    for L in Ls:
        J = 2
        while True:        # add dyadic rings while the ring stays well inside the box and N/2^J stays divisible by 8
            V, T, circ, outer, layer = build_mesh(N, H, J + 1, L) if (N // 2**(J + 1)) % 8 == 0 else (None,) * 5
            if V is None:
                break
            Rj = np.max(np.linalg.norm(V[list(outer)], axis=1)) if False else None
            ring_r = np.linalg.norm(V[T[layer == J + 1].ravel()], axis=1).max()
            if ring_r > 0.6 * L:
                break
            J += 1
        V, T, circ, outer, layer = build_mesh(N, H, J, L)
        minang, Th, nsing = mesh_checks(V, T, circ)
        S = svn.Space(V, T, circ, outer)
        ell = 2 * math.sin(math.pi / N)
        print(f"L={L:5.1f} J={J} ntri={len(T)}  h={S.h:.4f}  h_G={S.hGamma:.5f}  rho={S.h/S.hGamma:7.2f}  "
              f"min angle={minang:.2f} deg  Theta_0={Th:.4f}  singular vertices={nsing}", flush=True)
        if L == Ls[0] and mode != "g":
            # patches: boundary star at vertex 0, windows of W consecutive boundary stars, first layer
            def star_mask(vs):
                return np.array([layer[k] == 0 and any(v in T[k] for v in vs) for k in range(len(T))])
            for W in (1, 2, 3, 4, 8, 16):
                vs = [(i - W // 2) % N for i in range(W)]
                lam = patch_lambda(V, T, circ, outer, star_mask(vs))
                print(f"   {W:2d} consecutive boundary star(s): lambda = {('<= 1e-3 (no witness)' if lam is None else f'{lam:.4f}')}"
                      f"   ({time.time()-t0:.0f}s)", flush=True)
            lam = patch_lambda(V, T, circ, outer, layer == 0)
            print(f"   whole first layer: lambda = {lam:.4f}   ({time.time()-t0:.0f}s)", flush=True)
        if mode in ("all", "g"):
            penalty.negpivots.__defaults__[-1].clear()
            mu = penalty.threshold(S, lo=1e-2, hi=1e6, tol=1e-4)
            print(f"   GLOBAL: mu* = {mu:.4f}   mu*/rho = gamma* = mu* h_G/h = {mu*S.hGamma/S.h:.4f}"
                  f"   ({time.time()-t0:.0f}s)", flush=True)
