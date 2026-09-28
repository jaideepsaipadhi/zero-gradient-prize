"""
General unstructured meshes for Section 10 (general smooth domains): a polygonal outer box Sigma
(strong Dirichlet, or partly natural for the channel) and m >= 1 smooth closed curves Gamma_i, each
approximated by an INSCRIBED polygon Gamma_{h,i} whose N_i vertices lie on Gamma_i (equal arclength
spacing, so (G1) holds).  Triangulated with Shewchuk's Triangle (pip package `triangle`), quality
switch q30, no Steiner points on any segment (switch Y), so the Gamma_h vertices are exactly the
curve points and the Gamma_h edges are exactly the chords.  Local size

    s(x) = min(R_cap * s_G,  s_G * (1 + grad * dist(x, Gamma)))       (s_G = min chord length)

so h/hG stays roughly constant under refinement (like the structured meshes, where h/hG ~ 3.3).

After meshing, the paper's mesh assumptions are checked and repaired:
  (M1)  no singular vertex: every boundary vertex that is not a corner of Sigma has >= 3 triangles
        (a straight-side vertex with 2 triangles has its edges on two lines = singular), corners of
        Sigma have >= 2; interior vertices must not be "X" vertices;
  (M1') every Gamma_h vertex has >= 3 triangles;
  (M2)  Theta(z) = max_j |sin(theta_j + theta_{j+1})| >= Theta_0 (wrap-around omitted at boundary
        vertices).
Repair = barycentric (centroid) split of the incident triangle with the largest angle at the bad
vertex (raises its triangle count by one, creates a non-singular 3-valent interior vertex).  Every
repair is counted and reported.

API
    curve objects:  Circle(cx, cy, r), Ellipse(cx, cy, a, b), PolarCurve(cx, cy, s, eps, m)
                    (r(t) = s (1 + eps cos(m t)); eps = 0.3, m = 3 is the nonconvex test curve)
    make_mesh(box, curves, N, R_cap=3.5, grad=1.0, extra_lines=(), s_G=None) -> Mesh
    Mesh fields: verts, tris, gamma_edges [(a, b, comp)], sigma_edges [(a, b, side)], stats (dict)
"""
import numpy as np
import triangle as _tr
from scipy.spatial import cKDTree


# ---------------------------------------------------------------- curves
class Curve:
    """closed C^inf curve x(t), t in [0, 2 pi), positively oriented; subclasses give pos, dpos, ddpos."""
    centre = (0.0, 0.0)

    def length(self, n=20000):
        t = np.linspace(0, 2 * np.pi, n + 1)[:-1]
        return np.linalg.norm(self.dpos(t), axis=1).sum() * 2 * np.pi / n

    def arclength_params(self, N, n=20000):
        """N parameters equally spaced in arclength (spectrally accurate trapezoid + interpolation)."""
        t = np.linspace(0, 2 * np.pi, n + 1)
        sp = np.linalg.norm(self.dpos(t), axis=1)
        s = np.concatenate([[0], np.cumsum((sp[1:] + sp[:-1]) / 2 * (t[1] - t[0]))])
        target = np.arange(N) * s[-1] / N
        tt = np.interp(target, s, t)
        for _ in range(3):                           # Newton polish on the arclength equation
            si = np.interp(tt, t, s)
            tt = tt - (si - target) / np.linalg.norm(self.dpos(tt), axis=1)
        return tt

    def curvature(self, t):
        d = self.dpos(t); dd = self.ddpos(t)
        return (d[:, 0] * dd[:, 1] - d[:, 1] * dd[:, 0]) / np.linalg.norm(d, axis=1) ** 3


class Circle(Curve):
    def __init__(self, cx, cy, r):
        self.cx, self.cy, self.r = cx, cy, r; self.centre = (cx, cy)
        self.name = f"circle(r={r})"

    def pos(self, t):
        return np.stack([self.cx + self.r * np.cos(t), self.cy + self.r * np.sin(t)], -1)

    def dpos(self, t):
        return np.stack([-self.r * np.sin(t), self.r * np.cos(t)], -1)

    def ddpos(self, t):
        return np.stack([-self.r * np.cos(t), -self.r * np.sin(t)], -1)


class Ellipse(Curve):
    def __init__(self, cx, cy, a, b):
        self.cx, self.cy, self.a, self.b = cx, cy, a, b; self.centre = (cx, cy)
        self.name = f"ellipse(a={a},b={b})"

    def pos(self, t):
        return np.stack([self.cx + self.a * np.cos(t), self.cy + self.b * np.sin(t)], -1)

    def dpos(self, t):
        return np.stack([-self.a * np.sin(t), self.b * np.cos(t)], -1)

    def ddpos(self, t):
        return np.stack([-self.a * np.cos(t), -self.b * np.sin(t)], -1)


class PolarCurve(Curve):
    """r(t) = s (1 + eps cos(m t)); for eps = 0.3, m = 3 the curvature changes sign (nonconvex)."""
    def __init__(self, cx, cy, s, eps, m):
        self.cx, self.cy, self.s, self.eps, self.m = cx, cy, s, eps, m; self.centre = (cx, cy)
        self.name = f"polar(1+{eps}cos{m}t)"

    def _r(self, t):
        return (self.s * (1 + self.eps * np.cos(self.m * t)),
                -self.s * self.eps * self.m * np.sin(self.m * t),
                -self.s * self.eps * self.m ** 2 * np.cos(self.m * t))

    def pos(self, t):
        r, _, _ = self._r(t)
        return np.stack([self.cx + r * np.cos(t), self.cy + r * np.sin(t)], -1)

    def dpos(self, t):
        r, r1, _ = self._r(t)
        return np.stack([r1 * np.cos(t) - r * np.sin(t), r1 * np.sin(t) + r * np.cos(t)], -1)

    def ddpos(self, t):
        r, r1, r2 = self._r(t)
        return np.stack([r2 * np.cos(t) - 2 * r1 * np.sin(t) - r * np.cos(t),
                         r2 * np.sin(t) + 2 * r1 * np.cos(t) - r * np.sin(t)], -1)


# ---------------------------------------------------------------- mesh
_DEBUG = False


class Mesh:
    pass


def _angles(P):
    """interior angles (nt,3) of triangles P (nt,3,2); angle i at vertex i."""
    out = np.zeros(P.shape[:2])
    for i in range(3):
        a = P[:, (i + 1) % 3] - P[:, i]; b = P[:, (i + 2) % 3] - P[:, i]
        c = (a * b).sum(1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1))
        out[:, i] = np.arccos(np.clip(c, -1, 1))
    return out


def _orient(V, T):
    P = V[T]
    ar = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    T = T.copy(); T[ar < 0] = T[ar < 0][:, [0, 2, 1]]
    return T


def _boundary_edges(T):
    from collections import Counter
    E = np.sort(np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]]), axis=1)
    c = Counter(map(tuple, E))
    return {e for e, n in c.items() if n == 1}


def _theta(V, T, bverts, corner):
    """(M2) Theta(z) for every vertex; bverts = set of boundary vertices."""
    nv = len(V)
    ang = _angles(V[T])
    inc = [[] for _ in range(nv)]
    for t in range(len(T)):
        for i in range(3):
            inc[T[t, i]].append((t, i))
    Th = np.full(nv, np.nan); cnt = np.zeros(nv, int)
    for z in range(nv):
        L = inc[z]; cnt[z] = len(L)
        if not L:
            continue
        cen = np.array([V[T[t]].mean(0) for t, _ in L]) - V[z]
        phi = np.arctan2(cen[:, 1], cen[:, 0]); o = np.argsort(phi); phi = phi[o]
        th = np.array([ang[L[j][0], L[j][1]] for j in o])
        if z in bverts:
            gaps = np.diff(np.concatenate([phi, [phi[0] + 2 * np.pi]]))
            s = (np.argmax(gaps) + 1) % len(phi)
            th = np.roll(th, -s)
            pairs = [th[j] + th[j + 1] for j in range(len(th) - 1)]
        else:
            pairs = [th[j] + th[(j + 1) % len(th)] for j in range(len(th))]
        Th[z] = max(abs(np.sin(p)) for p in pairs) if pairs else 0.0
    return Th, cnt, inc, ang


def _edge_split(V, T, t, z, frac=2.0 / 3.0):
    """split triangle t = (z, x, y) and its neighbour across x-y at m = x + frac (y - x) (x = the vertex
    after z in t).  Raises the triangle count at z by one without creating a new boundary triangle
    that is thinner than 2/3 of the old one (the centroid split would make it 1/3 as high)."""
    tri = list(T[t]); i = tri.index(z); x, y = tri[(i + 1) % 3], tri[(i + 2) % 3]
    m = V[x] + frac * (V[y] - V[x])
    V = np.vstack([V, m]); g = len(V) - 1
    nb = [s for s in np.where((T == x).any(1) & (T == y).any(1))[0] if s != t]
    T = T.copy(); T[t] = [z, x, g]; T = np.vstack([T, [[z, g, y]]])
    for s in nb:
        o = [v for v in T[s] if v != x and v != y][0]
        T[s] = [o, y, g]; T = np.vstack([T, [[o, g, x]]])
    return V, T, len(nb) == 0


def _centroid_split(V, T, t):
    c = V[T[t]].mean(0)
    V = np.vstack([V, c]); g = len(V) - 1
    a, b, d = T[t]
    T = np.vstack([T, [[b, d, g], [d, a, g]]]); T[t] = [a, b, g]
    return V, T


def make_mesh(box, curves, N, R_cap=3.5, grad=1.0, extra_lines=(), s_G=None, minangle=30,
              theta_min=0.15, maxfix=10000, area_fac=2.0, layer=np.sqrt(3) / 2):
    """box = (x0, x1, y0, y1); curves = list of Curve; N = chords per curve (int or list);
    extra_lines = list of ((xa, ya), (xb, yb)) interior constraint segments whose endpoints lie on
    the box (used to align a piecewise-polynomial pressure with the mesh)."""
    x0, x1, y0, y1 = box
    Ns = [N] * len(curves) if np.isscalar(N) else list(N)
    cv_pts, cv_par = [], []
    for c, n in zip(curves, Ns):
        tt = c.arclength_params(n); cv_par.append(tt); cv_pts.append(c.pos(tt))
    chords = [np.linalg.norm(np.roll(P, -1, 0) - P, axis=1) for P in cv_pts]
    cG = [ch.mean() for ch in chords]              # per-curve chord length
    sG = s_G or min(cG)
    trees = [cKDTree(c.pos(np.linspace(0, 2 * np.pi, 4000, endpoint=False))) for c in curves]

    def size(X):
        s = np.full(len(X), R_cap * sG)
        for c, tr in zip(cG, trees):
            d, _ = tr.query(X)
            s = np.minimum(s, c * (1 + grad * d))
        return s

    # ---- PSLG
    pts, segs, smark = [], [], []

    def addpt(p):
        pts.append(np.asarray(p, float)); return len(pts) - 1

    corners = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
    must = {0: [], 1: [], 2: [], 3: []}          # side 0 bottom, 1 right, 2 top, 3 left
    for (A, B) in extra_lines:
        for P in (A, B):
            if abs(P[1] - y0) < 1e-12: must[0].append(P[0])
            elif abs(P[0] - x1) < 1e-12: must[1].append(P[1])
            elif abs(P[1] - y1) < 1e-12: must[2].append(P[0])
            elif abs(P[0] - x0) < 1e-12: must[3].append(P[1])
            else: raise ValueError("extra line endpoint not on box")
    cid = [addpt(c) for c in corners]
    ptid = {}
    for side in range(4):
        A = np.array(corners[side]); B = np.array(corners[(side + 1) % 4])
        Ls = np.linalg.norm(B - A); tvec = (B - A) / Ls
        knots = sorted(set([0.0, Ls] + [abs(m - (A[0] if side in (0, 2) else A[1])) for m in must[side]]))
        ids = [cid[side]]
        for ka, kb in zip(knots[:-1], knots[1:]):
            samp = A + np.linspace(ka, kb, 200)[:, None] * tvec
            n = max(1, int(np.ceil((kb - ka) / size(samp).min())))
            for j in range(1, n + 1):
                s = ka + (kb - ka) * j / n
                if j == n and abs(s - Ls) < 1e-12:
                    ids.append(cid[(side + 1) % 4])
                else:
                    q = addpt(A + s * tvec); ids.append(q)
                    ptid[(side, round(s, 12))] = q
        for a, b in zip(ids[:-1], ids[1:]):
            segs.append((a, b)); smark.append(side + 1)
    for (A, B) in extra_lines:                     # interior constraint lines (marker 50)
        A = np.array(A, float); B = np.array(B, float)
        def endpoint_id(P):
            for q in range(len(pts)):
                if np.linalg.norm(pts[q] - P) < 1e-10:
                    return q
            raise ValueError
        ia, ib = endpoint_id(A), endpoint_id(B)
        samp = A + np.linspace(0, 1, 400)[:, None] * (B - A)
        n = max(1, int(np.ceil(np.linalg.norm(B - A) / size(samp).min())))
        ids = [ia] + [addpt(A + (B - A) * j / n) for j in range(1, n)] + [ib]
        for a, b in zip(ids[:-1], ids[1:]):
            segs.append((a, b)); smark.append(50)
    curve_vid = []
    for i, P in enumerate(cv_pts):
        ids = [addpt(p) for p in P]; curve_vid.append(ids)
        for j in range(len(ids)):
            segs.append((ids[j], ids[(j + 1) % len(ids)])); smark.append(100 + i)
    if layer:
        # boundary layer: one guard point over each chord midpoint, at height layer*|e| on the Omega side
        # (the curves are positively oriented and Omega lies outside them, so the Omega side of the chord
        # A->B is its right-hand side).  With layer = sqrt(3)/2 every chord carries an equilateral
        # triangle and every Gamma_h vertex gets 3 triangles.
        for i, P in enumerate(cv_pts):
            Q = np.roll(P, -1, 0); e = Q - P; L_ = np.linalg.norm(e, axis=1)
            nrm = np.stack([e[:, 1], -e[:, 0]], 1) / L_[:, None]
            for g in (P + Q) / 2 + layer * L_[:, None] * nrm:
                addpt(g)
    nin = len(pts)
    holes = []
    for c in curves:
        holes.append(c.centre)
    A = dict(vertices=np.array(pts), segments=np.array(segs), segment_markers=np.array(smark)[:, None],
             holes=np.array(holes))
    Tm = _tr.triangulate(A, f"pq{minangle}Ya{np.sqrt(3) / 4 * (R_cap * sG) ** 2:.12f}")
    for it in range(12):
        V = Tm["vertices"]; T = Tm["triangles"]
        P = V[T]; area = 0.5 * np.abs(np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]))
        sc = size(P.mean(1)); tgt = np.sqrt(3) / 4 * sc ** 2
        dm = np.max(np.stack([np.linalg.norm(P[:, i] - P[:, (i + 1) % 3], axis=1) for i in range(3)], 1), 1)
        long_ = dm > 1.5 * R_cap * sG               # cap the global diameter, so h/hG is stable
        if (area <= area_fac * tgt).all() and not long_.any():
            break
        if _DEBUG: print("refine", it, len(T), long_.sum(), (area > area_fac * tgt).sum(), flush=True)
        tgt = np.where(long_, 0.5 * area, area_fac * tgt)
        B = dict(vertices=V, triangles=T, segments=Tm["segments"],
                 segment_markers=Tm["segment_markers"], triangle_max_area=np.where(area > tgt, tgt, 2 * area + tgt))
        Tm = _tr.triangulate(B, f"rpq{minangle}Ya")
    V = np.array(Tm["vertices"], float); T = _orient(V, np.array(Tm["triangles"]))
    assert np.allclose(V[:nin], np.array(pts)), "Triangle renumbered the input vertices"
    # ---- classify boundary edges (no Steiner points on segments => boundary = input segments)
    gamma_set = {}
    for i, ids in enumerate(curve_vid):
        for j in range(len(ids)):
            a, b = ids[j], ids[(j + 1) % len(ids)]
            gamma_set[tuple(sorted((a, b)))] = i
    seg_side = {tuple(sorted(s)): m for s, m in zip(segs, smark) if m < 50}
    # the box sides may have been split? (Y forbids it) -> verify
    bE = _boundary_edges(T)
    assert all(e in gamma_set or e in seg_side for e in bE), "boundary edge not an input segment"
    assert len(bE) == len(gamma_set) + len(seg_side)
    gverts = set(v for ids in curve_vid for v in ids)
    # ---- repair (M1), (M1'), (M2)
    nfix = dict(gamma_lt3=0, sigma_lt3=0, corner_lt2=0, theta=0)
    corner_set = set(cid)
    for it in range(maxfix):
        bverts = set(v for e in bE for v in e)
        Th, cnt, inc, ang = _theta(V, T, bverts, corner_set)
        bads = []
        for z in range(len(V)):
            if z in corner_set:
                if cnt[z] < 2: bads.append((z, "corner_lt2"))
            elif z in gverts:
                if cnt[z] < 3: bads.append((z, "gamma_lt3"))
            elif z in bverts:
                if cnt[z] < 3: bads.append((z, "sigma_lt3"))
            if z not in corner_set and cnt[z] >= 2 and Th[z] < theta_min and not any(b[0] == z for b in bads[-1:]):
                bads.append((z, "theta"))
        if not bads:
            break
        touched = set()                            # batch: split each triangle at most once per sweep

        def on_bdry(t):
            return any(tuple(sorted((T[t, a_], T[t, (a_ + 1) % 3]))) in bE for a_ in range(3))
        for z, why in bads:
            cand = sorted(inc[z], key=lambda ti: -ang[ti[0], ti[1]])
            if why == "theta":
                # interior vertex: barycentric split of the widest incident triangle with no boundary edge
                # (splitting a triangle that sits on Gamma_h makes it 3x thinner and raises the Nitsche
                # threshold, which is set by the worst boundary star)
                inner = [ti for ti in cand if not on_bdry(ti[0])]
                if inner:
                    t = inner[0][0]
                    if t in touched:
                        continue
                    nt0 = len(T); touched.add(t); nfix[why] += 1
                    V, T = _centroid_split(V, T, t)
                    touched.update([nt0, nt0 + 1])
                    continue
            # boundary vertex (or interior vertex whose triangles all touch the boundary): split the edge
            # opposite z in the widest incident triangle whose opposite edge is interior
            done = False
            for (t, _) in cand:
                tri = list(T[t]); j = tri.index(z); x, y = tri[(j + 1) % 3], tri[(j + 2) % 3]
                if tuple(sorted((x, y))) in bE:
                    continue
                nbs = set(np.where((T == x).any(1) & (T == y).any(1))[0].tolist())
                if nbs & touched:
                    done = True; break
                nt0 = len(T); touched.update(nbs); nfix[why] += 1
                V, T, _ = _edge_split(V, T, t, z)
                touched.update(range(nt0, len(T))); done = True
                break
            if not done:                           # last resort
                t = cand[0][0]
                if t not in touched:
                    nt0 = len(T); touched.add(t); nfix[why] += 1
                    V, T = _centroid_split(V, T, t); touched.update([nt0, nt0 + 1])
        T = _orient(V, T)
    else:
        raise RuntimeError("mesh repair did not terminate")
    # ---- stats
    M = Mesh()
    M.verts, M.tris = V, T
    M.gamma_edges = [(a, b, i) for (a, b), i in gamma_set.items()]
    M.sigma_edges = [(a, b, m) for (a, b), m in seg_side.items()]
    M.corners = cid
    M.curves = curves; M.curve_vid = curve_vid; M.curve_par = cv_par; M.box = box
    P = V[T]
    ang = _angles(P) * 180 / np.pi
    diam = np.max(np.stack([np.linalg.norm(P[:, i] - P[:, (i + 1) % 3], axis=1) for i in range(3)], 1), 1)
    glen = np.concatenate(chords)
    bverts = set(v for e in bE for v in e)
    Th, cnt, _, _ = _theta(V, T, bverts, corner_set)
    gcnt = np.array([cnt[v] for v in gverts])
    xcnt = [cnt[v] for v in cid]
    scnt = [cnt[v] for v in bverts if v not in gverts and v not in corner_set]
    M.stats = dict(ntri=len(T), nvert=len(V), h=float(diam.max()), hG=float(glen.max()),
                   hG_min=float(glen.min()), h_over_hG=float(diam.max() / glen.max()),
                   rho=float(diam.max() / glen.min()), minangle=float(ang.min()), maxangle=float(ang.max()),
                   theta_min=float(np.nanmin(Th)), gamma_min_tris=int(gcnt.min()),
                   sigma_min_tris=int(min(scnt)), corner_min_tris=int(min(xcnt)),
                   fixes=nfix, N=Ns, curves=[c.name for c in curves])
    return M


if __name__ == "__main__":
    import sys, time
    for name, box, curves in [
            ("twodisks", (-3.5, 3.5, -2.5, 2.5), [Circle(-1.5, 0, 0.5), Circle(1.5, 0, 0.8)]),
            ("ellipse", (-2.5, 2.5, -2.5, 2.5), [Ellipse(0, 0, 1.2, 0.7)]),
            ("star", (-2.5, 2.5, -2.5, 2.5), [PolarCurve(0, 0, 1.0, 0.3, 3)])]:
        for N in (16, 32, 64):
            t0 = time.time()
            M = make_mesh(box, curves, N)
            print(name, N, f"{time.time()-t0:.1f}s", M.stats)
