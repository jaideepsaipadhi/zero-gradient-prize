"""Meshes with deliberately nearly singular wall stars (m = 3, one needle at z).

Base meshes: 'one' (notes/v6/strong/one_lock.py: one locked wall vertex z_1 with apex w = d_1) and
'alt' (strong_bc.make_mesh_alt: every odd wall vertex locked, apex d_i).  At every locked z we then either
  II  ('split'): split the apex w into two vertices x_L, x_R = w -+ (delta/2) tau, tau _|_ (w - z),
                 delta = 2|w-z| tan(eps/2).  The star of z becomes (z_-, z, x_L), (z, x_R, x_L) [needle,
                 angle eps at z], (z, z_+, x_R); a second needle (x_L, x_R, f) closes the slit (f = the
                 neighbour of w on the ray z->w).  eps -> 0 recovers the lock; eps ~ 1 an (M2) star.
  I   ('wall'):  insert y = z_+ + eta * u (u = inward bisector of the angle of (z, z_+, w) at z_+), chosen so
                 that the angle between z->z_+ and z->y is eps; (z,z_+,w) -> (z, z_+, y) [needle at the wall],
                 (z, y, w), (y, z_+, w).
Mesh names: '{one|alt}{split|wall}' ; eps is given separately.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
import numpy as np
import svn, strong_bc


def make_mesh_one(N, L=2.5):          # identical to notes/v6/strong/one_lock.py
    verts, tris, circ, outer = svn.make_mesh(N, L=L)
    tris = tris.copy(); m = N // 4
    vid = np.arange((m + 1) * N).reshape(m + 1, N)
    i = 1; i1 = 2
    a, b, c, d = vid[0, i], vid[0, i1], vid[1, i1], vid[1, i]
    tris[2 * i] = [a, b, d]; tris[2 * i + 1] = [b, c, d]
    return _orient(verts, tris), circ, outer


def _orient(verts, tris):
    tris = np.array(tris)
    P = verts[tris]
    ar = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    tris[ar < 0] = tris[ar < 0][:, [0, 2, 1]]
    return verts, tris


def locked(tris, circ):
    inc = {}
    for t, tri in enumerate(tris):
        for v in tri:
            if int(v) in circ:
                inc.setdefault(int(v), []).append(t)
    return {z: T for z, T in inc.items() if len(T) == 2}


def split_apex(verts, tris, z, eps, diag=False, w=None):
    """split vertex w (default: the apex of the locked wall vertex z) along the ray z->w into a slit."""
    if w is None:
        T1, T2 = [t for t in range(len(tris)) if z in tris[t]]
        w = (set(tris[T1]) & set(tris[T2])) - {z}
        w = w.pop()
    d = verts[w] - verts[z]; r = np.linalg.norm(d); tau = np.array([-d[1], d[0]]) / r
    delta = 2 * r * np.tan(eps / 2)
    # neighbour f of w on the ray z->w (beyond w)
    nb = set()
    for t in np.where((tris == w).any(1))[0]:
        nb |= set(tris[t].tolist())
    nb -= {w}
    f = min(nb, key=lambda v: abs(np.cross(d / r, verts[v] - verts[w])) - 1e-9 * np.dot(d, verts[v] - verts[w]))
    assert np.dot(verts[f] - verts[w], d) > 0 and abs(np.cross(d / r, verts[f] - verts[w])) < 1e-9 * r, "no ray neighbour"
    if diag:        # use instead the neighbour beyond w making the largest positive angle < 90 deg with d
        cand = [v for v in nb if np.dot(verts[v] - verts[w], d) > 0 and v != f]
        f = max(cand, key=lambda v: abs(np.cross(d / r, (verts[v] - verts[w]) / np.linalg.norm(verts[v] - verts[w]))))
    az = np.arctan2(*(verts[z] - verts[w])[::-1]); af = np.arctan2(*(verts[f] - verts[w])[::-1])
    at = np.arctan2(tau[1], tau[0])
    def arc(a):                          # True if angle a lies on the arc from az to af (ccw) that contains tau
        ccw = lambda x: (x - az) % (2 * np.pi)
        inside = ccw(a) < ccw(af)
        return inside == (ccw(at) < ccw(af))
    xL = len(verts); xR = w                   # left copy is new, right copy reuses index w
    V = np.vstack([verts, verts[w]])
    newt = []
    for t in range(len(tris)):
        tri = tris[t]
        if w in tri:
            cen = verts[tri].mean(0) - verts[w]
            if arc(np.arctan2(cen[1], cen[0])):
                tri = np.where(tri == w, xL, tri)
        newt.append(list(tri))
    V[xL] = verts[w] + 0.5 * delta * tau
    V[xR] = verts[w] - 0.5 * delta * tau
    newt += [[z, xR, xL], [xL, xR, f]]
    return V, newt


def wall_needle(verts, tris, z, eps):
    T = [t for t in range(len(tris)) if z in tris[t]]
    # triangle (z, z_+, w): take the one whose other vertices include the ccw wall neighbour; choose the
    # triangle with the larger polar angle of its wall neighbour
    tri_info = []
    for t in T:
        o = [v for v in tris[t] if v != z]
        wall = [v for v in o if abs(np.linalg.norm(verts[v]) - 1) < 1e-12]
        apex = [v for v in o if v not in wall]
        tri_info.append((t, wall[0], apex[0]))
    t, zp, w = tri_info[0]
    Z, ZP, Wv = verts[z], verts[zp], verts[w]
    u1 = (Z - ZP) / np.linalg.norm(Z - ZP); u2 = (Wv - ZP) / np.linalg.norm(Wv - ZP)
    u = (u1 + u2) / np.linalg.norm(u1 + u2)
    # choose eta with angle(z->z_+, z->y) = eps  (bisection)
    ez = (ZP - Z) / np.linalg.norm(ZP - Z)
    ang = lambda eta: np.arccos(np.clip(np.dot(ez, (ZP + eta * u - Z) / np.linalg.norm(ZP + eta * u - Z)), -1, 1))
    lo, hi = 0.0, 0.5 * np.linalg.norm(ZP - Z)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if ang(mid) < eps else (lo, mid)
    y = len(verts)
    V = np.vstack([verts, ZP + 0.5 * (lo + hi) * u])
    newt = [list(tr) for k, tr in enumerate(tris) if k != t]
    newt += [[z, zp, y], [z, y, w], [y, zp, w]]
    return V, newt


def make(name, N, eps):
    base = name[:3]; kind = name[3:]
    if base == "one":
        (verts, tris), circ, outer = make_mesh_one(N)
    else:
        verts, tris, circ, outer = strong_bc.make_mesh_alt(N)
    lk = sorted(locked(tris, circ))
    tris = [list(t) for t in tris]
    for z in lk:
        tr = np.array(tris)
        if kind in ("split", "splitd"):
            verts, tris = split_apex(verts, tr, z, eps, diag=(kind == "splitd"))
        elif kind == "wall":
            verts, tris = wall_needle(verts, tr, z, eps)
        elif kind == "":
            pass
        else:
            raise ValueError(kind)
    verts, tris = _orient(verts, np.array(tris))
    return verts, tris, circ, outer, lk


def angles(verts, tris):
    P = verts[tris]; out = np.zeros((len(tris), 3))
    for a in range(3):
        u = P[:, (a + 1) % 3] - P[:, a]; v = P[:, (a + 2) % 3] - P[:, a]
        out[:, a] = np.arccos(np.clip(np.sum(u * v, 1) / np.linalg.norm(u, axis=1) / np.linalg.norm(v, axis=1), -1, 1))
    return out


def check(verts, tris, circ):
    P = verts[tris]
    ar = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    assert (ar > 0).all()
    # total area = square - polygon
    A = ar.sum() / 2
    ang = angles(verts, tris)
    stars = {}
    for t, tri in enumerate(tris):
        for a, v in enumerate(tri):
            if int(v) in circ:
                stars.setdefault(int(v), []).append(ang[t, a])
    ms = sorted(set(len(s) for s in stars.values()))
    return dict(area=A, minang=float(np.degrees(ang.min())), maxang=float(np.degrees(ang.max())), wall_m=ms)


if __name__ == "__main__":
    for name in ("onesplit", "onewall", "altsplit", "altwall"):
        for N, eps in ((16, 0.3), (32, 0.01)):
            v, t, c, o, lk = make(name, N, eps)
            print(name, N, eps, len(lk), check(v, t, c))
