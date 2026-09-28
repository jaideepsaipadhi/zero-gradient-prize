"""Floating-point cross-check of the closed-form witness against logs/lower_bound_stars.log. Run: python3 crosscheck_meshes.py"""
from witness_lib import *
ref = {int(l.split()[0]): float(l.split()[3]) for l in open('/home/claude/zero-gradient-prize/logs/lower_bound_stars.log') if l.split() and l.split()[0].isdigit()}
for N in (8, 12, 16, 32, 64, 256, 1024):
    v, t, c, o = svn.make_mesh(N)
    ks = []
    for i in range(N):
        tri = [2 * i, 2 * i + 1, 2 * ((i - 1) % N)]
        z = i   # circle vertex i is vid[0,i] = i
        nb = set()
        for tt in tri: nb |= set(t[tt].tolist())
        nb.discard(z)
        P = v[list(nb)] - v[z]
        ang = np.arctan2(P[:, 1], P[:, 0])
        chords = [q for q in nb if q in c]
        # start the ccw fan at the chord such that the fan (interior of Omega_h) is swept ccw
        best = None
        for q0 in chords:
            a0 = math.atan2(*(v[q0] - v[z])[::-1])
            rel = sorted((0.0 if q == q0 else (a - a0) % (2 * math.pi), q) for a, q in zip(ang, nb))
            if rel[-1][1] in c:
                best = [q for _, q in rel]
        pts = v[best] - v[z]
        rot = math.atan2(pts[0][1], pts[0][0])
        R = np.array([[math.cos(-rot), -math.sin(-rot)], [math.sin(-rot), math.cos(-rot)]])
        pts = pts @ R.T
        M, nr = kappa(pts)
        # e_z as in lower_bound_tests: the chord of triangle 2i
        ez = [q for q in t[2 * i] if q in c and q != z][0]
        ks.append(M / nr / (np.linalg.norm(v[ez] - v[z]) / np.linalg.norm(pts[0]))**2)
    ks = np.abs(ks)
    print(f"N={N:4d}: witness kappa min {ks.min():.5f} mean {ks.mean():.5f} max {ks.max():.5f}   (log, optimum over numerically computed Y^1(star): min {ref[N]:.5f})")
