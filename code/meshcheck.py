import svn, numpy as np
from collections import defaultdict
def theta_min(v, t, circ, outer):
    inc = defaultdict(list)
    for T in t:
        for k in range(3): inc[T[k]].append(T)
    sing = []; th = []
    for z, Ts in inc.items():
        dirs = set()
        for T in Ts:
            for w in T:
                if w != z:
                    d = v[w] - v[z]; d = d / np.linalg.norm(d)
                    if d[0] < -1e-12 or (abs(d[0]) < 1e-12 and d[1] < 0): d = -d
                    dirs.add((round(d[0], 7), round(d[1], 7)))
        if len(dirs) <= 2: sing.append(z)
    return sing
for N in (16, 24, 32, 48, 64, 96, 128):
    v, t, c, o = svn.make_mesh(N)
    s = theta_min(v, t, c, o)
    print(N, len(t), 'singular vertices:', len(s), 'circle', sum(z in c for z in s), 'outer', sum(z in o for z in s),
          [tuple(np.round(v[z], 3)) for z in s[:6]])
