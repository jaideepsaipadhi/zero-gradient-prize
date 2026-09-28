"""per-star detail of the production-star seminorm (robust4 REPORT Sec. 3.4): absolute singular values, kappa_KerT,
kappa_K, and the polar angle of z, to locate small third singular values.  python prod_detail.py N"""
import sys, json
import numpy as np
import stars as S
import svn

N = int(sys.argv[1]); v, t, circ, outer = svn.make_mesh(N, L=2.5)
rows = []
for i in range(N):
    z = i; fan = [tr for tr in t if z in tr]
    loc = sorted({int(a) for tr in fan for a in tr}); mp = {z: 0}
    for a in loc:
        if a != z: mp[a] = len(mp)
    P = np.zeros((len(mp), 2))
    for a, b in mp.items(): P[b] = v[a]
    g1 = mp[(i + 1) % N]; g2 = mp[(i - 1) % N]; S._ns["GAMMA"] = (g1, g2)
    fl = [[mp[int(a)] for a in tr] for tr in fan]
    M, info = S._star_forms_M(P, fl, 4)
    Te = [tr for tr in fl if 0 in tr and g1 in tr][0]; apex = [a for a in Te if a not in (0, g1)][0]
    T = P[[0, g1, apex]]
    ev = np.sqrt(np.clip(np.linalg.eigvalsh(M), 0, None))
    ang = np.degrees(np.arctan2(P[0][1], P[0][0]))
    rows.append(dict(i=i, theta=float(ang), sig=[float(x) for x in ev], kKer=S.kappa(M, S.KerT(T)), kK=S.kappa(M, S.KT(T)),
                     ntri=len(fl)))
    print(json.dumps(rows[-1]))
