"""Signed apex offset of the first-layer triangles T_e (projection of the apex onto the chord, minus the
chord midpoint, over |e|; positive = counterclockwise) on the production meshes.  usage: python3 apex_offsets.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../code"))
import numpy as np, svn
for N in (16, 32, 64, 128):
    v, t, c, o = svn.make_mesh(N); S = svn.Space(v, t, c, o)
    offs = []
    for (tt, a, b) in S.bedges:
        A, B = v[t[tt, a]], v[t[tt, b]]; C = v[t[tt, 3 - a - b]]
        L = np.linalg.norm(B - A); tv = (B - A) / L; m = (A + B) / 2
        sgn = 1 if A[0] * B[1] - A[1] * B[0] > 0 else -1
        offs.append(np.dot(C - m, tv) / L * sgn)
    offs = np.array(offs)
    print(f"N={N:4d}: apex offset/|e|: min {offs.min():+.3f} max {offs.max():+.3f} mean|.| {np.abs(offs).mean():.3f}")
