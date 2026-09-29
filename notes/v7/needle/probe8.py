"""Bent slit: split w along the polyline z - w - f' where f' is the DIAGONAL neighbour of w (the two needle axes then
meet at w at 135 deg instead of 180 deg), via notes/v6/nearsing/ns_mesh.split_apex(diag=True).  Compare with the straight
slit (diag=False).  Prediction (REPORT Prop. B / Remark): the cross-difference functional has ratio ~ eps only for a
straight kite; the eps-modes disappear for a bent kite.  Usage: python3 probe8.py n"""
import os, sys, json, time, numpy as np
from slit_lab import *
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "v6", "nearsing"))
import ns_mesh

def run(n, eps, diag):
    t0 = time.time()
    verts, tris, vid = square_mesh(n)
    i0 = n // 2
    w = int(vid[i0, n // 2]); z = int(vid[i0 - 1, n // 2])
    V, T = ns_mesh.split_apex(verts, np.array(tris), z, eps, diag=diag, w=w)
    V, T = ns_mesh._orient(V, np.array(T))
    S = svn.Space(V, T, {int(vid[0, 0]), int(vid[1, 0])}, None)
    P = V[T]; ar = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]) / 2
    assert abs(ar.sum() - 4.0) < 1e-10 and (ar > 0).all()
    Sd, M, one, sysm, _ = schur(S)
    E = interior_basis(S)
    b, _ = beta_sub(Sd, M, E, k=4)
    ang = ns_mesh.angles(V, T)
    print(json.dumps(dict(n=n, eps=eps, diag=diag, min_angle_deg=round(float(np.degrees(ang.min())), 4),
                          max_angle_deg=round(float(np.degrees(ang.max())), 2),
                          beta=[float(f"{x:.4g}") for x in b], beta_over_eps=[float(f"{x/eps:.4g}") for x in b],
                          beta_over_sqrteps=[float(f"{x/np.sqrt(eps):.4g}") for x in b], secs=round(time.time() - t0, 1))), flush=True)

if __name__ == "__main__":
    n = int(sys.argv[1])
    for diag in (False, True):
        for eps in (1e-2, 1e-3, 1e-4):
            run(n, eps, diag)
