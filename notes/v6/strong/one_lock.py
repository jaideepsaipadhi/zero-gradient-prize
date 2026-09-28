"""Counting test for Theorem T1/T4: mesh (s) with the diagonal of ONE first-layer quad flipped, so that
exactly one wall vertex (z_1) lies in two triangles (locked) and its neighbour z_2 in four.
Prediction (REPORT.md, Thm T4):  E ~ sqrt( A^2 hG^3 + B^2 hG^2 )  -> rate 1 (not 1/2, not 3/2).
Usage: PYPARDISO_MKL_RT=... SVN_OUT=... python3 one_lock.py N [lu|pardiso]
Writes $SVN_OUT/strong/N{N}_A_zero_one.jsonl (via strong_bc._write).
"""
import os, sys, json
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "..", "..", "..", "code"))
import numpy as np
import svn, strong_bc

def make_mesh_one(N, L=2.5):
    verts, tris, circ, outer = svn.make_mesh(N, L=L)
    tris = tris.copy(); m = N // 4
    vid = np.arange((m + 1) * N).reshape(m + 1, N)
    i = 1; i1 = 2
    a, b, c, d = vid[0, i], vid[0, i1], vid[1, i1], vid[1, i]
    tris[2 * i] = [a, b, d]; tris[2 * i + 1] = [b, c, d]
    P = verts[tris]
    ar = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    tris[ar < 0] = tris[ar < 0][:, [0, 2, 1]]
    return verts, tris, circ, outer

strong_bc.make_mesh_alt = make_mesh_one          # build(N, mesh != "std") now uses the one-flip mesh

if __name__ == "__main__":
    N = int(sys.argv[1]); solver = sys.argv[2] if len(sys.argv) > 2 else "lu"
    strong_bc.OUT = os.environ.get("SVN_OUT", os.path.join(ROOT, "runs"))
    for r in strong_bc.run(N, "A", 100.0, "zero", solver, "one"):
        print(json.dumps({k: (float(f"{v:.5g}") if isinstance(v, float) else v) for k, v in r.items()}))
