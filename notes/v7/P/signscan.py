"""notes/v7/P: sign/zero scan of the cell coefficient A1(offset, a, lam) on sheared lattices L_j = (j*delta, j*a),
apex offset of the wall triangle o = -(delta + 1/2) in cell units... reported as off = 0.5 + delta (cell.py convention,
off = 0 is the mirror-symmetric isosceles lattice).  lam = rho * c, c = k(k+1)/a = 20/a.
Checks: A1(-off) = -A1(off) (reflection), A1 != 0 for off != 0, sign, and the residue (lam - c) A1 at lam -> c.
also min over offsets of ||P_cell||_{L2} and ||T_cell||_{L2(e)} (both must be > 0).
usage: python3 signscan.py   (output signscan.log)"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cell import cell_from_lattice
for a in (0.5, 1.0, 2.0):
    c = 20.0 / a
    for rho in (1.05, 1.5, 3.0, 10.0, 100.0):
        lam = rho * c; row = []; pl = []; tr = []
        for off in (-0.75, -0.5, -0.25, -0.1, 0.0, 0.1, 0.25, 0.5, 0.75):
            cp = cell_from_lattice(a, off - 0.5, lam, J=int(np.ceil(10 / a)) + 2)
            o, _, _ = cp.run("bump"); row.append(o["A"]); pl.append(o["PL2"]); tr.append(o["trac_L2"])
        row = np.array(row)
        print(json.dumps(dict(a=a, rho=rho, lam=lam, A1=[float(f"{v:.6g}") for v in row],
                              odd_defect=float(np.max(np.abs(row + row[::-1]))),
                              min_abs_offzero=float(np.min(np.abs(np.delete(row, 4)))),
                              min_PL2=float(min(pl)), min_trac=float(min(tr)), signs="".join("+" if v > 1e-10 else ("-" if v < -1e-10 else "0") for v in row))), flush=True)
