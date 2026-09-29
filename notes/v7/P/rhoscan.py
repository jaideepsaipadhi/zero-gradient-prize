"""notes/v7/P: fine scan of A1 in lam > c = 20/a at the production-type offset (cell delta = -1, off = -1/2) and at
off = -0.25, a in {0.5,1,2,3}: does A1 cross zero above the resonance?  Also (lam - c) A1 -> residue at lam -> c+.
usage: python3 rhoscan.py   (output rhoscan.log)"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cell import cell_from_lattice
rhos = np.concatenate([[1.001, 1.01], np.geomspace(1.05, 1e4, 26)])
for a in (0.5, 1.0, 2.0, 3.0):
    c = 20.0 / a
    for off in (-0.5, -0.25):
        vals = []
        for rho in rhos:
            cp = cell_from_lattice(a, off - 0.5, rho * c, J=int(np.ceil(10 / a)) + 2)
            vals.append(cp.run("bump")[0]["A"])
        vals = np.array(vals)
        print(json.dumps(dict(a=a, off=off, rho=[float(f"{r:.4g}") for r in rhos], A1=[float(f"{v:.5g}") for v in vals],
                              min_A1=float(vals.min()), argmin_rho=float(rhos[vals.argmin()]),
                              residue=float((rhos[0] - 1) * c * vals[0]), all_positive=bool((vals > 0).all()))), flush=True)
