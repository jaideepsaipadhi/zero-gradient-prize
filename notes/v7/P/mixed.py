"""notes/v7/P: does A1 depend only on the first-layer triangle T_e?  Mixed columns (a = 1, lam = 30, split 'f'):
  X: T_e isosceles (mirror-symmetric first triangle), rows >= 1 rectangular (production type);
  Y: T_e right-angled (production type), rows >= 1 isosceles;
  R: all rectangular; S: all isosceles (reference).
usage: python3 mixed.py   (output mixed.log)"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cell import CellProblem
J, a, lam = 12, 1.0, 30.0
cols = {"R": [-j for j in range(J + 1)], "S": [-j / 2 for j in range(J + 1)],
        "X": [0.0] + [-0.5 - (j - 1) for j in range(1, J + 1)], "Y": [0.0] + [-1.0 - (j - 1) / 2 for j in range(1, J + 1)]}
for k, xs in cols.items():
    rows = np.array([[x, j * a] for j, x in enumerate(xs)])
    c = CellProblem(rows, 'f', lam)
    o, _, _ = c.run("bump"); oc, _, _ = c.run("trac")
    print(json.dumps(dict(col=k, A1=o["A"], Ac=oc["A"], PL2=o["PL2"], trac=o["trac_L2"])), flush=True)
