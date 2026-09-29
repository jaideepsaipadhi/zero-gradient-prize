"""notes/v7/P: cell prediction of lim P_h(phi)/h_G on the self-similar annulus meshes of annulus.py.
Every wall edge has the same cell (aspect a, shear delta_cell, lamhat), so
  P_h(phi)/h_G -> A1 * int_0^{2pi} W(theta) phi(theta) dtheta + t0-term (zero for mean-free phi, Ac constant),
W = d_r(u . tau_c) = 2(1+cos th) for S1 (tau_c = -e_theta).  Annulus delta = 0 (apex over the ccw endpoint)
is the cell lattice delta_cell = -1; annulus delta = -1/2 is delta_cell = -1/2 (mirror symmetric).
usage: python3 annpred.py        (output annpred.log)"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cell import cell_from_lattice
for (a, lam, dcell, lab) in [(1.0, 30.0, -1.0, "rect"), (2.0, 30.0, -1.0, "rect"), (1.0, 100.0, -1.0, "rect"),
                              (1.0, 25.0, -1.0, "rect"), (1.0, 30.0, -0.5, "sym")]:
    c = cell_from_lattice(a, dcell, lam, J=int(np.ceil(10 / a)) + 2)
    o, _, _ = c.run("bump"); oc, _, _ = c.run("trac")
    print(json.dumps(dict(a=a, lamhat=lam, lattice=lab, c_res=20.0 / a, A1=o["A"], Ac=oc["A"],
                          Ac_formula=(20.0 / a) / (20.0 / a - lam), pred_S1_cos1_over_hG=2 * np.pi * o["A"],
                          pred_S1_sin1_over_hG=0.0, PL2=o["PL2"], Pw_odd=o["Pw_odd"], trac=o["trac_L2"])), flush=True)
