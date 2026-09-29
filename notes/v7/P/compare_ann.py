"""notes/v7/P: annulus (exactly locally periodic) measurements vs the cell predictions (Theorem 3.2 of REPORT.tex):
  P_h(cos)/hG      -> 2 pi A1                        (S1: W = 2(1 + cos th))
  ||e_p - mean||^2 ~  hG^3 int || W P1 + bar_t0 Pc ||^2 dtheta   (P1 bump cell, Pc traction cell)
  ||T_h - sigma n||^2 ~ hG^2 int || W T1 + bar_t0 Tc ||^2 dtheta
Richardson column: 2 f(N) - f(N/2) (first-order correction).  usage: python3 compare_ann.py  (output compare_ann.log)"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from cell import cell_from_lattice
S = np.sqrt(12 * np.pi)
cases = [("ann_rect_S1_a1_l30.log", 1.0, 30.0, -1.0), ("annS_rect_a1_l30.log", 1.0, 30.0, -1.0),
         ("ann_rect_S1_a2.log", 2.0, 30.0, -1.0), ("ann_rect_S1_l100.log", 1.0, 100.0, -1.0),
         ("ann_sym_a1_l30.log", 1.0, 30.0, -0.5), ("annS_sym_a1_l30.log", 1.0, 30.0, -0.5)]
for f, a, lam, dc in cases:
    if not os.path.exists(os.path.join(HERE, f)):
        continue
    rs = [json.loads(l) for l in open(os.path.join(HERE, f)) if l.startswith("{")]
    if not rs:
        continue
    cp = cell_from_lattice(a, dc, lam, J=int(np.ceil(10 / a)) + 2)
    o = cp.run("bump")[0]; oc = cp.run("trac")[0]
    tb = 2 * o["A"] / (1 - oc["A"])            # bar t0 = int W A1 / int (1 - Ac),  W = 2(1+cos)
    def norms(s):                              # cell norms of the load F1 + s Fc
        g1, t1 = cp.unit_data("bump"); gc, tc = cp.unit_data("trac")
        U, P, _ = cp.solve(cp.rhs_edge(g1 + s * gc, t1 + s * tc)); r = cp.analyse(U, P)
        return r["PL2"] ** 2, r["trac_L2"] ** 2
    n0 = norms(0.0); n1 = norms(1.0); nc = (oc["PL2"] ** 2, oc["trac_L2"] ** 2)
    ip = [(n1[i] - n0[i] - nc[i]) / 2 for i in range(2)]
    # int_0^{2pi} || W X1 + tb Xc ||^2 dtheta,  int W^2 = 12 pi, int W = 4 pi
    full = [12 * np.pi * n0[i] + 2 * tb * 4 * np.pi * ip[i] + 2 * np.pi * tb ** 2 * nc[i] for i in range(2)]
    pred = dict(P=2 * np.pi * o["A"], L2=np.sqrt(full[0]), T=np.sqrt(full[1]))
    print(f"  [t0 term: bar_t0={tb:.4f}; without it L2/hG^1.5={S * o['PL2']:.4f}, T/hG={S * o['trac_L2']:.4f}]")
    print(f"{f}: a={a} lamhat={lam} rout={rs[0].get('rout', 2.5)}  predicted P/hG={pred['P']:.4f}  "
          f"L2/hG^1.5={pred['L2']:.4f}  T/hG={pred['T']:.4f}")
    prev = None
    for r in rs:
        hG = r["hG"]; m = dict(P=r["P_cos1"] / hG, L2=r["ep_L2"] / hG ** 1.5, T=r["trac_L2"] / hG,
                               Psin2=r["P_sin1"] / hG ** 2)
        s = f"  N={r['N']:4d} hG={hG:.4f}  P/hG={m['P']:+.4f}  L2/hG^1.5={m['L2']:.4f}  T/hG={m['T']:.4f}  Psin/hG^2={m['Psin2']:+.4f}"
        if prev is not None and abs(prev[0] / hG - 2) < 0.01:
            s += "  | Richardson: " + "  ".join(f"{k}={2 * m[k] - prev[1][k]:+.4f} ({(2 * m[k] - prev[1][k]) / pred[k]:.3f})"
                                              if abs(pred[k]) > 1e-9 else f"{k}={2 * m[k] - prev[1][k]:+.4f}" for k in ("P", "L2", "T"))
        print(s, flush=True); prev = (hG, m)
