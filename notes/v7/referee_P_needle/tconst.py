"""Referee check: does the traction prediction of Thm 3.2(iii) need the global pressure constant?
Theta = W T1 + tb (Tc + s n): s=0 report formula; s=+-1 add a constant normal traction tb*n."""
import os, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "P"))
from cell import cell_from_lattice
import svn
for a, lam, meas in [(1.0, 30.0, 3.6087), (2.0, 30.0, 3.8082), (1.0, 100.0, 6.7941)]:
    cp = cell_from_lattice(a, -1.0, lam, J=int(np.ceil(10 / a)) + 2)
    o = cp.run("bump")[0]; oc = cp.run("trac")[0]
    tb = 2 * o["A"] / (1 - oc["A"])
    (t, ref, Xe, we, n) = cp.wall
    def tracfield(s):
        g1, t1 = cp.unit_data("bump"); U, P, _ = cp.solve(cp.rhs_edge(g1, t1))
        gc, tc = cp.unit_data("trac"); Uc, Pc, _ = cp.solve(cp.rhs_edge(gc, tc))
        def Tr(U, P):
            S = cp.S; q3 = svn._mono(svn.MON3, ref[:, 0], ref[:, 1]); Pw = q3 @ P[10*t:10*t+10]
            dX, dY = svn.dbasis(ref[:, 0], ref[:, 1]); Ji = S.Jinv[t]
            gx = (Ji[0,0]*dX + Ji[1,0]*dY) @ U[svn.dofs(S.ids[t])]; gy = (Ji[0,1]*dX + Ji[1,1]*dY) @ U[svn.dofs(S.ids[t])]
            G = np.stack([gx, gy], -1); return np.einsum('gcd,d->gc', G + G.transpose(0,2,1), n) - Pw[:, None]*n[None, :]
        return Tr(U, P), Tr(Uc, Pc) + s * n[None, :]
    res = []
    for s in (0.0, 1.0, -1.0):
        T1, Tc = tracfield(s)
        th = np.linspace(0, 2*np.pi, 721)[:-1]; W = 2*(1+np.cos(th))
        tot = sum(we @ ((Wk*T1 + tb*Tc)**2).sum(-1) for Wk in W) * (2*np.pi/len(th))
        res.append(np.sqrt(tot))
    print(f"a={a} lam={lam} tb={tb:.4f} Ac={oc['A']:.4f}  pred s=0:{res[0]:.4f} s=+1:{res[1]:.4f} s=-1:{res[2]:.4f}  meas(Richardson)={meas}")
