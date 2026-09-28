"""Sanity checks for notes/v4/pressure_drag.tex (pressure layer, leak drag, CNS* geometric drag term).
Uses code/svn.py (k=4, L=2.5, meshes of [P, Section 13]); direct solves with SuperLU.
Usage (from repo root):
  python3 notes/v4/pressure_drag_check.py main  KIND MU N1,N2,...   # pressure errors + J^N, J^beta, Theta_f   (KIND in P,B1,A,C)
  python3 notes/v4/pressure_drag_check.py xi    KIND MU               # Xi_h = int xi_b n (pointwise-traction defect of GS)
  python3 notes/v4/pressure_drag_check.py K     KIND MU               # K_e = int (p-pbar) R_e via CNS* dual + volume formula
  python3 notes/v4/pressure_drag_check.py lam   KIND MU ND            # Lambda_h(e) from a CNS* dual on mesh ND
  python3 notes/v4/pressure_drag_check.py chi   KIND MU               # Gjerde-Scott volume functional omega_h(chi_h)
  python3 notes/v4/pressure_drag_check.py torque KIND MU              # torque (psi = x^perp): J^N error vs -Lambda_h
"""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "code"))
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
import svn
from svn import *

def solve_gen(S, sysm, g, theta, extra_rhs=None):
    A, B, F = sysm["A"], sysm["B"], sysm["F"].copy()
    if extra_rhs is not None: F = F + extra_rhs
    Bt = coupling(sysm, theta)
    ndof = A.shape[0]
    dD = dofs(S.dnodes).ravel(); free = np.setdiff1d(np.arange(ndof), dD)
    ug = np.zeros(ndof); ug[dD] = g(S.xy[S.dnodes, 0], S.xy[S.dnodes, 1]).ravel()
    K = sp.bmat([[A[free][:, free], -Bt[:, free].T], [B[:, free], None]], format="csc")
    rhs = np.concatenate([F[free] - A[free][:, dD] @ ug[dD], -B[:, dD] @ ug[dD]])
    lu = spl.splu(K, permc_spec="COLAMD", diag_pivot_thresh=1.0)
    x = lu.solve(rhs)
    for _ in range(3): x = x + lu.solve(rhs - K @ x)
    u = ug.copy(); u[free] = x[:len(free)]; p = x[len(free):]
    return u, p

def edge_eval(S, u, p, t, ref):
    ph = basis(ref[:, 0], ref[:, 1]); dX, dY = dbasis(ref[:, 0], ref[:, 1])
    Ji = S.Jinv[t]
    gxe = Ji[0, 0] * dX + Ji[1, 0] * dY; gye = Ji[0, 1] * dX + Ji[1, 1] * dY
    Ut = u[dofs(S.ids[t])]
    uh = ph @ Ut; gx = gxe @ Ut; gy = gye @ Ut      # gx[:,c] = d_x u_c
    q3 = svn._mono(MON3, ref[:, 0], ref[:, 1]); ph_ = q3 @ p[10*t:10*t+10]
    return uh, gx, gy, ph_

def dirichlet_rhs(S, mu, vec):
    ndof = 2*S.nn; R = np.zeros(ndof)
    for (t, ref, Xe, we, n) in S.edge_data():
        ph = basis(ref[:, 0], ref[:, 1]); dX, dY = dbasis(ref[:, 0], ref[:, 1])
        Ji = S.Jinv[t]
        gxe = Ji[0, 0] * dX + Ji[1, 0] * dY; gye = Ji[0, 1] * dX + Ji[1, 1] * dY
        dn = n[0]*gxe + n[1]*gye
        d = dofs(S.ids[t])
        for c in range(2):
            R[d[:, c]] += we @ (-vec[c]*dn + (mu/S.h)*vec[c]*ph)
    return R

def analyse(N, kind, mu, theta, R_dual=None):
    ex = make_exact(kind)
    verts, tris, circ, outer = make_mesh(N)
    S = Space(verts, tris, circ, outer)
    sysm = assemble(S, mu, ex["f"], theta)
    u, p = solve_gen(S, sysm, ex["u"], theta)
    nt = len(tris)
    # pressure errors
    X = S.phys(QX, QY); wq = QW[None, :]*np.abs(S.detJ)[:, None]
    ph = np.einsum('qj,tj->tq', Q3_Q, p.reshape(nt, 10))
    pex = ex["p"](X[..., 0], X[..., 1])
    # Gamma_h mean of exact p
    Ed = S.edge_data()
    num = sum((we*ex["p"](Xe[:,0],Xe[:,1])).sum() for (t,ref,Xe,we,n) in Ed); per = sum(we.sum() for (t,ref,Xe,we,n) in Ed)
    pbar = num/per
    err = pex - pbar - ph
    layer = np.zeros(nt, bool); layer[[t for (t,_,_) in S.bedges]] = True
    tot = np.sqrt((wq*err**2).sum()); mean = (wq*err).sum()/wq.sum()
    modc = np.sqrt((wq*(err-mean)**2).sum())
    lay = np.sqrt((wq[layer]*(err[layer]-mean)**2).sum()); rest = np.sqrt((wq[~layer]*(err[~layer]-mean)**2).sum())
    # xi_b prediction: representer of q -> int_e g q, g = (h/mu)(p - pbar)(x*)
    xib = np.zeros((nt, len(QW)))
    xinorm2 = 0.0
    for (t, ref, Xe, we, n) in Ed:
        xs = Xe/np.linalg.norm(Xe, axis=1)[:, None]
        g = (S.h/mu)*(ex["p"](xs[:,0], xs[:,1]) - pbar)
        q3 = svn._mono(MON3, ref[:, 0], ref[:, 1])
        b = q3.T @ (we*g); c = np.linalg.solve(sysm["Ml"][t], b)
        xib[t] = Q3_Q @ c; xinorm2 += b @ c
    rho = err - xib; rhom = (wq*rho).sum()/wq.sum()
    rhon = np.sqrt((wq*(rho-rhom)**2).sum())
    linf_layer = np.abs(err[layer]-mean).max()
    # drags
    out = dict(N=N, kind=kind, mu=mu, th=theta, h=S.h, hG=S.hGamma, pL2=tot, pL2modc=modc, player=lay, prest=rest,
               xib_pred=np.sqrt(xinorm2), rho=rhon, pLinf_layer=linf_layer, pmean_err=mean)
    for e in ([1.0,0.0],[0.0,1.0]):
        e = np.array(e)
        JN = 0.0; Jb = 0.0; Theta = 0.0
        for (t, ref, Xe, we, n) in Ed:
            uh, gx, gy, pht = edge_eval(S, u, p, t, ref)
            dnu = gx*n[0] + gy*n[1]                  # (ng,2): d_n u_c
            # D(u)n: (grad u + grad u^T) n ; grad u[c,d] = d_d u_c ; (grad u) n = dnu; (grad u^T n)_c = sum_d d_c u_d n_d
            gT = np.stack([gx[:,0]*n[0]+gx[:,1]*n[1], gy[:,0]*n[0]+gy[:,1]*n[1]], -1)
            pmean_h = 0.0
            JN += (we*((dnu - (mu/S.h)*uh) @ e)).sum() - theta*(we*pht).sum()*(n@e)
            Jb += (we*(((dnu+gT) @ e) - pht*(n@e))).sum()
        # Theta_f = int_{S_h} f . e  over segments between chords and arcs
        th_edges = []
        for (t, u_, w_) in S.bedges:
            A_ = verts[tris[t,u_]]; B_ = verts[tris[t,w_]]
            a0 = np.arctan2(A_[1],A_[0]); a1 = np.arctan2(B_[1],B_[0])
            da = (a1-a0+np.pi)%(2*np.pi)-np.pi
            gth, gw = np.polynomial.legendre.leggauss(12)
            for s_, ws in zip((gth+1)/2, gw/2):
                ang = a0 + s_*da
                dirv = np.array([np.cos(ang), np.sin(ang)])
                # chord radius along this ray
                # solve for r: point r*dirv on segment AB
                M = np.array([dirv, A_-B_]).T; r_, lam = np.linalg.solve(M, A_)
                rr, rw = np.polynomial.legendre.leggauss(6)
                rs = r_ + (1-r_)*(rr+1)/2; rws = (1-r_)*rw/2
                P_ = rs[:,None]*dirv[None,:]
                fv = ex["f"](P_[:,0], P_[:,1])
                Theta += abs(da)*ws*(rws*rs*(fv@e)).sum()
        # exact J_e on the unit circle
        th = np.linspace(0, 2*np.pi, 4001)[:-1]; xc = np.cos(th); yc = np.sin(th)
        nn = -np.stack([xc, yc], 1)
        Gm = ex["gu"](xc, yc); pv = ex["p"](xc, yc)
        Du = Gm + np.transpose(Gm, (0,2,1))
        tr = np.einsum('kij,kj->ki', Du, nn) - pv[:,None]*nn
        Je = (tr @ e).mean()*2*np.pi
        out[f"J{int(e[1])}"] = Je; out[f"JN{int(e[1])}"] = JN - Je; out[f"Jb{int(e[1])}"] = Jb - Je; out[f"Th{int(e[1])}"] = Theta
        if R_dual is not None:
            out[f"pred{int(e[1])}"] = -(S.h/mu)*R_dual(e)
    return out, S, u, p

import sys, time

def dual_pressure_corr(Nd, mu, kind):
    """compute I_e = int_Gamma (p - pbar) R  for e = e_x, e_y via CNS* dual on mesh Nd"""
    ex = make_exact(kind)
    verts, tris, circ, outer = make_mesh(Nd)
    S = Space(verts, tris, circ, outer)
    zero = lambda X, Y: np.zeros(np.shape(X) + (2,))
    sysm = assemble(S, mu, zero, 1)
    res = {}
    Ed = S.edge_data()
    num = sum((we*ex["p"](Xe[:,0],Xe[:,1])).sum() for (t,ref,Xe,we,n) in Ed); per = sum(we.sum() for (t,ref,Xe,we,n) in Ed)
    pbar = num/per
    for j, e in enumerate(([1.0,0.0],[0.0,1.0])):
        u, p = solve_gen(S, sysm, zero, 1, extra_rhs=dirichlet_rhs(S, mu, e))
        I = 0.0; Rm = 0.0
        for (t, ref, Xe, we, n) in Ed:
            uh, gx, gy, pht = edge_eval(S, u, p, t, ref)
            xs = Xe/np.linalg.norm(Xe,axis=1)[:,None]
            I += (we*(ex["p"](xs[:,0],xs[:,1]) - pbar)*pht).sum()
        res[j] = I
    return res

def xi_pred(N, kind, mu):
    ex = make_exact(kind)
    verts, tris, circ, outer = make_mesh(N)
    S = Space(verts, tris, circ, outer)
    Ed = S.edge_data()
    num = sum((we*ex["p"](Xe[:,0],Xe[:,1])).sum() for (t,ref,Xe,we,n) in Ed); per = sum(we.sum() for (t,ref,Xe,we,n) in Ed)
    pbar = num/per
    Xi = np.zeros(2)
    for (t, ref, Xe, we, n) in Ed:
        A_ = S.verts[S.tris[t]]
        # c_T = int_e r_T  = b^T M^{-1} b with q-> int_e q ; weighted by g: int_e xi_b = int_e g r_T
        xs = Xe/np.linalg.norm(Xe,axis=1)[:,None]
        g = (S.h/mu)*(ex["p"](xs[:,0],xs[:,1]) - pbar)
        q3 = svn._mono(MON3, ref[:,0], ref[:,1])
        Ml = np.einsum('q,qi,qj->ij', QW*abs(S.detJ[t]), Q3_Q, Q3_Q)
        c = np.linalg.solve(Ml, q3.T @ (we*g))
        intxi = we @ (q3 @ c)
        Xi += intxi*n
    return Xi



def dual_data(Nd, mu, kind):
    ex = make_exact(kind)
    verts, tris, circ, outer = make_mesh(Nd)
    S = Space(verts, tris, circ, outer)
    zero = lambda X, Y: np.zeros(np.shape(X) + (2,))
    sysm = assemble(S, mu, zero, 1)
    Ed = S.edge_data()
    num = sum((we*ex["p"](Xe[:,0],Xe[:,1])).sum() for (t,ref,Xe,we,n) in Ed); per = sum(we.sum() for (t,ref,Xe,we,n) in Ed)
    pbar = num/per
    out = {}
    for j, e in enumerate(([1.0,0.0],[0.0,1.0])):
        e = np.array(e)
        u, p = solve_gen(S, sysm, zero, 1, extra_rhs=dirichlet_rhs(S, mu, e))
        I = 0.0; th_all = []; wz_all = []
        for (t, ref, Xe, we, n) in Ed:
            uh, gx, gy, pht = edge_eval(S, u, p, t, ref)
            xs = Xe/np.linalg.norm(Xe,axis=1)[:,None]
            I += (we*(ex["p"](xs[:,0],xs[:,1]) - pbar)*pht).sum()
            dnu = gx*n[0] + gy*n[1]
            tn = dnu - (mu/S.h)*(uh - e[None,:])        # tangential part of Nitsche flux
            th = np.arctan2(xs[:,1], xs[:,0]); tau = np.stack([-np.sin(th), np.cos(th)],1)
            th_all.append(th); wz_all.append((tn*tau).sum(1))
        th_all = np.concatenate(th_all); wz_all = np.concatenate(wz_all)
        # Fourier fit
        K = 12
        Amat = np.concatenate([np.ones((len(th_all),1))] + [np.stack([np.cos(k*th_all), np.sin(k*th_all)],1) for k in range(1,K+1)], 1)
        coef, *_ = np.linalg.lstsq(Amat, wz_all, rcond=None)
        out[j] = (I, coef)
    return out

def WZ(coef, th):
    K = (len(coef)-1)//2
    v = coef[0]*np.ones_like(th)
    for k in range(1, K+1):
        v += coef[2*k-1]*np.cos(k*th) + coef[2*k]*np.sin(k*th)
    return v

def Lambda(N, kind, coef):
    ex = make_exact(kind)
    verts, tris, circ, outer = make_mesh(N)
    S = Space(verts, tris, circ, outer)
    L = 0.0
    for (t, ref, Xe, we, n) in S.edge_data():
        r = np.linalg.norm(Xe,axis=1); xs = Xe/r[:,None]; d = 1 - r
        th = np.arctan2(xs[:,1], xs[:,0]); tau = np.stack([-np.sin(th), np.cos(th)],1); nn = -xs
        G = ex["gu"](xs[:,0], xs[:,1]); dnu = np.einsum('kij,kj->ki', G, nn)
        W = (dnu*tau).sum(1)
        L += (we*d*W*WZ(coef, th)).sum()
    return L



def I_volume(Nd, mu, kind):
    ex = make_exact(kind)
    verts, tris, circ, outer = make_mesh(Nd)
    S = Space(verts, tris, circ, outer)
    zero = lambda X, Y: np.zeros(np.shape(X) + (2,))
    sysm = assemble(S, mu, zero, 1)
    Ed = S.edge_data()
    th = np.linspace(0,2*np.pi,20001)[:-1]
    pbar = ex["p"](np.cos(th), np.sin(th)).mean()      # mean over Gamma
    def chi(r):
        s = np.clip((r-1.2)/0.8, 0, 1)
        return np.where(s<=0,1.0,np.where(s>=1,0.0, np.exp(-1/np.maximum(1-s,1e-300))/(np.exp(-1/np.maximum(1-s,1e-300))+np.exp(-1/np.maximum(s,1e-300)))))
    X = S.xy; r = np.linalg.norm(X,axis=1); xs = X/r[:,None]
    W = (chi(r)*(ex["p"](xs[:,0],xs[:,1]) - pbar))[:,None]*(-xs)
    Wv = np.zeros(2*S.nn); Wv[0::2] = W[:,0]; Wv[1::2] = W[:,1]
    res = {}
    for j, e in enumerate(([1.0,0.0],[0.0,1.0])):
        u, p = solve_gen(S, sysm, zero, 1, extra_rhs=dirichlet_rhs(S, mu, np.array(e)))
        val = Wv @ (sysm["Avol"] @ u) - p @ (sysm["B"] @ Wv)
        # + <chi_h - e, d_n W_I> correction  (flux identity) -> subtract it
        corr = 0.0
        for (t, ref, Xe, we, n) in Ed:
            ph = basis(ref[:,0],ref[:,1]); dX, dY = dbasis(ref[:,0],ref[:,1]); Ji = S.Jinv[t]
            gxe = Ji[0,0]*dX+Ji[1,0]*dY; gye = Ji[0,1]*dX+Ji[1,1]*dY
            Wt = Wv[dofs(S.ids[t])]; Ut = u[dofs(S.ids[t])]
            dnW = (n[0]*gxe+n[1]*gye) @ Wt; uh = ph @ Ut
            corr += (we*((uh-np.array(e)[None,:])*dnW).sum(1)).sum()
        res[j] = (-val, -(val-corr))
    return res


KX, KY = 17.16, -25.36   # int (p-pbar) R, extrapolated (check4)

def Jchi(N, kind, mu, theta):
    ex = make_exact(kind)
    verts, tris, circ, outer = make_mesh(N)
    S = Space(verts, tris, circ, outer)
    sysm = assemble(S, mu, ex["f"], theta)
    u, p = solve_gen(S, sysm, ex["u"], theta)
    zero = lambda X, Y: np.zeros(np.shape(X) + (2,))
    sd = assemble(S, mu, zero, theta)
    out = []
    for j, e in enumerate(([1.0,0.0],[0.0,1.0])):
        chi, pi = solve_gen(S, sd, zero, theta, extra_rhs=dirichlet_rhs(S, mu, np.array(e)))
        w = chi.copy()
        om = w @ (sysm["Avol"] @ u) - p @ (sysm["B"] @ w) - sysm["F"] @ w
        out.append(om)
    return out, S


# torque (psi = x^perp) : J^N error vs -Lambda_h(psi) - Theta_f(psi), test A

def rot_rhs(S, mu):
    """l^hat_psi(w) = -<d_n psi, w> - <psi, d_n w> + (mu/h)<psi, w>,  psi = (-y, x)"""
    R = np.zeros(2*S.nn)
    for (t, ref, Xe, we, n) in S.edge_data():
        ph = basis(ref[:,0],ref[:,1]); dX, dY = dbasis(ref[:,0],ref[:,1]); Ji = S.Jinv[t]
        gxe = Ji[0,0]*dX+Ji[1,0]*dY; gye = Ji[0,1]*dX+Ji[1,1]*dY; dn = n[0]*gxe+n[1]*gye
        psi = np.stack([-Xe[:,1], Xe[:,0]],1); dnpsi = np.array([-n[1], n[0]])
        d = dofs(S.ids[t])
        for c in range(2):
            R[d[:,c]] += we @ (-psi[:,c][:,None]*dn + (mu/S.h)*psi[:,c][:,None]*ph - dnpsi[c]*ph)
    return R

def run_torque(N, kind, mu):
    ex = make_exact(kind)
    verts, tris, circ, outer = make_mesh(N)
    S = Space(verts, tris, circ, outer)
    res = {}
    for th in (0,1):
        sysm = assemble(S, mu, ex["f"], th)
        u, p = solve_gen(S, sysm, ex["u"], th)
        JN = 0.0
        for (t, ref, Xe, we, n) in S.edge_data():
            uh, gx, gy, pht = edge_eval(S, u, p, t, ref)
            dnu = gx*n[0]+gy*n[1]
            psi = np.stack([-Xe[:,1], Xe[:,0]],1); dnpsi = np.array([-n[1], n[0]])
            ph = basis(ref[:,0],ref[:,1])
            # omega_h(v_psi) = <t_h, psi> + <u_h, d_n psi>
            JN += (we*(((dnu-(mu/S.h)*uh)*psi).sum(1) - th*pht*(psi@n) + (uh*dnpsi[None,:]).sum(1))).sum()
        # theta-part: pbar(p_h) * int psi.n = 0 for rigid psi (check numerically skipped)
        th_ = np.linspace(0,2*np.pi,4001)[:-1]; xc=np.cos(th_); yc=np.sin(th_); nn=-np.stack([xc,yc],1)
        Gm = ex["gu"](xc,yc); pv = ex["p"](xc,yc); Du = Gm+np.transpose(Gm,(0,2,1))
        tr = np.einsum('kij,kj->ki',Du,nn) - pv[:,None]*nn
        Jex = ((tr*np.stack([-yc,xc],1)).sum(1)).mean()*2*np.pi
        res[th] = JN - Jex
    # dual Z_psi on same mesh (CNS*) -> Lambda_h
    zero = lambda X, Y: np.zeros(np.shape(X)+(2,))
    sd = assemble(S, mu, zero, 1)
    zeta_rhs = rot_rhs(S, mu)
    chi, pi = solve_gen(S, sd, zero, 1, extra_rhs=zeta_rhs)
    return res, S, chi, pi, Jex

def Lam_from(Sd, chi, mu, N, kind):
    # tangential Nitsche flux of dual minus d_n psi : W_Z - (d_n psi . tau)
    ths=[]; vals=[]
    for (t, ref, Xe, we, n) in Sd.edge_data():
        ph = basis(ref[:,0],ref[:,1]); dX, dY = dbasis(ref[:,0],ref[:,1]); Ji = Sd.Jinv[t]
        gxe = Ji[0,0]*dX+Ji[1,0]*dY; gye = Ji[0,1]*dX+Ji[1,1]*dY
        Ut = chi[dofs(Sd.ids[t])]; uh = ph@Ut; dn = (n[0]*gxe+n[1]*gye)@Ut
        psi = np.stack([-Xe[:,1], Xe[:,0]],1); dnpsi = np.array([-n[1], n[0]])
        tn = dn - (mu/Sd.h)*(uh-psi) - dnpsi[None,:]
        xs = Xe/np.linalg.norm(Xe,axis=1)[:,None]; th=np.arctan2(xs[:,1],xs[:,0]); tau=np.stack([-np.sin(th),np.cos(th)],1)
        ths.append(th); vals.append((tn*tau).sum(1))
    th=np.concatenate(ths); v=np.concatenate(vals); K=12
    A=np.concatenate([np.ones((len(th),1))]+[np.stack([np.cos(k*th),np.sin(k*th)],1) for k in range(1,K+1)],1)
    coef,*_=np.linalg.lstsq(A,v,rcond=None)
    ex = make_exact(kind); verts,tris,circ,outer = make_mesh(N); S = Space(verts,tris,circ,outer)
    L=0.0
    for (t, ref, Xe, we, n) in S.edge_data():
        r=np.linalg.norm(Xe,axis=1); xs=Xe/r[:,None]; d=1-r; th=np.arctan2(xs[:,1],xs[:,0])
        tau=np.stack([-np.sin(th),np.cos(th)],1); G=ex["gu"](xs[:,0],xs[:,1]); dnu=np.einsum('kij,kj->ki',G,-xs)
        W=(dnu*tau).sum(1)
        Az=np.concatenate([np.ones((len(th),1))]+[np.stack([np.cos(k*th),np.sin(k*th)],1) for k in range(1,K+1)],1)
        L+=(we*d*W*(Az@coef)).sum()
    return L


KX, KY = 17.16, -25.36   # int (p - pbar) R_e for p = x^2 y - y + x/2, L = 2.5 (extrapolated from the 'K' mode)

if __name__ == "__main__":
    mode = sys.argv[1]; kind = sys.argv[2]; mu = float(sys.argv[3])
    if mode == "main":
        for N in [int(a) for a in sys.argv[4].split(",")]:
            for th in (0, 1):
                t0 = time.time(); o, *_ = analyse(N, kind, mu, th)
                print({k: (float(f"{v:.5g}") if isinstance(v, float) else v) for k, v in o.items()}, f"{time.time()-t0:.1f}s", flush=True)
    elif mode == "xi":
        for N in (16, 32, 64): print(N, "Xi_h (e_x, e_y):", xi_pred(N, kind, mu), flush=True)
    elif mode == "K":
        for Nd in (32, 64, 96): print(Nd, "K (volume, volume+flux corr):", I_volume(Nd, mu, kind), flush=True)
    elif mode == "lam":
        D = dual_data(int(sys.argv[4]), mu, kind)
        for j in (0, 1):
            for N in (16, 32, 64, 96): print("e", j, "N", N, "Lambda_h", Lambda(N, kind, D[j][1]), flush=True)
    elif mode == "chi":
        Jex = {"P": (1.5708, -2.3562), "B1": (1.5708, 0.15708)}[kind]
        for N in (16, 32, 64):
            for th in (0, 1):
                (o0, o1), S = Jchi(N, kind, mu, th)
                print(N, th, "J^chi err", o0 - Jex[0], o1 - Jex[1], " (h/mu)K:", S.h/mu*KX, S.h/mu*KY, flush=True)
    elif mode == "torque":
        _, Sd, chi, pi, Jex = run_torque(96, kind, mu)
        print("J torque exact", Jex)
        for N in (16, 32, 64):
            res, *_ = run_torque(N, kind, mu)
            print(N, "torque err GS, CNS*:", res[0], res[1], " -Lambda:", -Lam_from(Sd, chi, mu, N, kind), flush=True)
