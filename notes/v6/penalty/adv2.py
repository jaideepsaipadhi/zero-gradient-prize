"""Adversarial search variants. mode 'eq': equal boundary edges (r_m = r_0 = 1)."""
import sys, numpy as np
from scipy.optimize import minimize
from patch import fan_polar
m = int(sys.argv[1]); thmin = float(sys.argv[2])*np.pi/180; ns=int(sys.argv[3]); seed=int(sys.argv[4]); mode=sys.argv[5]
rng=np.random.default_rng(seed)
def decode(x):
    w=np.exp(x[:m]); al=np.pi*w/w.sum(); ang=np.concatenate([[0],np.cumsum(al)]); ang[-1]=np.pi
    r=np.concatenate([[1.0],np.exp(x[m:2*m-1]),[1.0]]) if mode=='eq' else np.concatenate([[1.0],np.exp(x[m:2*m])])
    return ang,r,al
def minangle(ang,r):
    out=[]
    for i in range(m):
        O=np.zeros(2);A=r[i]*np.array([np.cos(ang[i]),np.sin(ang[i])]);B=r[i+1]*np.array([np.cos(ang[i+1]),np.sin(ang[i+1])])
        for p,q,s in ((O,A,B),(A,B,O),(B,O,A)):
            u,v=q-p,s-p; out.append(np.arccos(np.clip(u@v/np.linalg.norm(u)/np.linalg.norm(v),-1,1)))
    return min(out)
def obj(x):
    ang,r,al=decode(x); ma=minangle(ang,r)
    if ma<thmin: return 100+1e3*(thmin-ma)
    return fan_polar(ang,r)[0]
best=(np.inf,None)
nx = 2*m-1 if mode=='eq' else 2*m
for s in range(ns):
    x0=rng.normal(0,.5,nx)
    if obj(x0)>90: continue
    res=minimize(obj,x0,method='Nelder-Mead',options=dict(maxiter=4000,xatol=1e-7,fatol=1e-9))
    ang,r,al=decode(res.x)
    print(f"{s}: lam={res.fun:.5f} al={np.round(al*180/np.pi,2)} r={np.round(r,3)}",flush=True)
    if res.fun<best[0]: best=(res.fun,res.x)
ang,r,al=decode(best[1]); print("BEST",best[0],list(np.round(al*180/np.pi,4)),list(np.round(r,5)),fan_polar(ang,r)[1])
