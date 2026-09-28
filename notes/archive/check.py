import numpy as np, warnings; warnings.filterwarnings("ignore")
import lamstar as L
# rebuild matrices for uniform m=3 fan and evaluate N(v,v)=A-2B+t*C on the maximizer directly
import types
src=open('lamstar.py').read()
src=src.replace("    lam = np.linalg.eigvalsh(Li @ Seff @ Li.T).max()",
"""    ew,ev_=np.linalg.eigh(Li @ Seff @ Li.T); lam=ew.max(); y=Li.T@ev_[:,-1]
    k=-np.linalg.solve(Mkk,Mky@y) if Mkk.size else 0
    v=Z@(ranC@y+(kerC@k if Mkk.size else 0))
    vA=v@A@v; vB=v@((B+B.T)/2)@v; vC=v@Cm@v; divres=np.abs(D@v).max()
    return dict(lam=lam, quotient=(2*vB-vA)/vC, A=vA,B=vB,C=vC,divres=divres,
                N_at_0_9lam=vA-2*vB+0.9*lam*vC, N_at_1_1lam=vA-2*vB+1.1*lam*vC,
                A_psd=np.linalg.eigvalsh(Ar).min(), C_psd=np.linalg.eigvalsh(Cr).min())""")
ns={}; exec(compile(src,'x','exec'),ns)
print(ns['run'](np.linspace(0,np.pi,4)))
