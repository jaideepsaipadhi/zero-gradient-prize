"""Split the CNS* hydrostatic slip W|Gamma_h into normal/tangential and edge-mean/oscillatory parts (REPORT.tex Sec. 5, item N4).
usage: python slip_split.py N1,N2 mu phi"""
import sys, json
import numpy as np
import leak as Lm
import svn_k as K
Ns = [int(s) for s in sys.argv[1].split(",")]; mu = float(sys.argv[2]); name = sys.argv[3] if len(sys.argv) > 3 else "sin"
for N in Ns:
    phi, gphi = Lm.PHIS[name]
    S = K.build(N, 4); E = Lm.edge_info(S, phi)
    sysm = K.assemble(S, mu, gphi)
    u, p, rel = Lm.solve(S, sysm["A"], K.coupling(S, sysm, 1), sysm["B"], sysm["F"])
    eta = np.sqrt(sum(e["w"] @ e["eta"]**2 for e in E)); eps = S.h / mu
    acc = dict(n=0, t=0, nm=0, tm=0)
    for e in E:
        U = e["ph"] @ u[e["d"]]; n = e["n"]; tau = np.array([-n[1], n[0]])
        un = U @ n; ut = U @ tau; w = e["w"]; L = w.sum()
        acc["n"] += w @ un**2; acc["t"] += w @ ut**2
        acc["nm"] += (w @ un)**2 / L; acc["tm"] += (w @ ut)**2 / L
    r = {k: float(np.sqrt(v) / (eps * eta)) for k, v in acc.items()}
    print(json.dumps(dict(N=N, mu=mu, phi=name, normal=r["n"], tangential=r["t"], normal_mean=r["nm"], tangential_mean=r["tm"])), flush=True)
