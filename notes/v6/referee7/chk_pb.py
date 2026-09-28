import sys, os, numpy as np
sys.path.insert(0, os.path.abspath("../robust3"))
import pb
import leak as LK
phi, gphi = LK.PHIS["sin"]
for N in (8, 16):
  for mu in (1.0, 100.0):
    S, sysm, E, ed, Gt, n2t = pb.setup(N, 4, phi, gphi, mu)
    st = pb.stars(S)
    worst_c = 0; worst_pb = 0; dims=set(); worst_sym=0
    Bt = pb.K.coupling(S, sysm, 1)
    for z, tris in sorted(st.items()):
        dof, alld, ns = pb.patch_space(S, sysm, ed, n2t, tris, "sum")
        dims.add(ns.shape[1])
        Nm = sysm["A"][alld][:, dof].toarray() @ ns
        pos = {d:i for i,d in enumerate(alld)}
        for c in range(2):
            w = np.zeros(len(alld)); w[c::2] = 1.0  # dof layout 2g, 2g+1
            idx = np.array([pos[d] for d in alld if d % 2 == c]); w = np.zeros(len(alld)); w[idx] = 1
            worst_c = max(worst_c, np.abs(w @ Nm).max() / np.abs(Nm).max())
        v = np.zeros(2*S.nn); 
        for j in range(ns.shape[1]):
            v[:] = 0; v[dof] = ns[:, j]
            worst_pb = max(worst_pb, np.abs(Bt @ v).max())
    print(N, mu, "dimV", dims, "const-annihilation rel", worst_c, "B* resid", worst_pb, flush=True)
