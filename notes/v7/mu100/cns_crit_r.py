"""Memory-light CNS* critical mu (the saddle LU of cns_crit.py OOMs at N=128 in 7 GB).
Perturbed saddle [A_mu, -Bt^T; B, -M/r] eliminates p = r M^{-1} B u exactly, giving the velocity-only
nonsymmetric matrix  K_r(mu) = A0 + mu P + r Bt^T M^{-1} B,  Bt = B - Cm  (CNS*)  or Bt = B (GS).
Its singular mu converge to those of the saddle as r -> infinity (O(1/r)); validated against cns_crit.py at N <= 64.
usage: python cns_crit_r.py N1,N2 [r=1e6] [which=CNS]
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, scipy.sparse.linalg as spl
from cns_crit import Prob


def eigmu_r(pb, which, sigma=100.0, r=1e6, nev=8):
    Bt = pb.B if which == "GS" else (pb.B - pb.Cm).tocsr()
    K = (pb.A0 + sigma * pb.P + r * (Bt.T @ pb.Minv @ pb.B)).tocsc()
    lu = spl.splu(K, permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.1)
    op = spl.LinearOperator(K.shape, matvec=lambda x: lu.solve(pb.P @ x), dtype=float)
    nu = spl.eigs(op, k=nev, which="LM", return_eigenvectors=False, tol=1e-10, maxiter=5000)
    del lu
    return sigma - 1 / nu


if __name__ == "__main__":
    Ns = [int(s) for s in sys.argv[1].split(",")]
    r = float(sys.argv[2]) if len(sys.argv) > 2 else 1e6
    which = sys.argv[3] if len(sys.argv) > 3 else "CNS"
    for N in Ns:
        t0 = time.time(); pb = Prob(N)
        ms = eigmu_r(pb, which, 100.0, r)
        real = sorted([float(z.real) for z in ms if abs(z.imag) < 1e-4 * abs(z)], reverse=True)
        print(json.dumps(dict(N=N, which=which, r=r, muc=real[0] if real else None, real=real,
                              all=[(float(z.real), float(z.imag)) for z in ms], secs=round(time.time() - t0, 1))), flush=True)
