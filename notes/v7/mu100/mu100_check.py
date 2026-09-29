"""notes/v7/mu100: referee_robust item (4) -- is mu = 100 above the CNS* stability threshold on production meshes (k = 4)?

Reuses cns_crit.Prob (production mesh svn.make_mesh(N), P4/P3-disc, free velocity dofs, Cm = mean-corrected C exactly as
svn.coupling / solve_direct build it).  Every factorisation is fresh (no penalty.py cache is used).
Modes:
  thr N1,N2,..    mu*_GS, mu*_EFF (EFF = A_mu - C^T M^-1 C on Z_h, the effective-penalty form of notes/v7/robust)
                  by inertia bisection; negpiv at mu=100; CNS* critical mu (largest real mu where the saddle is singular)
                  via the memory-light r-form of cns_crit_r.py (sigma=100) plus a far probe sigma=300.
  ldl N           dense LDL^T inertia of the SYMMETRIC saddle [[K, B^T],[B, 0]] for K = GS and EFF at several mu:
                  expected #neg = #pressure dofs exactly when K is PD on Z_h (and B has full row rank).
  cns N           CNS* r-form critical mu only (for large N).
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, scipy.linalg as sla
from cns_crit import Prob
from cns_crit_r import eigmu_r


def real_sorted(ms, tol=1e-4):
    return sorted([float(z.real) for z in ms if abs(z.imag) < tol * abs(z)], reverse=True)


def cns(pb, rec):
    for sig in (100.0, 300.0):
        ms = eigmu_r(pb, "CNS", sig, 1e6)
        rs = real_sorted(ms)
        rec[f"CNS_real_sigma{int(sig)}"] = rs[:4]
        rec[f"CNS_nearest_sigma{int(sig)}"] = [(round(float(z.real), 4), round(float(z.imag), 4)) for z in
                                               sorted(ms, key=lambda z: abs(z - sig))[:4]]
    rec["muc_CNS"] = rec["CNS_real_sigma100"][0] if rec["CNS_real_sigma100"] else None
    rec["dist100_nearest_complex"] = float(min(abs(z - 100.0) for z in eigmu_r(pb, "CNS", 100.0, 1e6)))


def main():
    mode = sys.argv[1]; Ns = [int(s) for s in sys.argv[2].split(",")]
    for N in Ns:
        t0 = time.time(); pb = Prob(N); rho = pb.S.h / pb.S.hGamma
        rec = dict(mode=mode, N=N, nf=pb.nf, npr=pb.npr, rho=rho, gamma100=100 / rho)
        if mode == "thr":
            for w in ("GS", "EFF"):
                mu = pb.threshold(w); rec[f"mustar_{w}"] = mu; rec[f"gammastar_{w}"] = mu / rho
                rec[f"negpiv100_{w}"] = pb.negpiv(100.0, w)
            cns(pb, rec)
            rec["gammac_CNS"] = rec["muc_CNS"] / rho
        elif mode == "cns":
            cns(pb, rec); rec["gammac_CNS"] = rec["muc_CNS"] / rho
        elif mode == "ldl":
            Bd = pb.B.toarray()
            for w in ("GS", "EFF"):
                for mu in (40.0, 50.0, 60.0, 70.0, 80.0, 100.0):
                    K = pb.K(mu, w).toarray(); K = (K + K.T) / 2
                    Mfull = np.block([[K, Bd.T], [Bd, np.zeros((pb.npr, pb.npr))]])
                    ev = np.linalg.eigvalsh(Mfull)
                    tol = 1e-9 * np.abs(ev).max()
                    rec[f"{w}_mu{int(mu)}"] = dict(neg=int((ev < -tol).sum()), zero=int((np.abs(ev) <= tol).sum()),
                                                   pos=int((ev > tol).sum()))
            rec["rankB"] = int(np.linalg.matrix_rank(Bd))
        rec["secs"] = round(time.time() - t0, 1)
        print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
