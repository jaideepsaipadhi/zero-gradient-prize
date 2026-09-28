"""Collects run_*.log (unifmu.py output) into tables and fits d = a hG^{1/2} (1 + b hG + c hG^2).
Usage: python fit.py > fit.log"""
import json, glob, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
geo, fld, gs = {}, {}, {}
for fn in sorted(glob.glob(os.path.join(HERE, "run_*.log"))):
    for line in open(fn):
        if not line.startswith("{"):
            continue
        r = json.loads(line)
        if "field" not in r:
            geo[r["N"]] = r
        elif r["field"] == "GS":
            gs[(r["N"], r["mu"])] = r
        else:
            fld[(r["N"], r["field"])] = r
Ns = sorted(geo)
K, G = 2.7565157, np.sqrt(7 * np.pi / 8)
NU = K * G
print("== geometry / trace quantities (norms in L2(Gamma_h); H^{+-1/2} by Fourier in arclength)")
print(f"{'N':>4} {'hG':>8} {'|zeta|/hG':>10} {'|t*-g|/hG':>10} {'|PiTz|_1/2':>11} {'/hG^.5':>7} {'cert':>8} {'/hG^.5':>7}")
for N in Ns:
    g = geo[N]; s = np.sqrt(g["hG"])
    print(f"{N:4d} {g['hG']:8.5f} {g['zeta']/g['hG']:10.5f} {g['tstar_minus_g']/g['hG']:10.5f} "
          f"{g['PiTzeta_H12']:11.5f} {g['PiTzeta_H12']/s:7.4f} {g['cert_H12']:8.5f} {g['cert_H12']/s:7.4f}")
print("\n== discrete Dirichlet fields: d = |grad(. - U~)| / |grad U|;  ext = |grad X^saw| / |Pi_T zeta|_{H^1/2}")
print(f"{'N':>4} {'d(X^D)':>9} {'/hG^.5':>7} {'d(X^U)':>10} {'ext':>7}")
for N in Ns:
    if (N, "XD") not in fld:
        continue
    g = geo[N]; dD = fld[(N, "XD")]["d"]; dU = fld[(N, "XU")]["d"]
    print(f"{N:4d} {dD:9.5f} {dD/np.sqrt(g['hG']):7.4f} {dU:10.3e} {fld[(N,'XD-XU')]['d']*NU/g['PiTzeta_H12']:7.4f}")
print("\n== GS: d_h, distance to X^D, capture phi = <X-U~, q_h>/<g-U~, q_h>, |X - t*|_{L2(Gamma_h)}")
mus = sorted({m for (_, m) in gs})
for mu in mus:
    print(f"-- mu = {mu:.0e}")
    print(f"{'N':>4} {'gamma':>9} {'gam*hG':>9} {'d_h':>9} {'/hG^.5':>7} {'d(X,X^D)':>10} {'phi':>8} {'|X-t*|':>10}")
    for N in Ns:
        if (N, mu) not in gs:
            continue
        r = gs[(N, mu)]; g = geo[N]
        print(f"{N:4d} {r['gamma']:9.3g} {r['gammahG']:9.3g} {r['d']:9.5f} {r['d']/np.sqrt(g['hG']):7.4f} "
              f"{r['d_XD']:10.3e} {r['phi']:8.5f} {r['X_minus_tstar']:10.3e}")


def fit(h, d, nb):
    A = np.stack([h ** 0.5 * h ** j for j in range(nb)], 1)
    c, *_ = np.linalg.lstsq(A, d, rcond=None)
    res = np.abs(A @ c - d).max() / d.min()
    return c, res


for Nmin in (24, 32):
  print(f"\n== fits d = a hG^{{1/2}} (1 + c1 hG + c2 hG^2)  (N >= {Nmin})")
  for label, get in (("X^D (= GS, mu=inf)", lambda N: fld[(N, "XD")]["d"]),
                   ("GS mu=1e8", lambda N: gs[(N, 1e8)]["d"]),
                   ("GS mu=1e2", lambda N: gs[(N, 1e2)]["d"])):
    try:
        NN = []
        for N in Ns:
            try:
                get(N); NN.append(N)
            except KeyError:
                pass
        NN = [N for N in NN if N >= Nmin]
        h = np.array([geo[N]["hG"] for N in NN]); d = np.array([get(N) for N in NN])
    except KeyError:
        continue
    for nb in (2, 3):
        c, res = fit(h, d, nb)
        print(f"{label:20s} nb={nb}: a={c[0]:.4f} " + " ".join(f"c{j}={c[j]/c[0]:+.3f}" for j in range(1, nb))
              + f"  max rel resid {res:.1e}")
    lr = [np.log(d[i] / d[i + 1]) / np.log(h[i] / h[i + 1]) for i in range(len(h) - 1)]
    print(f"{label:20s} local rates: " + " ".join(f"{x:.3f}" for x in lr))
