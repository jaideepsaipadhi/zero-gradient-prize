"""NUMERICAL sanity check of the constants: sample the compact families F3(td,dd), Fs(td,dd) of
notes/v4/m3_proof.tex (Def. m3:fam), evaluate the witness constant kappa_w by quadrature, and locally
minimise.  Compares with the certified lower bounds (v4 log) and with kappa_el.  Not a proof.
Run: python3 sample_families.py 30 10 4000     (a few minutes)"""
import sys, math, numpy as np
from scipy.optimize import minimize
from witness_lib import kappa as kappa_parts
td, dd, n = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
T, Dl = math.radians(td), math.radians(dd)
rng = np.random.default_rng(1)

def fan(z):
    al, be = z[:3], z[3:]
    ga = math.pi - al - be
    r = [1.0]; th = [0.0]
    for j in range(3):
        r.append(r[-1] * math.sin(be[j]) / math.sin(ga[j])); th.append(th[-1] + al[j])
    return np.array([[r[i] * math.cos(th[i]), r[i] * math.sin(th[i])] for i in range(4)]), ga

def feasible(z, fam):
    al, be = z[:3], z[3:]; ga = math.pi - al - be
    if min(al.min(), be.min(), ga.min()) < T - 1e-12: return False
    S = al.sum()
    if fam == 'F3': return abs(S - math.pi) <= Dl + 1e-12
    return S <= math.pi + Dl - T + 1e-12

def kap(z, fam):
    p, _ = fan(z)
    cr = lambda u, v: u[0] * v[1] - u[1] * v[0]
    A = cr(p[1], p[3]) * cr(p[2], p[3]); A3 = cr(p[0], p[1]) * cr(p[0], p[2])
    _, nr = kappa_parts(p)                       # ||grad w|| for |p0| = 1
    r3 = np.linalg.norm(p[3])
    if fam == 'F3':
        return (A + r3**5 * A3) / 60 / nr / max(1.0, r3)**2
    return A / 60 / nr

def sample(fam):
    while True:
        al = rng.uniform(T, math.pi - 2 * T, 3); be = rng.uniform(T, math.pi - 2 * T, 3)
        if fam == 'F3': al[2] = math.pi + rng.uniform(-Dl, Dl) - al[0] - al[1]
        z = np.concatenate([al, be])
        if feasible(z, fam): return z

for fam in ('F3', 'Fs'):
    Z = [sample(fam) for _ in range(n)]
    K = np.array([kap(z, fam) for z in Z])
    order = np.argsort(K)
    best = K[order[0]]; bz = Z[order[0]]
    for i in order[:8]:
        f = lambda z: kap(z, fam) if feasible(z, fam) else 1e3
        r = minimize(f, Z[i], method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-10, maxiter=4000))
        if r.fun < best: best, bz = r.fun, r.x
    print(f"{fam}({td:g},{dd:g}): {n} samples, sample min {K.min():.6f}, all positive: {bool((K>0).all())}, "
          f"local-min estimate {best:.6f} at angles(deg) {np.round(np.degrees(bz),2).tolist()}")
