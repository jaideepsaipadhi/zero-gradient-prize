"""Exact verification of the explicit (P^T) fields (pt_fields.py) against the independent constraint assembly fan.py:
continuity, boundary conditions, pressure-blindness, (M), (D) (full vectors), for random rational shapes."""
import random, sys
from fractions import Fraction as Fr
import fan, pt_fields as PF

def to_vec(Fd, Fn):
    v = [Fr(0)] * Fn.nu
    for t in (1, 2):
        for comp in range(2):
            for al in PF.I4:
                v[Fn.u(t, comp, al)] = Fr(Fd[t][comp][al])
    return v

def check(x1, y1, x2, y2):
    N = PF.NumExact
    G, V, info = PF.star_fields(x1, y1, x2, y2, N)
    assert info["mC"] == 0 and info["mA"] != 0 and info["dC"] != 0
    z = (0, 0); a = (1, 0)
    Fn = fan.Fan([z, a, G["c"], G["b"]], k=4, last_is_gamma=False)
    rows = Fn.constraints()
    worst = 0; traces = []
    for F in V:
        v = to_vec(F, Fn)
        for r in rows:
            val = sum(cf * v[key] for key, cf in r.items())
            worst = max(worst, abs(val))
        traces.append([v[Fn.u(1, 0, (4 - i, i, 0))] * Fn.n_e[0] + v[Fn.u(1, 1, (4 - i, i, 0))] * Fn.n_e[1] for i in range(5)])
    return worst, traces, info

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
for it in range(int(sys.argv[2]) if len(sys.argv) > 2 else 20):
    p = [Fr(random.randint(-5, 20), 10) for _ in range(4)]
    if p[0] + p[1] <= 0 or p[2] + p[3] <= 0:
        continue
    w, tr, info = check(*p)
    print([str(q) for q in p], "max |constraint| =", w, " traces:", [[str(q) for q in t] for t in tr],
          " mA,dC:", str(info["mA"]), str(info["dC"]))
