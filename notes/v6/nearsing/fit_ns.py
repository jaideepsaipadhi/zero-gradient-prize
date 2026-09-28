"""Two-term fits E^2 = A hG^3 + B hG^q (q = 2 for eps ~ hG, q = 1 for eps ~ hG^2; Theorem N3 prediction), and the
free-exponent alternative E = a hG^p, on the alt-based sequences.  Prints coefficients and max relative residual."""
import json, numpy as np
from scipy.optimize import curve_fit
for name, q in (("altsplit_h1.0", 2), ("altwall_h1.0", 2), ("altsplit_q1.0", 1), ("altwall_q1.0", 1), ("alt_c0", 1)):
    R = sorted([json.loads(l) for l in open(f"runs/{name}.jsonl")], key=lambda r: r["N"])
    h = np.array([r["hG"] for r in R]); E = np.array([r["H1"] for r in R])
    M = np.stack([h ** 3, h ** q], 1); c, *_ = np.linalg.lstsq(M / E[:, None] ** 2, np.ones_like(E), rcond=None)
    res = np.sqrt(M @ c) / E - 1
    p, _ = curve_fit(lambda x, a, pp: a * x ** pp, h, E, p0=(1, 1))
    res2 = p[0] * h ** p[1] / E - 1
    print(f"{name:14s} A={c[0]:.3f} B={c[1]:.4f} maxrelres={abs(res).max():.1e}  | free: p={p[1]:.3f} maxrelres={abs(res2).max():.1e}"
          f"  | B-term share at finest N: {c[1] * h[-1] ** q / E[-1] ** 2:.2f}")
