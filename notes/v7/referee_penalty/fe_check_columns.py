"""Referee check of Prop. cols with the production FE assembly (code/svn.py):
v = curl(Phi(x) G(y)) sampled at the P4 Lagrange nodes (exact: curl F is P4 on every triangle),
Q = h_G (2B - A)/C computed from the assembled Nitsche matrices, and ||div v|| from the divergence matrix."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "code"))
import svn
from fractions import Fraction

zero_f = lambda X, Y: np.zeros(np.shape(X) + (2,))


def bs(t, k):  # uniform cubic B-spline in t (support [-2,2]) and derivatives
    t = np.asarray(t, float); a = np.abs(t); s = np.sign(t)
    out = np.zeros_like(t)
    m1 = a <= 1; m2 = (a > 1) & (a < 2)
    if k == 0:
        out[m1] = (4 - 6 * a[m1]**2 + 3 * a[m1]**3) / 6; out[m2] = (2 - a[m2])**3 / 6
    elif k == 1:
        out[m1] = s[m1] * (-12 * a[m1] + 9 * a[m1]**2) / 6; out[m2] = -s[m2] * 3 * (2 - a[m2])**2 / 6
    return out


def G(y, H, k):
    c = 1 / (8 * H); y = np.asarray(y, float); out = np.zeros_like(y)
    m1 = y <= H; m2 = (y > H) & (y <= 3 * H); m3 = (y > 3 * H) & (y < 5 * H)
    if k == 0:
        out[m1] = y[m1] - y[m1]**2 / (2 * H); out[m2] = H / 2 - c * (y[m2] - H)**2 / 2
        out[m3] = c * (y[m3] - 5 * H)**2 / 2
    else:
        out[m1] = 1 - y[m1] / H; out[m2] = -c * (y[m2] - H); out[m3] = c * (y[m3] - 5 * H)
    return out


def run(H, d, pad=2):
    nx = 4 * d + 2 * pad; ny = 5 + pad
    xs = np.arange(nx + 1) - nx / 2; ys = np.arange(ny + 1) * H
    V = np.array([(x, y) for y in ys for x in xs]); idx = lambda i, j: j * (nx + 1) + i
    T = []
    for j in range(ny):
        for i in range(nx):
            a, b, c, e = idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)
            T += [(a, b, c), (a, c, e)]          # forward diagonal
    T = np.array(T); circ = set(range(nx + 1)); outer = set(idx(i, ny) for i in range(nx + 1))
    S = svn.Space(V, T, circ, outer)
    s0 = svn.assemble(S, 0.0, zero_f, None); s1 = svn.assemble(S, 1.0, zero_f, None)
    X, Y = S.xy[:, 0], S.xy[:, 1]
    vx = bs(X / d, 0) * G(Y, H, 1); vy = -(bs(X / d, 1) / d) * G(Y, H, 0)
    u = np.stack([vx, vy], 1).ravel()
    AmB = u @ (s0["A"] @ u); Cm = u @ ((s1["A"] - s0["A"]) @ u) * S.h
    Avol = u @ (s0["Avol"] @ u)
    div = s0["B"] @ u; divn = np.sqrt(div @ (s0["Minv"] @ div))
    Q = -AmB / Cm * S.hGamma
    Qex = (2265 * (d / H)**4 - 2800 * (d / H)**2 - 6944) / (2416 * (d / H)**4) / H
    print(f"H={H} d={d}: hG={S.hGamma} FE Q={Q:.12f}  exact={Qex:.12f}  diff={Q-Qex:.2e}  ||div||={divn:.1e}"
          f"  A={Avol:.6f} 2B={Avol-AmB:.6f} C={Cm:.6f}", flush=True)


for H, d in ((5, 20), (1, 4), (8, 32), (2.5, 10)):
    run(H, d)
