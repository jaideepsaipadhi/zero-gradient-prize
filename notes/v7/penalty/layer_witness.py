"""NUMERICAL companion to Theorem 1 of REPORT.tex (penalty necessity on every (M0) mesh).

Witness: v = curl I_A F, I_A = Argyris (C^1-P5) interpolant, F(x,y) = Phi(x/X) chi(y/Xc) Psi(y) with
  Psi(0)=0, Psi'=q, q(0)=1, q' = -s*omega(y/y1),
  omega = 1 on [0,1]; 1-(1-beta)S(t-1) on [1,2]; beta on [2,T1]; beta(1-S((t-T1)/T1)) on [T1,2T1]; 0 after,
  S = C^5 smoothstep (degree 11), s fixed by q(inf)=0.  Psi is EXACTLY quadratic for y <= y1.
Half-plane meshes (boundary {y=0}, boundary edges of length 1):
  cols(H)   : columns 1 x H split by the forward diagonal, stacked (min angle atan(1/H)); the star of every
              boundary vertex is the counterexample star S^(H) of Prop. pb:counter
  dyadic    : row 0 unit squares, row j side 2^j (hanging-midpoint split)  (min angle 26.57 deg)
  colgrad(H): row 0 columns 1 x H, then dyadic rows of side 2^j (j>=1)
Reports Q = h_G (2B - A)/C with h_G = 1 (a lower bound for (mu/h)^* h_G), together with the 1D/continuum value.
"""
import numpy as np, math, sys, time
from numpy.polynomial import polynomial as P

gx, gw = np.polynomial.legendre.leggauss(8)
MONS = [(i, d - i) for d in range(6) for i in range(d + 1)]

# ---------------------------------------------------------------- profile
S_c = np.zeros(12); S_c[6:] = [462, -1980, 3465, -3080, 1386, -252]      # C^5 smoothstep on [0,1]


def S(t, k=0):
    t = np.clip(t, 0, 1)
    return P.polyval(t, P.polyder(S_c, k) if k else S_c)


class Profile:
    """piecewise-polynomial q' on knots; q, Psi integrated exactly (float)"""

    def __init__(self, y1, beta, T1):
        self.y1, self.beta, self.T1 = y1, beta, T1
        # pieces of omega(t), t=y/y1, as polynomials in t on [a,b]
        one = np.array([1.0])
        sh = lambda c, a, sc: P.polyval(np.poly1d([0]), [0]) if False else None
        pieces = [(0, 1, one)]
        # 1-(1-beta) S(t-1): substitute u=t-1
        Su = self._subst(S_c, 1.0, -1.0)                 # S(t-1)
        pieces.append((1, 2, P.polysub(one, (1 - beta) * Su)))
        pieces.append((2, T1, np.array([beta])))
        Sv = self._subst(S_c, 1.0 / T1, -1.0)            # S(t/T1 - 1)
        pieces.append((T1, 2 * T1, beta * P.polysub(one, Sv)))
        self.pieces = pieces
        # int omega
        Iom = sum(P.polyval(b, P.polyint(c)) - P.polyval(a, P.polyint(c)) for a, b, c in pieces)
        self.s = 1.0 / (y1 * Iom)
        # build q and Psi pieces in y: q'(y) = -s omega(y/y1)
        self.ypieces = []
        q0, Psi0 = 1.0, 0.0
        for a, b, c in pieces:
            cy = self._subst(c, 1.0 / y1, 0.0)            # omega(y/y1) in y
            dq = -self.s * cy
            qpoly = P.polyint(dq, lbnd=a * y1, k=q0)       # q(y) with q(a y1)=q0
            Ppoly = P.polyint(qpoly, lbnd=a * y1, k=Psi0)
            self.ypieces.append((a * y1, b * y1, Ppoly))
            q0 = P.polyval(b * y1, qpoly); Psi0 = P.polyval(b * y1, Ppoly)
        self.qinf, self.Psiinf = q0, Psi0
        self.A1d = sum(P.polyval(b, P.polyint(P.polymul(P.polyder(Pp, 2), P.polyder(Pp, 2))))
                       - P.polyval(a, P.polyint(P.polymul(P.polyder(Pp, 2), P.polyder(Pp, 2))))
                       for a, b, Pp in self.ypieces)

    @staticmethod
    def _subst(c, a, b):
        """coeffs of c(a*t+b) in t"""
        out = np.zeros(1); base = np.array([b, a]); pw = np.array([1.0])
        for k, ck in enumerate(c):
            out = P.polyadd(out, ck * pw); pw = P.polymul(pw, base)
        return out

    def psi(self, y, k=0):
        y = np.asarray(y, float); out = np.zeros_like(y)
        # y below 0: quadratic continuation of first piece
        a0, b0, P0 = self.ypieces[0]
        m = y <= b0
        out[m] = P.polyval(y[m], P.polyder(P0, k) if k else P0)
        for a, b, Pp in self.ypieces[1:]:
            m = (y > a) & (y <= b)
            out[m] = P.polyval(y[m], P.polyder(Pp, k) if k else Pp)
        m = y > self.ypieces[-1][1]
        out[m] = self.Psiinf if k == 0 else 0.0
        return out


def phi(t, k=0):
    """(1-t^2)^6 on [-1,1] and its derivatives"""
    c = P.polypow([1, 0, -1], 6)
    t = np.asarray(t, float)
    out = P.polyval(t, P.polyder(c, k) if k else c)
    return np.where(np.abs(t) < 1, out, 0.0)


def chi(t, k=0):
    """1 on t<=1, 1-S(t-1) on [1,2], 0 after"""
    t = np.asarray(t, float)
    if k == 0:
        return np.where(t <= 1, 1.0, np.where(t >= 2, 0.0, 1 - S(t - 1)))
    return np.where((t > 1) & (t < 2), -S(t - 1, k), 0.0)


class Witness:
    def __init__(self, prof, X, Xc):
        self.p, self.X, self.Xc = prof, X, Xc

    def G(self, y, k):
        p, Xc = self.p, self.Xc
        if k == 0:
            return chi(y / Xc) * p.psi(y)
        if k == 1:
            return chi(y / Xc, 1) / Xc * p.psi(y) + chi(y / Xc) * p.psi(y, 1)
        return (chi(y / Xc, 2) / Xc**2 * p.psi(y) + 2 * chi(y / Xc, 1) / Xc * p.psi(y, 1)
                + chi(y / Xc) * p.psi(y, 2))

    def D(self, x, y, a, b):
        X = self.X
        return phi(x / X, a) / X**a * self.G(y, b)

    def continuum(self):
        """exact separable forms on the half plane: A = 2 int Phi'^2 int G'^2 + int Phi^2 int G''^2 + int Phi''^2 int G^2"""
        X = self.X
        t = np.linspace(-1, 1, 20001); w = np.full_like(t, t[1] - t[0]); w[0] = w[-1] = w[0] / 2
        I0 = X * np.sum(w * phi(t)**2); I1 = np.sum(w * phi(t, 1)**2) / X; I2 = np.sum(w * phi(t, 2)**2) / X**3
        ymax = 2 * self.Xc; y = np.linspace(0, ymax, 400001); wy = np.full_like(y, y[1] - y[0]); wy[0] = wy[-1] = wy[0] / 2
        g0, g1, g2 = self.G(y, 0), self.G(y, 1), self.G(y, 2)
        A = 2 * I1 * np.sum(wy * g1**2) + I0 * np.sum(wy * g2**2) + I2 * np.sum(wy * g0**2)
        B = self.p.s * I0; C = I0
        return (2 * B - A) / C


# ---------------------------------------------------------------- Argyris, batched over triangles
def mon_d(c, p, a, b, s):
    """c (n,2), p (n,m,2), s (n,) -> (n,m,21): d^a_x d^b_y of ((x-cx)/s)^i ((y-cy)/s)^j"""
    X = (p[..., 0] - c[:, None, 0]) / s[:, None]; Y = (p[..., 1] - c[:, None, 1]) / s[:, None]
    out = np.zeros(p.shape[:2] + (len(MONS),))
    for k, (i, j) in enumerate(MONS):
        if i < a or j < b:
            continue
        ci = math.factorial(i) / math.factorial(i - a); cj = math.factorial(j) / math.factorial(j - b)
        out[..., k] = ci * cj * X**(i - a) * Y**(j - b) / s[:, None]**(a + b)
    return out


_a = (gx + 1) / 2; _wa = gw / 2
QU = np.repeat(_a, 8); QV = np.tile(_a, 8) * (1 - QU); QW = np.repeat(_wa, 8) * np.tile(_wa, 8) * (1 - QU)


def quotient(T, W, ytop, chunk=4000):
    T = np.array(T)
    keep = ~((T[:, :, 1].min(1) >= 2 * W.Xc) | (T[:, :, 0].min(1) >= W.X) | (T[:, :, 0].max(1) <= -W.X))
    T = T[keep]
    assert T[:, :, 1].max() < ytop, "support reaches the top of the mesh"
    A = B = C = 0.0
    for k0 in range(0, len(T), chunk):
        Pt = T[k0:k0 + chunk]; n = len(Pt)
        c = Pt.mean(1)
        s = np.max(np.stack([np.linalg.norm(Pt[:, i] - Pt[:, (i + 1) % 3], axis=1) for i in range(3)], 1), 1)
        rows, rhs = [], []
        for v in range(3):
            V = Pt[:, v:v + 1]
            for (a, b) in [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]:
                rows.append(mon_d(c, V, a, b, s)[:, 0]); rhs.append(W.D(V[:, 0, 0], V[:, 0, 1], a, b))
        for i in range(3):
            p, q = Pt[:, (i + 1) % 3], Pt[:, (i + 2) % 3]
            mid = ((p + q) / 2)[:, None]; tv = q - p; nn = np.stack([tv[:, 1], -tv[:, 0]], 1)
            nn /= np.linalg.norm(nn, axis=1)[:, None]
            rows.append(nn[:, :1] * mon_d(c, mid, 1, 0, s)[:, 0] + nn[:, 1:] * mon_d(c, mid, 0, 1, s)[:, 0])
            rhs.append(nn[:, 0] * W.D(mid[:, 0, 0], mid[:, 0, 1], 1, 0) + nn[:, 1] * W.D(mid[:, 0, 0], mid[:, 0, 1], 0, 1))
        M = np.stack(rows, 1); r = np.stack(rhs, 1)
        coef = np.linalg.solve(M, r[..., None])[..., 0]                    # (n,21)
        e1, e2 = Pt[:, 1] - Pt[:, 0], Pt[:, 2] - Pt[:, 0]
        Q = Pt[:, None, 0] + QU[None, :, None] * e1[:, None] + QV[None, :, None] * e2[:, None]
        area = np.abs(e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0])
        Fxx = np.einsum('nmk,nk->nm', mon_d(c, Q, 2, 0, s), coef)
        Fxy = np.einsum('nmk,nk->nm', mon_d(c, Q, 1, 1, s), coef)
        Fyy = np.einsum('nmk,nk->nm', mon_d(c, Q, 0, 2, s), coef)
        A += np.sum(area[:, None] * QW[None] * (4 * Fxy**2 + (Fyy - Fxx)**2))
        onb = (Pt[:, :, 1] == 0)
        idx = np.where(onb.sum(1) == 2)[0]
        if len(idx):
            o = np.argsort(~onb[idx], axis=1, kind='stable')[:, :2]
            p = Pt[idx, o[:, 0]]; q = Pt[idx, o[:, 1]]
            u = (gx + 1) / 2
            Qb = p[:, None] + u[None, :, None] * (q - p)[:, None]
            Wb = (gw / 2)[None] * np.linalg.norm(q - p, axis=1)[:, None]
            cc, ss, cf = c[idx], s[idx], coef[idx]
            Fx = np.einsum('nmk,nk->nm', mon_d(cc, Qb, 1, 0, ss), cf)
            Fy = np.einsum('nmk,nk->nm', mon_d(cc, Qb, 0, 1, ss), cf)
            Hyy = np.einsum('nmk,nk->nm', mon_d(cc, Qb, 0, 2, ss), cf)
            Hxy = np.einsum('nmk,nk->nm', mon_d(cc, Qb, 1, 1, ss), cf)
            B += np.sum(Wb * (-Hyy * Fy - Hxy * Fx)); C += np.sum(Wb * (Fy**2 + Fx**2))
    return (2 * B - A) / C, len(T)


# ---------------------------------------------------------------- meshes
def cols(H, W, rows):
    T = []
    for j in range(rows):
        for i in range(-W, W):
            a, b, c, d = (i, j * H), (i + 1, j * H), (i + 1, (j + 1) * H), (i, (j + 1) * H)
            T += [(a, b, c), (a, c, d)]
    return [np.array(t, float) for t in T], rows * H


def graded(H0, J, W):
    """row 0: columns 1 x H0 (H0=1: unit squares); row j>=1: squares of side 2^j (hanging midpoint)"""
    T = []
    for i in range(-W, W):
        a, b, c, d = (i, 0), (i + 1, 0), (i + 1, H0), (i, H0)
        T += [(a, b, c), (a, c, d)]
    y0 = H0
    for j in range(1, J + 1):
        s = 2**j
        for i in range(-W // s, W // s):
            x0 = i * s
            bl, br, tr, tl, mid = (x0, y0), (x0 + s, y0), (x0 + s, y0 + s), (x0, y0 + s), (x0 + s / 2, y0)
            T += [(bl, mid, tl), (mid, br, tr), (mid, tr, tl)]
        y0 += s
    return [np.array(t, float) for t in T], y0


def minangle(T):
    m = 180
    for Pt in T:
        for i in range(3):
            u, v = Pt[(i + 1) % 3] - Pt[i], Pt[(i + 2) % 3] - Pt[i]
            m = min(m, math.degrees(math.acos(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v))))
    return m


if __name__ == "__main__":
    t0 = time.time()
    # (y1, beta, T1, Xc); X is swept.  y1 >= max height of the boundary triangles => Psi exactly quadratic there.
    cases = [
        ("cols H=5 ", lambda W: cols(5, W, 26), [(5, 0.3, 4, 60), (5, 1.0, 1, 60), (2.5, 0.3, 4, 60)]),
        ("cols H=8 ", lambda W: cols(8, W, 17), [(8, 0.3, 4, 64), (8, 1.0, 1, 64)]),
        ("dyadic   ", lambda W: graded(1, 7, W), [(1, 0.3, 4, 50), (2, 0.3, 4, 50), (1, 1.0, 1, 50)]),
        ("colgrad 5", lambda W: graded(5, 7, W), [(5, 0.3, 4, 60), (5, 1.0, 1, 60)]),
    ]
    print("Q = h_G(2B-A)/C of v = curl I_A F (h_G = 1).  beta=1,T1=1 is the one-slope profile (no two-slope spreading).")
    for name, mk, pars in cases:
        T, ytop = mk(64)
        print(f"{name}: min angle {minangle(T):.2f} deg", flush=True)
        for (y1, beta, T1, Xc) in pars:
            pr = Profile(y1, beta, T1)
            print(f"   y1={y1:4.1f} beta={beta:4.2f} T1={T1:3.0f} Xc={Xc:3.0f}: s={pr.s:.5f}, 1D limit 2s-int q'^2 = "
                  f"{2*pr.s-pr.A1d:.5f}", flush=True)
            for X in (200, 800):
                W = Witness(pr, X, Xc)
                Wm = 2**int(math.ceil(math.log2(X + 2)))
                T, ytop = mk(Wm)
                q, n = quotient(T, W, ytop)
                print(f"      X={X:5d}: continuum {W.continuum():9.5f}   DISCRETE Q = {q:9.5f}   [{n} triangles, "
                      f"{time.time()-t0:.0f}s]", flush=True)
