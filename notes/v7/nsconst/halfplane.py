"""
notes/v7/nsconst: the flat-boundary (half-plane) part of the Stokes regularity constant (REPORT.tex, Lemma 2.3).

Scale-invariant problem: -Delta phi + grad pi = g in {y > 0}, phi = 0 on y = 0, one Fourier mode e^{ix} (by scaling, every
mode k != 0 gives the same value).  Grisvard's identity on a flat boundary gives ||D^2 phi|| = ||Delta phi|| = ||grad pi - g||, so
    q_flat := sup_g ( ||D^2 phi||^2 + ||grad pi||^2 ) / ||g||^2 = sup_g ( ||grad pi - g||^2 + ||grad pi||^2 ) / ||g||^2 ,
the top of the flat-boundary essential spectrum of the quadratic form of C_R (lower-order terms scale out).
1D Taylor--Hood in y (P2 velocity, P1 pressure) on [0, Y] (phi = 0 at y = Y too), g continuous P1 on a coarser grid.

python halfplane.py
"""
import numpy as np, scipy.linalg as sla


def mesh(n, Y, grade=2.0):
    s = np.linspace(0, 1, n + 1)
    return Y * s ** grade


def run(n, Y, ng):
    y = mesh(n, Y); ne = n
    # P2 nodes: vertices + midpoints ; dof index: vertices 0..n, midpoints n+1..2n
    nv = 2 * n + 1
    gp, gw = np.polynomial.legendre.leggauss(6); gp = (gp + 1) / 2; gw = gw / 2
    # reference P2 basis on [0,1]: nodes 0, 1, 1/2
    P2 = lambda t: np.stack([2 * (t - .5) * (t - 1), 2 * t * (t - .5), 4 * t * (1 - t)], -1)
    dP2 = lambda t: np.stack([4 * t - 3, 4 * t - 1, 4 - 8 * t], -1)
    P1 = lambda t: np.stack([1 - t, t], -1); dP1 = lambda t: np.stack([-np.ones_like(t), np.ones_like(t)], -1)
    # g: continuous P1 on a coarse grid yg (ng intervals), evaluated at fine quadrature points
    yg = mesh(ng, Y)
    nU = 2 * nv; nP = n + 1; ngd = 2 * (ng + 1)
    K = np.zeros((nU + nP, nU + nP), complex); Fg = np.zeros((nU + nP, ngd), complex)
    Gm = np.zeros((ngd, ngd), complex)
    # quadrature storage for Q
    rows = []
    for e in range(ne):
        a, b = y[e], y[e + 1]; hlen = b - a
        vd = [e, e + 1, n + 1 + e]; pd_ = [e, e + 1]
        t = gp; w = gw * hlen; yq = a + t * hlen
        ph = P2(t); dph = dP2(t) / hlen; pp = P1(t); dpp = dP1(t) / hlen
        # g basis at yq: hat functions on yg
        G = np.zeros((len(t), ng + 1))
        for j in range(ng + 1):
            G[:, j] = np.interp(yq, yg, np.eye(ng + 1)[j])
        for c in range(2):
            idx = [2 * d + c for d in vd]
            K[np.ix_(idx, idx)] += np.einsum('q,qa,qb->ab', w, dph, dph) + np.einsum('q,qa,qb->ab', w, ph, ph)
            # load (g_c, v_c): g component c basis columns c*(ng+1)+j
            Fg[np.ix_(idx, [c * (ng + 1) + j for j in range(ng + 1)])] += np.einsum('q,qa,qj->aj', w, ph, G)
            Gm[np.ix_([c * (ng + 1) + j for j in range(ng + 1)], [c * (ng + 1) + j for j in range(ng + 1)])] += \
                np.einsum('q,qi,qj->ij', w, G, G)
        # pressure: - int pi conj(div v), div v = i v1 + v2' ; constraint int conj(q) (i phi1 + phi2') = 0
        i1 = [2 * d for d in vd]; i2 = [2 * d + 1 for d in vd]; ip = [nU + d for d in pd_]
        D1 = 1j * np.einsum('q,qp,qa->pa', w, pp, ph)          # int q * (i phi1)
        D2 = np.einsum('q,qp,qa->pa', w, pp, dph)              # int q * phi2'
        K[np.ix_(ip, i1)] += D1; K[np.ix_(ip, i2)] += D2
        K[np.ix_(i1, ip)] += -D1.conj().T; K[np.ix_(i2, ip)] += -D2.conj().T
        rows.append((w, pp, dpp, G, pd_))
    # Dirichlet phi = 0 at y = 0 and y = Y
    fix = [0, 1, 2 * n, 2 * n + 1]
    keep = np.setdiff1d(np.arange(nU + nP), fix)
    Kk = K[np.ix_(keep, keep)]; Fk = Fg[keep]
    X = np.zeros((nU + nP, ngd), complex); X[keep] = np.linalg.solve(Kk, Fk)
    Pi = X[nU:]
    Q = np.zeros((ngd, ngd), complex)
    for (w, pp, dpp, G, pd_) in rows:
        pq = pp @ Pi[pd_]; dpq = dpp @ Pi[pd_]                   # pi and pi' at q-points, per g-basis column
        gx = np.zeros((len(w), ngd), complex); gy = np.zeros((len(w), ngd), complex)
        gx[:, :ng + 1] = G; gy[:, ng + 1:] = G
        e1 = 1j * pq - gx; e2 = dpq - gy                         # grad pi - g
        f1 = 1j * pq; f2 = dpq
        for A_ in (e1, e2, f1, f2):
            Q += np.einsum('q,qi,qj->ij', w, A_.conj(), A_)
    ev = sla.eigh((Q + Q.conj().T) / 2, (Gm + Gm.conj().T) / 2, eigvals_only=True)
    return ev[-1]


if __name__ == "__main__":
    for (n, Y, ng) in ((200, 12, 50), (400, 12, 100), (800, 16, 200)):
        print(dict(n=n, Y=Y, ng=ng, q_flat=round(float(run(n, Y, ng)), 6)), flush=True)
