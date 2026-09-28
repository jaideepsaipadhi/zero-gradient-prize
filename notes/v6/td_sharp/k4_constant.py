"""Reference constant for the k=4 Alfeld witness (NUMERICAL value; existence is exact).
kappa_hat = min singular value of the joint edge-moment map E: Y^1 -> R^6, with Y^1 carrying
the H^1-seminorm ||grad v||_{L2(T_hat)} (right inverse norm ||R|| = 1/kappa_hat)."""
import pickle, numpy as np, sympy as sp
import alfeld_witness as aw
k = 4
subs, polys, nd, nunk, A = aw.build(k)
d = pickle.load(open('alfeld_k4.pkl', 'rb'))
L = np.array(d['L'].tolist(), dtype=float); N = np.array(d['N'].tolist(), dtype=float)
# collapsed-coordinate Gauss quadrature on a tetrahedron (exact to high degree)
g, w = np.polynomial.legendre.leggauss(6); g = (g + 1) / 2; w = w / 2
pts = []; wts = []
for i, a in enumerate(g):
    for j, b in enumerate(g):
        for l, c in enumerate(g):
            u = a; v = b * (1 - a); t = c * (1 - a) * (1 - b)
            pts.append((u, v, t)); wts.append(w[i] * w[j] * w[l] * (1 - a) ** 2 * (1 - b))
pts = np.array(pts); wts = np.array(wts)
x, y, z = aw.X
G = np.zeros((nunk, nunk))
for K, bl in zip(subs, polys):
    V = np.array([[float(c) for c in v] for v in K])
    J = np.array([V[1] - V[0], V[2] - V[0], V[3] - V[0]]).T
    X = V[0][None, :] + pts @ J.T; vol = abs(np.linalg.det(J))
    grads = []
    for gi, p in bl:
        f = [sp.lambdify((x, y, z), p.diff(X_).as_expr(), 'numpy') for X_ in aw.X]
        gr = np.stack([np.broadcast_to(fi(X[:, 0], X[:, 1], X[:, 2]), (len(X),)) for fi in f], 1)
        grads.append((gi, gr))
    for gi, gr in grads:
        for gj, gr2 in grads:
            val = vol * np.sum(wts * np.sum(gr * gr2, 1))
            for c in range(3):
                G[3 * gi + c, 3 * gj + c] += val
GY = N @ G @ N.T                      # Gram of |grad v|^2 on Y basis
Cm = L[0:2]; E = L[2:]
_, sv, Vt = np.linalg.svd(Cm); K1 = Vt[2:].T   # kernel of mean map in Y coords
GY1 = K1.T @ GY @ K1
ev, U = np.linalg.eigh(GY1); Winv = U @ np.diag(ev ** -0.5) @ U.T
S = np.linalg.svd(E @ K1 @ Winv, compute_uv=False)
print('dim Y1 =', K1.shape[1], ' singular values of E (H1-orthonormal):', S)
print('kappa_hat = sigma_min =', S.min(), '   ||R|| = 1/kappa_hat =', 1 / S.min())
