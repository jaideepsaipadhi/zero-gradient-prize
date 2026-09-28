"""Independent float check of alfeld_tet.py: rebuild the optimal witness as piecewise polynomials, then
(i) check continuity across internal faces, zero trace on the 3 apex faces, div=0 at random points;
(ii) recompute A (Duffy/Gauss-Legendre quadrature on each subtet, gradients by central finite differences
     of the polynomial -- independent of the exact monomial tables), B and C on F (Duffy quadrature on the triangle);
(iii) compare the Rayleigh quotient with the exact value."""
import numpy as np, sys, alfeld_tet as a
k = int(sys.argv[1]) if len(sys.argv) > 1 else 3
nm = sys.argv[2] if len(sys.argv) > 2 else "near-regular"
d, nr, A, B, C, dF = a.build(a.TETS[nm], k)
Kl, M, sub, P, N = a.build.last
L, vec = a.lam_float(A, B, C)
coef = np.array([[float(x) for x in row] for row in Kl]) @ vec          # full coefficient vector
Pf = [np.array([float(c) for c in p]) for p in P]
subf = [[np.array([float(c) for c in p]) for p in S] for S in sub]
E = np.array(M)
def v(s, X):  # X: (n,3)
    mon = np.prod(X[:, None, :] ** E[None, :, :], axis=2)            # n x N
    return np.stack([mon @ coef[(s * 3 + c) * N:(s * 3 + c + 1) * N] for c in range(3)], axis=1)
def grad(s, X, h=1e-5):
    G = np.zeros((X.shape[0], 3, 3))
    for j in range(3):
        e = np.zeros(3); e[j] = h
        G[:, :, j] = (v(s, X + e) - v(s, X - e)) / (2 * h)
    return G
rng = np.random.default_rng(0)
def rand_simplex(V, n):
    w = rng.dirichlet(np.ones(len(V)), n); return w @ np.array(V)
err = 0
for i in range(4):
    for j in range(i + 1, 4):
        kl = [t for t in range(4) if t not in (i, j)]
        X = rand_simplex([subf[0][0], Pf[i], Pf[j]], 50)
        err = max(err, np.abs(v(kl[0], X) - v(kl[1], X)).max())
print("max continuity jump on internal faces:", err)
z = 0
for s in (1, 2, 3):
    f = [Pf[j] for j in range(4) if j != s]; z = max(z, np.abs(v(s, rand_simplex(f, 50))).max())
print("max |v| on apex faces:", z)
dv = max(np.abs(np.trace(grad(s, rand_simplex(subf[s], 50)), axis1=1, axis2=2)).max() for s in range(4))
print("max |div v| (FD):", dv)
g, w = np.polynomial.legendre.leggauss(10); g = (g + 1) / 2; w = w / 2
U, V_, W = np.meshgrid(g, g, g, indexing="ij"); WW = (w[:, None, None] * w[None, :, None] * w[None, None, :]).ravel()
U, V_, W = U.ravel(), V_.ravel(), W.ravel()
xi = np.stack([U, V_ * (1 - U), W * (1 - U) * (1 - V_)], 1); J = (1 - U) ** 2 * (1 - V_)
Aq = 0
for s in range(4):
    o = subf[s][0]; Mt = np.stack([subf[s][t] - o for t in (1, 2, 3)], 1)
    X = o + xi @ Mt.T; G = grad(s, X); D = G + np.transpose(G, (0, 2, 1))
    Aq += abs(np.linalg.det(Mt)) * np.sum(WW * J * 0.5 * np.sum(D * D, axis=(1, 2)))
U2, V2 = np.meshgrid(g, g, indexing="ij"); W2 = (w[:, None] * w[None, :]).ravel(); U2, V2 = U2.ravel(), V2.ravel()
xi2 = np.stack([U2, V2 * (1 - U2)], 1); J2 = (1 - U2)
o = Pf[1]; Mt = np.stack([Pf[2] - o, Pf[3] - o], 1); X = o + xi2 @ Mt.T
area2 = abs(np.linalg.det(Mt[:2, :]))
vv = v(0, X); dnv = -grad(0, X)[:, :, 2]
Bq = area2 * np.sum(W2 * J2 * np.sum(dnv * vv, 1)); Cq = area2 * np.sum(W2 * J2 * np.sum(vv * vv, 1))
print(f"quadrature: A={Aq:.8g} B={Bq:.8g} C={Cq:.8g}  (2B-A)/C = {(2*Bq-Aq)/Cq:.8f}   exact-pipeline lam* = {L:.8f}")
