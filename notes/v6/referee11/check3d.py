"""referee11: 3D spot checks (no writes to num3d/): Gamma_h normal orientation; P-test (u=0, f=grad p) GS slip / leak
ratios recomputed with an independent face quadrature (Dunavant-free: collapsed Gauss on each face, own normals)."""
import sys, numpy as np
sys.path.insert(0, '../../../code3d')
import svn3d, exact3d
k, N, mu = 3, 2, 400.0
S = svn3d.build("sph", N, k, L=2.0)
# orientation: fluid outward normal on Gamma_h must point INTO P_h (towards origin); on box faces away from origin
for name, faces, sgn in (("Gamma", S.fGamma, -1), ("box", S.fBox, +1)):
    ok = True
    for fd in S.face_data(faces):
        ctr = fd["X"].mean(1); ok &= np.all(sgn * np.einsum('ti,ti->t', fd["n"], ctr) > 0)
    print(f"{name}: normals oriented as claimed: {ok}")
ex = exact3d.make("P")
sysm = svn3d.assemble(S, mu, ex["f"], [(S.fGamma, None, "grad")], [])
u, p, hist, rr = svn3d.solve_ipm(S, sysm, ex["u"])
E = svn3d.errors(S, u, ex, mu, nq=k + 3)
# independent boundary quadrature: own normals from vertex coordinates and own reference points
g, w = np.polynomial.legendre.leggauss(10); g = (g + 1) / 2; w = w / 2
U_, V_ = np.meshgrid(g, g, indexing='ij'); W_ = np.outer(w, w) * (1 - U_)
bt, bi = S.fGamma; FL = [(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)]
REF = svn3d.Ref.REFV
num = den = un2 = pint = area = 0.0; vals = []
for t, i in zip(bt, bi):
    Vt = S.V[S.T[t]]; A, B, C = (Vt[j] for j in FL[i])
    nv = np.cross(B - A, C - A); ar2 = np.linalg.norm(nv); n = nv / ar2
    if np.dot(n, A) > 0: n = -n          # into P_h (A on unit sphere, P_h convex containing origin)
    s_ = U_.ravel(); r_ = (V_ * (1 - U_)).ravel(); wt = W_.ravel() * ar2
    X = A + s_[:, None] * (B - A) + r_[:, None] * (C - A)
    ra, rb, rc = (REF[j] for j in FL[i])
    ref = ra + s_[:, None] * (rb - ra) + r_[:, None] * (rc - ra)
    uh = S.R.basis(ref) @ u[S.dofs[t]]
    pe = ex["p"](X[:, 0], X[:, 1], X[:, 2])
    vals.append((wt, uh, n, pe)); pint += wt @ pe; area += wt.sum()
pbar = pint / area
G2 = sum(wt @ (pe - pbar) ** 2 for wt, uh, n, pe in vals)
b2 = sum(wt @ (uh ** 2).sum(1) for wt, uh, n, pe in vals)
leak = sum(wt @ ((uh @ n) * (pe - pbar)) for wt, uh, n, pe in vals)
lam = S.h / mu; G = np.sqrt(G2)
print(f"P test k=3 N=2 mu=400: own slip={np.sqrt(b2)/(lam*G):.4f} own leak={leak/(lam*G2):.4f} | solver slip={E['slip_ratio']:.4f} leak={E['leak_corr']:.4f} | h={S.h:.4f} |div|={hist[-1]:.1e}")
