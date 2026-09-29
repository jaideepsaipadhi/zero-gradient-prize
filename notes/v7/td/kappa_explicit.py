"""Explicit (non-compactness) constants for the Piola reductions of Theorem sh:alfeld (kappa_0, k>=4) and of the
3D leak lower bound (kappa_L, k>=4).  Uses the CERTIFIED reference norms of kappa_cert.py (k=4):
   ||R^||   = 91.678155469  (quadratic edge moments),   ||R^_1|| = 3.645838299  (first moments, leak witness).
Per macro-tetrahedron T with boundary face F (diam F = l), for the reference map A (P1->p1, P2->p2, P3->p3, P0->apex):
   sharpness:  t.M_F(v)/|grad v| >= |w| (det A)^{1/2} / (2 h_T^2 ||A|| ||A^-1||^2 ||R^||)       (Prop. sh:joint)
               kappa_0(T) := c_w(F) (det A)^{1/2} / (2 h_T^2 ||A|| ||A^-1||^2 ||R^|| l^{1/2}),  c_w(F):=sigma_min(H->w(H), Frobenius)/l^2 <= min_{|H|_op=1}|w(H)|/l^2
   leak:       S:M_F(v)/|grad v| >= sigma_min(A_F)^2 (det A)^{1/2} |S|_Fro / (h_T^2 ||A|| ||A^-1|| ||R^_1||)
               kappa_L(T) := sigma_min(A_F)^2 (det A)^{1/2} / (h_T^2 ||A|| ||A^-1|| ||R^_1|| l^{3/2})
All three labelings of F's vertices are tried (the bound holds for each; we report the best).  Float evaluation of
explicit formulas (NUMERICAL evaluation of a PROVED formula); also the crude closed-form worst-case bounds in (c_s, theta_0).
"""
import itertools, numpy as np
R4, R1 = 91.678155469, 3.645838299


def cw(F):
    z = F; E = [z[1] - z[0], z[2] - z[1], z[0] - z[2]]
    n = np.cross(E[0], -E[2]); n /= np.linalg.norm(n)
    u = E[0] / np.linalg.norm(E[0]); v = np.cross(n, u)
    e2 = [np.array([e @ u, e @ v]) for e in E]
    l = max(np.linalg.norm(e) for e in E)
    # rigorous lower bound: |w(H)| >= sigma_min(L)|H|_Fro >= sigma_min(L)|H|_op, L in an orthonormal basis of Sym2
    L = np.array([[e[0] ** 2, e[1] ** 2, np.sqrt(2) * e[0] * e[1]] for e in e2])
    return np.linalg.svd(L, compute_uv=False)[-1] / l ** 2


def factors(p0, F):
    out = []
    for perm in itertools.permutations(range(3)):
        p1, p2, p3 = [F[i] for i in perm]
        A = np.column_stack([p2 - p1, p3 - p1, p0 - p1])
        if np.linalg.det(A) < 0:
            continue
        det = np.linalg.det(A); sv = np.linalg.svd(A, compute_uv=False)
        nA, nAi = sv[0], 1 / sv[-1]
        nrm = np.cross(p2 - p1, p3 - p1); hT = abs((p0 - p1) @ nrm) / np.linalg.norm(nrm)
        AF = np.column_stack([p2 - p1, p3 - p1]); sF = np.linalg.svd(AF, compute_uv=False)[-1]
        l = max(np.linalg.norm(F[i] - F[j]) for i in range(3) for j in range(3))
        k0 = det ** 0.5 / (2 * hT ** 2 * nA * nAi ** 2 * R4 * l ** 0.5)
        kL = sF ** 2 * det ** 0.5 / (hT ** 2 * nA * nAi * R1 * l ** 1.5)
        out.append((k0, kL))
    return max(o[0] for o in out), max(o[1] for o in out)


def inradius_diam(P):
    vol = abs(np.linalg.det(np.column_stack([P[1] - P[0], P[2] - P[0], P[3] - P[0]]))) / 6
    S = sum(np.linalg.norm(np.cross(P[b] - P[a], P[c] - P[a])) / 2 for a, b, c in itertools.combinations(range(4), 3))
    D = max(np.linalg.norm(P[i] - P[j]) for i in range(4) for j in range(4))
    return 3 * vol / S, D


shapes = {
    'regular': (np.array([.5, np.sqrt(3) / 6, np.sqrt(2 / 3)]), [np.array([0., 0, 0]), np.array([1., 0, 0]), np.array([.5, np.sqrt(3) / 2, 0])]),
    'near-regular (pen table)': (np.array([.5, 7 / 24, 13 / 16]), [np.array([0., 0, 0]), np.array([1., 0, 0]), np.array([.5, 7 / 8, 0])]),
    'right corner': (np.array([0., 0, 1]), [np.array([0., 0, 0]), np.array([1., 0, 0]), np.array([0., 1, 0])]),
    'flat h=1/2': (np.array([.5, 7 / 24, .5]), [np.array([0., 0, 0]), np.array([1., 0, 0]), np.array([.5, 7 / 8, 0])]),
    'skew': (np.array([1.5, 1 / 3, .75]), [np.array([0., 0, 0]), np.array([1., 0, 0]), np.array([.5, 7 / 8, 0])]),
}
print(f"{'shape':28s} {'c_w(F)':>8s} {'kappa0(T)':>11s} {'kappaL(T)':>11s} {'c_s=r/D':>8s} {'theta0':>7s} {'crude k0':>10s} {'crude kL':>10s}")
for name, (p0, F) in shapes.items():
    c = cw(F); k0, kL = factors(p0, F)
    r, D = inradius_diam([p0] + F)
    cs = r / D
    ang = []
    for i in range(3):
        a, b = F[(i + 1) % 3] - F[i], F[(i + 2) % 3] - F[i]
        ang.append(np.arccos(a @ b / np.linalg.norm(a) / np.linalg.norm(b)))
    th = min(ang); s = np.sin(th)
    crude0 = s ** 5 * cs ** 2.75 * (16 * np.pi / np.sqrt(3)) ** 0.25 / (np.sqrt(30) * R4)
    # leak crude: sigma_min(A_F)^2 >= l^2 s^4/2 ; det A = 2|F| h_T >= l s * ... ; see REPORT
    crudeL = s ** 5 * cs ** 3.25 * (16 * np.pi / np.sqrt(3)) ** 0.75 / (np.sqrt(6) * R1)
    print(f"{name:28s} {c:8.4f} {c*k0:11.4e} {kL:11.4e} {cs:8.4f} {np.degrees(th):7.2f} {crude0:10.3e} {crudeL:10.3e}")

# ---- explicit Lambda_0 = 2 + C_PT C_I on the Alfeld sub-tetrahedron T_F = conv(B, F) of a REGULAR macro-tet (k=4),
#      C_I^2 <= k(k+2) diam F / h_sub  (Warburton-Hesthaven trace-inverse, degree k-1, d=3),
#      C_PT^2 <= (D_sub/h_sub)(3/pi^2 + 2/pi)  (divergence identity with (x-p) + Payne-Weinberger)
k = 4
p0, F = shapes['regular']
Bc = (p0 + sum(F)) / 4
hsub = Bc[2]; Dsub = max(1.0, max(np.linalg.norm(Bc - f) for f in F)); l = 1.0
CI = np.sqrt(k * (k + 2) * l / hsub); CPT = np.sqrt(Dsub / hsub * (3 / np.pi ** 2 + 2 / np.pi))
Lam0 = 2 + CPT * CI
kL = factors(p0, F)[1]; k0 = cw(F) * factors(p0, F)[0]
cA = np.sqrt(3) / 4                         # |F| >= c_A hG^2 for equilateral faces of diameter hG (c0 = 1)
cL = np.sqrt(4 / np.sqrt(3)) * kL / Lam0     # leak lower-bound constant, c0 = 1
csh = k0 * np.sqrt(cA) / Lam0                # sharpness constant before the factor (1/2)(3/4)^{1/2} of h_1
print(f"\nregular Alfeld macro, k=4: h_sub={hsub:.4f}, C_I<={CI:.3f}, C_PT<={CPT:.3f}, Lambda_0<={Lam0:.3f}")
print(f"  leak lower-bound constant c_L = (4/sqrt3)^(1/2) kappa_L / Lambda_0 = {cL:.4e}   (||grad(X-U)|| >= c_L G_S hG^(1/2) - C hG^(3/2))")
print(f"  sharpness constant kappa_0 c_A^(1/2)/Lambda_0 = {csh:.4e}; with the h_1 factor (1/2)(3/4)^(1/2): {csh*0.5*np.sqrt(0.75):.4e}")
