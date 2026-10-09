"""L^2_0 variant of the naive (full-traction) Nitsche form: archived reproduction.

Letter, Section 3, "Naive and L^2_0 forms"; extended paper Section 7 (07_consistent_method.tex).
Pressures and continuity tests in Pi_h cap L^2_0 are realized with a scalar multiplier alpha:
    N_h(u,v) - (p, div v) + <p, v.n> = (f, v)     for v in V_R
    (q, div u) + alpha (q, 1)        = 0          for q in Pi_h
    (p, 1)                           = 0
so that div u_h = -alpha is constant and ||div u_h||_{L^2(Omega_h)} = |alpha| |Omega_h|^{1/2}.
Outer data: nodal P_k interpolant g_I of the exact velocity (not flux-corrected), as in the tests.
For comparison the CNS* velocity (mean-free traction) is solved on the same mesh.

Run:  python3 code/l20check.py [N] [mu] [kind]      (defaults 16 100 A)
"""
import sys, json, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import svn


def l20_solve(S, sysm, g):
    A, B, C, F, Ml = sysm["A"], sysm["B"], sysm["C"], sysm["F"], sysm["Ml"]
    Bt = (B - C).tocsr()                      # naive: -(p,div v) + <p, v.n>  =  -p^T (B - C) v
    ndof = A.shape[0]; npr = B.shape[0]
    dD = svn.dofs(S.dnodes).ravel(); free = np.setdiff1d(np.arange(ndof), dD)
    ug = np.zeros(ndof); ug[dD] = g(S.xy[S.dnodes, 0], S.xy[S.dnodes, 1]).ravel()
    m = np.concatenate([Ml[t][:, 0] for t in range(len(Ml))])   # m_j = int q_j (monomial 0 = 1)
    mcol = sp.csr_matrix(m[:, None])
    K = sp.bmat([[A[free][:, free], -Bt[:, free].T, None],
                 [B[:, free], None, mcol],
                 [None, mcol.T, None]], format="csc")
    rhs = np.concatenate([F[free] - A[free][:, dD] @ ug[dD], -B[:, dD] @ ug[dD], [0.0]])
    lu = spl.splu(K, permc_spec="COLAMD", diag_pivot_thresh=1.0)
    x = lu.solve(rhs)
    for _ in range(3):
        x = x + lu.solve(rhs - K @ x)
    relres = np.linalg.norm(K @ x - rhs) / np.linalg.norm(rhs)
    nf = len(free)
    u = ug.copy(); u[free] = x[:nf]; p = x[nf:nf + npr]; alpha = x[-1]
    Bu = B @ u
    divL2 = np.sqrt(abs(Bu @ (sysm["Minv"] @ Bu)))   # ||P_{Pi_h} div u_h|| = ||div u_h|| (div V_h = Pi_h)
    area = float(sum(Ml[t][0, 0] for t in range(len(Ml))))   # |Omega_h| = sum_T int_T 1
    return u, p, alpha, divL2, area, relres


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    mu = float(sys.argv[2]) if len(sys.argv) > 2 else 100.0
    kind = sys.argv[3] if len(sys.argv) > 3 else "A"
    t0 = time.time()
    ex = svn.make_exact(kind)
    verts, tris, circ, outer = svn.make_mesh(N)
    S = svn.Space(verts, tris, circ, outer)
    sysm = svn.assemble(S, mu, ex["f"], 1)
    u, p, alpha, divL2, area, relres = l20_solve(S, sysm, ex["u"])
    E = svn.errors(S, u, ex, mu)
    # CNS* on the same mesh
    uc, pc, rc, relc = svn.solve_lu(S, sysm, ex["u"], 1)
    Ec = svn.errors(S, uc, ex, mu)
    out = dict(N=N, mu=mu, kind=kind, k=4, h=S.h, hGamma=S.hGamma, area_Omega_h=area,
               L20_alpha=alpha, L20_c=-alpha, L20_divL2=divL2,
               L20_abs_c_times_sqrt_area=abs(alpha) * np.sqrt(area),
               L20_relres=relres, L20_H1=E["H1"],
               CNS_divL2=rc[0], CNS_relres=relc, CNS_H1=Ec["H1"], seconds=time.time() - t0)
    print(json.dumps({k: (float(f"{v:.6g}") if isinstance(v, float) else v) for k, v in out.items()}, indent=1))


if __name__ == "__main__":
    main()
