"""Figure: the GS leak field on test P (pure pressure, u = 0), N = 32, mu = 100.

    python3 plot_leak.py            -> ../paper/figures/leak_field.pdf

Test P has u = 0, so the GS velocity u_h is exactly minus the error; Corollary B' predicts
u_h ~ (h/mu)(p - pbar) n on Gamma_h (n = outward normal of Omega_h, pointing into the polygon).
Left: radial velocity u_h . e_r near the wall (colour) with u_h arrows; black arrows = prediction on Gamma_h.
Right: normal and tangential traces of u_h along Gamma_h against the prediction.
One SuperLU solve (about 1 s, < 1 GB).
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import svn

ROOT = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(ROOT, "..", "paper", "figures"); os.makedirs(FIG, exist_ok=True)
N, MU = 32, 100.0
plt.rcParams.update({"font.size": 8, "savefig.bbox": "tight"})

ex = svn.make_exact("P1")
v, t, c, o = svn.make_mesh(N); S = svn.Space(v, t, c, o)
sysm = svn.assemble(S, MU, ex["f"], 0)
u, p, dres, relres = svn.solve_lu(S, sysm, ex["u"], 0)
E = svn.errors(S, u, ex, MU)
h = S.h

# ---- traces on Gamma_h and the prediction
th, un, ut, pred, pbar_num, per = [], [], [], [], 0.0, 0.0
for (tt, ref, Xe, we, n) in S.edge_data():
    pbar_num += we @ ex["p"](Xe[:, 0], Xe[:, 1]); per += we.sum()
pbar = pbar_num / per
bx, by, bu, bv, px_, py_, pu, pv = [], [], [], [], [], [], [], []
for (tt, ref, Xe, we, n) in S.edge_data():
    ss = np.linspace(0, 1, 9)
    R = np.array([[0, 0], [1, 0], [0, 1]], float)
    uu, ww = [x for x in S.bedges if x[0] == tt][0][1:]
    refp = R[uu][None, :] + ss[:, None] * (R[ww] - R[uu])[None, :]
    X = S.P0[tt] + refp @ S.J[tt].T
    U = svn.basis(refp[:, 0], refp[:, 1]) @ u[svn.dofs(S.ids[tt])]
    tau = np.array([-n[1], n[0]])
    th += list(np.arctan2(X[:, 1], X[:, 0])); un += list(U @ n); ut += list(U @ tau)
    pred += list((h / MU) * (ex["p"](X[:, 0], X[:, 1]) - pbar))
    m = 4                                                    # chord midpoint for the arrows
    bx.append(X[m, 0]); by.append(X[m, 1]); pu.append(pred[-9 + m] * n[0]); pv.append(pred[-9 + m] * n[1])
th, un, ut, pred = map(np.array, (th, un, ut, pred))
order = np.argsort(th)
err_rms = np.sqrt(np.mean((un - pred) ** 2)) / np.sqrt(np.mean(pred ** 2))   # RMS over 9 equispaced points per chord
# relative L2(Gamma_h) misfit of the normal trace, with the edge quadrature of svn.edge_data
num = den = 0.0
for (tt, ref, Xe, we, n) in S.edge_data():
    Ue = svn.basis(ref[:, 0], ref[:, 1]) @ u[svn.dofs(S.ids[tt])]
    pe = (h / MU) * (ex["p"](Xe[:, 0], Xe[:, 1]) - pbar)
    num += we @ (Ue @ n - pe) ** 2; den += we @ pe ** 2
err_rel = np.sqrt(num / den)
print(f"relative L2(Gamma_h) misfit of u_h.n vs (h/mu)(p-pbar): {err_rel:.4f}   (RMS over 9 points/chord: {err_rms:.4f})")

# ---- nodal field near the wall
xy = S.xy; r = np.hypot(xy[:, 0], xy[:, 1]); keep = r < 1.55
U2 = u.reshape(-1, 2)
er = xy / r[:, None]
ur = (U2 * er).sum(1)
tri = mtri.Triangulation(xy[keep, 0], xy[keep, 1])
cen = np.stack([tri.x[tri.triangles].mean(1), tri.y[tri.triangles].mean(1)], 1)
rc = np.hypot(cen[:, 0], cen[:, 1])
edge_len = np.max(np.hypot(np.diff(tri.x[tri.triangles][:, [0, 1, 2, 0]], axis=1),
                           np.diff(tri.y[tri.triangles][:, [0, 1, 2, 0]], axis=1)), 1)
tri.set_mask((rc < np.cos(np.pi / N)) | (edge_len > 0.25))

fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 3.1), gridspec_kw=dict(width_ratios=[1.0, 1.25], wspace=0.55))
vmax = np.abs(ur[keep]).max()
cf = a1.tricontourf(tri, ur[keep], levels=np.linspace(-vmax, vmax, 21), cmap="RdBu_r")
poly = np.vstack([v[:N], v[:1]])
a1.plot(poly[:, 0], poly[:, 1], color="0.2", lw=0.6)
circ = np.linspace(0, 2 * np.pi, 400); a1.plot(np.cos(circ), np.sin(circ), color="0.5", lw=0.4, ls=":")
sel = keep & (r < 1.4)
idx = np.where(sel & (np.arange(len(r)) % 7 == 0))[0]
sc = 0.35 / np.abs(U2[sel]).max()
a1.quiver(xy[idx, 0], xy[idx, 1], U2[idx, 0], U2[idx, 1], color="0.25", angles="xy", scale_units="xy", scale=1 / sc,
          width=0.003, headwidth=3, alpha=0.8)
a1.quiver(bx, by, pu, pv, color="k", angles="xy", scale_units="xy", scale=1 / sc, width=0.008, headwidth=3.5,
          label=r"prediction $(h/\mu)(p-\bar p)n$")
a1.set_aspect("equal"); a1.set_xlim(-1.55, 1.55); a1.set_ylim(-1.55, 1.55); a1.set_xticks([]); a1.set_yticks([])
cb = fig.colorbar(cf, ax=a1, fraction=0.046, pad=0.02, format="%.0e"); cb.ax.set_title(r"$u_h\cdot e_r$", fontsize=7)
cb.ax.tick_params(labelsize=6)
a1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.01), fontsize=6.5, frameon=False)
a1.set_title(r"GS velocity $u_h$, test P", fontsize=8, loc="left")

a2.plot(th[order], pred[order], color="k", lw=1.2, label=r"prediction $(h/\mu)(p-\bar p)$")
a2.plot(th[order], un[order], color="#2a78d6", lw=0, marker="o", ms=1.6, label=r"computed $u_h\cdot n$ on $\Gamma_h$")
a2.plot(th[order], ut[order], color="#c2410c", lw=0.8, label=r"computed $u_h\cdot\tau$ on $\Gamma_h$")
a2.axhline(0, color="0.6", lw=0.5)
a2.set_xlim(-np.pi, np.pi); a2.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
a2.set_xticklabels([r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
a2.set_xlabel(r"polar angle $\theta$"); a2.set_ylabel("boundary velocity")
a2.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
a2.legend(fontsize=6.5, frameon=False, loc="best")
a2.set_title(rf"$N={N}$, $\mu={MU:g}$, $h/\mu={h/MU:.2e}$; rel. $L^2$ misfit {err_rel:.3f}", fontsize=7.5)
out = os.path.join(FIG, "leak_field.pdf")
fig.savefig(out); fig.savefig(os.path.join(ROOT, "..", "notes", "v3", "leak_field_preview.png"), dpi=130)
print(f"wrote {out}; relres={relres:.1e}; slip*mu/(hG)={E['bL2']*MU/(h*np.sqrt(7*np.pi/8)):.4f}; "
      f"normal-trace misfit {err_rel:.4f}; max|u.tau|/max|u.n| = {np.abs(ut).max()/np.abs(un).max():.3f}")
