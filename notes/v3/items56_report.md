# Items 5–6 report: the H¹ leak constant and broader penalty necessity

New files only. Nothing in `paper/` or in the existing `code/` files was edited.

- `notes/v3/h1_constant.tex` (labels `hc:`)
- `notes/v3/penalty_broad.tex` (labels `pb:`)
- `code/h1_limit.py`
- `code/lamstar_family.py`
- `logs/h1_limit_series.log`, `logs/h1_limit_fem.log`
- `logs/lamstar_family_{bubble,certify,mesh,limit,limitstar}.log`

Both notes are fragments in the macros of `paper/main.tex`. Compiled together with sections 02–09, the build has no errors and no undefined labels (the only warnings are the expected missing citations, since no .bib was run).

## Item 5 — the H¹ constant (`h1_constant.tex`)

### The limit
The limit of (μ/h)(ũ − u_h) is the **leak flow** U: the Stokes solution in Ω with Dirichlet data U = (p − p̄)e_r = −(p − p̄)n on Γ and U = 0 on ∂R.

- The data are compatible, since p̄ is the mean.
- The trace of U equals the trace of the test field w, but U is Stokes and w is not.
- The penalty term of N_h cancels to leading order against Rc. What remains is the chord/arc mismatch ζ, which is odd in s, plus O(h/μ) terms.
- The normal-load cancellation lemma (Lemma 6.4) controls the ⟨Ũ, ∂_n v⟩ term.

### Theorem hc:thm (proved)
With ε = δ_R − (h/μ)z_U,

‖∇((μ/h)δ_R − Ũ)‖ ≤ C [ min(γh_Γ^{1/2}, (γh_Γ)^{1/2}) + h_Γ^{1/2} + (h/μ)^{1/2} + E_U ],  E_U = O(h).

- **Corollary hc:cor.** Under the hypotheses of Cor B′, plus γ(1+γ^{1/2})h_Γ^{1/2} → 0 and (μ/h)ε_h → 0, the H¹ constant converges to K = ‖∇U‖_{L²(Ω)}/G.
- **Test P.** ũ − u_h = δ_R exactly, so no extra hypothesis is needed.
- **Theorem C's remainder.** (μ/h)r_h → U − w in H¹. This remainder is the Stokes correction of the cut-off field.

### Numbers for test P
- **Rigorous bracket** (domain monotonicity of the Dirichlet minimum principle, exact annulus mode energies): **2.0746004 ≤ K ≤ 3.0525769**.
- **Series solution** (Michell expansion; r = 1 conditions exact, least squares on ∂R; stable to 1e-9 in both truncation and quadrature): **K = 2.7565157**, ‖∇U‖ = 4.5702445.
- **FEM check** (paper's solver, SuperLU, N ≤ 96, μ = 1e2 to 1e8, compared pointwise with the series U):
  - It reproduces Table 3 exactly: 2.464, 2.612, 2.684, 2.708, 2.770, 2.790.
  - At μ = 100, K_h rises to K and the direct distance d_h = ‖∇((μ/h)e − U)‖/‖∇U‖ decays like h^0.6 (0.193 → 0.067).
  - At large μ (γh_Γ ≫ 1, outside the theorem), K_h falls to K from above, and d_h decays like h_Γ^0.8 (0.72 → 0.175).
  - The excess at large μ is an almost orthogonal element-scale oscillation: K_h² ≈ K²(1 + d_h²).
  - So the measured 2.71–2.79 bracket the limit 2.7565. The μ-trend at fixed N is a crossover between the two regimes.
- **Dependence on the box.** K is not universal:
  - E_3(R) → 45π (exact).
  - E_1(R) → π, split into a dipole plus a translation. The translation's energy vanishes like 1/log R (Stokes paradox).
  - Hence K(L) → 3/√7 ≈ 1.134 as L → ∞, but only logarithmically: annulus values are 1.42, 1.24 and 1.18 at R = 10, 1e2, 1e4.

### What remains open for item 5
Convergence uniformly in μ ≥ μ0, i.e. in the Dirichlet-dominated regime γh_Γ → ∞, is observed but not proved (Remark hc:regime).

## Item 6 — broader penalty necessity (`penalty_broad.tex`)

### Theorem pb:thm (proved)
Take **any** boundary triangle T with edge e on Γ_h. The field v = curl(λ₁²λ₂²(c₀λ₀ + c₁λ₁ + c₂λ₂)) is a P4 field in Z_h for every k ≥ 4, whatever the rest of the mesh is. It needs no nonsingularity, straightness or perturbation hypothesis.

- The forms A_T, B_T, C_T are explicit 3×3 rational matrices in the apex (a, b), printed in the note.
- N_h is indefinite on Z_h if μ|e|/h < λ_T.
- The forms were cross-checked by independent quadrature (4.5e-9) and against the paper's penalty code on single-triangle patches (10.48912 vs 10.48912).

### Proposition pb:cert (exact certificates)
Proved by exact rational interval arithmetic on Φ_q = 2520b³ qᵀ(2B − A − cC)q over dyadic boxes, with rational q per box:

- λ_T > 9.3 on all triangles with apex angle ≥ 60° and base angles ≥ 5°. The true minimum is 9.31381, at the equilateral triangle.
- λ_T > 4 on apex ≥ 45°.
- λ_T > 1 on apex ≥ 35°.

The negative controls (c = 9.4 and 4.5) fail as they should.

### Corollary pb:rho
If the triangle on a shortest boundary edge has apex angle ≥ 35°, then μ > ρ is necessary (and μ > 9.3ρ if the apex angle is ≥ 60°). So μ0 ≍ ρ is necessary and sufficient on that whole class.

**Limitation.** λ_T < 0 for tall triangles (apex ≤ ~30°). A uniform c(Θ0) > 0 for every nonsingular star is still open.

### The gap between 9 and γ* ≈ 18–21
This was computed with the paper's own threshold code on the Section 9 meshes, then in the N → ∞ limit geometry.

| N | max triangle | max star | 4-edge window | full first layer | 2 layers ≈ global |
|---|---|---|---|---|---|
| 64 | 10.49 | 17.70 | 20.29 | 20.69 | 20.80 |

In the limit geometry at θ = 0, the layer has height 3/4 per chord:

- The star S_{3/4} has λ* = 17.7064, certified exactly ≥ 17.
- The W = 4 window has λ* = 20.4872, certified exactly ≥ 20.
- The windows converge: W = 32 gives 21.376 (1 layer) and 21.397 (2 layers).
- **Predicted γ*_∞ ≈ 21.4.** The paper says "a limit near 22"; it should say ≈ 21.4.

The gap has four sources, in order of size:

1. **Flatter stars.** The actual stars have height 3/4, not 1, at the longest chords, which is where γ* is decided.
2. **Collective modes.** The optimal witness spreads over 4–8 boundary edges, which adds about 3.7.
3. **Depth.** The second layer adds only ~0.1, so the threshold is a boundary-layer quantity.
4. **Finite N.**

The perturbation transfer to the actual meshes (γ*(N) > 20 for N ≥ N0) is sketched, with a non-explicit N0.

## Errors found in the paper (Section 9)
1. **The chords are not equal.** The mesh "has N equal boundary chords" is false. Square points are equispaced and projected radially, so max|e|/min|e| = 1.43 (N = 16), 1.72 (32), 1.87 (64), 1.97 (256). Only N = 8 is uniform.
2. **ρ is mislabelled.** "ρ ≈ 3.3" (Section 9 and Table 4) is h/max|e|, as printed by `penalty.py`. It is not ρ = h/min|e| as defined in Section 2; the true ρ is 4.70, 5.73, 6.35 at N = 16, 32, 64. The reported γ* = μ*h_Γ/h is correct as defined.
3. **The limit near 22.** It should read ≈ 21.4 (see above).

## Operational note
- **The pypardiso path gives garbage for N ≥ 48.** `svn.run(..., method="direct")` returned values around 1e60 on test P, and was already off in the 3rd digit at N = 32. SuperLU (`solve_lu`, as the paper states) is fine, so `h1_limit.py` uses SuperLU. It is worth checking that no production number came from `solve_direct`.
- **An N = 128 FEM run was dropped** because it exceeded the memory budget. At about the same time the cgroup OOM-killer killed another python3 process (pid 2285, 3.5 GB RSS), possibly the other agent's job. All reported FEM data are N ≤ 96.
