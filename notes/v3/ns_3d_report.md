# Navier–Stokes and 3D extensions: what is proved, under what, and what is open

Files:
- `notes/v3/navier_stokes.tex` (labels `ns:`)
- `notes/v3/three_d.tex` (labels `td:`)
- `notes/v3/face_mean_check.py`, the numerical check behind Lemma td:faces(d).

Both `.tex` files are fragments meant to be `\input` after the paper's sections, in the order navier_stokes, then three_d. They use the paper's macros and labels. I compiled them together with the full paper in a scratch wrapper: no errors, and no undefined references or citations.

Citations that are not in `refs.bib` (Zhang 2005 and 2011, Neilan 2015, Guzmán–Neilan 2018, Farrell–Mitchell–Scott) are written inline in the text, because existing files were not to be edited. They need bib entries before these fragments go into the paper.

---

## Part A: steady Navier–Stokes (2D)

### What Gjerde–Scott actually print
I checked the arXiv HTML of 2306.12362.
- The continuous problem is −νΔu + u·∇u + ∇p = 0.
- The drag functional (their (7)) uses the standard term (u·∇u)·v.
- The Nitsche system (27)–(28) is printed only for the Stokes operator.
- No discrete convective form is printed: no skew-symmetrisation and no inflow/outflow term. The nonlinear solve is "iterated penalty".

So the notes analyse both
- c₁(w;u,v) = ((w·∇)u, v), and
- c_s = ½[((w·∇)u,v) − ((w·∇)v,u)].

### Is b(u_h;u_h,u_h) = 0?
For div-free w and v|∂R = 0 (Lemma ns:bdry):
- c₁(w;v,v) = ½⟨w·n, |v|²⟩_{Γ_h}
- c₁ − c_s = ½⟨(w·n)u, v⟩_{Γ_h}
- c_s(w;v,v) = 0

So with c₁ the answer is no: the term is nonzero whenever u_h·n ≠ 0 on Γ_h. Nitsche never enforces u_h·n = 0, and for GS the normal trace is exactly the leak.

The answer to "is an inflow/outflow correction needed?" is **no**, for the reasons below.
- **What the proofs use.** Stability is assumed (and, under small data, proved) for the linearisation at ũ, and ũ·n = O(h_Γ²) on Γ_h. Switching c₁ ↔ c_s perturbs the linearised operator 𝔄 by O(h_Γ² h/μ). The discrete leak enters only through the quadratic term c(e;e,v). The contraction controls that term by continuity alone (Lemma ns:tri), so no sign is needed.
- **Quantitative criterion (Lemma ns:pe).** νN_h(v,v) + c₁(w;v,v) ≥ (νc₁ − ½ Pe_Γ(w)) |||v|||², where Pe_Γ(w) = (h/μ)‖(w·n)₋‖_∞ is a wall Péclet number. For the GS solution this is formally O((h/μ)²), but that estimate assumes an L^∞ version of the leak asymptotics, which is not proved. A correction (skew form, or an upwind term −⟨(w·n)₋u, v⟩) matters only for iterates far from ũ: coarse meshes, high Re at fixed μ, or the early sweeps of a Picard or iterated-penalty solve.
- **Physical picture.** GS's discrete wall is porous, with Darcy law (νμ/h) u_h·n ≈ p − p̄. Fluid leaves Ω_h where p > p̄ (the front stagnation region) and re-enters where p < p̄.

### Assumptions
- **(NS1)** A smooth solution (u,p) on Ω̄, a smooth extension f̃, and (M0)–(M2) and (D) from the paper.
- **(NS2)** Div-free extension ũ = curl ψ̃. The defect F = f̃ + νΔũ − (ũ·∇)ũ − ∇p̃ is supported in the strip.
- **Penalty convention.** Wall penalty νμ/h, i.e. νN_h. With an unscaled μ̂, replace νμ by μ̂ throughout (Remark ns:conv).
- **Stability**, in one of two forms:
  - **(L)** for GS: a uniform inf-sup of 𝔄 = νN_h + c(·;ũ,·) + c(ũ;·,·) on Z_h in |||·|||.
  - **(L\*)** for CNS\*: uniform stability of the mixed linearised CNS\* operator on Z_h × Π_h, tested on V_h^R.
- **(SD)** λ = N C₁ ‖ũ‖_{H¹}/(νc₁) ≤ 1/4. It implies (L) with β_L = νc₁/2, and implies (L\*) for γ ≥ γ₂′ (Lemma ns:SD). γ₂′ is independent of ν because ν cancels in the pressure absorption.
  - Caveat: (SD) uses the discrete Nitsche coercivity constant c₁, so it is stronger than the continuous small-data condition.

### Proved (complete proofs in navier_stokes.tex)
The Stokes theorems are used as black boxes, with the Navier–Stokes pressure p and G = ‖p − p̄‖_{L²(Γ)}. Write X = h/μ + (1+γ^{1/2})h_Γ^{3/2} + ε_h.

1. **Theorem ns:thmGS (GS), under (L) and γ ≥ γ₀, for either convective form, when X ≤ c\*:**
   - existence and local uniqueness;
   - ‖∇(ũ−u_h)‖ ≤ C[h/μ + (1+γ^{1/2})h_Γ^{3/2} + ε_h];
   - |||ũ−u_h||| ≤ ν⁻¹(h/μ)^{1/2}G + CX;
   - the NS error equals the linearised (Oseen–Nitsche) error up to O(X²) in energy.
   - Method of proof: Banach fixed point for δ = z − u_h around δ₀ = 𝔄⁻¹ℓ₀. The linear analysis of δ_R reuses the Theorem C ansatz δ_R = (h/νμ)v_φ + r. The only new term, c(v_φ;ũ,·) + c(ũ;v_φ,·), is bounded in the H¹-dual norm.
2. **Corollary ns:Bp (leading term).** Under the hypotheses of Corollary B′:
   - ũ − u_h = (h/νμ) v_φ + o((h/μ)^{1/2}) in energy and o(h/μ) in L²(Γ_h);
   - slip and energy ratios → 1, with the universal constant unchanged;
   - ‖∇(ũ−u_h)‖ ≥ c(h/νμ)G, so the H¹ rate 1 is sharp.
   - The H¹ constant is not universal: it now involves the Oseen response.
3. **Theorem ns:thmD (CNS\*), under (L\*):** |||ũ−u_h||| ≤ C[(1+γ^{1/2})h_Γ^{3/2} + ε_h]. Hence ‖ũ−u_h‖_{H¹} ≤ C(h_Γ^{3/2} + h^k) under the same provisos as Theorem D (bounded γ, h_Γ ≥ h^{2(σ_k−k)}).
4. **Remark ns:drag.** For cylinder flow, G ≥ |pressure drag|/√(2π) > 0. So in Gjerde–Scott's own application the printed method always leaks, and at fixed μ its H¹ rate is 1. With μ = 10⁶ the leak is small in absolute terms.

### Open in Part A
- **Nonsingular branch without small data (Remark ns:branch).** Deriving (L)/(L\*) from continuous nonsingularity needs a Schatz/Gårding argument. It requires convergence, with no rate, of the discrete Nitsche–Stokes solution operators on the moving domains Ω_h → Ω, for rough data (L^{3/2} / H^{−1+s}), where Lemma 3.5 (transposed-gradient term) must become an o(1) statement. I expect this to work, but it is not written.
- **Existence for arbitrary (large) data.** Not attempted. It needs the skew form plus a discrete Hopf extension of g_I inside W_h.
- **The drag functional.** The leak's effect on the computed drag is not analysed.
- **L^∞ leak asymptotics.** An L^∞ version of the leak asymptotics (which would make the Péclet estimate for u_h rigorous) is not proved.
- **Outflow boundary.** The outer boundary is kept Dirichlet, as in the paper. A do-nothing outflow boundary, as in channel benchmarks, would bring its own convective-outflow issue and is not treated.

---

## Part B: three dimensions

### Setting
- R = (−L,L)³, D the unit ball.
- P_h a convex inscribed polyhedron with triangular faces F.
- (M0³): shape-regular mesh, c₀h_Γ ≤ diam F ≤ h_Γ.
- **(IS3)**: k ≥ 3 and a uniform inf-sup for (V_h⁰, Π_h ∩ L²₀).

### Findings on the specific questions
- **Sagitta.** Still O(h_Γ²) (r_F²/2). The area Jacobian is 1 + O(h_Γ²).
- **Odd-in-s cancellation.** It **fails** on non-equilateral faces.
  - n(x\*) − n_F = −(x − o_F) + O(h_Γ²), where o_F is the **circumcentre**, not the centroid.
  - So the face mean is −(x_c − o_F) + O(h_Γ²), which is nonzero unless F is equilateral.
  - Checked numerically in `face_mean_check.py`: 0, 0.105 h, and 0.122√2 h for three faces.
  - In 2D the circumcentre of a chord is its midpoint, which is why the cancellation exists there.
- **Transposed-gradient term.**
  - In the ‖∇v‖-dual norm it is generically only **h_Γ** in 3D, not h_Γ^{3/2} (Proposition td:transposed(iii), with an explicit nondegeneracy condition).
  - In the penalty-weighted (energy) dual norm it is C h_Γ (h/μ)^{1/2} = C γ^{−1/2} h_Γ^{3/2}, which is all the error analysis needs. The geometric consistency is therefore still (1+γ^{1/2})h_Γ^{3/2} (Lemma td:G).
- **Side observation, also true in 2D (Remark td:nocancel2d).** Lemma 3.5's cancellation is not needed for Lemma 5.2, Theorem C or Theorem D. It only sharpens a dual-norm bound. This could simplify the paper, but I did not edit it.
- **Normal-load cancellation (the analogue of Lemma 6.3).** It holds in 3D (Lemma td:cancel).
  - div v = 0 gives ∂_n v·n_F = −div_F v_F.
  - Integrating by parts on each face leaves edge terms v·(σ_F ν_{E,F} + σ_{F′} ν_{E,F′}).
  - Their coefficient is O(h_Γ), because the conormals satisfy ν_{E,F} + ν_{E,F′} = t_E × (n_F − n_{F′}) = O(h_Γ) and σ_F − σ_{F′} = O(h_Γ).
  - The edge L¹ norms are converted to face L¹ norms by a polynomial inverse inequality.
- **v_φ without stream functions.** There is no C¹/Morgan–Scott analogue in 3D, so v_φ is built differently.
  - An explicit div-free field w with w|_Γ = −(p−p̄)n (Lemma td:w): w = χ r⁻² q e_r − χ′ r⁻¹ ∇_S φ₀ with Δ_S φ₀ = q.
  - Then v_φ = I_h w − a·(face-bubble sum) − (a V_h⁰ inf-sup correction) (Lemma td:vphi).
  - The bounds it provides replace the W^{1,∞} bound at every place it is used. The construction also works in 2D.
- **Uniform Bogovskii constant over Ω_h.** Proved in 3D (Lemma td:bog), using cones and the fact that face planes separate the star-balls from P_h.
- **Net-flux term.** The 3D flux exponent is k+1, since there is no Newton–Cotes-type extra exactness on triangles. So the Theorem D proviso becomes h_Γ ≥ h².

### Scott–Vogelius stability in 3D (literature, checked via web)
- **Alfeld/barycentric refinements.** k ≥ 3 (Zhang, Math. Comp. 74 (2005)). k ≥ d in general dimension (Guzmán–Neilan, SINUM 2018). Farrell–Mitchell–Scott (arXiv:2211.05494 / SISC) summarise this as "k ≥ d for Alfeld splits".
- **Worsey–Farin splits.** Lower degrees are listed in arXiv:2211.05494 Table 2.1, with a robustness caveat. I did not pin down the exact degree from the primary source.
- **Freudenthal meshes.** k ≥ 6 (Zhang, Math. Comp. 80 (2011)), conjectured k ≥ 4. These meshes cannot fit an inscribed polyhedron, so this is irrelevant here.
- **General tetrahedral meshes.** Only conditional: k ≥ 6 in Neilan (Math. Comp. 84 (2015), Prop. 6.5), under a mesh condition he states as a conjecture.
- The task's guess "Neilan 2015 general meshes k ≥ 6" is therefore correct, but conditional.
- The paper's EickmannGNST2026 reference is 2D only.

### Proved in 3D
- **Theorem td:main.** Under (M0³), (IS3) and (D), the following hold verbatim, with ε_h replaced by ε_h⁽³⁾ = (1+(γh_Γ)^{1/2})h^k + h^{k+1}(h_Γ^{−1/2} + (μ/h)^{1/2}):
  - Theorems A, B, C, D, Corollaries A′ and B′, Proposition 5.1, and Proposition 8.1(a);
  - that is, the GS leak (h/μ rate in H¹, the exact leading term with constant 1) and the CNS\* rate h_Γ^{3/2} + h^k (for h_Γ ≥ h²).
  - The proof is an item-by-item replacement list. The two non-mechanical changes are: Theorem C step 4 now uses an energy-dual bound (the result is unchanged, in fact better: |||r_h||| ≲ h/μ + h_Γ²), and Theorem B step 2 has an O(h_Γ) rather than O(h_Γ^{3/2}) term, which is absorbed by h_Γ ≤ h₀(G).
- **Corollary td:ns.** All of Part A in 3D, since H¹ ⊂ L⁴.
- **Unconditional case.** On Alfeld-refined meshes with k ≥ 3, subject to the caveat on uniformity below.

### Open in 3D (exact reasons)
1. **Penalty necessity (Proposition 8.1(b)).** The scaling carries over (N_h(v,v) = ℓ[A − 2B + (μℓ/h)C]). But no 3D witness v ∈ Z_h on a boundary star with 2B > A has been computed or certified. That needs a 3D analogue of `lamstar_exact.py`.
2. **(IS3) on unsplit tetrahedral meshes.** This is an open problem in the literature.
3. **Uniformity of the split-mesh inf-sup constant in h.** The domain part is proved (Lemma td:bog). That the discrete part of Zhang's and Guzmán–Neilan's proofs depends only on shape regularity is taken from their statements; I did not re-audit their proofs.
4. **Sharpness of h_Γ^{3/2} for CNS\*.** Not proved, in 2D or 3D. In 2D it is numerical evidence only.
5. **Non-convex D.** Not treated, in either dimension, because Ω ⊂ Ω_h is used.
6. **No 3D computations were run.**
