# Pressure, boundary stress and drag for GS and CNS*: report

This report accompanies `notes/v4/pressure_drag.tex`, which compiles standalone to 14 pp. with no warnings. The sanity-check script is `notes/v4/pressure_drag_check.py`; it imports `code/svn.py` and uses SuperLU. Result numbers such as "[P, Lemma 4.5]" refer to the compiled `paper/main.pdf` on `v3-draft`. No existing file was edited, and nothing was committed.

The setting is that of [P, §2]: the disk in the square with L = 2.5, k ≥ 4, and assumptions (M0)–(M2) and (D). Throughout, γ ≥ γ₀ for GS and γ ≥ γ₂ for CNS*. Notation:

- p° := p̃ − mean_{Γ_h}(p̃), which is the p* of Theorem D.
- q := p − p̄ on Γ, and G = ‖q‖_{L²(Γ)}.
- L₁ is the **first layer**: the triangles that have an edge on Γ_h. It has area of order h_Γ.
- H_e is the height of T_e over the edge e.
- Y := (1+γ^{1/2})h_Γ^{3/2} + ε_h, and X₁ := h/μ + Y.
- m_R := ∫_{∂R} g_I·ν. It is zero for flux-corrected data and for all test data in [P].

## 1. Gjerde–Scott's drag (checked in the arXiv HTML of 2306.12362)

GS compute drag in two ways and compare them.

- **β (boundary integral).** β_h(χ_h) = ∮_{Γ_h} ((νD(u_h) − p_h I) χ_h)·n, using pointwise derivatives and the pointwise discrete pressure.
- **ω (volume, Babuška–Miller).** ω_h(χ_h) = ∫_{Ω_h} (ν/2) D(u_h):D(χ_h) + (u_h·∇u_h)·χ_h. They drop the term −p div χ because div χ_h = 0.
- **χ_h** is computed "with the same scheme": it is the SV–Nitsche Stokes solution with χ = (1,0) on Γ_h (imposed by Nitsche) and zero on the outer boundary.
- **Inconsistency parameter.** They report ε = (ω − β)/ω (their eq. 10, Table 1b). They use μ = 10⁶, h = the maximum mesh size, and k = 4.

## 2. Results

### 2.1 Tools (§2 of the note)

- **Pressure identity on V_h⁰ (Lemma 2.1), for both methods:** (p° − p_h, div v) = a_h(e_u, v) + ⟨u_h, ∂_n v⟩ + (F, v).
  - The Nitsche symmetric term applied to the discrete trace (the "Nitsche coupling") is the only boundary source in the pressure.
  - For v ∈ V⁰, div v|_{T_e} = ∂_n v·n_e on e. So the **normal trace u_h·n acts on the pressure exactly like an edge source**.
- **Edge representers (Lemma 2.2).** On a triangle T_e, the P_{k−1} representer of q ↦ ∫_e q has
  - ‖·‖²_{T_e} = ĉ_k |e|/H_e, with **ĉ_k = k(k+1)** (exact rational computation for k ≤ 8; ĉ₄ = 20);
  - ‖·‖²_{L²(e)} = d̂_k |e|/H_e².
- **Normal-load cancellation without div = 0 (Lemma 2.3).** This extends [P, Lemma 6.3] by the extra term ‖div δ|_{T_e}‖_{L¹(e)}.

### 2.2 Pressure

**Theorem 3.1 (CNS*), proved.**

- ‖(p° − p_h) − mean‖ ≤ C Y.
- |mean| ≤ C[(1 + (γh_Γ)^{1/2}) Y + (μ/h)|m_R|].
- Hence ‖p° − p_h‖_{L²(Ω_h)} = O(h_Γ^{3/2} + ε_h) for bounded γ and flux-corrected data.
- The CNS* discrete pressure is unique (Thm D). Its Γ_h-mean λ_h = p̄_{Γ_h}(p_h) → 0.
- **Not proved:** a lower bound.

**Proposition 4.1 (GS pressure is fully determined), proved.**

- The GS pressure is unique, with **no free constant**. The constant is fixed by the flux-carrying test functions, because div V_h^R = Π_h.
- It approximates p − p̄_{Γ_h}(p).

**Theorem 4.3 (GS pressure), proved.** Define the explicit first-layer field ξ_b ∈ Π_h, supported in L₁, by

  (ξ_b, r) = (h/μ) ⟨q∘π, r|_{T_e}⟩_{Γ_h}  for all r ∈ Π_h.

Put ρ_h := p° − p_h − ξ_b. Then:

- (a) ‖ρ_h − mean‖ ≤ C_p X₁, and ‖ρ_h‖ ≤ C_p[(1 + (γh_Γ)^{1/2}) X₁ + (μ/h)|m_R|].
- (b) ‖ξ_b‖² = (h/μ)² [ĉ_k Σ_e (|e|/H_e) q(m_e*)² + O(1)], which is ≍ (h/μ)² h_Γ^{−1} G².
- (c) For every constant c, ‖p° − p_h − c‖_{L²} ≥ c₁(h/μ) h_Γ^{−1/2} G − C_p Q. So the **L² rate is exactly 1/2** at fixed γ, and sharp.
- (d) Off the first layer, the error is O(h/μ): rate 1.
- (e) For every constant c, ‖p° − p_h − c‖_{L^∞(L₁)} ≥ c₂ G/γ − o(1). So there is **no pointwise convergence** in the first layer at fixed γ.
- If G = 0, then ξ_b = 0 and GS behaves like CNS*.

Mechanism: the leak u_h·n ≈ (h/μ)(p − p̄) feeds the pressure through the Nitsche coupling, as an edge source that the discontinuous pressure resolves in the first layer, at pointwise size (h/μ)/H_e ≍ 1/γ. The leak flow (U, P) itself only contributes (h/μ)P to ρ_h.

### 2.3 Drag, lift and torque (ψ = a + b x^⊥, rigid motions)

Definitions:

- J_ψ = ∫_Γ σ(u,p) n·ψ, where n is the outward normal of Ω.
- v_ψ := curl Π_{k+1}(χ₀ φ_ψ) ∈ Z_h, which equals ψ on the first layer.
- ω_h(w) := a_h(u_h, w) − (p_h, div w) − (f̃, w).
- Θ_f(ψ) := ∫_{S_h} f̃·ψ = Σ_e (|e|³/12)(f·ψ)(m_e*) + O(h_Γ³). This is the body force in the chord–arc strip; it is **zero when f = 0 near Γ**, as for cylinder drag and test A.

**Lemma 5.2, proved.**

- (a) The **variational (Babuška–Miller) drag with any discrete test field equal to ψ on L₁ is algebraically the Nitsche-flux drag**: J^N_h := ω_h(v_ψ) = ∫_{Γ_h} t_h·ψ + ⟨u_h, ∂_nψ⟩. The Nitsche flux is t_h = ∂_n u_h − θ(p_h − p̄_{Γ_h}(p_h)) n − (μ/h) u_h.
  - For GS, J^N_h contains no pressure at all. The penalty term (μ/h)u_h ≈ −(p − p̄)n reproduces the missing pressure traction.
- (b) J^N_h − J_ψ = −a_h(e_u, v_ψ) − Θ_f(ψ), for both methods.

**Lemma 5.3 (reciprocity), proved.** Let (Z_ψ, R_ψ) be the dual Stokes problem, with Z = ψ on Γ and 0 on ∂R. Then

  J_ψ(U) = K_ψ := ∫_Γ (p − p̄)(R_ψ − R̄_ψ).

Here U is the leak flow of [P, Def. 6.8]. The proof uses σ(Z, R)n·n = −R on Γ, by the wall structure and the skewness of ∇ψ.

**Theorem 5.4 (GS, J^N), proved.**

- J^N_h − J_ψ = −(h/μ) K_ψ − Θ_f + R_h, with |R_h| ≤ C(h/μ) B_U + C[(1+γ^{1/2}) h_Γ^{3/2} + ε_h]. Here B_U is the bound of [P, Thm 6.10], which is O(h^{1/2}) at fixed γ.
- So the rate is exactly 1 when K_ψ ≠ 0.
- The unconditional bound |J^N − J| ≤ C_p h/μ + C Y holds without Thm 6.10.

**Theorem 5.6 (CNS*, J^N), proved.** Assume γ ≥ max(1, γ₂) and γ h_Γ ≤ 1.

- (a) |J^N − J + Θ_f| ≤ C Y (rate 3/2).
- (b) |J^N − J + Θ_f + Λ_h| ≤ C[γ^{−1/2} h_Γ² + γ h_Γ³ + γ^{1/2} h_Γ^{3/2} E_Z + ε_h]. This gives **rate 2**.
  - The geometric term is Λ_h := ⟨ũ, ∂_n(Z̃_ψ − ψ)⟩_{Γ_h} = Σ_e (|e|³/12)(∂_n u·∂_n(Z_ψ − ψ))(m_e*) + O(h_Γ³).
  - −Θ_f − Λ_h is the shape derivative of J_ψ under the polygonal displacement ½(ℓ²/4 − s²). It is method-independent.
  - The leading term is identified exactly in the limit γ → ∞ with γh_Γ → 0. At fixed γ it is identified only up to O(γ^{−1/2} h_Γ²), which comes from the weighted normal slip (Lemma 5.7).

**Proposition 5.8 (GS's own volume formula J^χ = ω_h(χ_h), translations).**

- GS: J^χ = J − **2**(h/μ) K_e − Θ_f + O(h^{3/2}), for fixed γ, bounded ρ and E_Z + E_U ≤ Ch. The test field χ_h leaks as well, which doubles the constant.
- CNS*: |J^χ − J^N| = O(h_Γ² + (h/μ)^{1/2} E_Z).

**Proposition 5.9 (tractions on Γ_h), proved (upper bounds, plus a two-sided GS pointwise bound).**

- CNS*:
  - ‖t_h − σ°n‖_{L²(Γ_h)} and ‖σ(u_h,p_h)n − σ°n‖_{L²(Γ_h)} are both O(h_Γ).
  - |J^β − J| = O(h_Γ^{−1/2} Y), i.e. rate 1.
- GS:
  - ‖t_h − σ°n‖_{L²(Γ_h)} = O(h_Γ^{1/2}).
  - The pointwise traction satisfies **‖σ(u_h,p_h)n − σ°n‖_{L²(Γ_h)} ≍ G/γ, which does not converge**.
  - J^β = J + Ξ_h(ψ) + O(h_Γ^{−1/2} X₁), where
    - Ξ_h(ψ) := ⟨ξ_b, ψ·n⟩_{Γ_h}, and for translations Ξ_h(e) = (h/μ) ĉ_k Σ_e (|e|/H_e) q(m_e*) n_e·e + O(h/μ);
    - Ξ_h = O(1/γ) and **does not vanish** as h → 0 at fixed γ;
    - a *lower* bound on |Ξ_h| is not proved (the weighted sum is sign-indefinite).
- Remark 5.10: the GS inconsistency parameter ε → −Ξ_h/J (plus O(h/μ)). It plateaus at order ĉ_k (h/(μh_Γ))(|e|/H_e) |J^p|/|J|, which is about 10⁻⁵–10⁻⁴ for μ = 10⁶. That is the order of GS's reported 8·10⁻⁵. This is a remark only: we have not checked it against their meshes.

### 2.4 Navier–Stokes (§6 of the note)

Assume (L) for GS and (L*) for CNS*, and use the scheme's convective form in ω_h.

Proved:

- Both pressure theorems transfer. The GS layer is the same ξ_b with the Navier–Stokes pressure. In the νN_h convention it does not depend on ν; with the unscaled penalty μ̂ it is g = ν(h/μ̂)q.
- The drag error identity acquires the linearised convective terms.
- |J^N − J| = O(X) (rate 1) for GS and O(Y) (rate 3/2) for CNS*.

**Not proved:**

- The GS leading constant. It needs an Oseen analogue of [P, Thm 6.10]. The conjectured constant is K† = ∫(p − p̄)R†, with R† the pressure of the adjoint Oseen dual.
- The CNS* rate 2. It needs an adjoint (L*).

### 2.5 General domains (Remark 6.2)

- The pressure theorems and the CNS* traction bounds hold on the general domains of [P, §10].
- S_h becomes Ω_h △ Ω (with signs), and Λ_h carries the curvature through d = ½ϰ(ℓ²/4 − t²).
- The GS leading drag term is only proved on the disk (it uses Thm 6.10). On general domains only the O(h/μ) bound is proved.

### 2.6 Not proved (§7)

- CNS* pressure lower bound.
- Sharpness of the CNS* rate 2 at fixed γ.
- A lower bound on |Ξ_h|.
- J^χ for torque.
- The Navier–Stokes leading terms.
- ĉ_k = k(k+1) for k > 8.

## 3. Sanity checks run here (k = 4, meshes of [P, §13], μ = 100 unless stated)

These were run with `python3 notes/v4/pressure_drag_check.py <mode> ...`.

**Test P, GS (G = 1.658).**

| N | ‖p°−p_h‖_{L²(L₁)} | ‖ξ_b‖ (pred.) | ‖·‖ off L₁ | ‖·‖_{L∞(L₁)} | J^N−J (e_x) | pred. −(h/μ)K−Θ_f | J^β−J (e_x) | Ξ_h (e_x) |
|---|---|---|---|---|---|---|---|---|
| 16 | 0.174 | 0.169 | 0.214 | 0.63 | −0.2669 | −0.3022 | −1.107 | −1.018 |
| 32 | 0.130 | 0.130 | 0.127 | 0.73 | −0.1417 | −0.1507 | −1.170 | −1.132 |
| 64 | 0.0958 | 0.0961 | 0.0689 | 0.76 | −0.0731 | −0.0753 | −1.211 | −1.192 |
| 96 | 0.0793 | 0.0796 | 0.0472 | 0.77 | −0.0492 | −0.0502 | −1.225 | — |

- **Layer.** The first-layer rate goes 0.45 → 0.47, approaching 1/2. The rest has rate 0.90 → 0.93, approaching 1. The L^∞ norm in L₁ does not decrease.
- **Leak drag.** The ratio of the actual to the predicted J^N defect is 0.883 → 0.940 → 0.970 → 0.981. For e_y at N = 64 and 96 it is 0.970 and 0.980.
- **GS's J^χ.** Its ratio to −2(h/μ)K − Θ_f is 0.905 (N = 32) and 0.952 (N = 64).
- **K values.** K_x ≈ 17.16 and K_y ≈ −25.36 for p = x²y − y + x/2 and L = 2.5 (mode `K`, CNS* dual with a volume formula; converged to about 0.1%). Scale by a for B_a.
- **CNS* on test P is exact:** its J^N error equals −Θ_f to machine precision.

**Test B₁, CNS*.**

| N | ‖p°−p_h‖_{L²} | first layer | off L₁ | L^∞(L₁) | lift J^N−J | −Θ_f−Λ_h | drag J^N−J | drag J^β−J |
|---|---|---|---|---|---|---|---|---|
| 16 | 0.497 | 0.0958 | 0.488 | 0.54 | −0.4169 | −0.4676 | −0.0813 | −0.296 |
| 32 | 0.0767 | 0.0354 | 0.0681 | 0.34 | −0.1153 | −0.1218 | −0.0258 | −0.187 |
| 64 | 0.0162 | 0.0123 | 0.0106 | 0.19 | −0.0298 | −0.0308 | −0.00727 | −0.104 |
| 96 | 0.0080 | 0.0066 | 0.0045 | 0.13 | −0.01338 | −0.01369 | −0.00334 | −0.0708 |

- **Pressure.** The first-layer rate is 1.55; the L^∞(L₁) rate is about 0.95.
- **Lift.** The ratio of (J^N − J + Θ_f) to (−Λ_h) is 0.788 → 0.895 → 0.941 → 0.955.
- **Drag component.** Here Λ_h ≈ 10⁻⁶, and the residual J^N − J + Θ_f = −0.0145, −0.0044, −0.0021 has rate approaching 2. It shrinks with μ: at N = 64 it is −0.0044 (μ = 100), −0.0010 (μ = 10³) and −0.00015 (μ = 10⁴). This is the O(γ^{−1/2} h_Γ²) remainder.
- **Pointwise traction.** J^β for the drag component converges at rate 0.73 → 0.86 → 0.95, approaching 1.
- **Penalty dependence.** At fixed N the CNS* pressure error *grows* with μ, consistent with the (1+γ^{1/2}) factor. At N = 64 the L² error is 0.016, 0.039 and 0.059 for μ = 10², 10³, 10⁴; the L^∞(L₁) error is 0.19, 0.66 and 1.10.

**Test A, torque (ψ = x^⊥, J = −4π).**

| N | CNS* error | GS error | −Λ_h(x^⊥) | ratio (CNS*) |
|---|---|---|---|---|
| 16 | 0.348 | 0.370 | 0.413 | 0.84 |
| 32 | 0.0974 | 0.100 | 0.106 | 0.92 |
| 64 | 0.0254 | 0.0258 | 0.0267 | 0.95 |

- Drag and lift of test A are 0 by symmetry, for both methods, to round-off. Use the **torque**.

## 4. Predicted rates for the numerics agent

Assume fixed μ, k = 4, and the meshes of [P] (h/h_Γ ≈ 3.3), so rates in h and h_Γ coincide. Report every functional for e_x, e_y and x^⊥.

**Pressure.** Measure ‖p_h − (p − mean_{Γ_h} p)‖ in full L², mod constants, on L₁, off L₁, and in L^∞(L₁).

| | GS (G > 0) | CNS* (and GS with G = 0) |
|---|---|---|
| L²(Ω_h) | **1/2** (asymptotically); ratio to ‖ξ_b‖ → 1 on L₁ | **3/2** (preasymptotic 2–3 at N ≤ 64; the L₁ part is already at 1.5) |
| L²(Ω_h \ L₁) | **1** | ≥ 3/2 (observed ~2 up to N = 96) |
| L^∞(L₁) | **0** (≈ const·G/γ; 0.77 on test P at μ = 100) | ≈ **1** (observed) |
| μ-dependence at fixed N | error on L₁ ∝ 1/μ | grows like γ^{1/2} preasymptotically |

- Test P with CNS* is exact to round-off.
- Test C exercises η = p° − π_h p°, which is higher order.

**Drag, lift and torque (Stokes tests).**

| functional | GS | CNS* |
|---|---|---|
| J^N = Nitsche flux = volume formula with v_ψ | **1**, error → −(h/μ)K_ψ − Θ_f | **2**, error → −Θ_f − Λ_h (+ O(γ^{−1/2}h_Γ²)) |
| J^χ = ω_h(χ_h) (GS's formula) | **1**, error → −2(h/μ)K_e − Θ_f | **2** (≈ J^N) |
| J^β = pointwise σ(u_h,p_h)n | **0**: tends to Ξ_h = O(1/γ) (≈ −1.2 for e_x on P/B₁ at μ = 100) | **1** (one component observed ~2) |
| t_h in L²(Γ_h) | ≥ 1/2 | 1 |

- **Remove Θ_f** (or integrate f over Ω ∩ Ω_h instead of Ω_h) to see the pure discretisation error. Θ_f is O(h_Γ²), with coefficient (1/12)Σ|e|³ f·ψ. It is nonzero for P, B and C, and zero for A.
- Λ_h and K_ψ need the dual Stokes solution (Z_ψ, R_ψ); the modes `lam` and `K` compute them with CNS*.
  - Λ_h uses the tangential Nitsche flux of the dual minus ∂_nψ.
  - K uses the volume identity K = −[a(Z, W) − (R, div W)] with W = (p − p̄)n extended. The pressure-trace version converges only at rate 1.

**Cylinder drag (GS's setting: f = 0, so Θ_f = 0; μ = 10⁶).**

- **Both methods:** the drag error is dominated by the polygonal geometric term **−Λ_h ≈ −(1/12) Σ_e |e|³ (∂_n u·∂_n Z_e)(m_e*)**, i.e. **rate 2 in h_Γ**. Λ_h is computable from the primal and dual wall shears.
  - For CNS* this is proved as the limit γ → ∞, γh_Γ → 0; μ = 10⁶ is deep in that regime once γh_Γ ≤ 1.
  - For GS the leak adds −(h/μ)K_e (J^N) or −2(h/μ)K_e (J^χ). Relative to the drag this is ~10⁻⁶h: invisible.
- **Pointwise-traction drag of GS:** it plateaus at Ξ_h/J ~ ĉ_k (h/(μ H)) |J^p|/|J|, i.e. 10⁻⁵–10⁻⁴. Prediction: GS's ε plateaus there instead of → 0.
  - A μ-sweep at fixed mesh separates the pieces: Ξ_h ∝ 1/μ and K-terms ∝ 1/μ, while Λ_h is μ-independent.
- **Navier–Stokes:** proved rates are 1 (GS) and 3/2 (CNS*). The expected rates are as for Stokes, with adjoint-Oseen constants; these are not proved.

## 5. Caveats

- K, Λ_h and E_Z involve the dual problem. The corners of the square limit Z to H²-type regularity in the proofs (E_Z ≤ Ch is all that is used). Numerically this is not limiting.
- The GS leading-term remainders are proved only as O(h^{3/2}) (relative O(h^{1/2})). They are observed as relative O(h).
- Everything is for rigid motions ψ, except J^χ, which is for translations only.
