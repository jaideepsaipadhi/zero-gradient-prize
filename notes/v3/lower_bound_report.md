# H¹ lower bound for CNS\* and GS: report

Files:
- `notes/v3/lower_bound.tex`: the theorem and its proof, as a drop-in section. It uses the paper's macros and labels `lb:`. A copy of the paper with `\input{sections/lower_bound}` after §7 compiles with no undefined references.
- `code/lower_bound_tests.py`: the numerics. Logs are `logs/lower_bound_{stars,field,global}.log`.
- `code/plot_leak.py`: writes `paper/figures/leak_field.pdf`, plus a PNG preview at `notes/v3/leak_field_preview.png`.

## 1. What is proved

**Theorem (lb:thm).** Assume (M0)–(M2), (D), and a new local condition (M3). Let u_h be the velocity of CNS\* (γ ≥ γ₂), or of GS (γ ≥ γ₀, **any** G). Let W = ∂_n u·τ be the wall shear, with W ≢ 0. Then

  ‖∇(ũ − u_h)‖ ≥ c_LB κ₀ ‖W‖_{L²(Γ)} h_Γ^{3/2} − C h_Γ².

- c_LB depends only on the mesh constants.
- The bound is independent of γ, μ, ρ, h_Ω, the pressure, and ε_h.

Together with Theorem D, and with Theorem A at G = 0, this gives c h_Γ^{3/2} ≤ ‖∇e‖ ≤ |||e||| ≤ C h_Γ^{3/2}. So the H¹ rate 3/2 is exact for CNS\* and for GS with constant boundary pressure, which includes Scott's shear-flow benchmark.

**Status of (M3).**
- **k ≥ 7:** proved (Lemma lb:k7). The test field is an explicit single-triangle stream-function bubble, φ = b_T²((λ₁−λ₂)² − 1/7).
- **k = 4 (the case in the paper):** not proved for general meshes. It is verified numerically for every boundary vertex of the production meshes, N = 16…1024 (below).
- **k = 5, 6:** not checked.

(M3) is a scale-invariant, finite-dimensional condition on each boundary star. It asks for a field in Y_h¹ supported on the star with a nonzero second-moment functional M.

### Structure of the proof

1. **Test space.** Test the error equation only with Y_h = Z_h ∩ V_h⁰, the divergence-free fields that vanish on Γ_h. On Y_h the following all vanish identically:
   - the pressure error (both methods);
   - GS's traction defect ℛ;
   - the transposed-gradient term;
   - the penalty term.

   What is left is the identity (lb:identity):

   a_h(e,v) − ⟨e, ∂_n v⟩ = −⟨ũ, ∂_n v⟩ − (F, v).

2. **Controlling ⟨e, ∂_n v⟩ by ‖∇e‖ (lb:poinc).** For v ∈ Y_h, ∂_n v·n = −∂_τ(v·τ) = 0 on each chord. If also ∫_e ∂_n v = 0 on every chord (the space Y_h¹), then a Poincaré–trace inequality on T_e bounds ⟨e, ∂_n v⟩ by ‖∇e‖‖∇v‖. No slip bound, and hence no upper-bound machinery, is needed. This gives

   K₀‖∇e‖ ≥ sup over v ∈ Y_h¹ of |⟨ũ, ∂_n v⟩ + (F, v)| / ‖∇v‖.

3. **Leading term (lb:expand).** On a chord, ũ·τ_e = ½(ℓ²/4 − s²)W_e + O(h³). The mean part drops out against the mean-zero ∂_n v·τ, leaving

   ℓ_W(v) = −½ Σ_e W_e ∫_e s² ∂_n v·τ_e,

   with error O(h^{5/2})‖∇v‖.

4. **Test field (lb:assemble).** Take v = Σ_z W(z) v_z, built from the (M3) star fields. This gives ℓ_W(v) ≳ κ₀ ‖W‖ h^{3/2} ‖∇v‖.

## 2. What failed, and why: the suggested transposed-term route

1. **The trace must be normal, not tangential.** The transposed term's leading part s·a_e points along n, since a_e = n(m_e\*)W_e. It therefore pairs with the **normal** trace v·n. A test field with an odd *tangential* trace gives zero at leading order.
2. **Divergence-free does not obstruct an odd normal trace.** For v = curl φ we have ∫_e v·n = 0 exactly when φ takes equal values at the chord ends. The relation ∂_n v·n = −∂_τ(v·τ) constrains the normal derivative, not the trace.
3. **The penalty term, computed as asked.** On an odd normal trace, ũ·n_e = d(s)·W·(τ\*·n_e) with τ\*·n_e = ±s. So the penalty multiplies the transposed term pointwise by 1 + (μ/h)d(s), where 0 ≤ (μ/h)d ≤ γh_Γ/8.
   - For fixed γ this is lower order, and it has the same sign, so it reinforces rather than cancels.
   - It becomes leading once γ ≳ h_Γ⁻¹.
   - On a tangential trace the penalty pairs the even O(h²) bump with v·τ. That is not lower order unless v·τ is zero or odd on the chord.
4. **The Nitsche symmetric term is the same order as the transposed term.** −⟨ũ, ∂_n v⟩ pairs the O(h²) chordal bump with ∂_n v·τ ~ A/h. This is the term the proof above uses.
5. **The fatal problem for CNS\* is the pressure error.** Any test field with v·n ≠ 0 brings in ⟨ξ − ξ̄, v·n⟩. Theorem D's proof bounds this only by Cγ^{-1/2}β^{-1}(|||e||| + h^{3/2})|||v|||, which is the same order as the target. So the lower bound cannot be closed without uncontrolled constants.
6. **For GS with G = 0 the route works, but only in the energy norm.** It gives |||e||| ≳ h^{3/2}/(1+γ^{1/2}). That bound can be carried entirely by the slip part. For CNS\*, test A, N = 64, μ = 100, the slip part (μ/h)^{1/2}‖e‖_Γ is 0.077, against ‖∇e‖ = 0.021. So it says nothing about H¹.

**Why k = 4 needs stars.**
- A single triangle, or one first-layer quad (two triangles), carries no nonzero field in Y_h¹. For the production meshes the computed dimension is 0.
- An Argyris (C²-vertex) stream function cannot work either. Its Hessian at a boundary vertex is killed by the two boundary lines, which forces ∂_nnφ = 0 at the chord ends. With ∂_nnφ|_e ∈ P₃, a mean-zero profile vanishing at both ends is then necessarily odd, and its pairing with the even bump is zero.
- The admissible fields use the Hessian *jumps* across interior edges that general C¹ quintics allow at a boundary vertex with at least 3 triangles, which is (M1′). The boundary star is the smallest patch where these exist.
- Proving that the star space is nontrivial with M ≠ 0 for every admissible geometry is the open step. It needs a Bernstein–Bézier dimension argument that I did not complete.

## 3. Numerics (k = 4, μ = 100; all runs ≤ 3 GB, one core)

**(M3) on the production meshes** (`stars`). Here κ_z := M(v_z)/|e_z|² with ‖∇v_z‖ = 1, and v_z is supported on the three triangles at z.

| N | dim Y¹(star) | min κ_z | mean | max |
|---|---|---|---|---|
| 16 | 1 | 0.00555 | 0.0109 | 0.0179 |
| 64 | 1 | 0.00366 | 0.0145 | 0.0252 |
| 256 | 1 | 0.00303 | 0.0149 | 0.0260 |
| 1024 | 1 | 0.00288 | 0.0149 | 0.0261 |

- The minimum converges to about 0.0028.
- It occurs in the diagonal directions of the square, where the first-layer triangles are most elongated.
- A scratch check on 15 random 3- and 4-triangle boundary stars gave dim Y¹ = 2m − 5 for m triangles, and κ ∈ [0.013, 0.028] in the same normalisation, never zero.

**Assembled test field** v = Σ W(z)v_z (`field`). The table gives ⟨ũ, ∂_n v⟩ / (‖∇v‖ h_Γ^{3/2}):

| N | test A | test B₁ |
|---|---|---|
| 32 | 0.0606 | 0.0184 |
| 128 | 0.0665 | 0.0209 |
| 512 | 0.0670 | 0.0212 |

- The rates reach 1.497 (A) and 1.496 (B₁).
- The leading term ℓ_W(v) agrees to three digits from N = 64 on.

**Exact dual norms on the full spaces** (`global`, test A, N = 16→64; rates in h_Γ):

| quantity | rates |
|---|---|
| sup over Y_h¹ of ⟨ũ, ∂_n v⟩/‖∇v‖ | 1.12, 1.27, 1.35, 1.40 |
| Z_h, energy-dual norm of 𝒢 | 1.70, 1.61, 1.56, 1.54 |
| Z_h, energy-dual norm of the transposed part | 1.61, 1.55, 1.52, 1.51 |
| Z_h, ‖∇·‖-dual norm of 𝒢 | 1.23 → 1.06 |
| Z_h, ‖∇·‖-dual norm of the transposed part | 1.78 → 1.58 |
| CNS\* H¹ error | 1.82 → 1.64 |
| CNS\* slip | 2.18 → 2.02 |

- **Y_h¹ dual norm.** Its ratio to h^{3/2} rises from 0.086 to 0.115, approaching the limit from below. It is 14–24% of the actual CNS\* ‖∇e‖.
- **Z_h dual norms, which had not been measured before.**
  - In the energy norm, both 𝒢 and its transposed part tend to rate 3/2. So Lemma 3.9 is sharp on Z_h too.
  - In the ‖∇·‖ norm, 𝒢 decays only like h. The penalty part dominates there, at size γh. So Theorem D's h^{3/2} genuinely uses the mesh-dependent norm.

## 4. Leak figure (`paper/figures/leak_field.pdf`)

The figure shows GS on test P at N = 32, μ = 100.
- **Left panel:** the radial velocity u_h·e_r near the wall (colour), with u_h arrows. Black arrows show the prediction (h/μ)(p − p̄)n at the chord midpoints.
- **Right panel:** u_h·n and u_h·τ along Γ_h against (h/μ)(p − p̄).

Results:
- The relative L² misfit of the normal trace is 0.053.
- The tangential trace is about 12% of the normal trace and oscillates at the element scale.
- The slip ratio is 0.945, matching Table 3.

## 5. Open items

1. Prove (M3) for k = 4 on general meshes. The route would be to show that the C¹-quintic space on a 3-triangle boundary star that vanishes to second order on the star boundary has dimension 3, and that M does not vanish on it. The numerics say dim = 2m − 3 before the two mean constraints.
2. The theorem's constant c_LB κ₀ is pessimistic. The actual Y_h¹ dual norm is about 0.115 h^{3/2} (test A, N = 64), against ‖∇e‖ ≈ 0.48 h^{3/2}.
