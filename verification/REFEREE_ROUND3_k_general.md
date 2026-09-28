# Third referee pass: generalisation to every k ≥ 4 (independent agent), and responses

## Report (summary of findings)

No fatal defects. The new approximation lemma (interpolate, then correct the divergence) and the uniform inf-sup lemma are correct and non-circular.

- **Lemma 4.6 (approximation):** interpolation on boundary triangles, the net-flux bound, the edge-bubble field and its norms, the divergence correction via Lemma 4.5, and membership z ∈ W_h were all verified. The claimed ε_h is a valid upper bound.
- **Lemma 4.5 (uniform inf-sup):** Guzmán–Scott Theorem 1 needs no quasi-uniformity, no small h and no interior vertex per triangle. The domain enters only through the continuous (Bernardi–Raugel) inf-sup constant and a Poincaré step, and their boundary Θ matches (M2). So the Galdi cover argument closes the gap, and the argument is not circular.
- **Degree k elsewhere:** checked and correct. This covers v_φ (C¹ P_{k+1} interpolant), the O(h^k) substitutions in Theorems A, B and C, Corollary B′ and Theorem D (the η term, and absorption with γ ≳ β⁻²), the P4 certificate implying necessity for all k ≥ 4, and the Scott–Vogelius range for k ≥ 4.

Defects raised:

1. **Major.** The abstract and introduction stated Theorem D without its boundary-edge proviso. The flux term is a genuine floor: every z ∈ W_h carries the flux mismatch of g_I through Γ_h.
2. **Major.** The flux bound |m| ≤ Ch^{k+1} was cruder than v1's. For even k, closed Newton–Cotes gives h^{k+2}, which is Boole's h⁶ at k = 4. So v2 as drafted weakened k = 4.
3. **Minor.** ε_h contained redundant or artefact terms, including ρ^{1/2}.
4. **Minor.** The draft remark claiming that v1's Argyris proof lost a power of h could not be reproduced. ∇ curl δs is a second derivative of δs, so v1's h⁴ was right.
5. **Minor.** Lemma 4.5 citation details were missing: Guzmán–Scott Theorem 1, the Poincaré constant, arcsin(r_B/r_b), α < π/2, and how the wedges are truncated.
6. **Minor.** Lemma 4.3 step 1 justified div V_h⁰ = Π_h ∩ L²₀ on the doubly connected domain only loosely.
7. **Minor.** The dependence of the Morgan–Scott interpolant's constants on the angles was not stated. v_φ ∈ Z_h also needs h small, not only h_Γ.
8. **Minor, carried over from v1.** The Theorem A lower bound dropped h^k terms in steps 2 and 4.
9. **Cosmetic.** Several small fixes.

## Response (v2)

| # | Action |
|---|---|
| 1 | The abstract, introduction and Theorem D now state the boundary-edge proviso (h_Γ ≥ h_Ω⁴ for even k, h_Ω² for odd k). The lower bound showing the flux term is genuine is reinstated after Lemma 4.6. |
| 2 | Flux bound sharpened to h^{σ_k}, with σ_k = k+2 for even k and k+1 for odd k, using closed Newton–Cotes exactness. k = 4 is now as strong as in v1. |
| 3 | ε_h is now (1+(γh_Γ)^{1/2})h^k + h^{σ_k}(h_Γ^{-1/2} + (μ/h)^{1/2}); the ρ term is gone. |
| 4 | Accepted. The claim was withdrawn before release; the remark now says only that v1 used Argyris matching specific to k = 4. |
| 5 | All citation details fixed as listed. |
| 6 | Lemma 4.3 step 1 now derives the divergence range from Lemma 4.5, with no circularity. |
| 7 | The angle dependence is stated, and (D) now assumes h ≤ h_*. |
| 8 | The h^k terms are kept, and the hypothesis is now h ≤ h_0(G). |
| 9 | Fixed. |
