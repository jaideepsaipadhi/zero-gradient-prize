# Second-round hostile referee report (independent agent), and responses

Second-round referee report on `/home/claude/svn115/proof.md`. I cross-checked it against `/home/claude/podres/results/report.txt`. No files were edited.

The main defects are:
- the regime condition in Theorem B(i), which does not match what its proof needs;
- the claim in Section 7 that D isolates the traction response, which is false;
- a misreading of the μ = 10⁴ data;
- the perturbation claim for the reference star 𝒮₄ in Proposition 6.1.

Nothing is fatal. The proofs of Theorem D and Corollary B′ are correct.

## Defects

**1. Major. Theorem B, regime hypothesis and the proof of (i), step 5.**
- **Problem:** Write X := C_p h_Γ/γ + (1+γ^{1/2})h_Γ^{3/2} + ε_h. Step 5 bounds the term by C·h_Γ^{-1/2}·X·‖v‖_{Γ_h}, with ‖v‖_{Γ_h} ≈ G. For that to be o(G²), you need X ≤ δ·h_Γ^{1/2}·G. The stated regime only gives X ≤ δ(h_Γ/γ)^{1/2}G², which implies what is needed only when G ≲ γ^{1/2}. Test B₁₀₀ has G = 166 and γ^{1/2} ≈ 5.5, so it falls outside.
- The hypothesis also scales wrongly with pressure amplitude a: X is linear in a but the right-hand side is quadratic, so for large a it allows h/μ that is not small.
- Step 4 has the same issue: it needs h/μ ≤ δG, which does not follow.
- There is also a bookkeeping problem: Section 0 lets C depend on the norms of p, so tracking G explicitly is inconsistent. For example, ‖∂_n v_φ‖ is proportional to C_p, not to C.
- **Fix:** Replace the regime with X ≤ δ·h_Γ^{1/2}·G and h_Γ/γ ≤ δG. Or, more simply, derive the lower bound in B(i) from Corollary B′, whose slip ratio tends to 1, and delete step 5.

**2. Major. Section 7: "D := u_GS − u_CNS\* isolates the response to the missing traction" is false.**
- **Problem:** Split the data into the ψ part and the pressure part. Then D = D₀ + a·D₁, where D₀ = (GS − CNS\*) applied to the ψ part. D₀ is nonzero because p_h enters the extra term of CNS\*.
- Test A has G = 0, yet its D is nonzero: at μ = 100, N = 256, ‖∇D_A‖ = 1.346e−3. That is 27% of ‖∇D_{B₁}‖ = 4.930e−3 (a different ψ, but the same order of magnitude). The contamination is O(h^{3/2}) against O(h).
- This likely explains why the normalized H¹ constant still drifts (2.655 → 2.735) and depends on μ (2.94 and 3.10 at N = 192).
- **Fix:** Run GS on a pure-pressure problem: u = 0, f = ∇p. CNS\* reproduces it exactly (u_h = 0, p_h = p\*, since p is cubic and so lies in Π_h), so the GS output is exactly the traction response. Also report the ratios of Corollary B′ for GS itself. From the report, at B₁, μ = 100, N = 256, the GS slip ratio is 0.995 and the energy ratio 1.001, which support the corollary better than D does.

**3. Major. Section 7, Table 3 paragraph: "with μ = 10⁴ … reach 0.987 by N = 192, as the regime condition γh_Γ ≪ G requires."**
- **Problem:** γh_Γ = μh_Γ²/h. At N = 192 this is about 74, and at N = 16 about 1000, against G = 1.66. So the regime condition is badly violated throughout.
- The 0.987 is D's ratio. The quantity the corollary is about, ũ − u_GS, has an energy ratio of 4.914e−2·(10⁴/0.1445)^{1/2}/1.658 ≈ 7.8 and a slip ratio of ≈ 7.6 there, far from 1.
- At μ = 100, N = 256, γh_Γ ≈ 0.55 against G = 1.66, so even the main series is only marginally in the regime.
- **Fix:** State that D cancels most of the geometric error, so D's ratios do not test the regime condition. Report the GS ratios and say that the corollary's hypothesis fails for μ = 10⁴ at these mesh levels.

**4. Major. Proposition 6.1(b), "same holds, by continuity, … small perturbation of 𝒮_m" for 𝒮₄.**
- **Problem:** The rim vertex (0,1) of 𝒮₄ is singular: its rim edges are collinear along y = 1, and its spoke is vertical. That is why the kernel has dimension 17 rather than 16. The count is 28 free scalar nodes × 2 = 56, minus 40 divergence conditions gives 16, and the singular vertex removes one divergence constraint. For 𝒮₃ the same count gives 44 − 30 = 14, which matches.
- A generic perturbation drops the dimension to 16, so the paper's own premise, that "the dimension of the divergence-free subspace is locally constant", fails for 𝒮₄. The witness need not persist.
- **Fix:** Restrict the continuity claim to 𝒮₃. Or move (0,1) off the line, for example to (0, 6/5), and redo the certificate.

**5. Minor. Lemma 2.4, the flux floor.**
- ∫_{Γ_h} z·n = −F_g, not +F_g (sign).
- "Genuine floor" needs |F_g| ≳ h⁶, which is plausible but not shown.
- For test B, F_g = 0 exactly, because g·ν is a polynomial of degree ≤ 5 on each side of ∂R and Boole's rule is exact for it. Test A presumably also has F_g = 0. So the floor is never exercised numerically.
- **Fix:** Use flux-corrected boundary data (∫_{∂R} g_I·ν = 0). This removes both the μ^{1/2}h^{5.5} term and the proviso h_Γ ≥ h⁴ in Theorem D.

**6. Minor. Lemma 2.4, step 4.**
- "Values and gradients need no correction" is false once the carrier is subtracted: g̃ − g_I = O(h⁶) at the vertices.
- The midpoint normal-derivative mismatch is 0 for g_I, not O(h⁵), because the midpoint is both a P4 Lagrange node and an Argyris degree of freedom. It is O(h⁶) only through the carrier.
- Neither affects the result. Say that every carrier-induced mismatch is O(h⁶).

**7. Minor. Theorem C, steps 4–5.** Two terms are dropped. Since h_Γ ≪ h is allowed, neither is automatically controlled.
- In m: an O(h⁴)·‖v‖_{L¹(Γ_h)} term, which is not ≤ h_Γ² in general.
- In ℓ: a (h/μ)·Ch⁴·h_Γ^{-1/2} term.
- **Fix:** Write |m| ≤ C_p(h_Γ^{3/2} + h⁴)‖∇v‖ and |ℓ| ≤ (h/μ)(C_p + Ch⁴h_Γ^{-1/2})‖∇v‖. Corollary B′ is unaffected under its ε_h hypothesis.

**8. Cosmetic. Theorem C, step 1.** δ_G should solve N_h(δ_G, v) = 𝒢(v) − N_h(ũ − z, v), not + N_h(ũ − z, v).

**9. Cosmetic. Theorem A lower bound, step 1.** The displayed integrand, with its factor (−n(x\*)·n_e), equals −G², while the text says G². Drop the minus sign, since ℛ = −⟨p̃ − p̄, v·n⟩ and v·n ≈ −(p − p̄).

**10. Minor. Corollary A′.** The hypothesis should read ε_h ≪ (h_Γ/γ)^{1/2}G, not h⁴ ≪ …. The terms μ^{1/2}h^{5.5} and (γh_Γ)^{1/2}h⁴ must also be small.

**11. Minor. Theorem B(i), step 2.** ũ·v also contains an O(d·h⁴) piece from v_φ − w, so the penalty term is ≤ C(γh_Γ³ + γh_Γh⁴).

**12. Minor. Theorem D, the Galdi cover.**
- The pieces are not "fixed": they are wedge ∩ Ω_h, which depends on h.
- Star-shapedness near the chords needs the ball on the outer side of every chord line that meets the wedge. That requires a quantified margin, such as (r_b − r_B)·cos(α + r_B/r_b + h_Γ) > 1, not just r_b cos α > 1.
- Galdi's constant also depends on the overlap measures between pieces. These are uniform, but the paper should say so.
- Cite the exact step of Guzmán–Scott where the domain enters only through the continuous inf-sup constant, and confirm that their argument needs shape regularity only, not quasi-uniformity.
- With these added, the argument is sound.

**13. Minor. Section 4 closing remark and Section 7, the H¹ constant.**
- "≈ 2.73–2.77 for test B (Table 3)": Table 3 only shows 2.655–2.735.
- "Settles" overstates it. For B₁₀₀ at μ = 100 the constant is still rising (+0.012 per level); for B₁₀₀ at μ = 10⁴ it is still falling (−0.008 per level). Together with item 2, the data do not yet establish a μ-independent constant.

**14. Minor. Table 1: "independently of the pressure."**
- This holds identically, for a structural reason: p is cubic, so p ∈ Π_h, F ≡ 0, and CNS\* reproduces the pressure part exactly. The η term of Theorem D is therefore never tested.
- B₁₀₀ was not run at N = 256, so "agree for every N" is false.
- **Fix:** Add a test with a non-polynomial pressure.

**15. Minor. Table 4 text.**
- "Increases by 0.66, 0.30 over the last two doublings": those are the steps 32 → 48 and 48 → 64 (factors 1.5 and 1.33).
- Per doubling, the increases are +1.76 (16 → 32) and +0.95 (32 → 64).
- **Fix:** Soften "saturating" accordingly.

**16. Minor. Section 7, "222 solves".** The report has 222 *records*: 74 cases × (GS, CNS, D). Only 148 of those are solves.

**17. Minor. Numbers quoted in the paper but absent from `report.txt`.** These cannot be verified from the raw report:
- the 10⁻⁵ code-to-code agreement;
- the "typically 10⁻¹⁵ and 10⁻¹¹" residuals;
- the discrete divergence 4.5·10⁻³ in Proposition 5.1;
- the whole lemma-test table (≈ 1.75, +0.05, −0.52, 1.57 / 1.56), which is also the only support for Lemma 1.5's "1.56–1.57";
- the floating-point λ\* values (≈ 15; 14.8–38.5);
- the Rayleigh quotients 9.97854 and 9.98584;
- the constraint counts 80 / 110.

My node counts do reproduce the kernel dimensions 14 and 17. Add the rest to the report or an appendix.

**18. Cosmetic. Table 1 rate column.** Each rate is taken against the previous run, whose rows (N = 24, 48) are omitted from the table. For example, the 1.69 at N = 64 is the 48 → 64 rate; for A, 32 → 64 gives 1.72. Say so in the caption.

**19. Minor. Lemma 1.5 sharpness.** It is tested on V_h^R, but the error equation only tests on Z_h. Sharpness for the method should be argued from Table 1 (CNS\* H¹ rate → 1.5), not from the V_h^R dual norm.

**20. Minor. Proposition 6.1 and Table 4(b).**
- Necessity is shown only for meshes that contain a star near 𝒮₃ (or 𝒮₄). The main-results summary and the "necessary and sufficient" line should say so.
- Table 4(b) varies ρ by coarsening the far field. This shows the threshold is set locally (μ/h ≳ c/ℓ_z), and that the ρ-scaling comes from Scott's use of the global h. The text's "boundary resolved more finely than the bulk" should say this.

**21. Cosmetic. Theorem D.** "‖u − u_h‖_{H¹(Ω_h)}" should be ‖ũ − u_h‖, since u is only defined on Ω.

## Checked and correct
- **Lemma 2.4:**
  - The carrier's curl is continuous across the cut: the Argyris traces jump by the constant 1, and the normal-derivative traces agree.
  - Its flux is ±1, and its size is C(1 + (μ/h)^{1/2})|F_g|.
  - The minors sin φ₁ sin φ₂ sin(φ₂ − φ₁), the requirement m ≥ 3 on straight edges, and the corner case m ≥ 2 are all right, as is the bound on δs.
- **Theorem D:** the consistency identity, the pressure estimate on V_h^0, the vanishing of the ⟨p̄_Γ(e_p), e_h·n⟩ and ξ̄ terms, the absorption with γ₂ ≍ max(γ₀, β⁻²), and the uniqueness argument through the edge bubble all hold.
- **Corollary B′:**
  - The relative errors are O((h/μ)^{1/2}/G + (γ^{1/2} + γ)h_Γ/G + ε_h(γ/h_Γ)^{1/2}/G), which tends to 0 under the stated hypotheses.
  - The constant 1 in both the slip norm and the energy norm follows.
- **Theorem B(ii)** is correct.
- **Report cross-check:** every entry of Tables 1–4, the numbers in the Table 2 and Table 3 text, G = a(7π/8)^{1/2}, the defect ≈ 0.105·h, the ρ and γ\* ranges, and μ\* = 60.343 all match the report.

---

## Response (all items addressed in `paper/` and `notes/proof.md`)

| # | Action |
|---|---|
| 1 | Theorem B regime restated as X ≤ δ·G·min(G, h_Γ^{1/2}), h_Γ/γ ≤ δG; steps 3–5 rewritten against it. |
| 2 | New test P (u = 0, f = ∇p): CNS* returns u_h = 0 to 1e-13, so GS output is exactly the traction response; Table 3 now uses P. |
| 3 | Out-of-regime μ = 10⁴ case stated as such, with the GS ratio ≈ 7.8. |
| 4 | 𝒮₄ replaced by nonsingular 𝒮₄′ (apex (0, 6/5)); exact certificate redone: kernel dim 16, λ ≥ 9 (Rayleigh 9.26747). |
| 5 | Flux sign fixed (−F_g); F_g = 0 for the test data noted; flux-corrected option noted. |
| 6 | Lemma 2.4 step 4 corrected (midpoint datum exact; carrier mismatches O(h⁶)). |
| 7 | Theorem C: h⁴ and (h/μ)h⁴h_Γ^{-1/2} terms kept and shown ≤ Cε_h. |
| 8, 9, 21 | Sign / notation fixes. |
| 10 | Corollary A′ hypothesis uses ε_h. |
| 11 | 𝒢(v) bound includes γh_Γh⁴. |
| 12 | Galdi cover quantified: (r_b − r_B)cos(α + r_B/r_b + h_Γ) > 1, overlap measures uniform; Guzmán–Scott needs shape regularity + (M2) only. |
| 13 | H¹ constant claim weakened: consistent with a common limit ≈ 2.7–2.8, not established. |
| 14 | Non-polynomial pressure test C added: CNS* differs from B₁ by 8e-7 → 1e-8 relative. "Every N" corrected to N ≤ 192. |
| 15 | γ* increments stated per doubling (3.4, 1.8, 0.95). |
| 16 | "148 solves (74 cases × 2 methods)". |
| 17 | All cited numbers now regenerated into `logs/` (lemma tests, exact certificates, float λ*) or `results/`. |
| 18 | Table captions say rates are between consecutive runs of the full sequence. |
| 19 | Lemma 1.5 sharpness for the method argued from the CNS* rates. |
| 20 | Necessity restricted to meshes with near-reference boundary stars; locality of the threshold stated. |
