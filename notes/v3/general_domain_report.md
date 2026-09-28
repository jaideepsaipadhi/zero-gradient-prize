# General curved domains: report

Companion to `notes/v3/general_domain.tex` (compiles standalone to 12 pp., no warnings). Result numbers such as "Lemma 3.2" refer to the compiled `paper/main.pdf` on `v3-draft`. The task's "Corollary 3.3" is **Corollary 5.3** (`cor:GSerr`) in the compiled paper. No existing file was edited, and nothing was committed.

## 1. Setting covered

- Ω ⊂ R² is bounded, open and **connected**, with ∂Ω = Γ ⊔ Σ.
- Γ = Γ₁ ⊔ … ⊔ Γ_m (m ≥ 1). Each Γᵢ is a whole boundary component and a closed embedded **C³** curve. Convex, nonconvex, inner and outer components are all allowed.
- **Σ ≠ ∅** is a union of polygonal components carrying strong Dirichlet data. It plays the role of ∂R.
- Γ_h is made of inscribed polygons with vertices on Γ. The curvature sign is arbitrary, so Ω_h \ Ω and Ω \ Ω_h can both be nonempty.
- Mesh assumptions (M0)–(M2) and (M1′), and data assumption (D), are unchanged, with ∂R → Σ.

## 2. New or changed assumptions (precise)

| Label | Content | Needed for |
|---|---|---|
| Ω connected, Σ ≠ ∅ | Korn and Friedrichs with zero trace on Σ; div V_h^0 onto mean-zero pressures; kernel argument | everything |
| Γ ∈ C³ | Tubular coordinates are C², κ ∈ C¹, and the O(ℓ²) remainders hold | Lemmas 3.x, 4.x, 5.x, Thm A (upper bound), Thm D |
| (G1) | Vertices are in cyclic order on each Γᵢ; the arc Γ_e of each chord has \|Γ_e\| ≤ 2\|e\| (no wrap-around chord) | Chord lemma, normal-graph lemma |
| (R1) | u ∈ C^{k+1}(Ω̄), p ∈ C^k(Ω̄) (the paper said "smooth") | as in the paper |
| (R2) | **Γ ∈ C^{k+3}** and p\|_Γ ∈ C^{k+1}(Γ) | Only the results that use the test field v_φ: lower bound of Thm A, Thms B and C, Cor B′. Needed so that w ∈ W^{k+1,∞}, which gives ‖v_φ − w‖_{W^{1,∞}} ≤ Ch^k. With Γ only C³, these results are **not proved** (the O(h^k) would become a weaker rate, which I have not tracked). |
| C depends on | ‖κ‖_{C¹}, tube radius r, graph scale a₀, the fixed star-shaped cover of Ω_* | |

Curvature is written ϰ in the note, because κ is the Korn constant of the paper.

## 3. What holds, section by section

**§3 Geometry**

- **Chord lemma (curvature).** The signed distance is d(x) = ½ϰ(x*)(ℓ²/4 − t²) + O(ℓ³), and |d| ≤ ‖ϰ‖∞ℓ²/8 + Cℓ⁴. The chord lies in Ω̄ where ϰ ≥ 0 and outside Ω where ϰ ≤ 0.
  - Normal rotation: n(x*) − n_e = ϰ_e t τ_e + O(ℓ²). For ϰ ≡ −1 this is the paper's −sτ_e.
  - The Jacobian of π: e → Γ_e is 1 + O(ℓ²), and the turning angle at a vertex is O(ĥ).
  - The proof works in a graph chart at the point where the tangent is parallel to the chord. It uses the midpoint divided-difference identity f′(ξ_m) = O(ℓ²), which holds uniformly even where ϰ → 0.
- **Normal graph.** Γ_{h,i} = {γ(s) + η_h(s)N(s)}, with ‖η_h‖∞ ≤ Cĥ² and ‖η_h′‖∞ ≤ Cĥ. The note defines Ω_h this way.
- **Lemma 3.2.**
  - Φ_h(Ψ(s,d)) = Ψ(s, d + χ(d)η_h(s)) in normal coordinates, with ‖DΦ_h − I‖∞ ≤ Cĥ.
  - The map moves points in both normal directions, so it is valid for nonconvex Γ.
  - It also shows Ω_h is connected, and gives a uniform Poincaré constant for mean-zero functions.
- **Extension.**
  - A global stream function need not exist (with several Σ components, g can carry flux through each). It is replaced by **local stream functions on the collars of the Γᵢ**. Each is single-valued because the flux through Γᵢ is 0. They are Whitney-extended across Γᵢ.
  - F is supported in the closure of Ω_h \ Ω, and |(F,v)| ≤ Cĥ²‖∇v‖ by integrating along normals.
  - Nothing needs extending into Ω \ Ω_h: all integrations by parts are on Ω_h.
  - The error on Ω_h \ Ω depends on the extension only at O(ĥ³) in H¹.
- **Lemma 3.5 (transposed-gradient term).** Holds with a_e = −ϰ_e n(∂_n u·τ_e).

**§4 Discrete tools**

- Lemmas 4.1 and 4.2 are mesh-only. Lemmas 4.3 and 4.6 hold verbatim with ∂R → Σ.
- **Lemma 4.5 (uniform inf-sup) is re-proved** without wedges and without P_h ⊂ D:
  - A self-contained decomposition lemma for the divergence equation over an ordered cover (explicit constants).
  - Bogovskiĭ's theorem on star-shaped pieces, with a constant depending only on diam/radius.
  - An h-uniform cover made of:
    - N₁ chart pieces {|ξ|<a, f_h(ξ)<η<a}, where f_h = I₁f is Lipschitz with constant ≤ 1/5 (slopes of the interpolant are values of f′). An explicit computation shows each piece is star-shaped with respect to B((0, 3a/4), a/8).
    - A fixed Galdi cover of the h-independent, connected Lipschitz set Ω_* = {dist(·,Γ) > a/4}.
  - The chart pieces overlap Ω_* in area ≥ a²/2.
  - Still inherited from the paper without re-proof: the claim that the Guzmán–Scott discrete constant depends on the domain only through the continuous inf-sup and Poincaré constants.

**§5 Error equation and fluxes**

- Lemmas 5.1 and 5.2 hold verbatim.
- **New flux lemma.** For v ∈ Z_h, only the *total* flux vanishes: Σᵢ ∫_{Γ_{h,i}} v·n = 0. For m ≥ 2 the flux vector maps Z_h *onto* the hyperplane Σxᵢ = 0.
- So **Corollary 5.3 allows exactly one global constant c**. Piecewise constants cᵢ add Σ(cᵢ − c)Φᵢ(v), which is nonzero.

**§6 Printed method GS**

- **G must be global**: G = ‖p − p̄‖_{L²(Γ)}, where p̄ = |Γ|⁻¹∫_Γ p is the mean over all curved components. It is not the per-component G_per.
  - The two are related by G² = G_per² + Σ|Γᵢ|(p̄ᵢ − p̄)².
  - Proposition (proved): with u = 0 and p ≡ cᵢ on Γᵢ (not all equal), G_per = 0, but |||u_h||| ≥ (2M)⁻¹(h/μ)^{1/2}G > 0. So no estimate in terms of G_per holds.
  - Physically, GS leaks net volume −(h/μ)|Γᵢ|(p̄ᵢ − p̄) through each obstacle.
- **New test field.**
  - The prescribed recipe φ = χ(d)∫(p − p̄)ds is single-valued only with the per-component mean. It works for m = 1 and is kept as the special case.
  - For the global mean, w must carry net flux through each Γᵢ, so it has no global stream function. It is built as follows:
    - an explicit divergence-free-by-construction field in normal coordinates, w₀ = (χq/J)N − χ′Q_i T;
    - a Bogovskiĭ correction supported in a fixed connected union of balls inside Ω, which removes the residual divergence χ′q̄ᵢ/J (it has total mean 0 because ∫_Γ q = 0);
    - the result satisfies w = (p − p̄)N exactly on Γ.
  - **New discrete-curl interpolant**: v|_T = curl Π_{k+1}(φ_T), where φ_T is a local stream function of w on a disk around T.
    - It is branch-independent and continuous.
    - It lies in Z_h and has the Morgan–Scott local error h^k.
    - It reduces to the paper's v_φ when a global stream function exists.
- Theorems A, B and C and Corollary B′ then hold with the global G and C_p = C(1 + ‖p − p̄‖_{W^{1,∞}(Γ)}). All geometric inputs are listed as explicit chord identities; the leading term of m(v) is b_e t with b_e = −q(m_e*)ϰ_e τ_e.

**§7 Consistent method CNS\***

- **One global mean suffices** for any number of components, both for well-posedness (full kernel analysis) and for the rate ĥ^{3/2} + h^k. The proof only ever uses the *total* flux of Z_h.
- A stronger fact is proved: for **any** normalisation ℓ with ℓ(1) = 1 (global boundary mean, domain mean, …), the method CNS\*_ℓ has the **same discrete velocity**, and pressures differing by a constant. The correction's only role is to kill the kernel vector (0,1) of Proposition 7.1.
- A **per-component** correction is *not* of this form and is **inconsistent**. Its residual is Σ(p̄ᵢ − p̄)Φᵢ(v), with dual norm ≥ c(h/μ)^{1/2}(Σ|Γᵢ|(p̄ᵢ − p̄)²)^{1/2}. It is still uniquely solvable for γ large.

**§8 Penalty threshold.** Proposition 8.1 holds verbatim (vertex angles are π + O(ĥ)).

## 4. Not proved, or inherited

1. **Γ merely C³ with the test-field results.** Theorems B and C, Corollary B′ and the Theorem A lower bound are proved only under (R2) (Γ ∈ C^{k+3}, p|_Γ ∈ C^{k+1}). With less smoothness, ‖v_φ − w‖_{W^{1,∞}} = O(h^{j}) for some j < k. That changes the admissible remainder terms, and I have not redone the bookkeeping.
2. **Error size of the per-component CNS variant.** The inconsistency is proved to be of order (h/μ)^{1/2}. I did **not** prove that the discrete error is actually that large: its pressure does not drop out on Z_h, so the Theorem A lower-bound argument does not apply as is.
3. **Σ = ∅** (the whole boundary curved, e.g. a curved channel with no straight part) is **not covered**. It breaks Korn/Friedrichs with zero trace on Σ and the V^R = "vanish on Σ" structure. A Nitsche-only version would need Korn with a boundary L² term, uniform in h. That is plausible but not done.
4. **Components mixing curved and polygonal pieces** (corners where Γ meets Σ) are not covered.
5. **Inherited from the paper without re-proof**, unchanged in status:
   - the Guzmán–Scott discrete inf-sup constant depends on the domain only via the continuous inf-sup and Poincaré constants;
   - the Morgan–Scott C¹ interpolant has element-local nodal degrees of freedom and local W^{2,∞} error h^k for all k+1 ≥ 5 under (M0), (M2);
   - Galdi's explicit Bogovskiĭ constant c₀(diam/R)²(1 + diam/R) (Thm III.3.1) and the W^{s,p} mapping property (Thm III.3.3).
6. **Standard facts about the fixed curve** are stated without proof: the tubular neighbourhood, regularity of d, σ and π, and uniform local graphs (Facts 1.2–1.3). They are h-independent and classical.
7. **No numerics** have been run for general domains. Suggested sanity checks:
   - an ellipse or a nonconvex "bean" obstacle;
   - two disks with piecewise-constant boundary pressure, i.e. the pure-pressure test of the G-global proposition, where G_per = 0 but GS must show rate-½ energy error.
8. The paper's sharpness remark after Lemma 3.5 (the ĥ^{3/2} rate is attained) was not re-examined. On straight parts of Γ (ϰ = 0) the leading term vanishes, so sharpness needs ϰ ≠ 0 somewhere.
