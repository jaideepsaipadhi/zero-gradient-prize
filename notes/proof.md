# Scott–Vogelius–Nitsche on a polygonally approximated boundary: convergence, the missing pressure traction, and the penalty threshold

**Setting.** Two-dimensional Stokes flow. The no-slip curve Γ is replaced by an inscribed polygon Γ_h. The velocity is continuous P4 and exactly divergence-free (Scott–Vogelius). Dirichlet conditions are imposed weakly on Γ_h by Nitsche's method, as in Gjerde–Scott (Adv. Comput. Sci. Eng. 2 (2024), eqs. (27)–(28)).

**Main results.**
- **Theorem A.** The method as printed (GS) converges, but in the natural mesh-dependent energy norm only at the exact rate (h/μ)^{1/2}·G, where G = ‖p − p̄‖_{L²(Γ)}.
- **Theorem C.** In the H¹ seminorm, GS converges at the sharp rate h/μ.
- **Theorem B.** The discrete velocity slips on Γ_h by ≍ (h/μ)G.
- **Corollary B′.** To leading order the GS error is exactly (h/μ)·v_φ, a field with boundary values −(h/μ)(p − p̄)n: the slip and energy errors equal (h/μ)G and (h/μ)^{1/2}G with constant 1 in the limit.
- **The cause** is a missing pressure traction. The naive fix is singular (**Proposition 5.1**); the corrected method CNS\* is well posed and attains Scott's expected rate h_Γ^{3/2} + h_Ω^k for general Stokes data (**Theorem D**).
- **Proposition 6.1.** The penalty threshold scales like ρ = h_Ω/h_Γ. This is sufficient, and necessary (by an exact-arithmetic certificate) for meshes whose boundary stars are near the certified reference stars.

---

## 0. Setting, methods, assumptions

**Domain.**
- R = (−L, L)² with L > 1; D is the open unit disk; Ω = R \ D̄; Γ = ∂D.
- P_h is a convex polygon with vertices on Γ, Γ_h = ∂P_h, and Ω_h = R \ P_h.
- Since P_h ⊂ D, Ω ⊂ Ω_h. The strip S_h = Ω_h \ Ω̄ lies between each chord and its arc.

**Mesh.** T_h is a conforming triangulation of Ω_h, and E_h^Γ is its set of edges on Γ_h. Mesh parameters:
- h = h_Ω := max_T diam T. This is the h in Scott's penalty μ/h.
- h_Γ := max_{e∈E_h^Γ} |e|.
- ρ := h / min_{e∈E_h^Γ} |e|.
- γ := μ h_Γ / h, the effective boundary penalty. Then μ/h = γ/h_Γ.

**Assumptions.**
- **(M0)** Shape regularity, and c₀h_Γ ≤ |e| ≤ h_Γ for all e ∈ E_h^Γ. Hence ρ/h ≤ 1/(c₀h_Γ). No quasi-uniformity is assumed: h_Γ ≪ h is allowed.
- **(M1)** No vertex of T_h is singular. Singular means an interior vertex whose incident edges lie on two lines, or a boundary vertex with one incident triangle, or one whose incident edges lie on two lines.
- **(M1′)** Every vertex of Γ_h has at least three incident triangles. This is needed because the angle of Ω_h at a vertex of Γ_h is π + O(h_Γ): with two triangles, Θ(z) = O(h_Γ) → 0, contradicting (M2) below.
- **(M2)** Θ(z) := max_j |sin(θ_j + θ_{j+1})| ≥ Θ₀ > 0 at every vertex (Guzmán–Scott; boundary vertices omit the wrap-around term).
- **(D)** f̃ is smooth and bounded on a fixed neighbourhood of Ω̄, with f̃ = f on Ω. h_Γ ≤ h\* is small.

**Spaces.**
- V_h: continuous piecewise-P4 vector fields.
- V_h^R := {v ∈ V_h : v|_∂R = 0}.
- V_h^0 := {v ∈ V_h : v|_{∂Ω_h} = 0}.
- Π_h: discontinuous piecewise P3.
- Z_h := {v ∈ V_h^R : div v = 0}.
- W_h := {z ∈ V_h : z|_∂R = g_I, div z = 0}, with g_I the P4 Lagrange interpolant of g.

**Forms.** n is the outward unit normal of Ω_h on Γ_h (pointing into P_h), and ⟨·,·⟩ := ⟨·,·⟩_{L²(Γ_h)}.

a_h(u, v) = ½(D u, D v), D u = ∇u + (∇u)ᵀ,

N_h(u, v) = a_h(u, v) − ⟨∂_n u, v⟩ − ⟨u, ∂_n v⟩ + (μ/h)⟨u, v⟩, with ∂_n u := (∇u)n.

**Method GS** (Gjerde–Scott (28), as printed). Find u_h ∈ V_h with u_h|_∂R = g_I, and p_h ∈ Π_h, such that

N_h(u_h, v) − (p_h, div v) = (f̃, v) ∀v ∈ V_h^R,  (q, div u_h) = 0 ∀q ∈ Π_h.

**Method CNS\*.** As GS, with momentum row

N_h(u_h, v) − (p_h, div v) + ⟨p_h − p̄_Γ(p_h), v·n⟩ = (f̃, v), where p̄_Γ(q) := |Γ_h|⁻¹ ∫_{Γ_h} q.

**Exact solution.** (u, p) is smooth with −Δu + ∇p = f, div u = 0 in Ω, u|_Γ = 0, u|_∂R = g. Here −div D(u) = −Δu. The stream function ψ with u = curl ψ = (∂_yψ, −∂_xψ) is single-valued, because u·n = 0 on Γ.

**Extensions.**
- ψ̃ is a smooth extension of ψ to a neighbourhood of Ω̄, and ũ := curl ψ̃, so div ũ = 0.
- p̃ is a smooth extension of p.
- F := f̃ + Δũ − ∇p̃ vanishes on Ω and is supported in S_h.

**Norms and constants.**
- |||v|||² := ‖∇v‖²_{Ω_h} + (μ/h)‖v‖²_{Γ_h} + Σ_{e∈E_h^Γ} |e| ‖∂_n v‖²_{L²(e)}. This is used for discrete fields and for errors ũ − z; ũ is smooth, so every trace exists.
- G := ‖p − p̄‖_{L²(Γ)}, with p̄ the mean over Γ.
- C, c depend on (M0)–(M2), k = 4, L and norms of (u, p, ψ̃). They are independent of h, h_Γ, μ, γ and ρ unless displayed. C_p denotes a constant that may also depend on ‖p − p̄‖_{W^{1,∞}(Γ)}. Where G is tracked explicitly, "δ small" means δ ≤ δ₀(C_p).

---

## 1. Geometry

**Lemma 1.1 (chords).** For a chord e of length ℓ and x ∈ e, dist(x, Γ) ≤ ℓ²/8 + O(ℓ⁴) (in particular ≤ ℓ²/4). Let s be the signed arclength from the chord midpoint m_e, x\* the nearest point on Γ, and τ_e the unit tangent of e. Then n(x\*) − n_e = −s τ_e + O(ℓ²), up to a fixed orientation sign. The arclength Jacobian of the radial projection Γ_h → Γ is 1 + O(ℓ²).

*Proof.* The sagitta is 1 − √(1 − ℓ²/4) = ℓ²/8 + ℓ⁴/128 + …. The normal at x\* is the unit radial vector. Its angle differs from that of n_e by the central angle to x\*, which is s + O(s³). ∎

**Lemma 1.2 (uniform trace, Friedrichs, Korn).** There are C_tr, C_F, κ, independent of h, such that for v ∈ H¹(Ω_h)² with v|_∂R = 0:

‖v‖_{L²(Γ_h)} ≤ C_tr ‖∇v‖, ‖v‖ ≤ C_F ‖∇v‖, ‖∇v‖² ≤ κ a_h(v, v).

*Proof.*
1. Write Γ_h in polar form as r = ρ_h(θ). On a chord with endpoint angles θ_m ± α, ρ_h = cos α / cos(θ − θ_m). Hence ‖ρ_h − 1‖_∞ ≤ Ch_Γ² and ‖ρ_h′‖_∞ ≤ Ch_Γ.
2. With χ smooth, χ(1) = 1, χ ≡ 0 outside a fixed collar in Ω, define Φ_h(r, θ) = (r + χ(r)(ρ_h(θ) − 1), θ).
3. Φ_h maps Ω onto Ω_h and Γ onto Γ_h, is the identity near ∂R, and is bi-Lipschitz with ‖DΦ_h − I‖_∞ ≤ Ch_Γ almost everywhere.
4. For w = v∘Φ_h: D(w) = D(v)∘Φ_h + E with |E| ≤ Ch_Γ |∇v∘Φ_h|.
5. On the fixed domain Ω, w vanishes on ∂R, a set of positive measure, so the first Korn, trace and Friedrichs inequalities hold.
6. Changing variables back gives ‖∇v‖ ≤ C(‖D(v)‖ + h_Γ‖∇v‖). Absorb the last term for h_Γ ≤ h\*. The other two inequalities transfer the same way. ∎

**Lemma 1.3 (extension).** ‖ũ‖_{L^∞(Γ_h)} ≤ Ch_Γ², so ‖ũ‖_{L²(Γ_h)} ≤ Ch_Γ². For v ∈ V_h^R, |(F, v)| ≤ Ch_Γ² ‖∇v‖.

*Proof.* The first bound is Lemma 1.1 plus ũ|_Γ = 0. For the second, use a 1D Poincaré inequality across strips of width ≤ Ch_Γ², the bounds |S_h| ≤ Ch_Γ² and ‖F‖_∞ ≤ C (from (D)), and Lemma 1.2. ∎

**Lemma 1.4 (wall structure).** On Γ: ∇u = ∂_n u ⊗ n, ∂_n u · n = 0, and (∇u)ᵀn = 0.

*Proof.* The tangential derivatives vanish, and div u = n·∂_n u. ∎

**Lemma 1.5 (transposed-gradient term).** For v ∈ V_h^R:

|⟨(∇ũ)ᵀn, v⟩| ≤ C Σ_e ( |e|² ‖v‖_{L¹(e)} + |e|² ‖∂_τ v‖_{L¹(e)} ) ≤ C h_Γ^{3/2} ‖∇v‖.

The rate h_Γ^{3/2} is sharp on V_h^R: a first-layer field with odd trace on each chord attains it, and the dual-norm test of Section 7 gives 1.56–1.57. For the method itself (tested on Z_h), sharpness is supported by the CNS\* H¹ rates of Table 1, which decrease monotonically toward 3/2.

*Proof.*
1. Taylor-expand about x\*, then apply Lemmas 1.1 and 1.4:

   (∇ũ(x))ᵀn_e = n(x\*)(∂_n u(x\*) · (n_e − n(x\*))) + O(h_Γ²) = s·a_e + O(h_Γ²),

   where a_e := n(m_e\*)(∂_n u(m_e\*) · τ_e) is constant on e. The O(h_Γ²) remainder collects |x − x\*|, the O(ℓ²) term of Lemma 1.1, and s·O(s).
2. Since ∫_e s ds = 0: |∫_e s a_e·v| = |∫_e s a_e·(v − v(m_e))| ≤ (|e|²/4)‖∂_τ v‖_{L¹(e)}.
3. Use ‖·‖_{L¹(e)} ≤ |e|^{1/2}‖·‖_{L²(e)} and the inverse trace |e|‖∂_τ v‖²_{L²(e)} ≤ C_I²‖∇v‖²_{T_e}.
4. Summing, Σ_e |e|^{5/2}‖∂_τ v‖_{L²(e)} ≤ C h_Γ² (#E_h^Γ)^{1/2} ‖∇v‖ ≤ C h_Γ^{3/2}‖∇v‖. Finish with Lemma 1.2. ∎

---

## 2. Discrete tools

**Lemma 2.1 (inverse traces).**
- For v ∈ V_h: Σ_e |e|‖∂_n v‖²_{L²(e)} ≤ C_I²‖∇v‖², hence ‖∂_n v‖_{Γ_h} ≤ C_I(ρ/h)^{1/2}‖∇v‖.
- For q ∈ Π_h: ‖q‖²_{Γ_h} ≤ C_I² Σ_{T∈T_h^Γ} |e_T|⁻¹‖q‖²_T ≤ C_I²(ρ/h)‖q‖².

(A triangle with two edges on Γ_h only doubles C_I.) ∎

**Lemma 2.2 (coercivity, continuity).** Let μ₀ := 4κρC_I². A sufficient condition for μ ≥ μ₀ is γ ≥ γ₀ := 4κC_I²/c₀. For μ ≥ μ₀:
- N_h(v, v) ≥ c₁|||v|||² for all v ∈ V_h^R, with c₁ = 1/max(2κ(1 + C_I²), 2).
- |N_h(w, v)| ≤ M|||w||| |||v|||, with M = 3 + (κC_I²)^{-1/2}.

Both constants are independent of h, μ and ρ.

*Proof.*
1. By Lemma 2.1 and Young's inequality with ε = 1/(2κ):

   2|⟨∂_n v, v⟩| ≤ 2C_I‖∇v‖(ρ/h)^{1/2}‖v‖_{Γ_h} ≤ ½a_h(v, v) + 2κρC_I² h⁻¹‖v‖²_{Γ_h}.

2. Hence N_h(v, v) ≥ ½a_h(v, v) + (μ − 2κρC_I²)h⁻¹‖v‖², which is at least (2κ)⁻¹‖∇v‖² + (μ/2h)‖v‖² for μ ≥ μ₀.
3. Add Σ|e|‖∂_n v‖² ≤ C_I²‖∇v‖² to get the coercivity constant.
4. Continuity: a_h(w, v) ≤ 2‖∇w‖‖∇v‖, and each ∂_n-term is ≤ (ρ/μ)^{1/2}|||w||| |||v||| ≤ (4κC_I²)^{-1/2}|||w||| |||v|||. ∎

**Lemma 2.3 (divergence range).** Under (M1):
- div V_h^R = Π_h;
- W_h ≠ ∅;
- the velocities of GS and CNS\* are pointwise divergence-free.

*Proof.*
1. With no singular vertices, Scott–Vogelius (1985; k ≥ 4) gives div V_h^0 = Π_h ∩ L²_0. That result is a local construction and applies to the doubly connected Ω_h.
2. The edge bubble b = β_e n_e on one e ⊂ Γ_h lies in V_h^R and has ∫ div b = ∫_{Γ_h} b·n = |e|/6 ≠ 0. This adds the constants, so div V_h^R = Π_h.
3. For any 𝒢 ∈ V_h with trace g_I, div 𝒢 ∈ Π_h = div V_h^R. So 𝒢 − v ∈ W_h for a suitable v ∈ V_h^R.
4. For the discrete velocity, u_h − 𝒢 ∈ V_h^R, and the continuity row forces div u_h = 0.

*Remark.* Suppose a vertex of ∂R is singular, e.g. a corner with one triangle. Then A_z(div 𝒢) = 0 is a condition on g_I that generic data violate, and W_h = ∅. ∎

**Lemma 2.4 (approximation in W_h).**

E_h := inf_{z∈W_h} |||ũ − z||| ≤ Cε_h,  ε_h := (1 + (γh_Γ)^{1/2}) h⁴ + μ^{1/2}h^{5.5}.

The last term is a genuine lower-order floor, not an artefact: every z ∈ W_h has ∫_{Γ_h} z·n = −∫_{∂R} g_I·ν =: −F_g (ν the outward normal of R), while ∫_{Γ_h} ũ·n = 0, so E_h ≥ (μ/h)^{1/2}|F_g|/|Γ_h|^{1/2}. F_g = O(h⁶) is the composite Boole-rule error, and |F_g| ≍ h⁶ for generic g. The term is ≤ h⁴ whenever μ ≤ h⁻³ (e.g. bounded γ and h_Γ ≥ h⁴). It vanishes if g_I is replaced by a flux-corrected interpolant with ∫_{∂R} g_I·ν = 0, and it vanishes for the test data of Section 7, where g·ν is a polynomial of degree ≤ 5 on each side of ∂R and Boole's rule is exact.

*Proof.*

1. *Interior and Γ_h.* Let s₀ = Π_A ψ̃ be the Argyris interpolant (C¹, piecewise P5, C² at vertices), with |ψ̃ − s₀|_{W^{m,∞}(T)} ≤ Ch_T^{6−m}. Then curl s₀ ∈ V_h is divergence-free. Its bulk error is O(h⁴). On Γ_h the error terms are (μ/h)^{1/2}·Ch_Γ⁵ = Cγ^{1/2}h_Γ^{4.5} and Ch_Γ⁴.
2. *C¹ gluing fact.* A piecewise P5 function is C¹ iff, on each interior edge, the trace and the normal-derivative trace agree from both sides. The trace (P5) is Hermite in (s, ∂_τ s, ∂_τ² s) at both endpoints. The normal-derivative trace (P4) is determined by (∂_ν s, ∂_τ∂_ν s) at both endpoints plus the midpoint value; the kernel t²(1 − t)² is nonzero there.

   Consequently each triangle may carry its own vertex Hessian H_T, provided H_{T₊} − H_{T₋} = α_e ν_e⊗ν_e across each edge. A symmetric jump with τᵀJτ = τᵀJν = 0 is exactly of this form.
3. *Matching g_I on ∂R.* The target trace is ∂_τ s = g̃·n and ∂_n s = −g̃·τ, with g̃ := g_I − F_g·curl Π_A(θ/2π)|_∂R. Here F_g := ∫_∂R g_I·n is the net flux of the interpolated data (step 5). The carrier's stream function jumps by exactly 1 across the cut, so its flux through ∂R is exactly F_g and ∫_∂R g̃·n = 0; g̃ is therefore the trace of a single-valued stream function.

   At each ∂R vertex this prescribes one-sided τᵀHτ and τᵀHn.
   - At a straight-edge vertex with m triangles: two equations in the m − 1 coefficients α_j, with coefficient vectors sin φ_j (sin φ_j, −cos φ_j). The 2×2 minors are sin φ₁ sin φ₂ sin(φ₂ − φ₁), so rank 2 holds iff m ≥ 3.
   - At a corner: τ₁ᵀHn₁ and τ₂ᵀHn₂ are the same Hessian entry, which forces compatibility if m = 1. For m ≥ 2 the system is solvable.
   - Both cases are exactly the nonsingular ones, so (M1) suffices. The solutions are bounded in terms of the minimum angle (M0).
4. *Sizes.*
   - For g_I itself, Lagrange data are exact at vertices, so values and gradients need no correction, and the edge-midpoint normal derivative matches exactly (the midpoint is both a P4 Lagrange node and an Argyris degree of freedom). Subtracting the carrier adds mismatches of size O(|F_g|) = O(h⁶) in every datum, which are lower order.
   - The Hessian mismatch is the derivative of the P4 interpolation error, O(h⁴).
   - Vertex stream-function values are cumulative edge fluxes. Equispaced P4 Lagrange interpolation integrates exactly like the 5-point Newton–Cotes (Boole) rule, so each per-edge flux error is O(h⁷) and cumulative errors are O(h⁶). Against the value basis (|D²| ~ h⁻²) this contributes O(h⁴).

   The correction δs therefore has ‖∇ curl δs‖_∞ ≤ Ch⁴ on the one-element layer at ∂R (area O(h)), so ‖∇ curl δs‖ ≤ Ch^{4.5}.
5. *Net flux.* F_g = ∫_∂R (g_I − g)·n = O(h⁶), because ∫_∂R g·n = 0 and composite Boole integration of a smooth function has error O(h⁶). The carrier F_g·curl Π_A(θ/2π) is defined on the domain cut along a path of mesh edges. θ is smooth on Ω_h because 0 ∈ P_h. Argyris uses the true vertex Hessian on both branches, so the traces jump by the constant F_g and the normal-derivative traces agree. Its curl is therefore continuous P4 and divergence-free, with flux F_g through Γ_h. Its size is |||F_g curl Π_A(θ/2π)||| ≤ C(1 + (μ/h)^{1/2})|F_g| ≤ C(h⁶ + μ^{1/2}h^{5.5}).

   Set z := curl(s₀ + δs) + F_g curl Π_A(θ/2π). Then z|_∂R = g_I and div z = 0, so z ∈ W_h. ∎

---

## 3. The error equation

**Lemma 3.1.** For v ∈ V_h^R:

N_h(ũ, v) − (p̃, div v) + ⟨p̃, v·n⟩ = (f̃, v) + 𝒢(v),

with 𝒢(v) := ⟨(∇ũ)ᵀn, v⟩ − ⟨ũ, ∂_n v⟩ + (μ/h)⟨ũ, v⟩ − (F, v).

*Proof.* The steps, in order:
- Symmetry of D(ũ): a_h(ũ, v) = (D(ũ), ∇v) = −(div D(ũ), v) + ⟨D(ũ)n, v⟩, since v = 0 on ∂R.
- div D(ũ) = Δũ.
- −Δũ = f̃ − F − ∇p̃, and −(∇p̃, v) = (p̃, div v) − ⟨p̃, v·n⟩.
- D(ũ)n − ∂_n ũ = (∇ũ)ᵀn. ∎

**Lemma 3.2.**
- For v ∈ V_h^R: |𝒢(v)| ≤ C(1 + γ^{1/2}) h_Γ^{3/2} |||v|||.
- For v ∈ V_h^0: |𝒢(v)| ≤ C h_Γ^{3/2} ‖∇v‖.

*Proof.*
- Transposed term: ≤ C h_Γ^{3/2}‖∇v‖ (Lemma 1.5).
- |⟨ũ, ∂_n v⟩| ≤ Ch_Γ² · C_I(ρ/h)^{1/2}‖∇v‖ ≤ C h_Γ^{3/2}‖∇v‖.
- Penalty term: ≤ (μ/h)^{1/2}Ch_Γ²|||v||| = Cγ^{1/2}h_Γ^{3/2}|||v|||.
- F term: Lemma 1.3.
- For v ∈ V_h^0 the transposed and penalty terms vanish. ∎

**Corollary 3.3 (GS on Z_h).** For v ∈ Z_h, (·, div v) = 0 and ∫_{Γ_h} v·n = 0. Hence for every constant c:

N_h(ũ − u_h, v) = ℛ(v) + 𝒢(v), ℛ(v) := −⟨p̃ − c, v·n⟩.

---

## 4. The printed method GS (γ ≥ γ₀ throughout)

**Test field.**
- Let φ₀(θ) := ∫₀^θ (p − p̄)(θ′)dθ′, which is periodic, and φ := χ₁(r)φ₀(θ), with χ₁ ≡ 1 near r = 1 (so χ₁′(1) = 0) and χ₁ ≡ 0 near ∂R.
- w := curl φ satisfies w_r = ∂_θφ/r and w_θ = −∂_rφ. On Γ: w = (p − p̄)e_r = −(p − p̄)n.
- v_φ := curl Π_A φ ∈ Z_h, and ‖v_φ − w‖_{W^{1,∞}} ≤ Ch⁴.

**Theorem A (mesh-dependent energy norm: exact rate 1/2).**

Upper bound:
|||ũ − u_h||| ≤ C[(h_Γ/γ)^{1/2}G + (1 + γ^{1/2})h_Γ^{3/2} + ε_h].

Lower bound, if G > 0 and h_Γ ≤ h₀(G), with h₀ independent of γ ≥ γ₀:
|||ũ − u_h||| ≥ (2M)⁻¹(h_Γ/γ)^{1/2}G − C[(1 + γ^{1/2})h_Γ^{3/2} + ε_h].

(h_Γ/γ)^{1/2} = (h/μ)^{1/2} in Scott's variables.

*Proof of the upper bound.* For z ∈ W_h, e_h := z − u_h ∈ Z_h. By Corollary 3.3 and Lemma 2.2,

c₁|||e_h|||² ≤ N_h(z − ũ, e_h) + ℛ(e_h) + 𝒢(e_h).

Bound |ℛ(e_h)| ≤ (G + Ch_Γ²)(h/μ)^{1/2}|||e_h|||, using the Jacobian bound of Lemma 1.1. Then apply Lemmas 3.2 and 2.4.

*Proof of the lower bound.*
1. ℛ(v_φ) = ⟨(p̃ − p̄)(p − p̄)(x\*)·(n(x\*)·n_e)⟩ + O(h⁴) = G² + O(h_Γ² + h⁴), since v_φ·n_e = −(p − p̄)(x\*)n(x\*)·n_e + O(h⁴), ℛ = −⟨p̃ − p̄, v·n⟩ and n(x\*)·n_e = 1 + O(s²).
2. |||v_φ|||² ≤ (μ/h)(G² + Ch_Γ²) + C.
3. By continuity, M|||ũ − u_h||| ≥ (|ℛ(v_φ)| − |𝒢(v_φ)|)/|||v_φ|||, and Lemma 3.2 bounds 𝒢.
4. |ℛ|/|||v_φ||| ≥ (h/μ)^{1/2}G(1 − Ch_Γ²/G² − C(h_Γ/γ₀)^{1/2}/G). This is ≥ ½ for h_Γ ≤ h₀(G). ∎

The lower bound is carried by the boundary parts of |||·|||. By Theorem C, ‖∇(ũ − u_h)‖ is only O(h/μ).

**Corollary A′.** If G > 0, γh_Γ ≪ G and ε_h ≪ (h_Γ/γ)^{1/2}G, then |||ũ − u_h||| ≍ (h/μ)^{1/2}. The mesh-dependent energy rate is exactly 1/2.

**Lemma 4.1 (normal-load cancellation).** Let S be edgewise W^{1,∞} on Γ_h, continuous at the vertices of Γ_h, with |S·τ_e| ≤ C h_Γ on every e ∈ E_h^Γ. (For instance, S smooth with S·τ = 0 on Γ.) Then for v ∈ Z_h:

|⟨S, ∂_n v⟩| ≤ C(‖v‖_{Γ_h} + h_Γ^{1/2}‖∇v‖) ≤ C‖∇v‖.

*Proof.*
1. Tangential part: |⟨(S·τ)τ, ∂_n v⟩| ≤ Ch_Γ‖∂_n v‖_{Γ_h} ≤ Ch_Γ^{1/2}‖∇v‖ (Lemma 2.1 and (M0)).
2. Normal part: on e, v is polynomial on T_e with div v = 0, so ∂_n v·n_e = −∂_τ(v·τ_e). With σ_e := S·n_e,

   ∫_e σ_e ∂_n v·n_e = ∫_e ∂_τσ_e (v·τ_e) − [σ_e v·τ_e]_{∂e}.

3. At a polygon vertex x shared by e and e′, the endpoint terms combine to v(x)·(σ_{e′}τ_{e′} − σ_eτ_e), of size ≤ C(|e| + |e′|)|v(x)|.
4. Σ_e |e|‖v‖_{L^∞(e)} ≤ C Σ_e |e|^{1/2}‖v‖_{L²(e)} ≤ C|Γ_h|^{1/2}‖v‖_{Γ_h}.
5. Finish with Lemma 1.2. ∎

Numerically (Section 7), the dual norm of v ↦ ⟨S, ∂_n v⟩ behaves as follows:
- On Z_h with normal S: it stays bounded.
- On V_h^R (no divergence constraint), or with tangential S: it grows like h_Γ^{-1/2}.

Incompressibility is exactly what the proof uses.

**Theorem C (H¹ seminorm: rate 1, sharp).**

‖∇(ũ − u_h)‖ ≤ C_p (h_Γ/γ) + C[(1 + γ^{1/2})h_Γ^{3/2} + ε_h],

with C_p := C(1 + ‖p − p̄‖_{W^{1,∞}(Γ)}).

*Proof.*
1. *Split.* For z ∈ W_h: z − u_h = δ_R + δ_G. Here δ_R, δ_G ∈ Z_h solve N_h(δ_R, v) = ℛ(v) and N_h(δ_G, v) = 𝒢(v) − N_h(ũ − z, v) on Z_h. Lax–Milgram on Z_h (Lemma 2.2) gives existence and uniqueness.
2. *Geometric part.* ‖∇δ_G‖ ≤ |||δ_G||| ≤ c₁⁻¹(C(1 + γ^{1/2})h_Γ^{3/2} + ME_h).
3. *Ansatz.* Heuristically, δ_R ≈ −(h/μ)(p − p̄)n = (h/μ)w on Γ_h. Set S_h := v_φ ∈ Z_h and r_h := δ_R − (h/μ)S_h ∈ Z_h. Then

   N_h(r_h, v) = m(v) + ℓ(v),
   - m(v) := ℛ(v) − ⟨S_h, v⟩ = −⟨(p̃ − p̄)n + S_h, v⟩,
   - ℓ(v) := −(h/μ)[a_h(S_h, v) − ⟨∂_n S_h, v⟩ − ⟨S_h, ∂_n v⟩].
4. *Bound on m.* On e: (p̃ − p̄)n_e + S_h = (p − p̄)(x\*)(n_e − n(x\*)) + O(h_Γ²) + O(h⁴) = (p − p̄)(m_e\*) s τ_e + O(h_Γ²), which is odd in s with a constant coefficient. The argument of Lemma 1.5 handles the odd part; the O(h⁴) remainder contributes Ch⁴‖v‖_{L¹(Γ_h)}. Hence |m(v)| ≤ C_p(h_Γ^{3/2} + h⁴)‖∇v‖.
5. *Bound on ℓ.* |a_h(S_h, v)| + |⟨∂_n S_h, v⟩| ≤ C_p‖∇v‖. Also |⟨S_h, ∂_n v⟩| ≤ |⟨w, ∂_n v⟩| + |⟨S_h − w, ∂_n v⟩| ≤ C_p‖∇v‖ + Ch⁴·h_Γ^{-1/2}‖∇v‖, by Lemma 4.1 applied to S = w. So |ℓ(v)| ≤ (h/μ)(C_p + Ch⁴h_Γ^{-1/2})‖∇v‖.
6. *Close.* c₁|||r_h|||² ≤ [C_p(h/μ + h_Γ^{3/2} + h⁴) + C(h/μ)h⁴h_Γ^{-1/2}]‖∇r_h‖. The last two terms are ≤ Cε_h, because (h/μ)h⁴h_Γ^{-1/2} = γ^{-1}h⁴h_Γ^{1/2} ≤ Ch⁴. Hence ‖∇δ_R‖ ≤ (h/μ)‖∇S_h‖ + ‖∇r_h‖ ≤ C_p(h_Γ/γ) + C(h_Γ^{3/2} + ε_h). ∎

Without incompressibility, Lemma 4.1 would give only C h_Γ^{-1/2}, and the H¹ error would be (h/μ)h_Γ^{-1/2}: rate 1/2, with element-scale oscillation in the trace. The inconsistency −p n is purely normal, and div v = 0 converts ∂_n v·n into a tangential derivative. That is why GS still converges at rate 1 in H¹.

**Theorem B (boundary slip; H¹ lower bound).** Let γ ≥ γ₀, G > 0, and put X := C_p(h_Γ/γ) + (1 + γ^{1/2})h_Γ^{3/2} + ε_h (the right side of Theorem C). Suppose

X ≤ δ·G·min(G, h_Γ^{1/2}),  h_Γ/γ ≤ δG,  h_Γ ≤ h₀(G),

for a fixed small δ. (All three hold for fixed γ and bounded ρ as h_Γ → 0.) Then:
- (i) c(h_Γ/γ)G ≤ ‖ũ − u_h‖_{L²(Γ_h)} ≤ C(h_Γ/γ)^{1/2}·|||ũ − u_h|||;
- (ii) ‖∇(ũ − u_h)‖ ≥ c(h_Γ/γ)G − Ch^{4.5}.

Combined with Theorem C: ‖∇(ũ − u_h)‖ ≍ h/μ in this regime.

*Proof of (i).* Let e := ũ − u_h and v := v_φ. Corollary 3.3 rearranges to

(μ/h)⟨e, v⟩ = ℛ(v) + 𝒢(v) − a_h(e, v) + ⟨∂_n e, v⟩ + ⟨e, ∂_n v⟩.

Terms, using that v is smooth up to O(h⁴):
1. ℛ(v) = G² + O(h_Γ² + h⁴).
2. |𝒢(v)| ≤ C(h_Γ^{3/2} + γh_Γ³ + γh_Γh⁴). On the chord, ũ = ±d·∂_n u(x\*) + O(d²) is tangential while v ≈ −(p − p̄)n is normal, so ũ·v = O(h_Γ⁴ + h_Γ²h⁴), the second term from v_φ − w.
3. |a_h(e, v)| ≤ C_p‖∇e‖ ≤ C_pX ≤ C_pδG², by Theorem C.
4. |⟨e, ∂_n v⟩| ≤ C_p‖e‖_{Γ_h} ≤ C_p(h/μ)^{1/2}|||e||| ≤ C_p[(h/μ)G + X] ≤ 2C_pδG², by Theorem A and h/μ = h_Γ/γ ≤ δG.
5. |⟨∂_n e, v⟩| ≤ ‖∂_n(ũ − z)‖_{Γ_h}‖v‖ + C_I(ρ/h)^{1/2}(‖∇(z − ũ)‖ + ‖∇e‖)‖v‖_{Γ_h}. By Theorem C and (ρ/h)^{1/2} ≤ (c₀h_Γ)^{-1/2}, this is ≤ C h_Γ^{-1/2}·X·G ≤ CδG² in the regime.

For δ small, |⟨e, v⟩| ≥ ½(h/μ)G². Dividing by ‖v‖_{Γ_h} ≤ G + Ch_Γ² gives the lower bound. The upper bound is ‖e‖_{Γ_h} ≤ (h/μ)^{1/2}|||e|||.

*Proof of (ii).* Let E ∈ H¹(Ω_h) lift (g − g_I)|_∂R, supported near ∂R, with E|_{Γ_h} = 0 and ‖E‖_{H¹} ≤ C‖g − g_I‖_{H^{1/2}(∂R)} ≤ Ch^{4.5}. Then ‖e‖_{Γ_h} = ‖e − E‖_{Γ_h} ≤ C_tr(‖∇e‖ + ‖∇E‖) by Lemma 1.2. ∎

**Corollary B′ (exact leading term: why Scott's rate fails for GS).** Let G > 0 and let h_Γ → 0 with h/μ → 0, γh_Γ → 0 and ε_h/(h_Γ/γ)^{1/2} → 0. Then

ũ − u_h = (h/μ)·v_φ + ϑ_h,  |||ϑ_h||| = o((h/μ)^{1/2}),  ‖ϑ_h‖_{Γ_h} = o(h/μ),

and consequently

‖ũ − u_h‖_{L²(Γ_h)} / ((h/μ)G) → 1,  |||ũ − u_h||| / ((h/μ)^{1/2}G) → 1.

On Γ_h, v_φ ≈ −(p − p̄)n. So to leading order the printed method replaces the no-slip condition by the leak condition (μ/h)u_h·n ≈ (p − p̄): the penalty, not a traction, balances the pressure.

*Proof.* Write ũ − u_h = (ũ − z) + δ_R + δ_G with the notation of Theorem C, and δ_R = (h/μ)S_h + r_h with S_h = v_φ.
1. From step 6 of Theorem C, c₁|||r_h|||² ≤ [C_p(h/μ + h_Γ^{3/2}) + Cε_h]|||r_h|||, so |||r_h||| ≤ C_p(h/μ + h_Γ^{3/2}) + Cε_h and ‖r_h‖_{Γ_h} ≤ (h/μ)^{1/2}|||r_h|||.
2. |||δ_G||| + |||ũ − z||| ≤ C[(1 + γ^{1/2})h_Γ^{3/2} + ε_h], and the Γ_h norm carries the extra factor (h/μ)^{1/2}.
3. Relative to (h/μ)^{1/2}G (energy) and (h/μ)G (slip), steps 1–2 contribute O((h/μ)^{1/2}/G + (γ^{1/2} + γ)h_Γ/G + ε_h(γ/h_Γ)^{1/2}/G) → 0.
4. ‖v_φ‖_{Γ_h} = G(1 + O(h_Γ²)) + O(h⁴) (Lemma 1.1: Jacobian and displacement of Γ_h from Γ), and |||v_φ|||² = (μ/h)‖v_φ‖²_{Γ_h} + O(1). ∎

The corollary gives the sharp answer to Scott's question for GS: the h_Γ^{3/2} + h^k estimate fails exactly when p is not constant on Γ, and the error is then the explicit field (h/μ)v_φ with the universal constant 1 in both the slip and the energy norm. The H¹ constant ‖∇(ũ − u_h)‖μ/(hG) is not universal: it depends on the global Stokes response to the leak, as r_h has the same order as (h/μ)S_h in H¹. In the numerics it is ≈ 2.7–2.8 (Table 3, test P).

---

## 5. The pressure-consistent method

**Proposition 5.1 (the naive consistent form is singular).** Replace ⟨p_h − p̄_Γ(p_h), v·n⟩ by ⟨p_h, v·n⟩. The resulting linear system is square: V_h^R × Π_h tested by V_h^R × Π_h. It has the kernel vector (0, 1) and is therefore singular.

*Proof.* −(1, div v) + ⟨1, v·n⟩ = −∫_{∂Ω_h} v·n + ∫_{Γ_h} v·n = 0 for v ∈ V_h^R. The continuity row vanishes at u = 0. ∎

The system is solvable only for right-hand sides orthogonal to a left null vector, which the argument does not identify. In computations with generic data (Section 7), both an iterated-penalty solver and a bordered direct solve returned discrete velocities that were not divergence-free (discrete divergence ≈ 4.5·10⁻³), consistent with incompatibility.

**Theorem D (CNS\*: well-posed, expected rate).** Let β be the Guzmán–Scott inf-sup constant for (V_h^0, Π_h ∩ L²_0) (Math. Comp. 88 (2019); k ≥ 4, depending on Θ₀, k and shape regularity).

β is uniform in h. Its domain dependence enters only through the continuous (Bogovskii) constant. Guzmán–Scott build the discrete Fortin operator locally from vertex/edge patches using shape regularity and (M2) only (no quasi-uniformity), and combine it with the continuous inf-sup constant of the domain. For the domains Ω_h, that constant is uniform. Cover Ω_h by the h-dependent pieces W_j ∩ Ω_h, where W_j are fixed angular wedges of half-angle α, together with fixed pieces near ∂R. Each W_j ∩ Ω_h is star-shaped with respect to a ball B(c_j, r_B) with |c_j| = r_b, provided every chord line meeting W_j leaves the ball on its outer side. This holds uniformly when (r_b − r_B)cos(α + r_B/r_b + h_Γ) > 1, which can be arranged with fixed α, r_b, r_B since P_h ⊂ D. The pieces have diameters, ball radii and pairwise overlap measures bounded above and below independently of h. Galdi's Lemma III.3.4 then gives a Bogovskii constant, hence a continuous inf-sup constant, independent of h.

There is γ₂ ≍ max(γ₀, β⁻²) such that for γ ≥ γ₂, CNS\* has a unique solution and

|||ũ − u_h||| ≤ C[(1 + γ^{1/2})h_Γ^{3/2} + ε_h].

Consequently ‖ũ − u_h‖_{H¹(Ω_h)} ≤ C(h_Γ^{3/2} + h_Ω^4) for γ in a bounded range [γ₂, γ_max], provided h_Γ ≥ h_Ω⁴ (so that the flux floor μ^{1/2}h^{5.5} ≤ C h⁴; without it that single extra term is kept). This is the rate Scott conjectured, for general Stokes data. The L² part of the H¹ norm uses the lift E of Theorem B(ii).

*Proof.*

1. *Consistency.* Write B\*(q, v) := −(q, div v) + ⟨q − p̄_Γ(q), v·n⟩.
   - For constants c: B\*(q + c, v) = B\*(q, v) − c∫_{Γ_h} v·n.
   - B\*(p̃, v) equals the pressure form of Lemma 3.1 minus p̄_Γ(p̃)∫_{Γ_h} v·n.
   - Hence p\* := p̃ − p̄_Γ(p̃) satisfies N_h(ũ, v) + B\*(p\*, v) = (f̃, v) + 𝒢(v) for all v ∈ V_h^R.
   - With e_u := ũ − u_h and e_p := p\* − p_h: N_h(e_u, v) + B\*(e_p, v) = 𝒢(v). (5.1)
2. *Split.* e_p = η + ξ with η := p\* − π_h p\* and ξ ∈ Π_h. For v ∈ V_h^0, (η, div v) = 0 because div V_h^0 ⊂ Π_h.
3. *Pressure.* For v ∈ V_h^0, (5.1) gives

   (ξ − ξ̄, div v) = N_h(e_u, v) − 𝒢(v) = a_h(e_u, v) − ⟨e_u, ∂_n v⟩ − 𝒢(v).

   The right side is ≤ C(|||e_u||| + h_Γ^{3/2})‖∇v‖, using (ρ/μ)^{1/2} ≤ C and Lemma 3.2. Guzmán–Scott gives β‖ξ − ξ̄‖ ≤ C(|||e_u||| + h_Γ^{3/2}).
4. *Velocity.* For z ∈ W_h, e_h := z − u_h ∈ Z_h. Test (5.1) with v = e_h. Then (e_p, div e_h) = 0 and ⟨p̄_Γ(e_p), e_h·n⟩ = 0, so

   c₁|||e_h|||² ≤ ME_h|||e_h||| + |𝒢(e_h)| + |⟨η, e_h·n⟩| + |⟨ξ − ξ̄, e_h·n⟩|.

   - |⟨η, e_h·n⟩| ≤ Ch_Γ⁴(h/μ)^{1/2}|||e_h|||.
   - |⟨ξ − ξ̄, e_h·n⟩| ≤ C_I(ρ/h)^{1/2}‖ξ − ξ̄‖(h/μ)^{1/2}|||e_h||| ≤ C(c₀γ)^{-1/2}‖ξ − ξ̄‖|||e_h|||.
5. *Absorb.* Insert step 3 and |||e_u||| ≤ E_h + |||e_h|||. For γ ≥ γ₂ with C(c₀γ₂)^{-1/2}β⁻¹ ≤ c₁/2, the pressure term is absorbed.
6. *Uniqueness.* Apply steps 3–5 to the homogeneous problem (g_I = 0, so z = 0 ∈ W_h). This gives u_h = 0, then p_h = c constant. Then (5.1) reads −c∫_{Γ_h} v·n = 0 for all v; taking the edge bubble of Lemma 2.3 gives c = 0. The system is square, so uniqueness gives existence. ∎

---

## 6. The penalty threshold scales with ρ, and this is sharp

**Proposition 6.1.**
- (a) *Sufficiency.* μ ≥ μ₀ = 4κC_I²ρ implies coercivity (Lemma 2.2).
- (b) *Necessity, with an exact certificate.* Let z be a vertex of a straight part of the Nitsche boundary. Suppose its star, dilated by the factor ℓ_z (the length of its boundary edges), is congruent to one of the reference stars 𝒮₃ or 𝒮₄′ below. If μ/(h/ℓ_z) < 9, then N_h is indefinite on Z_h. The same holds, by continuity, at a vertex of Γ_h with interior angle π + O(h_Γ) whose star is a small perturbation of 𝒮₃ or 𝒮₄′. This uses that both reference stars are nonsingular (every vertex, including rim vertices, has non-collinear incident edges), so the dimension of the divergence-free subspace equals the generic count and is locally constant; it fails at the singular m = 2 configuration, which (M1′) excludes.

So for meshes containing such a boundary star at a vertex with ρ ≍ h/ℓ_z, μ₀ ≍ ρ is both necessary and sufficient. The threshold is local: what is required is μ/h ≳ c/ℓ_z, and the ρ-scaling is the price of Scott's normalisation by the global h.

*Proof of (b).*
1. In 2D, a(v, v) and ⟨∂_n v, v⟩ are dilation-invariant, and ‖v‖²_{L²(Γ_h)} scales linearly.
2. A witness v on the reference star, extended by zero, lies in Z_h: it is continuous, vanishes on the star's outer chords and on ∂R, and is divergence-free.
3. It satisfies N_h(v, v) = a − 2⟨∂_n v, v⟩ + (μℓ_z/h)‖v‖²_ref.

*Exact certificate.* The reference stars are fans at O = (0, 0) with rational rim points:
- 𝒮₃: (1,0), (1,1), (−1,1), (−1,0);
- 𝒮₄′: (1,0), (1,1), (0,6/5), (−1,1), (−1,0). (The apex is moved off the line y = 1: with apex (0,1) the rim vertex is singular, the kernel dimension jumps from 16 to 17, and the perturbation argument does not apply.)

The Nitsche boundary is {y = 0}, and the rim chords carry zero Dirichlet data. In exact rational arithmetic:
- Continuity, the rim condition and pointwise div v = 0 are imposed as polynomial identities: 80 and 110 constraints, giving kernel dimensions 14 and 16, the generic counts.
- The quadratic forms A = Σ_T ∫ ½Dv:Dv, B = ∫(∂_n v)·v and C = ∫|v|² are integrated exactly.
- A rational witness in the kernel satisfies 2B − A − 9C > 0 **exactly** on each star. Its exact Rayleigh quotients are 9.97854… (𝒮₃) and 9.26747… (𝒮₄′).

Floating-point values on other fans: λ\* ≈ 15 (uniform, m = 3–4) and 14.8–38.5 (random 4-fans). A global check: the measured threshold γ\* = μ\*/ρ stays within 18.1–20.1 as ρ varies over 3.3–44 on two mesh levels, and is 14.7–20.8 under refinement at fixed ρ, saturating (Section 7). ∎

---

## 7. Numerical verification

**Code and meshes.** P4 Scott–Vogelius with Nitsche on the inscribed polygon and strong Dirichlet data on the square R = (−2.5, 2.5)². Divergence is enforced exactly (all ten P3 moments per triangle). Saddle systems are solved by SuperLU with partial pivoting. Over the 148 production solves (74 cases × two methods), the relative residual is ≤ 3.8·10⁻¹⁵ and the discrete divergence ≤ 7·10⁻¹¹ (median 9·10⁻¹²). The code reproduces the published zero-pressure-gradient benchmark (Cavalcante) to 10⁻⁵ relative on four meshes, and gives μ\* = 60.343 at N = 20 against the published 60.338. The mesh family has N boundary chords, h_Ω = 1.505 → 0.109 for N = 16 → 256, and ρ ≈ 3.3, so rates in h_Ω and h_Γ coincide. Rates in the tables are between consecutive runs of the full sequence N = 16, 24, 32, 48, 64, 96, 128, 192, 256, some of which are omitted from the tables.

**Test problems.**
- **A** (Scott's benchmark): shear flow, p ≡ 0, so G = 0.
- **B_a**: ψ = (r² − 1)²(x + y²/2)/10 and p = a(x²y − y + x/2), with G = a(7π/8)^{1/2} (G = 1.658 for a = 1, 165.8 for a = 100).
- **P**: pure pressure, u ≡ 0, f = ∇p, with p as in B₁. CNS\* reproduces this exactly (u_h = 0, p_h = p − p̄_Γ(p), since p ∈ Π_h); the computed CNS\* velocity is ≤ 10⁻¹³ in H¹. So the GS output is *exactly* the response to the missing traction.
- **C**: ψ as in B, non-polynomial pressure p = e^{x/2}cos y. This exercises the term p\* − π_h p\* of Theorem D.

A and B were run with μ ∈ {100, 1000, 10000} and N ≤ 192, plus N = 256 for (A, 100) and (B₁, 100). P was run with all three μ and C with μ = 100, both for N ≤ 96.

**Table 1. CNS\* attains Scott's rate h_Γ^{3/2} (Theorem D), independently of the pressure.** H¹ error of CNS\*, μ = 100.

| N | h_Ω | A | rate | B₁ | rate |
|---|---|---|---|---|---|
| 32 | 0.813 | 6.41e−2 | 1.80 | 2.34e−2 | 2.43 |
| 64 | 0.422 | 2.08e−2 | 1.69 | 6.83e−3 | 1.76 |
| 96 | 0.285 | 1.09e−2 | 1.64 | 3.54e−3 | 1.68 |
| 128 | 0.215 | 6.96e−3 | 1.60 | 2.24e−3 | 1.63 |
| 192 | 0.145 | 3.71e−3 | 1.58 | 1.19e−3 | 1.60 |
| 256 | 0.109 | 2.39e−3 | 1.56 | 7.58e−4 | 1.57 |

At N = 256 the energy-norm rate is 1.52 and the slip rate 2.02 (the slip carries the extra (h/μ)^{1/2}). The rates decrease monotonically toward 3/2, the geometric term of Lemma 1.5. The h_Ω⁴ term is invisible at these levels because the geometric term dominates.

Pressure independence has two parts:
- For polynomial pressures it is structural (p ∈ Π_h), and B₁ and B₁₀₀ agree to all printed digits wherever both were run (N ≤ 192, every μ).
- With the non-polynomial pressure of test C, the CNS\* H¹ error differs from B₁'s by a relative 7.9·10⁻⁷, 3.9·10⁻⁸ and 9.7·10⁻⁹ at N = 16, 64 and 96. This is the size of the ⟨p\* − π_h p\*, v·n⟩ contribution, which is of higher order as Theorem D predicts.

**Table 2. GS loses Scott's rate as soon as G > 0 (Theorems A–C).** H¹ rate between N = 128 and N = 192 (B₁₀₀).

| μ | GS | CNS\* |
|---|---|---|
| 100 | 0.99 | 1.60 |
| 1000 | 1.00 | 1.56 |
| 10000 | 1.06 | 1.58 |

For the non-polynomial test C (μ = 100), the GS H¹ rate from N = 64 to 96 is already 1.12 against 1.68 for CNS\*. For the smaller pressure B₁ at μ = 1000, GS still shows 1.51 at N = 192: the (h/μ)G term dominates only once (h/μ)G ≳ h_Γ^{3/2}, which always happens eventually because h_Γ^{3/2} decays faster. Larger μ or smaller G delays the crossover but does not remove it.

**Table 3. The leading term of Corollary B′ (test P: exact traction response).** Normalised GS error for test P (G = 1.658).

| N | μ | ‖∇e‖·μ/(hG) | slip·μ/(hG) | energy·(μ/h)^{1/2}/G |
|---|---|---|---|---|
| 16 | 100 | 2.464 | 0.894 | 1.051 |
| 32 | 100 | 2.612 | 0.945 | 1.007 |
| 64 | 100 | 2.684 | 0.972 | 0.996 |
| 96 | 100 | 2.708 | 0.981 | 0.996 |
| 96 | 1000 | 2.770 | 0.998 | 0.999 |
| 96 | 10000 | 2.790 | 0.9993 | 0.9994 |

Findings from Table 3:
- **Slip and energy ratios.** Both converge to 1, as Corollary B′ predicts. The defect in the slip ratio decreases as h/μ decreases: at N = 96 it is 0.019, 0.0025 and 0.0007 for μ = 10², 10³ and 10⁴.
- **H¹ constant.** It is not covered by the corollary. The μ = 100 sequence rises (2.46 → 2.71), while the μ = 10³ and 10⁴ sequences fall (2.78 → 2.77 and 3.06 → 2.79). This is consistent with a common, μ-independent limit near 2.7–2.8, but the data do not establish one.
- **Test B.** GS itself at B₁, μ = 100, N = 256 gives slip and energy ratios 0.995 and 1.001.
- **Out of regime.** For μ = 10⁴ with the small pressure B₁, γh_Γ ≈ 74 ≫ G at N = 192, so the corollary's hypothesis fails. The GS energy ratio there is ≈ 7.8, dominated by the geometric error, not by the traction response.
- **The difference D.** D := u_GS − u_CNS\* is not a clean proxy for the traction response on B: it also contains the ψ-part difference, which is O(h_Γ^{3/2}) and nonzero even for test A. Test P removes this.

**Table 4. Penalty study (Proposition 6.1).** The coercivity threshold μ\* on Z_h is computed from the inertia of A_μ + rBᵀM⁻¹B by bisection.

*(a) Refinement at fixed ρ.*

| N | 8 | 12 | 16 | 20 | 24 | 32 | 48 | 64 |
|---|---|---|---|---|---|---|---|---|
| ρ | 3.41 | 3.02 | 3.28 | 3.17 | 3.30 | 3.33 | 3.37 | 3.40 |
| μ\* | 50.3 | 51.8 | 59.3 | 60.3 | 63.3 | 66.0 | 69.1 | 70.7 |
| γ\* = μ\*/ρ | 14.7 | 17.1 | 18.1 | 19.1 | 19.2 | 19.9 | 20.5 | 20.8 |

Per doubling of N, γ\* increases by 3.4 (8 → 16), 1.8 (16 → 32) and 0.95 (32 → 64). The increments roughly halve, which suggests (but does not prove) a limit near 22.

*(b) Varying ρ via an enlarged outer box.*

| L | 2.5 | 5 | 10 | 20 | 40 |
|---|---|---|---|---|---|
| N = 16: ρ | 3.28 | 5.90 | 11.3 | 22.2 | 43.9 |
| N = 16: μ\* | 59.3 | 112.1 | 208.0 | 411.5 | 809.3 |
| N = 16: γ\* | 18.1 | 19.0 | 18.4 | 18.5 | 18.4 |
| N = 32: ρ | 3.33 | 5.79 | 10.8 | 21.1 | 41.5 |
| N = 32: μ\* | 66.0 | 113.9 | 218.0 | 422.2 | 835.7 |
| N = 32: γ\* | 19.9 | 19.7 | 20.1 | 20.0 | 20.1 |

Here ρ is varied by coarsening the far field while the boundary layer is fixed. μ\* is linear in ρ over more than a decade, with a constant slope. The threshold is thus set locally at Γ_h (μ/h ≳ c/h_Γ), and the ρ-scaling is the cost of Scott's normalisation by the global h. This is the computational penalty study the prize asks for:
- a safe choice is μ ≥ 25ρ, i.e. γ ≥ 25;
- Scott's fixed μ fails once the far field is coarsened beyond about ρ ≈ μ/20.

The exact certificate λ ≥ 9 on 𝒮₃ and 𝒮₄′ is a local lower bound, consistent with the global γ\* ≈ 18–21.

**Lemma tests (exact discrete dual norms; lemma_tests.py, N = 16–64).**

| Functional | Space | Rate in h_Γ |
|---|---|---|
| ⟨S, ∂_n v⟩, normal S | Z_h | bounded (≈ 1.75, rate +0.05) |
| ⟨S, ∂_n v⟩, normal S | V_h^R | h_Γ^{-0.52} |
| ⟨S, ∂_n v⟩, tangential S | Z_h | h_Γ^{-0.52} |
| ⟨(∇ũ)ᵀn, v⟩ | V_h^R | h_Γ^{1.57}, h_Γ^{1.56} (Lemma 1.5: 3/2) |

**Proposition 5.1.** The naive consistent form is singular, and generic data are incompatible: both an iterated-penalty solver and a bordered direct solve return a discrete divergence of 4.5·10⁻³.

**Reproducibility.** The scripts behind each number:
- svn.py: solver.
- job.py and run_all.py: Tables 1–4.
- report.py: tabulation.
- penalty.py: μ\*.
- lemma_tests.py: dual norms.
- lamstar_exact.py and lamstar_exact_v2.py: the exact certificates for 𝒮₃ and 𝒮₄′, with constraint counts, kernel dimensions and Rayleigh quotients printed.
- lamstar.py: the floating-point λ\* on other fans.

---

*External inputs, all classical:* Scott–Vogelius (1985); Guzmán–Scott (2019); Galdi, *An Introduction to the Mathematical Theory of the Navier–Stokes Equations*, Lemma III.3.4; Argyris interpolation; the first Korn inequality.
