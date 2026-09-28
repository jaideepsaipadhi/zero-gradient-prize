# Blind verification of the SVN claims (P1, P2, L, A, C, B′, D, PR, E)

Referee: independent adversarial pass. Part I was written **before** opening
`paper/sections/*.tex`, `notes/v4/rescue.tex` and `notes/v4/pressure_drag.tex`.
Part II (comparison) was written after. Numerical checks used `code/svn.py`
unmodified; the scripts are in the session scratchpad and are reproduced inline in
condensed form below.

---

## Part I: own derivations. Written 2026-09-28T14:35Z, before reading any proof

### 0. Conventions and two identities used everywhere

* `a_h(u,v) = ½(Du,Dv) = (Du,∇v) = (∇u,∇v) + (∇uᵀ,∇v)`, because `Du` is symmetric.
* **Consistency identity.** For the extended exact solution (div ũ = 0 on Ω_h), `div Dũ = Δũ`, and every
  v ∈ V_h^R (v = 0 on ∂R), with n the outward normal of Ω_h on Γ_h:

  `N_h(ũ,v) − (p̃,div v) + ⟨p̃, v·n⟩ = (f̃,v) + ⟨(∇ũ)ᵀn, v⟩ + (μ/h)⟨ũ,v⟩ − ⟨ũ,∂_n v⟩.`   (0.1)

  (∫∂Ω_h = ∫Γ_h since v|∂R = 0.) The transposed-gradient term exists because N_h uses the full
  gradient `∂_n u = (∇u)n` and not the traction `(Du)n`.
* **Constant pressure.** For c ∈ ℝ and v ∈ V_h^R, `(c, div v) = c∫_{Γ_h} v·n`. So in GS a constant
  pressure shift turns the missing load `⟨p̃n,v⟩` into `⟨(p̃−c)n,v⟩`. The optimal c is the Γ_h mean,
  so G appears.
* **Geometry.** Chord of length ℓ ≤ h_Γ in the unit circle has sagitta ≈ ℓ²/8, so
  dist(Γ_h,Γ) ≤ h_Γ²/8. The chord normal n differs from the radial normal ν by an angle **O(h_Γ)**
  (zero at the midpoint, ±θ/2 at the endpoints, odd about the midpoint).
* On Γ, ũ = 0 and div ũ = 0 give `∇ũ = a⊗ν` with a = ωτ ⊥ ν. Hence `(∇ũ)ᵀν = ν(a·ν) = 0` on Γ.
  On Γ_h, however, `(∇ũ)ᵀn = ν(ωτ·(n−ν)) + O(h_Γ²) = O(h_Γ)` **pointwise**, with mean zero on each
  edge to leading order. It is **not** O(h_Γ²).
* No triangle has two edges on Γ_h: the Ω_h-angle at a vertex of the convex polygon exceeds π.
* **Div-free trace identity.** On a straight edge with (n,τ) constant, if div v = 0 pointwise on the
  adjacent triangle then `n·(∂_n v) = ∂_n(v·n) = −∂_τ(v·τ)`. This holds for discretely div-free SV
  fields: div V_h ⊂ Π_h, so they are pointwise div-free.

Consistency sizes in the dual of |||·||| (for any v ∈ V_h^R):

| term | pointwise size on Γ_h | bound by |
|---|---|---|
| `(μ/h)⟨ũ,v⟩` | ũ = O(h_Γ²) | `(μ/h)^{½}h_Γ²·|||v||| = γ^{½}h_Γ^{3/2}|||v|||` |
| `⟨ũ,∂_n v⟩` | O(h_Γ²) | `h_Γ²(min|e|)^{-½}|||v||| = C h_Γ^{3/2}|||v|||` |
| `⟨(∇ũ)ᵀn, v⟩` | **O(h_Γ)** | `h_Γ‖v‖_{Γ_h} ≤ h_Γ(h/μ)^{½}|||v||| = γ^{-½}h_Γ^{3/2}|||v|||` |
| `⟨(p̃−p̄)n, v⟩` (GS only) | O(1) | `G_h (h/μ)^{½}|||v|||` |

In H¹-duality, with no penalty available (‖v‖_{Γ_h} ≤ C‖v‖_{H¹}), the transposed-gradient term is
only O(h_Γ)‖v‖_{H¹}. The O(h_Γ^{3/2}) rate therefore **uses the penalty** (factor γ^{-½}), or else
the edgewise odd symmetry of τ·n plus one more integration by parts.

Note on γ: μ ≥ μ₀ ∝ ρ = h/min|e| and h_Γ/min|e| ≥ 1 give γ = μh_Γ/h ≥ const·h_Γ/min|e| ≥ const.
So γ is automatically bounded below, and γ^{-½} factors are harmless for any fixed μ₀.

### P1. Naive system: square and singular. Its kernel is exactly span(0,1), and it is generically unsolvable

* Square: the unknowns are (u ∈ g_I + V_h^R, p ∈ Π_h) and the equations are indexed by
  (v ∈ V_h^R, q ∈ Π_h).
* (0,1) is in the kernel: `−(1,div v) + ⟨1,v·n⟩ = −∫_{Γ_h}v·n + ∫_{Γ_h}v·n = 0`.
* Integrating by parts elementwise, the naive pressure operator is `b̃(p,v) = (∇_h p, v) − Σ_{E int}⟨[[p]], v·n_E⟩`.
  The Γ_h terms cancel, so it is a discrete gradient that ignores Γ_h. Its kernel on V_h^R equals the
  kernel of `(p,div v)` on V_h^0 plus the constants. Under GS (no singular vertices, k ≥ 4) that kernel
  is exactly the constants.
* Kernel of the full matrix: the system is non-symmetric (the pressure operator in the momentum row
  differs from the constraint `(q,div u)`). Testing with v = u gives `N_h(u,u) = −⟨p,u·n⟩`, which
  does not close. So "kernel = span(0,1)" is not a two-line consequence of coercivity. Numerically
  (below) it is exactly one-dimensional.
* **Missed consequence.** Square with a one-dimensional kernel means a one-dimensional cokernel. The left
  null vector (v*,q*) has v* ≠ 0, because div: V_h^R → Π_h is onto. The consistency condition
  `F(v*) + (continuity data)(q*) = 0` is **not** automatic. So the naive system is not merely
  non-unique: for generic (f,g) it has **no solution**. The two readings are also directly linked:
  the P2 constant c (below) is exactly the residual of the dropped q = 1 row.

Numerics (svn.py mesh, kind B, μ = 20, dense SVD):

| N | size | σ_min/σ_max | next σ/σ_max | right null vector | P2: c |
|---|---|---|---|---|---|
| 8 | 832 | 4e-16 | ≥3e-9 (cluster, scaling) | ‖u‖ = 3e-13, p ∝ 1 (err 1e-9) | −1.49e-3 |
| 12 | 1872 | 6e-20 | ≥6.5e-10 | ‖u‖ = 3e-13, p ∝ 1 (err 2e-9) | +7.62e-3 |

**Verdict P1:** TRUE, with two additions. (i) The claim "kernel is exactly span" needs an argument
(non-symmetric system), or it should be stated as "contains, and numerically equals". (ii) The
cokernel means generic insolvability; that is the stronger and more useful statement.

### P2. p ∈ Π_h ∩ L²_0(Ω_h): div u_h ≡ c

* Constants lie in the kernel of the naive pressure operator, so the momentum row is the same.
  `(q,div u_h) = 0 ∀q ∈ Π_h∩L²_0` and div u_h ∈ Π_h give div u_h ≡ c.
* `c|Ω_h| = ∫_{∂R} g_I·ν + ∫_{Γ_h} u_h·n`. This was checked numerically to 1e-12.
* **c vanishes if and only if the naive full system is solvable**, that is, iff the discrete leak
  flux matches: `∫_{Γ_h}u_h·n = −∫_{∂R}g_I·ν`. Explicitly, c·∫_{Ω_h}q* = N_h(g_ext,v*) + (q*,div g_ext) − (f̃,v*).
  This is a nonzero linear functional of the data.
* Sufficient conditions for c = 0:
  * the naive method reproduces the data exactly, e.g. u ≡ 0 with f = ∇p and p ∈ Π_h. The naive method
    is pressure-consistent, so (0,p) solves it.
  * trivial data.
* In general c ≠ 0. Numerically c = −1.5e-3 at N = 8 and 7.6e-3 at N = 12, with
  ‖Bu − c·w‖/‖Bu‖ ~ 1e-11. That confirms div u_h is constant and nonzero.

**Verdict P2:** TRUE. Precise vanishing condition: c = 0 iff ∫_{Γ_h}u_h·n + ∫_{∂R}g_I·ν = 0,
iff the naive (unconstrained) system is consistent.

### L. Coercivity

* `N_h(v,v) = ½‖Dv‖² − 2⟨∂_n v,v⟩ + (μ/h)‖v‖²_{Γ_h}`.
* Each T has at most one Γ_h edge, and the inverse trace estimate gives `|e|‖∇v‖²_e ≤ C_I²‖∇v‖²_T`. So
  `2|⟨∂_n v,v⟩| ≤ 2C_I‖∇v‖(ρ/h)^{½}‖v‖_{Γ_h} ≤ ε‖∇v‖² + (C_I²ρ/(εh))‖v‖²`.
* Korn in the form `‖∇v‖² ≤ κ a_h(v,v)` on V_h^R, uniform over Ω_h. With ε = 1/(2κ):
  `N_h(v,v) ≥ (2κ)^{-1}‖∇v‖² + ((μ−2κC_I²ρ)/h)‖v‖²`.
* μ ≥ 4κC_I²ρ gives ≥ (2κ)^{-1}‖∇v‖² + (μ/2h)‖v‖². Adding `Σ|e|‖∂_n v‖² ≤ C_I²‖∇v‖²` gives
  `N_h(v,v) ≥ c|||v|||²` with `c = min(1/(2κ),½)/(1+C_I²)`.
* Hypotheses needed:
  * Korn constant κ uniform over the domain family Ω_h, for fields vanishing on ∂R only. This holds
    because ∂R is fixed and Ω_h is uniformly Lipschitz, but it must be stated.
  * The normalisation of κ: with `‖∇v‖² ≤ κ‖ε(v)‖²` instead, the threshold constant changes by a
    factor of 2.
  * ρ uses min|e| globally. The local version μ ≥ C h/|e| per edge would suffice.

**Verdict L:** TRUE, as a sufficient condition.

### A. GS energy error, upper and lower bounds

Error equation. Subtract GS from (0.1), shifting the pressure by p̄. For every v ∈ V_h^R:

`N_h(e,v) − (π,div v) = R(v) := ⟨(∇ũ)ᵀn − (p̃−p̄)n, v⟩ + (μ/h)⟨ũ,v⟩ − ⟨ũ,∂_n v⟩`.

**Upper bound.**
1. Take w_h: a div-free approximation of ũ with w_h|_{∂R} = g_I. Use Argyris / C¹-P_{k+1} interpolation
   of the stream function (k+1 ≥ 5), which gives exactly div-free fields with no inf-sup needed, plus
   flux adjustment on Γ_h.
2. Test with e_h = w_h − u_h, which is div-free, so the pressure drops out.
3. With the table in §0:
   `|||ũ−u_h||| ≤ C[(h/μ)^{½}G + (1+γ^{½}+γ^{-½})h_Γ^{3/2} + μ^{½}h^k]`. The last term is really the
   local mesh size to the power k near Γ_h.
4. The γ^{-½} factor is absorbed because γ is bounded below (§0).

**Lower bound.**
1. Build a smooth div-free Z = curl(Φ(θ)η(r)), where Φ' = −(p−p̄) on Γ (periodic because the mean is
   zero), η is a fixed cutoff with η(1) = 1 and η'(1) = 0, and Z = 0 on ∂R.
2. Then on Γ_h, Z·n = (p̃−p̄) + O(h_Γ²) (sign convention aside) and Z·τ = O(h_Γ).
3. Take v_h = curl(I_Argyris(Φη)) ∈ V_h^R, which is exactly div-free, so the pressure again drops out.
4. `|||v_h||| ≤ C(μ/h)^{½}G` for μ/h ≥ 1. Then
   `G² ≤ C|||e|||(μ/h)^{½}G + C(h_Γ + γh_Γ)G + h.o.t.`, which gives
   `|||e||| ≥ c(h/μ)^{½}G − C(γ^{½}+γ^{-½})h_Γ^{3/2} − C h^k-terms`.
5. This is informative only while γh_Γ → 0, i.e. (h/μ)^{½}G dominates γ^{½}h_Γ^{3/2}.

**Verdict A:** TRUE, including the lower bound. The γ-dependence is γ^{½} (from the penalty) and
γ^{-½} (from (∇ũ)ᵀn). The latter is harmless only because γ ≥ const, which follows from μ ≥ μ₀ρ.

### C. GS H¹ error at rate 1

1. Subtract the leading layer. Let z = (h/μ)Z with Z as in A, oriented so that (μ/h)⟨z,v⟩ cancels
   −⟨(p̃−p̄)n,v⟩ up to O(h_Γ) pointwise.
2. Set ẽ = ũ − z − u_h and ẽ_h = I_F(ũ−z) − u_h, which is div-free.
3. The new residual is R(v) − N_h(z,v). Its terms:
   * `a_h(z,v) ≤ C(h/μ)|||v|||`.
   * `⟨∂_n z,v⟩ ≤ C(h/μ)^{3/2}|||v|||`.
   * `⟨Z+(p̃−p̄)n, v⟩`: normal part O(h_Γ²), tangential part Z·τ = O(h_Γ), total ≤ Cγ^{-½}h_Γ^{3/2}|||v|||.
   * **`⟨z,∂_n v⟩ = (h/μ)⟨Z,∂_n v⟩`** is the dangerous term, because ‖∂_n v‖_{Γ_h} ≤ h_Γ^{-½}|||v||| only.
     * Split it into the normal part `⟨Z·n, n·∂_n v⟩ = −⟨Z·n, ∂_τ(v·τ)⟩` (the div-free identity, which
       uses the pointwise div-free property of the SV test function v = ẽ_h) and the tangential part.
     * Integrate by parts edge by edge: `⟨∂_τ(Z·n), v·τ⟩ + Σ_vertices O(h_Γ)|Z||v(x_i)|`.
     * The vertex sum is at most C‖v‖_{Γ_h}, using `|v(x_i)| ≤ C|e|^{-½}‖v‖_e`. So the normal part is
       at most C(h/μ)‖v‖_{Γ_h}.
     * The tangential part `⟨Z·τ, τ·∂_n v⟩` is at most h_Γ·h_Γ^{-½}|||v|||, which is o(1).
4. Hence `|||ẽ_h||| ≤ C[h/μ + (1+γ^{½}+γ^{-½})h_Γ^{3/2} + h^k-terms]`, and
   `‖∇(ũ−u_h)‖ ≤ ‖∇z‖ + ‖∇ẽ‖ ≤ C_p h/μ + …`, with C_p ≈ ‖Z‖_{H¹} ∝ ‖p−p̄‖_{H^{1/2+}(Γ)}.

Numerics (kind P: u = 0, f = ∇p, so e = −u_h is exactly the traction response):

| N | μ | γ | ‖∇e‖/(h/μ) |
|---|---|---|---|
| 16 | 40 | 12.2 | 3.93 |
| 16 | 160 | 48.8 | 4.23 |
| 32 | 40 | 12.0 | 4.05 |
| 32 | 160 | 48.1 | 4.42 |
| 64 | 40 | 11.8 | 4.30 |

The ratio is stable, so the rate is 1.

**Verdict C:** TRUE. The mechanism is correct: the identity is what removes the h_Γ^{-½} loss, and it
needs pointwise div-free test functions, which SV provides. The vertex terms in the edgewise
integration by parts must be accounted for, since n and τ jump at vertices. They are O(h_Γ) per
vertex, and their sum is bounded by ‖v‖_{Γ_h}.

### B′. Leading term and the ratios going to 1

1. From C: `‖ẽ‖_{Γ_h} ≤ (h/μ)^{½}|||ẽ_h||| + interp ≤ C(h/μ)^{½}[h/μ + (1+γ^{½})h_Γ^{3/2}] + …`, and
   `‖z‖_{Γ_h} = (h/μ)(G + O(h_Γ))`.
2. Slip ratio − 1 = O((h/μ)^{½} + γ^{½}h_Γ + γh_Γ + h_Γ + μh^{k−½})/G.
3. Energy ratio:
   * `(μ/h)‖e‖²_Γ / ((h/μ)G²)` → 1.
   * `‖∇e‖²/((h/μ)G²) = O(h/μ)`.
   * `Σ|e|‖∂_n e‖²/((h/μ)G²) = O(h/μ + γ²h_Γ²)`.
4. Hence both ratios → 1 provided G > 0, h/μ → 0, γh_Γ → 0, and μh^{k−½} → 0. The last is
   automatic for μ ≲ h^{-1}; for fixed γ all of these hold.

Numerics, same runs (slip ratio, then energy ratio as defined in svn.errors, which uses h·‖∂_n‖²):

| N | μ | slip ratio | energy ratio |
|---|---|---|---|
| 16 | 40 | 0.8135 | 1.619 |
| 32 | 40 | 0.8807 | 1.038 |
| 64 | 40 | 0.9348 | 1.012 |
| 16 | 160 | 0.9231 | 1.024 |
| 32 | 160 | 0.9626 | 1.002 |

The slip deficit roughly halves per refinement, consistent with O(h/μ). The claimed v_φ ≈ −(p−p̄)n
sign also agrees with the heuristic `(μ/h)e|_{Γ_h} ≈ −(p̃−p̄)n`.

**Verdict B′:** TRUE under the stated limits, plus G > 0 and h^k negligibility. Convergence of the
slip ratio is slow, O((h/μ)^{½}) in the proof and ~O(h/μ) observed.

### D. CNS*: well-posedness and error

1. Constants are **not** in the kernel of `b*(p,v) = −(p,div v) + ⟨p−p̄(p), v·n⟩`: b*(1,v) = −∫_{Γ_h}v·n.
   So p_h is unique and approximates `p* = p̃ − p̄_Γ(p̃)`. Using (0.1), p* is consistent, with no
   pressure load left: the residual is only the geometric R_geo.
2. **Well-posedness.**
   1. Uniform inf-sup on (V_h^0, Π_h∩L²_0) gives
      `‖p−m(p)‖ ≤ Cβ^{-1}(‖∇u‖ + h_Γ^{-½}‖u‖_{Γ_h}) ≤ Cβ^{-1}(1+γ^{-½})|||u|||`.
      The `⟨u,∂_n v⟩` term survives for v ∈ V_h^0.
   2. With the discrete trace `‖q‖_{Γ_h} ≤ C h_Γ^{-½}‖q‖`, the non-symmetric term satisfies
      `|⟨p−p̄,u·n⟩| ≤ Cβ^{-1}γ^{-½}(1+γ^{-½})|||u|||²`.
   3. So u = 0 once `Cβ^{-1}γ^{-½}(1+γ^{-½}) < c_coerc`. That defines γ₂.
   4. Then p = m and m∫v·n = 0 for all v, so m = 0.
   5. γ₂ depends on β (uniform inf-sup), c₀, shape regularity and κ.
3. **Why the flux correction.** The q = 1 continuity row forces `∫_{Γ_h}u_h·n = −∫_{∂R}g_I·ν = −m_R`,
   while `∫_{Γ_h}ũ·n = 0`. The pressure mean is the multiplier of the Γ_h flux: testing with v₁
   (v₁·n ≡ 1 on Γ_h, achievable in V_h) brings in `(μ/h)∫e·n = (μ/h)m_R ~ μh^{k}`. With μ ∝ h/h_Γ
   this is `γ h^{k+1}/h_Γ`, which needs h/h_Γ bounded (a "mesh proviso"). With g̃_I, ∫_{Γ_h}e·n = 0
   and the term vanishes, so the correction is well motivated.
4. **Energy / H¹ for fixed γ ≥ γ₂.** Same as A, but without the pressure load:
   `|||e||| ≤ C[(1+γ^{½}+γ^{-½})h_Γ^{3/2} + h^k]`.
5. **γ-uniform H¹ for all γ ≥ γ₂.** The penalty consistency `(μ/h)⟨ũ,v⟩` grows like γ^{½}, and so do
   the interpolation terms. A γ-uniform H¹ bound needs a **different comparison function**: a
   div-free w with w = 0 on Γ_h, obtained by correcting ψ̃ by a boundary-layer lifting of its Cauchy
   data on Γ_h (sizes h_Γ⁴ and h_Γ²).
   * Then ‖∇(ũ−w)‖ ~ h_Γ^{3/2}, and w needs a discretely div-free interpolant that also vanishes on Γ_h.
   * The penalty term drops because (μ/h)⟨e,e⟩ ≥ 0 on the good side.
   * Plausible, but the corrector plus an interpolant preserving both div-free and zero trace is a
     genuine piece of work.
6. The claim "(∇ũ)ᵀn = O(h_Γ^{3/2}) because ũ = 0 on Γ and chords deviate O(h_Γ²)" is **imprecise**.
   (∇ũ)ᵀn is O(h_Γ) pointwise, driven by the O(h_Γ) **normal** deviation. The O(h_Γ^{3/2}) dual bound
   comes from the penalty (γ^{-½}), or from edgewise odd symmetry plus integration by parts. Neither
   is "ũ = 0 and chords deviate O(h_Γ²)". In pure H¹-duality, without the penalty, it is only O(h_Γ).

Numerics (svn.py, CNS*, kind C, **g_I not flux-corrected**):

| μ | γ | ‖∇e‖ at N = 16 / 32 / 64 | rate vs h_Γ |
|---|---|---|---|
| 40 | ≈12 | 0.61 / 0.109 / 0.041 | 2.7, 1.43 |
| 400 | ≈120 | 0.137 / 0.030 / 0.0095 | 2.4, 1.69 |

The rate is consistent with h_Γ^{3/2}. The error at γ ≈ 12 is ~4× the error at γ ≈ 120. That is
consistent with a stability constant that degrades as γ ↓ γ₂ (the 1/(1 − Cγ^{-½}) structure), not
with the γ^{½} upper-bound growth. So the constant of the fixed-γ bound is **not** uniform near γ₂.
Any "γ-uniform" statement must be for γ ≥ γ₂′ > γ₂, or must carry an explicit 1/(1−(γ₂/γ)^{½})
factor.

**Verdict D:** well-posedness TRUE for γ ≥ γ₂. The fixed-γ error bound is TRUE, and flux correction
is the right way to remove the mesh proviso. γ-uniform H¹ is plausible but needs the zero-trace
div-free corrector. The consistency claim for (∇ũ)ᵀn needs its justification rewritten.

### PR. CNS* pressure

1. `‖π−m(π)‖` follows from the V_h^0 inf-sup. For v ∈ V_h^0 the residual is only `⟨ũ,∂_n v⟩` =
   O(h_Γ^{3/2})‖∇v‖ (γ-free), plus `⟨e_h,∂_n v⟩ ≤ γ^{-½}|||e_h|||‖∇v‖`.
2. **Mean m(π):** test with v₁ ∈ V_h, v₁·n ≡ 1 on Γ_h (possible: at each vertex solve
   v·n⁻ = v·n⁺ = 1), with ‖∇v₁‖ = O(1) and v₁·τ = O(h_Γ).
   * `(μ/h)⟨e_h,v₁⟩ = (μ/h)[∫e_h·n + ⟨e_h·τ, v₁·τ⟩]`. The first term is **0 only with flux-corrected
     data** (and a zero-flux w_h). The second is at most γ^{½}h_Γ^{½}|||e_h|||.
   * **`⟨∂_n e_h, v₁⟩` naively costs h_Γ^{-½}|||e_h||| = O(h_Γ)**, which loses the rate. The normal part
     must again be handled by the div-free identity:
     `∫ n·∂_n e_h = −Σ_edges[e_h·τ]`, a vertex sum at most C‖e_h‖_{Γ_h}.
3. Result for fixed γ ≥ γ₂: `‖p*−p_h‖ ≤ C(h_Γ^{3/2}+h^k)`.

**Verdict PR:** TRUE for fixed γ with flux correction, **provided** the mean is estimated with the
div-free identity. A crude trace/inverse bound gives only O(h_Γ).

### E. Penalty necessity

1. **Choice.** Let e = [a,b] and c the opposite vertex. Take λ₁ = λ_a, λ₂ = λ_b, the barycentrics of
   the **two endpoints of e**. Each vanishes on the interior edge opposite its vertex.
2. Why this choice works:
   * φ = λ_a²λ_b²ℓ has φ = ∇φ = 0 on edges bc and ac, so v = curl φ vanishes on ∂T\e (including the
     vertices).
   * v extends by zero to V_h^R ∩ {div = 0} when ℓ ∈ P_{k−3}; for k = 4, ℓ is linear.
   * The choice λ_c would not work: it vanishes on e itself and kills the trace, making the penalty
     irrelevant only through ⟨∂_n v, v⟩ = 0.
3. Scaling: N_h(v,v) < 0 for some ℓ iff `μ|e|/h < λ_T := max_ℓ (2⟨∂_n v,v⟩_e − ½‖Dv‖²_T)|e|/‖v‖²_e`,
   which is scale-invariant.
4. Computed exactly (sympy assembly, generalized eigenproblem, n = outward of Ω_h = away from T):

| triangle (|e| = 1) | k = 4 | k = 5 | k = 6 |
|---|---|---|---|
| equilateral | **9.3138** | 9.644 | 10.573 |
| apex (0.5, 0.6) | 16.29 | | |
| apex (0.2, 0.9) | 8.39 | | |
| apex (0.5, 1.2), tall | 4.54 | | |

**Verdict E:** TRUE for near-equilateral T (λ_T ≈ 9.31 > 9 at k = 4). λ_T is **shape-dependent**: it
drops to ≈ 4.5 for tall triangles. So "μ₀ ∝ ρ is necessary" requires a uniform positive lower bound
on λ_T over the shape-regular class, taken at the triangle attaining min|e|. The necessary constant is
λ_min(shape class), not 9.

(End of Part I.)

---

## Part II: comparison with the paper. Written 2026-09-28T14:42Z, after reading the proofs

Sources read:
* `paper/sections/02–08`
* `notes/v4/rescue.tex`, all of it
* `notes/v4/pressure_drag.tex` §§1–3: setting, tools, CNS* pressure

### Errata to Part I (found on re-reading; the substance is unchanged)

1. **E, "choice" bullet.** The sentence about λ_c is garbled. Correct statement: the factor must be
   λ_a²λ_b², where a and b are the endpoints of e. Using λ_c² would make v vanish on e, so the test
   field would not see the penalty or the Nitsche terms at all.
2. **μ = 40 runs.** These have γ ≈ 12, which is **below** the measured coercivity threshold on these
   meshes (γ* ≈ 18–21, paper §8 table). So they are not tests of the theorems, which assume γ ≥ γ₀ or
   γ ≥ γ₂.
   * The μ = 160 and μ = 400 runs (γ ≈ 48 and ≈ 120) are in range, and they confirm C, B′ and the D rate.
   * The larger CNS* error at γ ≈ 12 is therefore sub-threshold behaviour. It is not evidence against
     γ-uniformity.
3. **P1 kernel, post-hoc check.** Redone with an element-wise L²-orthonormal pressure basis. The
   smallest singular values are:

   | N | σ_min | next σ |
   |---|---|---|
   | 8 | 8e-17 | 2.4e-2 |
   | 12 | 3e-16 | 8.1e-3 |

   So the kernel is exactly one-dimensional. The 1e-9 "cluster" in Part I was a basis-scaling
   artefact.

### Claim-by-claim comparison

| claim | paper location | agreement | discrepancies (severity) |
|---|---|---|---|
| P1 | Prop. 7.1 + the paragraph after it | Agrees: square, kernel contains (0,1), incompatible for generic data. | (minor) The paper states "has the kernel vector" and says the left null vector is **not identified**. Part I identifies incompatibility exactly (P2 below), and numerics show the kernel is exactly one-dimensional. |
| P2 | **Not in the paper.** Only the remark "discrete divergence about 4.5e-3 … consistent with incompatibility". | n/a | (minor, recommended addition) State the exact result: with p ∈ L²_0 the momentum row is unchanged, div u_h ≡ c, `c|Ω_h| = ∫_{∂R}g_I·ν + ∫_{Γ_h}u_h·n`, and **c = 0 iff the naive system is consistent**. Well-posedness of the L²_0 system follows from the Thm 7.2 absorption argument for γ ≥ γ₂. The observed 4.5e-3 divergence is this c. |
| L | Lemma 4.2 | Identical proof, κ defined by ‖∇v‖² ≤ κ a_h(v,v), same constants. Uniform Korn via the bi-Lipschitz map Φ_h (Lemma 3.2) is correct. | none |
| A | Thm 6.1 | Upper and lower bounds agree. The paper bounds (∇ũ)ᵀn in the **H¹-dual** at O(h_Γ^{3/2}) by the odd-in-s cancellation (Lemma 3.5). That is sharper than the γ^{-½} penalty route in Part I, and the two are equivalent for γ ≥ γ₀. Lower bound: same test field curl(Π_{k+1}χφ₀), same continuity argument. | none. (Cor 6.2 needs γ^{½}h_Γ ≪ G as well as γh_Γ ≪ G. This is implied when γ ≥ 1; worth one line if γ₀ < 1 is possible.) |
| C | Thm 6.4 + Lemma 6.3 | Same mechanism: subtract (h/μ)v_φ; the div-free identity n·∂_n v = −∂_τ(v·τ); edgewise integration by parts; vertex terms O(|e|)|v(x)| summed to C‖v‖_{Γ_h}. | none |
| B′ | Thm 6.5, Cor 6.6 | Same hypotheses (G > 0, h/μ → 0, γh_Γ → 0, ε/(h/μ)^{½} → 0) and the same error budget. | none |
| D | Thm 7.2; rescue Lemmas 9.x, Thms rs:unif and rs:Dfc | Well-posedness: same absorption argument, γ₂ ∝ β^{-2}. The paper's version is cleaner: testing with e_h ∈ Z_h, the constant pressure mode drops out exactly. Flux correction (rs:approx, rs:floor): same mechanism as Part I. **γ-uniform H¹ (rs:unif)**: the paper closes the gap Part I left open, and more economically. It compares with the *discrete* strongly-imposed Dirichlet solution u^D (Lemma rs:dir), uses the identity on Y_h = Z_h ∩ V_h^0 (Lemma 7.x) and a div-free discrete lift (rs:lift). The only coupling is h_Γ^{-½}‖u_h‖_{Γ_h}, and γ^{-½}·ε cancels the γ-growing parts of ε. Checked line by line; correct. | (a) (typo) rs:unif step 5: "ũ−u_h = (ũ−u^D) − w − δ_R" should be "+ δ_R". Harmless for norms. (b) (presentation) The claim as summarised in the task ("(∇ũ)ᵀn = O(h_Γ^{3/2}) because ũ = 0 on Γ and chords deviate O(h_Γ²)") is **not** the paper's argument, and it would be wrong as stated. Pointwise, (∇ũ)ᵀn = O(h_Γ) because the chord *normal* tilts by O(h_Γ). The paper correctly uses the odd-in-s structure (Lemma 3.5), or the penalty route (pressure_drag Lemma pd:lem:G). Any summary or abstract text with the "chords deviate O(h_Γ²)" wording should be corrected. |
| PR | pressure_drag Thm pd:thm:cnsp | The mean-free part agrees (V_h^0 inf-sup). **Mean:** the paper uses the normal test field w₁ = I_h(χ₀n̂) and handles ⟨∂_n e_u, w₁⟩ with the div-free cancellation (pd:lem:cancel, with div δ = O(h_Γ^k) on the first layer). This is exactly the step Part I flagged as necessary: a crude trace/inverse bound would give only O(h_Γ). The flux term (μ/h)m_R appears explicitly and vanishes with correction. | (c) (minor, statement gap) As written, pd:thm:cnsp concludes "C(h_Γ^{3/2}+ε), and C(h_Γ^{3/2}+h^k) under the proviso h_Γ ≥ h^{2(σ_k−k)}", **even for flux-corrected data**. With correction, Y must be built with ε̃ (rescue Lemma rs:approx and Thm rs:Dfc), and then ε̃ ≤ C(h^k + h_Γ^{k+½}) for bounded γ. So the claimed PR (no proviso) holds, but only by combining the two notes; one sentence is needed. (d) (scope) PR is **not γ-uniform**: the constant carries (1+(γh_Γ)^{½})(1+γ^{½}). The note says so (Remark pd:rem:cnsp). PR must be stated for γ ∈ [γ₂, γ_max]. |
| E | §8, Thm pb:thm, Prop pb:cert, Cor pb:rho | Same witness (λ₁, λ₂ = barycentrics of the two endpoints of e, ψ = λ₁²λ₂²ℓ, ℓ ∈ P₁). λ_T(equilateral, k = 4) = 9.31381 in both, to all printed digits. Shape dependence as in Part I (my tall-triangle value 4.54 at apex angle ≈ 45° matches the paper's level-set minimum 4.47). | (scope) "Hence μ₀ ∝ ρ is necessary" is true only under a shape hypothesis on the triangle of a *shortest* Γ_h edge (Cor pb:rho: apex ≥ 35°, base angles ≥ 5°). For tall triangles λ_T < 0. A uniform necessity statement under (M0)–(M2) alone is open (Remark pb:limits). The paper states this correctly; the claim as summarised needs the hypothesis. |

### Items I could not fully certify (accepted on citation)

* **Lemma 4.4** (uniform inf-sup): this rests on Guzmán–Scott 2019, Thm 1, applied on the doubly
  connected Ω_h with (M1)/(M2). The Bogovskii step (a cover by pieces star-shaped w.r.t. balls, then
  Galdi III.3.4) is standard. I did not re-derive GS's local Fortin construction. A referee should
  confirm that GS need no simple-connectivity or boundary-vertex hypothesis beyond (M1)/(M2).
* **Morgan–Scott / Argyris C¹ interpolation** of the stream function (k+1 ≥ 5): standard.

### Final verdicts

| claim | verdict | precise fix |
|---|---|---|
| P1 | **TRUE** | Optional: state that the kernel is exactly span(0,1) (numerically verified) and that consistency ⇔ c = 0 (P2). |
| P2 | **TRUE** (not stated in paper) | Add: div u_h ≡ c with c|Ω_h| = ∫_{∂R}g_I·ν + ∫_{Γ_h}u_h·n. c = 0 iff the naive system is consistent. Generic data give c ≠ 0 (numerically −1.5e-3 at N = 8 and 7.6e-3 at N = 12). |
| L | **TRUE** | none |
| A | **TRUE** | Cor 6.2: also require γ^{½}h_Γ ≪ G, or note that it is implied by γ ≥ 1. |
| C | **TRUE** | none |
| B′ | **TRUE** under G > 0, h/μ → 0, γh_Γ → 0, ε(γ/h_Γ)^{½} → 0 | none |
| D | **TRUE** (well-posedness for γ ≥ γ₂; flux-corrected H¹ bound with no proviso; γ-uniform in H¹ via rs:unif) | Fix the "−δ_R" sign typo in rs:unif step 5. Correct any text that attributes the O(h_Γ^{3/2}) transposed-gradient bound to "chords deviate O(h_Γ²)": it is O(h_Γ) pointwise and O(h_Γ^{3/2}) only through odd-in-s cancellation (Lemma 3.5) or the penalty (pd:lem:G). |
| PR | **TRUE for γ ∈ [γ₂, γ_max]** with flux correction; **NEEDS FIX** in statement | In pd:thm:cnsp, for flux-corrected data replace ε by ε̃ (rescue Lemma rs:approx) and drop the proviso h_Γ ≥ h^{2(σ_k−k)}. State the γ-dependence (not γ-uniform). |
| E | **TRUE with extra hypothesis** | Necessity of μ₀ ∝ ρ needs the shortest Γ_h edge's triangle to lie in the certified shape class (apex ≥ 35°, base ≥ 5°). λ_T ≈ 9.314 > 9 near-equilateral (k = 4). λ_T < 0 for tall triangles; uniform necessity under (M0)–(M2) is open. |
