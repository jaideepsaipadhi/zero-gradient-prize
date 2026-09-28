# Rescuing the printed method with a growing penalty, and flux-corrected outer data: report

Files (all new):
- `notes/v4/rescue.tex`: drop-in section. It uses the paper's macros and labels `rs:`. Insert it after §8 with `\input{…/rescue}`. A copy of the paper with this input compiles with no undefined references (76 pages).
- `code/rescue.py` (modes `study`, `musweep`, `flux`, `cond`, `fluxsize`), `code/rescue_run.sh` (batch driver) and `code/rescue_report.py` (tables → `results/v4/rescue/summary.md`).
- `results/v4/rescue/{study,musweep,flux,cond}/*.jsonl`.
- `logs/rescue_{cond,flux,musweep,study}.log`.

Setup: k = 4, the paper's meshes (h/h_Γ = 3.28–3.43, so γ ≈ 0.29 μ), one thread. N ≤ 64 used SuperLU with at most 2.2 GB. N = 96 used PARDISO (1.2 GB) because SuperLU needs 5.6 GB there. The largest relative residual was 1.8e-12 and the largest divergence 1.6e-9. The μ = 100 test-P run reproduces the paper's value at N = 96 (1.2810e-2).

## Task A: GS with μ = μ(h)

### 1. Theorems A and C alone are not enough

Take γ ≍ h_Γ^{-α}. Theorem C then gives ‖∇e‖ ≲ h_Γ^{1+α} + h_Γ^{3/2−α/2}. The two terms balance at α = 1/3, so the best provable H¹ rate from Theorem C is **4/3**. Beyond α = 1/3 the bound gets worse, because the γ^{1/2}h_Γ^{3/2} term grows.

### 2. New estimate: an H¹ bound uniform in γ (Theorem rs:unif)

For every γ ≥ γ₀, with no upper bound,

  ‖∇(ũ − u_h)‖ ≤ C_p h_Γ/γ + C(h_Γ^{3/2} + h^k + h^{σ_k} h_Γ^{-1/2}).

CNS\* satisfies the same bound without the first term, for every γ ≥ γ₂.

How the proof works:
1. Remove the traction defect δ_R. Theorem C's bound for it is already uniform in γ.
2. On Y_h (fields vanishing on Γ_h), the rest of the discrete solution satisfies exactly the discrete Dirichlet equation.
3. Compare it with the discrete Dirichlet solution u^D on the polygon (Lemma rs:dir). u^D carries a flux-bubble trace of size m_R, and its error is ≲ h_Γ^{3/2} + h^k + |m_R|h_Γ^{-1/2}.
4. The difference is controlled by the trace of the discrete solution on Γ_h, using a divergence-free first-layer lift (Lemma rs:lift): ‖∇‖ ≲ h_Γ^{-1/2}‖trace‖.
5. The energy estimate makes that trace small by the factor (h/μ)^{1/2}. This factor exactly cancels the γ-growth of γ^{1/2}h_Γ^{3/2}, of (γh_Γ)^{1/2}h^k and of h^σ(μ/h)^{1/2}.

Side results:
- The H¹ part of Theorem D does not need γ ≤ γ_max.
- The flux proviso h_Γ ≥ h^{2(σ_k−k)} does not depend on μ.

### 3. The rescue (Corollary rs:rescue)

If γ ≳ h_Γ^{-1/2}, i.e. μ ≳ h h_Γ^{-3/2} (for h ≍ h_Γ this is μ ∝ h^{-1/2}, μ/h ∝ h^{-3/2}), then

  ‖ũ − u_h‖_{H¹} ≤ C(h_Γ^{3/2} + h^k)

for general Stokes data. The rate is sharp by Theorem lb:thm, which holds for GS for every γ and every G.

### 4. The energy norm cannot be rescued (Theorem rs:elower, Corollary rs:rates)

**New lower bound**, for GS and CNS\*, whenever W ≢ 0 and γ ≥ γ₁:

  |||e||| ≥ (μ/h)^{1/2}‖e‖_{Γ_h} ≥ (c₀²/72)‖W‖ γ^{1/2} h_Γ^{3/2} − C(h_Γ/γ)^{1/2}.

The proof tests with the tangential field v_W = curl Π(χ(r−1)W).

Combined with Theorem A:

  |||e||| ≍ (h_Γ/γ)^{1/2}G + γ^{1/2}h_Γ^{3/2}‖W‖.

This holds outside a window γh_Γ ∈ (λ₀, Λ₀). Inside the window the upper bound is ≍ h_Γ, and a matching lower bound is not proved. The product of the two terms is h_Γ², so the energy error can never go below h_Γ.

Rates as a function of α, for G > 0 and W ≢ 0:

| α | H¹ | energy |
|---|---|---|
| α < 1/2 | 1+α (exact) | (1+α)/2 |
| 1/2 ≤ α < 1 | **3/2** (exact) | (1+α)/2 |
| α = 1 | 3/2 | ≤ 1 (energy-optimal) |
| α > 1 | 3/2 | (3−α)/2; no convergence for α ≥ 3 |

- H¹-optimal: any α ≥ 1/2. The cheapest is α = 1/2.
- Energy-optimal: α = 1, giving rate 1. **No scaling reaches energy rate 3/2**; only the traction (CNS\*) does.
- For CNS\*, and for GS with G = 0, any growth of γ only costs. The energy rate drops to (3−α)/2, and the H¹ constant grows by a bounded factor.

### 5. Coercivity

Coercivity is only helped by a larger μ (Remark rs:coercive). N_h^{μ'} ≥ min(c₁,1)|||·|||²_{μ'} for μ' ≥ μ ≥ μ₀. Also, M and the CNS\* absorption threshold γ₂ are monotone.

### 6. Conditioning (Proposition rs:cond)

Setting: Lagrange velocity basis, L²-orthonormal pressure basis, quasi-uniform mesh.

- **Velocity block:** λ_max(A_μ) ≍ 1+γ and λ_min(A_μ) ≍ h², so κ(A_μ) ≍ (1+γ)h^{-2}.
- **Schur complement:** λ_min(S_μ) ≍ (1+μ/h)^{-1}. This mode is the constant pressure, since the penalty controls ∫_{Γ_h} v·n. The rest of the spectrum is ≥ cβ² and ≤ 2/c₁. So κ(S_μ) ≍ 1 + μ/h = 1 + γ/h_Γ.
- **Saddle matrix:** κ(K_μ) ≤ C(1+γ)·max(h^{-2}, 1+μ/h). This is quadratic in the penalty once the Schur branch is active.
- **Price for α = 1/2:** a factor h^{-1/2} in both κ(A) (h^{-2} → h^{-5/2}) and κ(S) (h^{-1} → h^{-3/2}).

### 7. Numerics (N = 16…96)

**Rates between N = 64 and N = 96** (Table rs:tab:rates):

| test | method | (α, c) | H¹ rate | energy rate | theory (H¹, energy) |
|---|---|---|---|---|---|
| P | GS | (0, 100) | 0.95 | 0.49 | 1, ½ |
| P | GS | (½, 100) | 1.45 | 0.73 | 1.5, ¾ |
| P | GS | (½, 1000) | 1.47 | 0.73 | 1.5, ¾ |
| P | GS | (1, 100) | 1.94 | 0.97 | 2, 1 |
| B₁ | GS | (0, 100) | 0.96 | 0.56 | 1, ½ |
| B₁ | GS | (½, 100) | 1.42 | 0.86 | 1.5, ¾ |
| B₁ | GS | (½, 1000) | 1.52 | 1.28* | 1.5, ¾ |
| B₁ | GS | (1, 100) | 1.47 | 1.06 | 1.5, 1 |
| A | GS | (1, 100) | 1.34** | 1.16 | 1.5, 1 |

\* This is the geometric branch (3−α)/2 = 1.25; the leak branch ¾ takes over only later.
\*\* This is a transient (see the test A item below).

What the runs show:
- **Test P:** the normalised leak constants are unchanged. ‖∇e‖μ/(hG) = 2.71–2.78 and |||e|||(μ/h)^{1/2}/G = 0.996–0.999 at N = 96.
- **B₁ at N = 96:** the GS H¹ error halves, from 1.27e-2 to 6.4e-3. CNS\* at μ = 100 gives 3.5e-3.
- **B₁₀₀:** the error drops by a factor of 22 (c = 1000).
- **Test A, μ-sweep at fixed N** (Table rs:tab:mu):
  - The H¹ error saturates at the Dirichlet-limit value: 4.557e-2 at N = 64. That is about 2.5× the μ = 100 value. The Dirichlet-limit rates are 1.90, 1.80, 1.72, 1.66, heading towards 3/2.
  - This explains the α = 1 transient rate on A: the ratio of the GS error to the Dirichlet limit climbs from 0.37 to 0.56 between N = 16 and N = 64.
  - The energy grows like γ^{1/2}h_Γ^{3/2}. The constant is 0.339 for A and 0.105 for B₁, which equals ‖ũ‖_{Γ_h}/h_Γ² exactly.
- **Conditioning, N = 16…48** (Table rs:tab:cond):
  - λ_min(A) ≈ 0.03h², independent of μ.
  - λ_max(A) ≈ max(102–118, 0.47γ). So for γ ≲ 200 the P4 stiffness dominates and moderate growth is free.
  - The smallest Schur eigenvalue equals (|Γ_h|/|Ω_h|)h/μ to within 2.5% at every N and μ. The second one is 0.007–0.010 and independent of μ.
  - κ(K) goes from 4.9e4 at μ = 100 to 6.0e9 at μ = 1e5 (N = 32).
  - The α = ½, c = 100 study stays in the linear regime: κ(K) rises by a factor ≤ 2.3.

## Task B: flux-corrected outer data

**Definition.** g̃_I := g_I − (m_R/(8L²))·x on ∂R.
- x is linear, so it is reproduced exactly by the P_k interpolant and is continuous at the corners.
- ∫_{∂R} x·ν = 2|R| = 8L², so ∫ g̃_I·ν = 0 exactly (≤ 1.3e-15 numerically).
- The accuracy is unchanged: |m_R| = O(h^{σ_k}) ≤ h^{k+1}.

Results:
- **Floor, exact with constant 1 (Proposition rs:floor).** Uncorrected: |||e||| ≥ (μ/h)^{1/2}‖e‖_{Γ_h} ≥ (μ/h)^{1/2}|m_R|/|Γ_h|^{1/2}, because ∫_{Γ_h} e·n = m_R. Corrected: ∫_{Γ_h} e·n = 0.
- **Approximation (Lemma rs:approx).** ε̃_h = (1+(γh_Γ)^{1/2})h^k + (1+γ^{1/2})h_Γ^{σ_k−1/2} + |m_R|.
  - The remaining flux term is the flux of the **Γ_h-interpolation error**, ∫_{Γ_h}(I_hũ − ũ)·n = O(h_Γ^{σ_k}). Its cost is exactly the referee's O(h_Γ^{σ_k−1/2}(1+γ^{1/2})), confirmed.
  - It is below h_Γ^{k+1/2}(1+γ^{1/2}).
  - It is an artefact of the comparison field, not a floor: u_h and ũ both have zero flux through Γ_h.
- **Theorem D without the proviso (Theorem rs:Dfc).** For every h_Γ ≤ h and γ ≥ γ₂ (the H¹ part uniformly in γ): ‖ũ − u_h‖_{H¹} ≤ C(h_Γ^{3/2} + h^k). For bounded γ the energy norm obeys the same bound. Corollary lb:cor (sharpness) then holds with no flux condition.

**Numerics.** The paper's tests all have m_R = 0 up to round-off (|m_R| ≤ 2e-13): B and P are polynomial, and A is symmetric.

New test F:
- ψ = (r²−1)³ sin(2x+1) cos(3y/2 + 7/10)/400.
- The cube gives W ≡ 0, so the geometric slip is O(h_Γ⁴).
- The phase shift breaks the symmetries.
- m_R = 2.0e-3 … 3.7e-7 for N = 16…64. It decays like h^{6.6}, against σ₄ = 6.

Results for CNS\*:
- **Uncorrected slip error at μ ≥ 10⁴:** it equals √(m_R²/|Γ_h| + ‖ũ‖²_{Γ_h}) to 3 digits, for example 8.204e-4 against the floor 8.19e-4 at N = 16. This is the floor plus the orthogonal geometric slip.
- **Corrected slip error:** it equals ‖ũ‖_{Γ_h} to 4–5 digits, for example 2.396e-5 at N = 16.
- **Flux through Γ_h after correction:** ≤ 2.5e-14.
- **Energy norm at μ = 10⁸:** 6.75 → 0.955 (N = 16), 0.70 → 0.37 (N = 24), 0.187 → 0.149 (N = 32). The uncorrected value sits on the floor (μ/h)^{1/2}|m_R|/|Γ_h|^{1/2}.
- **H¹ errors:** identical to relative 3e-6 … 1e-9, as Theorem rs:unif predicts (the flux enters H¹ only as |m_R|h_Γ^{-1/2}).

Test F2 has W ≠ 0, i.e. (r²−1)². There the geometric slip exceeds the floor by 23× to 2×10⁴, so the floor cannot be seen at feasible N. This is the practical reason the paper's proviso is harmless.

## Caveats

- **Energy window.** At the energy crossover γh_Γ ∈ (λ₀, Λ₀) the matching energy lower bound is not proved. The normal leak and the tangential geometric slip should not cancel, but step 2 of Theorem rs:elower does not resolve the tangential discrete trace finely enough to rule it out.
- **Mesh assumptions.** Proposition rs:cond assumes quasi-uniformity. The H¹ lower bound for α ≥ ½ uses (M3), as in §7.1.
- **Two-sided constants.** The Schur-complement lower bound constant is not explicit. The numerics give s₁ = (|Γ_h|/|Ω_h|)h/μ·(0.975–0.996).
- **Pre-asymptotic rates.** Some reported rates are pre-asymptotic: B₁ at α = ½, c = 1000 is in the geometric energy branch, and A at α = 1 is in the H¹ transient. Both are explained by the μ-sweep; see the notes in section 7 above.
