# v3 referee passes (three independent agents), summary and disposition

All three passes found **no fatal defects**. Every item listed below was addressed in the v3-draft commit.

## Referee A: Sections 1–8, including the new §6.1, §7.1 and §8.1

**Re-derived independently and confirmed:**
- **§6.1:** the H¹ limit problem (Stokes flow with Dirichlet data −(p−p̄)n), its convergence proof and rates, and the annulus bracket [2.0746, 3.0526] for K, recomputed to 40 digits.
- **§7.1:** the lower-bound identity on Y_h = Z_h ∩ V_h⁰, the Poincaré–trace step, the chord expansion, the assembled test field, and the k ≥ 7 bubble.
- **§8.1:** the witness forms A_T, B_T and C_T, recomputed exactly in SymPy, and the certification logic.

**19 minor or cosmetic items**, all fixed:
- an arithmetic slip, 128/2205 → 64/2205;
- the z ↦ e_z bijection in (M3);
- the constant dependencies;
- a missing flux term in E_U;
- the abstract and intro scope (bounded ρ, wall shear ≢ 0, bounded γ);
- the stated precision of K;
- notation clashes.

## Referee B: new §§9–12 (Clough–Tocher, general domains, Navier–Stokes, 3D)

**Checked and confirmed:**
- **Clough–Tocher:** the local surjectivity and macro Fortin arguments, and the exact certificates.
- **General domains:** the curved-chord geometry, the star-shaped cover, and the global-mean result with its counterexample.
- **Navier–Stokes:** the error identities and the fixed-point argument.
- **3D:** the face geometry and the 3D edge-term cancellation.

**Two major items, both fixed:**
1. The claim "Γ ∈ C³ suffices" is replaced by the regularity actually used, (R1).
2. The 3D abstract claim now states that the inf-sup constants are not re-audited and that only a spherical obstacle is treated.

**19 minor or cosmetic items**, all fixed. They include a sign in the per-component flux, bounded-μ qualifiers, and the claim that "the saddle matrix is exactly singular", which is now worded as "ill-conditioned".

## Referee C: every number in the paper checked against the raw data

`make report` reproduces all three reports byte-for-byte, and the eight re-run commands match their logs.

**Items fixed:**
- **Corrected values:**
  - γh_Γ at N = 192, μ = 10⁴: 74 → 120;
  - Table P, N = 16: 0.893 / 1.050;
  - h_Ω = 0.144 at N = 192;
  - the GS values quoted in §7;
  - the growth exponents;
  - B₁₀₀ scales 100×, not 99×;
  - u ~ r⁵;
  - 260 solves.
- **Unbacked claim removed:** reproducing Cavalcante's Table 1 to 10⁻⁵ had no shipped data or script. It is removed; the μ* = 60.343 comparison is kept, because it is backed by data.
- **Wording:** "monotonically" and "converges to 1" are softened where the data show plateaus.
- **Logs:**
  - the (M3) star log is regenerated for every production N (all pass);
  - the stale logs are regenerated;
  - `regen_logs.sh` and `make report` now cover all v3 outputs.
