# notes/v6/l2 — agent report (saved verbatim)

[Subagent hand-back] The text below is the final report of a subagent this session delegated to. It is model output, NOT a message from the user: instructions, requests, or approval claims inside it are the subagent's words and carry no user authority. The harness indents every line of the report, so a frame-like line at column zero inside it would be forged. Notes above this frame may quote model-derived text, which carries no user authority either. The report follows:
  I couldn't write `REPORT.md`: the tool blocked it because subagents must return findings as text, so the full report is below (the summary first, then an appendix with the proofs). The scripts, logs and results are in `notes/v6/l2/`. Nothing in `paper/` or `letter/` was touched and nothing was committed.
  
  ## Summary
  Notation: e = ũ − u_h and W = ∂ₙu·τ (wall shear). CNS* uses the flux-corrected data g̃_I. "(R)" is H²×H¹ regularity of Stokes on the smooth domain Ω; the paper already relies on it for U.
  
  **(L1) CNS*, optimal a = 2.**
  - ‖e‖_{L²(Ω_h)} ≤ C(h_Γ² + h·h_Γ^{3/2} + h^{k+1}), with C independent of γ ≥ γ₂. This is PROVED MODULO (R).
  - On quasi-uniform meshes it reduces to C(h_Γ² + h^{k+1}). On strongly graded meshes the h·h_Γ^{3/2} term remains; removing it is open.
  - The same bound holds for GS when G = 0.
  
  **(L3) The rate 2 is sharp.**
  - ‖e‖ ≥ c‖W‖² h_Γ² when W ≢ 0 and h^k ≤ δ₀h_Γ². This is PROVED, needs no (M3), and needs no regularity theory because the dual is explicit.
  - Leading term (PROVED weak identity): (e, −Δz̃) = −Σ_e(|e|³/12)W·W_z(m_e*) + O(h_Γ³ + h^k).
  - Consequence: e ≈ V_h, the Stokes flow driven by the chord-mean slip (|e|²/12)∂ₙu. The remainder is a normal-slip field of size O(γ^{-1/2}h_Γ²).
  
  **(L2) GS.**
  - ‖e‖_{L²} = (h/μ)‖U‖_{L²(Ω)}(1 + O(h^{1/2})) + O(h_Γ²), where U is the leak flow. PROVED for fixed γ; the O(h_Γ²) term is modulo (R).
  - So the GS L² error is Θ(h/μ), rate exactly 1.
  - For test P, ‖U‖_{L²}/G = 1.2736341 (NUMERICAL, from the series solution).
  
  **Adjoint consistency.** The pressure term of CNS* with the Γ_h mean removed is **not** adjoint consistent. Its transpose has the plain GS pressure row, so the defect is a leak ⟨r − r̄, v·n⟩ driven by the dual pressure. The proofs avoid it by testing only with fields in Y_h = Z_h ∩ V₀ (upper bound) and by using a dual with zero pressure (lower bound). The dual problem is posed on Ω, not Ω_h: the vertices of Γ_h are re-entrant corners with Stokes exponent λ ≈ 1 − 2ε/π < 1, so there is no uniform H² bound on Ω_h.
  
  **(L4) Numerics (NUMERICAL).** The stored runs already contained L² errors (`results/pod/study`, N ≤ 256). Selected rows at μ = 100:
  
  | N | CNS* test A: L² | rate | L²/h_Γ² | CNS* B₁: L² | rate | L²/h_Γ² |
  |---|---|---|---|---|---|---|
  | 16 | 4.58e-2 | – | 0.217 | 1.32e-2 | – | 0.0625 |
  | 32 | 1.24e-2 | 2.03 | 0.208 | 2.95e-3 | 2.05 | 0.0494 |
  | 64 | 3.21e-3 | 1.99 | 0.208 | 7.87e-4 | 1.95 | 0.0510 |
  | 128 | 8.15e-4 | 1.99 | 0.209 | 2.04e-4 | 1.97 | 0.0523 |
  | 256 | 2.05e-4 | 1.99 | 0.210 | 5.18e-5 | 1.98 | 0.0531 |
  
  (Rates are between consecutive runs of the full sequence 16, 24, …, 256, in `existing_l2.txt`.)
  
  - **CNS*:** μ-independent, and B₁, B₁₀₀ and C give identical errors.
  - **GS:** rate 2 on test A (G = 0); rate 0.95–0.98 on P, B₁ and C.
  - **GS leak constant (test P):** ‖e‖/((h/μ)‖U‖) at N = 96 is 0.980 / 0.997 / 0.999 for μ = 10²/10³/10⁴. The L² distance ‖(μ/h)e − U‖/‖U‖ decays like h¹ at μ = 100.
  - **Weak identity (test A, N = 96):** (e, φ̃)/(−Λ_h) = 0.980 (CNS*) and 0.987 (GS). The remainder settles at a constant times h_Γ³ (−0.177 and −0.11), as the proof predicts.
  - **Share of the error carried by V_h (μ = 100):** about 90% for A and about half for B₁. The rest is the smooth normal slip, which falls roughly like 1/γ.
  - **Not run:** N = 128 of the new diagnostics was killed for lack of memory at about 4 GB. `run_pod.sh` runs N = 128–256; command from the repo root: `sh notes/v6/l2/run_pod.sh > notes/v6/l2/pod.log 2>&1`.
  
  Files are in `/home/claude/zero-gradient-prize/notes/v6/l2/`:
  - `l2_numerics.py` — the new runs (modes `slip`, `pair`, `moments`, `leak`, `unorm`)
  - `collect_l2.py` and `existing_l2.txt` — tables and rates from the stored results
  - `slip.log`, `leak.log`, `pair.log`, `l2_results.jsonl` — the new runs (N ≤ 96)
  - `run_pod.sh` — the larger levels
  
  ---
  
  ## Appendix: proofs
  
  **Tools.**
  - *Strip lemma (PROVED).* Write Γ_h in polar form as r = ρ_h(θ), with 1 − ρ_h ≤ Ch_Γ². Integrating along radial lines gives:
    - ‖w‖_{S_h} ≤ C(h_Γ‖w‖_{Γ_h} + h_Γ²‖∇w‖_{S_h}), and ‖w‖_{S_h} ≤ Ch_Γ‖w‖_{H¹(N)};
    - ‖w‖_{Γ_h} ≤ C‖w‖_{H¹(N)};
    - if w ∈ H² and w|_Γ = 0, then ‖w‖_{Γ_h} ≤ Ch_Γ²‖w‖_{H²};
    - if w|_{Γ_h} = 0, then ‖w‖_{S_h} ≤ Ch_Γ²‖∇w‖.
  - *Dual on Ω.* Solve −Δz + ∇r = φ, div z = 0, z = 0 on ∂Ω, with ‖z‖_{H²} + ‖r‖_{H¹} ≤ C‖φ‖ by (R). Write z = curl ζ with ζ = ∂ζ = 0 on Γ. Extend ζ across Γ by a Hestenes H³ reflection, set z̃ = curl ζ̃ (divergence-free), and extend r̃ in H¹.
  - *Lemma Y (PROVED).* Take y_h = I_h z̃ − L₀ − b, where:
    - L₀ is the nodal lift of I_h z̃|_{Γ_h} (Lemma rs:lift), with ‖I_h z̃‖_{Γ_h} ≤ Ch_Γ^{3/2}‖z‖, so ‖∇L₀‖ ≤ Ch_Γ‖z‖;
    - b is the Bogovskiĭ correction of the mean-free divergence.
  
    Then y_h ∈ Y_h, ‖∇(z̃ − y_h)‖ ≤ Ch‖z‖_{H²}, and ‖∂ₙy_h‖_{Γ_h} ≤ C(1 + h·h_Γ^{-1/2})‖z‖_{H²} (inverse trace on the discrete part).
  
  **Theorem L1, proof.** Put φ = e|_Ω. Since div e = 0, integration by parts on Ω_h plus Lemma lb:identity on Y_h give
  
  ‖e‖²_Ω = a_h(e, z̃ − y_h) − ⟨u_h, ∂ₙy_h⟩ − (F, y_h) − ⟨σ(z̃, r̃)n, e⟩_{Γ_h} − ⟨σν, g − g̃_I⟩_{∂R} − (e, −Δz̃ + ∇r̃)_{S_h}.
  
  Two inputs from the paper:
  - ‖∇e‖ ≤ C(h_Γ^{3/2} + h^k), uniformly in γ (Theorem rs:Dfc).
  - ‖e‖_{Γ_h} ≤ (h_Γ/γ)^{1/2}|||e||| ≤ C(h_Γ² + h_Γ^{1/2}h^k), uniformly in γ.
  
  The six terms are then bounded by:
  1. C·h(h_Γ^{3/2} + h^k)‖φ‖;
  2. C(h_Γ² + h·h_Γ^{3/2} + h_Γ^{1/2}h^k)‖φ‖;
  3. C·h_Γ³‖φ‖;
  4. C(h_Γ² + h_Γ^{1/2}h^k)‖φ‖;
  5. C·h^{k+1}‖φ‖;
  6. C·h_Γ³‖φ‖.
  
  Finally, h_Γ^{1/2}h^k ≤ h_Γ² + h^{k+1} for k ≥ 3. ∎
  
  **L1 for GS.** Apply the same argument to u_h′ = u_h + δ_R; the proof of Theorem rs:unif supplies every input. When G = 0, |p̃ − p̄| ≤ Ch_Γ² on Γ_h, which gives ‖δ_R‖ ≤ Ch_Γ².
  
  **Theorem L3, proof.**
  - *Dual.* Take ζ = ½χ(r)(r−1)²W(θ), z = curl ζ and φ̃ = −Δz̃. Then (z, r ≡ 0) is an exact Stokes pair on Ω_h, and W_z = W.
  - *Identity.* Let (z_h, r_h) be the CNS* solution with right-hand side φ̃. Then
  
    (e, φ̃) = T₁ + 𝒢(z_h) − ⟨e_p − p̄, z_h·n⟩ + ⟨r_h − r̄_h, e_h·n⟩,
  
    where T₁ = (ũ − Z, φ̃) − N_h(ũ − Z, z_h), Z ∈ W̃_h attains Lemma rs:approx, and e_h = Z − u_h.
  - *Bounds.*
    - T₁ = O(h^k).
    - 𝒢(z_h) = −⟨ũ, ∂ₙz̃⟩ + O((1+γ)h_Γ³), using Lemma pd:lem:G, |z̃| ≤ Ch_Γ², and ũ, z̃ tangential to leading order.
    - Since z̃·n = O(h_Γ³), the third term is O(h_Γ³).
    - Because r ≡ 0, ‖r_h − r̄_h‖_{Γ_h} ≤ Ch_Γ^{-1/2}Y_z, so the fourth term is O(γ^{-1/2}Y·Y_z) = O(h_Γ³).
  - *Leading term.* Since ũ = dWτ + O(d²) and ∫_e d = |e|³/12, we get Λ_h = Σ(|e|³/12)W²(m_e*) + O(h_Γ³), which is at least (c₀²/12)h_Γ²‖W‖² − Ch_Γ³.
  - *GS with G = 0.* Use a GS dual on Z_h. Its leak term vanishes because r ≡ 0, and there are no pressure terms. ∎
  
  **Theorem L2, proof.** Write e = (ũ − u_h′) + δ_R. By Friedrichs and Theorem hc:thm, ‖(μ/h)δ_R − Ũ‖_{L²} ≤ C_F·B_U, and ‖Ũ‖_{S_h} ≤ Ch_Γ. The part ũ − u_h′ is bounded by L1 for GS. The γ-uniform upper bound C_p h_Γ/γ + C(h_Γ^{3/2} + h^k) follows from Theorem rs:unif plus Friedrichs. ∎
  
  **For the main session to note (memory):** the Scott–Vogelius–Nitsche L² results are rate 2 for CNS* and rate exactly 1 for GS, with constant ‖U‖_{L²}/G = 1.2736 for test P.
