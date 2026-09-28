# v5: numerics for the Applied Mathematics Letters version

All runs use k = 4 (P4 Scott–Vogelius velocity, P3-disc pressure), μ = 100, and the meshes of `svn.make_mesh` (L = 2.5), except where mesh (a) is stated. No existing file was modified.

## Files

| file | contents |
|---|---|
| `code/pressure_err.py` | Pressure errors for GS (θ = 0) and CNS* (θ = 1). Both are solved on one assembly, and the velocity norms are recorded as well. Writes `pressure/N{N}_{KIND}_mu100.jsonl`. |
| `code/strong_bc.py` | Strongly imposed no-slip on Γ_h, with no Nitsche terms. It reuses svn's mesh, space, volume assembly and error routines. It also runs GS/CNS* on the alternative mesh (a). Writes `strong/*.jsonl`. |
| `code/report_v5.py` | `python3 report_v5.py thresh` recomputes h/h_Γ and the true ρ for every `results/pod/thresh` run and writes `thresh_rho.jsonl`. `python3 report_v5.py` regenerates `letter_tables.txt` from the JSONL alone. |
| `letter_tables.txt` | Health summary, plain-text tables T1, T1′, T2, T3 and T4, then LaTeX booktabs tables `tab:A`, `tab:B100`, `tab:thresh`, `tab:pB1`, `tab:pC` and `tab:pA`. |
| `jobs*.txt`, `run.log`, `err.log` | Job lists and the log (see "Run history"). |

The H¹ data in T1 (GS and CNS*) and T2 come from `results/pod/study` and were not rerun. Velocity cross-check: every v5 pressure run (A, B1 and B100 against pod; C against `results/local`) reproduces the stored H¹ to a relative difference of at most 2.1e-7.

## Pressure normalisation (question (1))

- **p̄_Γ normalisation.** p* = p̃ − p̄_{Γ_h}(p̃), with the Γ_h mean computed by 10-point Gauss quadrature on the chords. With this choice, `pL2` equals `pL2opt` (optimal constant removed) to 3 digits for both methods and all tests. The Ω_h mean of p* − p_h is 1e-5 to 1e-6 at N = 96.
- **So the discrete pressures are already p̄_Γ-normalised, up to a vanishing constant**, for GS as well as CNS*.
- **The discrete Γ_h mean is not zero.** λ_h = p̄_{Γ_h}(p_h) → 0 at rate ≈ 1 (CNS* B1: −4.6e-3 → −9.7e-4). Forcing p̄_{Γ_h}(p_h) = 0 (`pL2G`) is slightly worse, because it imports the first-layer error into the constant. Use `pL2`.
- **CNS* does not depend on the pressure data.** CNS* pressure errors for B1, B100 and C agree to 1e-6 at N = 16 and 3e-11 at N = 96. This is expected:
  - CNS* reproduces a polynomial pressure exactly (test P);
  - C differs from B1 only through p* − π_h p*, which is O(h⁴).

## Pressure results (T4, N = 16 … 96, rates in h)

- **GS on B1, B100 and C** (G > 0):
  - L²(Ω_h) rate 0.60, 0.62 and 0.63 at N = 96, decreasing towards 1/2;
  - first layer 0.49, 0.48 and 0.50;
  - off the first layer 0.89, 0.96 and 1.02, i.e. → 1;
  - L^∞ on the first layer does not converge (0.78 for B1, 77 for B100, 0.66 for C).
  - B100 has pressure errors of 9.2 at N = 96, which is 100 times B1's layer error.
  - This is as predicted in notes/v4/pressure_drag_report.md.
- **CNS\*** (same numbers for B1, B100 and C):
  - L² rates 3.0, 3.0, 2.6, 2.1 and 1.80, decreasing towards 3/2 (pre-asymptotic, as predicted);
  - first layer ≈ 1.6;
  - L^∞(L₁) rate → 1 (0.98).
- **Test A** (p ≡ 0, G = 0):
  - GS and CNS* both have L² rate 1.68 / 1.65 at N = 96 and are still decreasing. GS here behaves like CNS*, as predicted for G = 0.
  - The GS first layer converges (1.67), unlike the G > 0 tests.

## Strong imposition (T1, T1′): anomaly and resolution

- **Mesh (s): no h^{1/2}.** Strong u_h = 0 at all Γ_h nodes on the mesh of GS/CNS* (mesh (s)) does **not** show h^{1/2}.
  - H¹ rates are 1.93 → 1.61 (N = 128), tracking GS/CNS* at about 2.5× their error.
  - The vertex-gradient defect `vgrad` decays like h (1.28 → 0.16).
  - Reason: on `svn.make_mesh` every polygon vertex lies in exactly 3 triangles. A divergence-free, continuous, piecewise-P4 field vanishing on both Γ_h edges then has a 1-parameter family of vertex gradients: 12 gradient unknowns against 11 conditions (2 per boundary edge, 2 per interior edge, 3 traces). So the zero-gradient lock is not triggered.
- **Mesh (a): h^{1/2} reproduced.** Mesh (a) alternates the first-layer diagonals, so polygon vertices lie alternately in 2 and 4 triangles (`strong_bc.make_mesh_alt`). At a vertex in 2 triangles, the 8 conditions fix the 8 gradient entries to 0.
  - H¹ rate 0.75 → 0.56 (N = 128), decreasing to 1/2.
  - vgrad → 1.99, i.e. |∂_n u_τ| = 2 at r = 1: the gradient error at the vertex is O(1).
  - ‖p_h‖_{L²} grows like h^{−1/2} (14 → 36 with p ≡ 0).
  - GS and CNS* on the same mesh (a) are unaffected: H¹ 6.00e-3 / 7.17e-3 at N = 128, versus 6.11e-3 / 6.96e-3 on mesh (s).
- **Interpolated control.** Strong imposition of u_h = ũ at the Γ_h nodes (`interp`, mesh (s)) converges at rate ≈ 4. The strong defect is therefore entirely the u = 0 transfer.
- **T1 has 4 error columns, not 3:** strong on (a) and strong on (s), then GS and CNS* (pod data, mesh (s), N ≤ 256). For a 3-column version, keep the strong (a) column. The mesh-(a) GS/CNS* numbers in T1′ show that swapping the mesh barely changes those two columns (≤ 6% difference at N = 16 and ≤ 3% at N = 128).
- **Size limit.** The strong runs stop at N = 128 because N = 192 would exceed 3 GB.

## Penalty thresholds (T3)

- **Source.** μ* comes from `results/pod/thresh`.
- **Definitions.** h/h_Γ = h/max|e|; the true ρ = h/min|e| over the Γ_h edges, recomputed from `svn.make_mesh` (the recomputed h and h_Γ match the stored values exactly).
- **Part (a).** At fixed h/h_Γ ≈ 3.3, ρ grows 3.41 → 6.35, because max|e|/min|e| grows 1.00 → 1.87 with the radial projection of the square.
- **Part (b).** The ρ-sweep gives ρ = 4.70 … 63.0 (N = 16) and 5.73 … 71.5 (N = 32).
- **Scaling.** γ* = μ* h_Γ/h stays within 18.1–19.0 (N = 16) and 19.7–20.1 (N = 32). μ*/ρ is also flat in the sweep: 12.6–13.2 and 11.4–11.7.
  - At fixed ρ-profile, μ* ∝ h/h_Γ and μ* ∝ ρ cannot be separated.
  - Under refinement, μ*/ρ drifts down (14.7 → 11.1) while γ* drifts up (14.7 → 20.8).
- **Paper check.** The paper's caption values ρ = 4.70, 5.73 and 6.35 are reproduced. The h/h_Γ entry at N = 20 is 3.1648 and is printed 3.16 here; the paper prints 3.17.

## Solver health

- **Records.** 83 v5 records (48 pressure, 35 strong/mesh-(a)). None has relres > 1e-8 or divres > 1e-6. Max relres is 1.3e-13.
- **Largest divres.** Max divres is 9.6e-7 (strong, mesh (a), N = 128, PARDISO). The mesh-(a) strong system is nearly singular in the pressure (the lock), and the 1e-12·M regularisation shows up there. At N = 64, SuperLU gives the same H¹ and ‖p_h‖ to 5 digits, with divres 1.5e-10.
- **Solvers.** N ≤ 64 use SuperLU (`svn.solve_lu`). N = 96 pressure runs, strong runs with N ≥ 64, and mesh-(a) runs with N ≥ 64 use PARDISO, with the `svn_k.solve_pardiso` algorithm (pressure block + 1e-12·M, refinement, relres of the unregularised system) transcribed onto svn's assembly.
- **PARDISO agreement.** PARDISO matches SuperLU at N = 64 (B1) to ≤ 4e-9 relative in H¹ and the pressure norms.
- **Memory.** Peak RSS was ≤ 1.9 GB.
- **pod data.** The 34 pod records used (A and B100, μ = 100) are all clean: relres ≤ 3.8e-15 and divres ≤ 5.1e-11.

## Run history

- **`jobs.txt`: the N = 96 SuperLU pressure jobs were killed on purpose.** They were projected at about 2.9 GB each, running two at a time next to other agents. They were replaced by PARDISO jobs in `jobs2.txt`.
- **`strong/*` for N ≤ 48 were rerun** (`jobs2.txt`) after `strong_bc.py` switched from a dense mean-zero bordering row to pinning one pressure unknown. The bordering row caused heavy SuperLU fill. The results were identical.
- **The strong N = 128 job was OOM-killed once** (`err.log`, while running concurrently) and completed on the rerun in `jobs4.txt`.
