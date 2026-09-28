# v3: higher degree (k = 5, 6) and Clough–Tocher (k = 2, 3) runs

Code: `code/svn_k.py` (degree-general solver: P_k velocity / P_{k-1}-disc pressure, `std` or `ct` mesh),
`code/job_k.py` (one case per process), `code/penalty_k.py` (coercivity threshold), `code/report_v3.py`
(tables -> `report_v3.txt`), `code/validate_k.py` (k = 4 reproduction check -> `validation.txt`).
Jobs: `jobs.txt`, run by the memory-aware `sched.py` (log in `run.log`). All study runs used
`SVN_SOLVER=pardiso` (see "Solver" below).

The `ct` mesh is the barycentric (Alfeld / Clough–Tocher) split of the `svn.make_mesh` macro mesh; the
Γ_h edges, the square and h = h_Ω (the longest macro edge survives in a micro-triangle) are unchanged,
so the h columns of `std` and `ct` tables are identical.

## Validation
- k = 4: `svn_k` is bit-identical to `svn.py` run on this machine (B1 and P1, N = 16, 24, μ = 100,
  θ = 0, 1). Against the stored JSONL: identical for `results/local` (P1); ≤ 2.6e-10 relative for
  `results/pod` (B1; different machine/BLAS).
- `penalty_k` at k = 4, N = 16 gives μ* = 59.271 (pod: 59.27102); μ* is independent of r = 1e5…1e7.
- Solver: SuperLU (`solve_lu`, the reference) runs out of memory well before the target sizes, so the
  studies use MKL PARDISO (`solve_pardiso`). The saddle matrix is ill-conditioned (an earlier version of this
  note called it exactly singular; a referee check found B full row rank on the free dofs, smallest relative
  singular values 1.7e-2 … 1.3e-4); plain PARDISO + refinement diverged for θ = 0 (errors ~1e+60 at k = 6,
  N = 32). The pressure block is therefore regularised by 1e-12·M. It agrees with
  SuperLU to ≤ 8e-10 (k = 4) and ≤ 1e-9 (k = 6, N = 32) relative. Side effect: CNS* on test P returns
  |u_h|_1 ≈ 4e-11 instead of the SuperLU round-off 3e-14; relres (of the unregularised system) ≈ 1e-12
  on test P, ≈ 1e-15 otherwise.
- Health: 390 records, max divres 1.0e-9, max relres 2.2e-12; nothing flagged.

## Coercivity thresholds (the important caveat)
γ* = μ* h_Γ/h_Ω, with h_Ω/h_Γ ≈ 3.3–3.45 (the report prints this ratio as "rho"; it is not the paper's ρ = h/min|e|):

| mesh, k | N=8 | N=16 | N=32 | N=64 | N=128 |
|---|---|---|---|---|---|
| std k=4 (pod) | 14.7 | 18.1 | 19.9 | 20.8 | |
| std k=5 | 25.7 | 30.3 | 32.1 | 33.0 | |
| std k=6 | | 45.2 | 47.4 | 48.3 | |
| ct k=2 | 5.7 | 9.1 | 10.5 | 11.2 | 11.5 |
| ct k=3 | | 31.6 | 33.9 | 34.6 | |
| ct k=4 | | 61.9 | 65.5 | | |

μ* for k = 5 is 99.2 (N = 16) → 112.3 (N = 64); k = 6: 148 → 164; ct k = 3: 104 → 118. **So the
requested μ = 100 is below the threshold for std k = 5 (N ≥ 32), std k = 6 (all N) and ct k = 3 (all
N)**. Those μ = 100 runs are kept (flagged `<<<` in the report) and are erratic, as they should be:
e.g. std k = 5 B1 CNS* H1 = 0.050, 0.138, 0.449, 0.077, 3.69 at N = 24, 32, 48, 64, 96; std k = 6 B1 GS
jumps to 0.65 at N = 32. To test the theory those cases were rerun at μ = 300 and/or μ = 1000 (above
threshold). γ* grows roughly like k(k+1) on `std`; on `ct` it is 1.5–3x larger at equal k (the
micro-triangles on Γ_h are ~3x thinner, while h in μ/h is the macro h). γ* still drifts up with N
(+7% to +24% from N = 16 to 64) and has not visibly saturated.

## What matches the theory (finest level; all μ above μ*)
- **Test P, GS**: slip·μ/(hG) → 1 and energy·(μ/h)^{1/2}/G → 1 in every case: std k = 5 μ = 1000
  N = 128: 0.9983 / 0.9994; std k = 6 μ = 1000 N = 96: 0.9976 / 0.9993; ct k = 2 μ = 100 N = 256:
  0.9928 / 0.9974; ct k = 2 μ = 1000 N = 128: 0.9982 / 0.9992; ct k = 3 μ = 1000 N = 128:
  0.9983 / 0.9996. Rates H1 1.00, energy 0.50, slip 1.00. H1·μ/(hG) → ≈ 2.74–2.77 (same constant as
  k = 4).
- **Test B1, D = GS − CNS***: same normalised ratios, slip 0.986–1.003, energy 0.991–1.006; D H1 rate
  1.01–1.15 → 1, energy rate 0.49–0.52.
- **GS H1 rate → 1 when G > 0**: visible where the D term is not swamped: std k = 5 B1 μ = 300 GS
  H1 rate 1.63, 1.58, 1.48, 1.38, 1.28, 1.20 (N = 24…128) and GS energy rate 0.76; std k = 6 μ = 300
  similarly 1.27 at N = 96. At μ = 1000 the crossover is later (GS H1 rate 1.53 at N = 128, k = 5),
  exactly as for k = 4 in `results/pod`.
- **CNS* H1 rate → 3/2**: std k = 5 B1 μ = 300 → 1.70, 1.68, 1.64, 1.61, **1.59** (N = 128), energy
  1.53; k = 5 μ = 1000 1.58; k = 6 μ = 300/1000 1.61 at N = 96. The k = 5, 6 CNS* numbers are almost
  the same as k = 4 (geometric term dominates, as predicted; e.g. k = 4/5/6 CNS* H1 at N = 64, μ = 1000:
  1.05e-2 / 9.02e-3 / 8.47e-3). ct k = 3 test A (μ = 1000): H1 rate **1.56**, energy 1.55 at N = 128.
- **CNS* on P returns u_h ≈ 0** for std k = 5, 6 (|u_h|_1 ≈ 4e-11, the regularisation floor, at every
  N and μ ≥ μ*). Since test B_a has the same u as B1 and f_a = f_1 + (a−1)∇p = f_1 + (a−1)f_P, linearity makes
  this exactly the statement that CNS* is independent of the pressure amplitude (u_h(B_a) = u_h(B1)).

## What does not (or only partly) match
1. **μ = 100 is sub-threshold** for k = 5, 6 and ct k = 3 (above). Not a failure of the theory, but the
   requested parameter set is outside it.
2. **Clough–Tocher k = 2, 3: CNS* on test P is NOT zero.** p = x²y − y + x/2 is cubic, not in
   P_{k−1} for k ≤ 3, so the consistent coupling uses p_h ≠ p on Γ_h. The spurious velocity is small and
   converges at rate ≈ k + 1/2: ct k = 2 H1 1.4e-6 at N = 256 (rate 2.56), ct k = 3 H1 4.7e-9 at
   N = 128 (rate 3.53); slip rates 3.0 and 4.05. So pressure-independence of CNS* holds only up to an
   O(h^{k+1/2})·|p| term for k ≤ 3; it is subdominant to h^{3/2} but scales with the pressure amplitude
   (for B100 it would be 99x larger).
3. **Clough–Tocher k = 2, test B1: the 3/2 geometric rate is not observable.** The H1 error is dominated
   by P2 interpolation of this steep velocity (u ~ r^6 on the square of side 5): H1 = 0.41 at N = 256
   (vs 2e-3 for P4), rate 1.91 rising towards 2; the crossover to h^{3/2} is far outside reach. The D
   part (the Nitsche/pressure effect) is perfectly on theory (ratios 1.0009 / 1.0056 at N = 256). On
   test A (bounded velocity) ct k = 2 gives H1 rate 1.78 and energy rate 1.55 at N = 256 — the energy
   norm is at 3/2, H1 is still pre-asymptotic.
4. **ct k = 3 test B1**: likewise interpolation-dominated, CNS* H1 rate 2.83 (≈ h^3) at N = 128; only
   test A shows the 3/2 (above).
5. **H1·μ/(hG) for D on test B** is still drifting down at the finest levels (k = 5 μ = 1000: 3.03,
   k = 6: 3.05; μ = 300: 2.82) toward the test-P value ≈ 2.77; not converged, but consistent.
6. μ* has not saturated in N (ct k = 2: μ* +4.5%, γ* +2.9% from N = 64 to 128).

## Not run (memory)
The shared 5.8 GiB bash memory cgroup (shared with another agent's jobs) limited sizes: std k = 6
N = 128 (≈ 6 GB peak estimated) and ct k = 3 N = 192 (≈ 4.6 GB) were not run; std k = 6 N = 96 B1
μ = 1000 was OOM-killed twice (the other agent's usage) and completed on retry. Peak RSS per run is in
the `maxrss_gb` field (std k = 5 N = 128: 3.6 GB; std k = 6 N = 96: 3.5 GB; ct k = 2 N = 256: 2.7 GB).
