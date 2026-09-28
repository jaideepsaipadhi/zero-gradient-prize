# End-to-end validation

Checks run on the packaged repository, from `code/`, after moving everything into this layout.

| Check | How | Result |
|---|---|---|
| Solver health, production campaign | `results/pod/report.txt` header; max over 148 solves | relres ≤ 3.8e-15, discrete divergence ≤ 7e-11 (median 9e-12); 0 records over thresholds |
| Tables regenerate from raw results | `make report` | `report.txt`, `report_local.txt` and (v3) `results/v3/report_v3.txt` reproduced byte-for-byte |
| Solver reproduces from this layout | `make smoke` (N = 16, B1 and P1, μ = 100) | GS H¹ = 1.3423e-01, CNS* H¹ = 1.2882e-01: identical to the pod records |
| CNS* consistency (exactness on pure pressure) | test P, all μ, N ≤ 96 | CNS* velocity ≤ 1e-13 in H¹, i.e. u_h = 0 to round-off |
| Code-to-code | coercivity threshold vs Cavalcante's published value | μ* = 60.343 (`results/pod/thresh/N20_L2.5_m5.json`) vs published 60.338 at N = 20. (An earlier claim of agreement with Cavalcante's Table 1 to 1e-5 on four meshes was removed in v3: no archived run backs it.) |
| Exact certificates | `code/regen_logs.sh` → `logs/lamstar_exact_*.log` (the script now also regenerates all v3 logs) | 𝒮₃: 80 constraints, dim 14, 2B−A−(99/10)C > 0 exactly (quotient 9.978546); 𝒮₄′: 110 constraints, dim 16, 2B−A−9C > 0 exactly (quotient 9.267472). The singular 𝒮₄ (dim 17) run is kept in the log for the record but is not used. |
| Lemma tests | `logs/lemma_tests.log` | normal S on Z_h bounded (rate +0.05); on V_h^R and tangential: −0.52; (∇ũ)ᵀn: 1.57 / 1.56 |
| Independent referee | `verification/REFEREE_ROUND2.md` | 21 items, none fatal; all addressed (table at the end of that file) |

Not rerun from this layout (too heavy for the packaging machine): the full N ≤ 256 campaign (`make campaign`, 32 cores / 125 GB, 36 min).
Its raw outputs, scheduler log and per-job stdout are shipped in `results/pod/`.
