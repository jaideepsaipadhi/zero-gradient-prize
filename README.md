# Scott–Vogelius–Nitsche on a polygonally approximated boundary

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23044296.svg)](https://doi.org/10.5281/zenodo.23044296)

**Preprint:** J. S. Padhi, *Scott–Vogelius–Nitsche on a polygonally approximated boundary: convergence, the missing pressure traction, and the penalty threshold*, Zenodo (2026). Published version **v3** (extended 2D paper, 156 pp): [doi:10.5281/zenodo.23044296](https://doi.org/10.5281/zenodo.23044296). Earlier: **v2** (all results for every fixed k ≥ 4): [doi:10.5281/zenodo.23007081](https://doi.org/10.5281/zenodo.23007081). v1 (k = 4): [doi:10.5281/zenodo.23006372](https://doi.org/10.5281/zenodo.23006372).

> **Current state (unreleased drafts, Sept 2026).** The work is now split into three papers, all drafts and not yet on Zenodo (the published version is still v2 above):
>
> | Directory | Paper | Status |
> |---|---|---|
> | `letter/` | C. de Souza Cavalcante and J. S. Padhi, 6-page letter for *Applied Mathematics Letters*: the printed method's failure, the mean-free correction at Scott's rate, unconditional sharpness for every k ≥ 4, the penalty threshold | Joint draft; co-author's sections pending |
> | `paper2d/` | **Journal version of the 2D paper** (50 pp): readable core — history, related work and novelty, printed method, corrected method at Scott's rate with unconditional sharpness, penalty threshold, strong imposition and the zero-gradient lock, selected numerics; everything else is cited to the extended version | Draft for journal submission |
> | `paper/` | Extended 2D version (156 pp): all results with full proofs, the supplement to `paper2d/`; §14 lists what is proved, conditional, numerical or open | Draft |
> | `paper3d/` | 3D companion (59 pp): inscribed polyhedra on Alfeld splits, sharpness for every k ≥ 3, general smooth obstacles, penalty necessity, Navier–Stokes, first 3D numerics (`code3d/`) | Draft |
>
> Working notes, proofs, referee reports and scripts for the latest rounds are in `notes/v6/` (merged into the papers; `notes/v6/MERGE_MANIFEST.txt`, eleven referee reports) and `notes/v7/` (latest round, not yet merged or refereed; see `notes/v7/STATUS.txt`). **The results table below is a snapshot of the v4 draft and is partly superseded** — in particular Theorem D′ is now unconditional for every k ≥ 4, penalty necessity is proved much more broadly, and the 3D results (Theorem I) moved to `paper3d/`. The papers themselves are authoritative.

This repository holds a paper, the code behind it and all the raw results, answering Ridgway Scott's **zero-gradient prize question (PPL 115)**. The method in question uses exactly divergence-free P4 velocities, with the no-slip condition imposed weakly by Nitsche's method on an inscribed polygon, as in Gjerde–Scott (2024). The prize asks whether the error behaves like `h_Γ^{3/2} + h_Ω^k` for general Stokes data, *and if not, why not*.

The core theorems hold for every fixed polynomial degree k ≥ 4 on general meshes, and (v3) for every k ≥ 2 on Clough–Tocher refinements. Most computations use k = 4, as Gjerde–Scott do; v3 adds k = 5, 6 and Clough–Tocher k = 2, 3.

**Short answer:** for the method as printed, no. The printed Nitsche form has no `−pn` traction term, so whenever the pressure varies along the wall, the penalty balances the pressure instead of a traction. The discrete velocity then leaks through the wall with normal velocity `(h/μ)(p − p̄)`. Adding the traction with the wall mean of the discrete pressure removed (`CNS*`, the mean-free correction) restores Scott's rate. That device is **not new**: its velocity coincides with the zero-net-flux formulation of Frachon, Nilsson and Zahedi (BIT 64 (2024) 39). What is new is the analysis.

## Main results

Notation: `h = h_Ω` is the bulk mesh size, `h_Γ` the longest boundary edge, `ρ = h/min|e|`, `γ = μh_Γ/h` and `G = ‖p − p̄‖_{L²(Γ)}`. Section numbers refer to the current draft (a new §8 was inserted, so the penalty section is now §9, Clough–Tocher §10, general domains §11, Navier–Stokes §12, 3D §13, computations §14).

| Result | Statement |
|---|---|
| **Theorem A** | The printed method converges, with energy-norm error `≍ (h/μ)^{1/2} G + O(h_Γ^{3/2})`. Upper and lower bounds are proved. |
| **Theorem C** | The H¹ seminorm error is `≤ C_p h/μ + C(1+γ^{1/2})h_Γ^{3/2} + …`. The rate is 1 (not 1/2) because div v = 0 turns the normal load into a tangential derivative. |
| **Theorem B** | Slip `≍ (h/μ)G` and H¹ error `≍ h/μ`, so Theorem C is sharp. |
| **Corollary B′** | `ũ − u_h = (h/μ)·v_φ + o(·)` for an explicit field `v_φ` with boundary values `≈ −(p − p̄)n`. The slip and energy errors divided by `(h/μ)G` and `(h/μ)^{1/2}G` tend to **exactly 1**. |
| **Theorem B″** (v3, §6.1) | `(μ/h)(ũ − u_h) → U` in H¹, where U is the Stokes flow driven by the leak as Dirichlet data. So the H¹ constant is `K = ‖∇U‖/G`. It depends on the domain: `K = 2.7565157` for test P in the box L = 2.5 (series solution, inside a rigorous bracket [2.0746, 3.0526]), and `K → 3/√7` as L → ∞. Proved for fixed γ; the large-γh_Γ regime is observed only. |
| **Prop. 7.1** | The naive consistent form `+⟨p_h, v·n⟩` is singular, with kernel `(0, 1)` (numerically exactly one-dimensional). Its L²₀ variant gives only `div u_h ≡ c`, `c|Ω_h| = ∫_{∂R} g_I·ν + ∫_{Γ_h} u_h·n`, and `c = 0` iff the naive system is solvable (iff the CNS* pressure has zero Γ_h-mean; the second "iff" assumes unique solvability of the L²₀ variant). Numerically `c ≠ 0`. |
| **Theorem D** (§7, §7.2) | The mean-free correction `CNS*`, which uses `+⟨p_h − p̄_Γ(p_h), v·n⟩` in the momentum row only (cf. Frachon–Nilsson–Zahedi), is well posed and keeps `u_h` exactly divergence-free. With **flux-corrected outer data** `g̃_I = g_I − (m_R/8L²)x`: `‖ũ − u_h‖_{H¹} ≤ C(h_Γ^{3/2} + h_Ω^k)` for general data, **every fixed k ≥ 4** and every γ ≥ γ₂, with C independent of γ and no relation between h_Γ and h_Ω (Theorem 7.22). With the plain interpolant `g_I` the same holds under the proviso h_Γ ≥ h_Ω⁴ (even k) / h_Ω² (odd k). This is Scott's rate. |
| **Theorem D′** (v3, §7.1) | Matching lower bound `‖∇(ũ − u_h)‖ ≥ c‖W‖ h_Γ^{3/2} − Ch_Γ²` (W the wall shear) for CNS* and for GS, so the rate 3/2 is exact. **Unconditional for k ≥ 7. For 4 ≤ k ≤ 6 it is conditional on a local mesh hypothesis (M3)**, which is not proved; it was verified numerically at every boundary vertex of the production meshes (N ≤ 1024). |
| **Prop. E** (§9) | `μ ≥ 4κC_I²ρ` suffices on every mesh. Necessity: on every mesh, `μ ≥ λ_T h/|e|` for each boundary triangle, and in exact rational arithmetic `λ_T > 1` (apex angle ≥ 35°) and `λ_T > 9.3` (apex ≥ 60°), base angles ≥ 5°. So **`μ ≳ ρ` is necessary on every mesh whose shortest boundary edge carries a triangle with apex ≥ 35° and base angles ≥ 5°**; for tall triangles `λ_T < 0` and necessity under (M0)–(M2) alone is open. Also certified on stars near two reference stars. |
| **Theorem R** (v4, §7.2) | The GS H¹ error is bounded uniformly in the penalty: `≤ C_p h_Γ/γ + C(h_Γ^{3/2} + h^k + flux term)`. So `μ ≳ h h_Γ^{−3/2}` restores H¹ rate 3/2 for GS; the energy norm cannot be rescued (lower bound `c‖W‖γ^{1/2}h_Γ^{3/2}`; energy rate ≤ 1 under any power-law penalty). Conditioning grows by `1+γ` (velocity block) and `1+μ/h` (Schur complement). |
| **Theorem P** (v4, §8) | Pressure: CNS* `‖p° − p_h‖_{L²} ≤ C(h_Γ^{3/2} + h_Ω^k)` with flux-corrected data, γ ∈ [γ₂, γ_max] (C depends on γ_max). GS carries an explicit first-layer pressure layer: L² rate exactly 1/2 when G > 0, no pointwise convergence there. Drag/lift/torque (Nitsche flux = Babuška–Miller volume formula): GS error `−(h/μ)K_ψ + O(h^{3/2})` (reciprocity constant), CNS* `O(h_Γ² + ε_h)`. |
| **Theorem F** (v3, §10) | Clough–Tocher refinements: Theorems A–D and Corollary B′ for **every k ≥ 2**, with no angle or singular-vertex hypothesis; penalty necessity recertified exactly on split stars. |
| **Theorem G** (v3, §11) | General smooth curved boundaries: finitely many closed C³ curves (C^{k+3} for results using the leak field), nonconvex and multiply connected, with a nonempty polygonal remainder carrying strong data. The pressure constant must be the global boundary mean. |
| **Theorem H** (v3, §12) | Steady Navier–Stokes, assuming discrete stability of the linearisation (proved for small data): same leak and rates for GS, Scott's rate for CNS*. |
| **Theorem I** (v3, §13) | Three dimensions, assuming a uniform discrete inf-sup constant (known for Alfeld splits, k ≥ 3): Theorems A–D and Corollary B′ with the same rates. Penalty necessity, sharpness and unsplit meshes are open in 3D. |

The computations run up to N = 256 boundary chords.
- The `CNS*` H¹ rate falls monotonically toward 3/2 (1.57 at N = 256 for k = 4; 1.58–1.61 for k = 5, 6), independently of the pressure.
- GS has rate 1 whenever G > 0.
- A pure-pressure test isolates the traction response exactly, and its slip and energy constants reach 0.99–0.999 for every element tested (k = 4, 5, 6 and Clough–Tocher k = 2, 3).
- The penalty threshold is `γ* = μ*h_Γ/h ≈ 18–21` for k = 4, linear in `h/h_Γ` over more than a decade; it is larger for higher k (≈ 33 for k = 5, 48 for k = 6) and on Clough–Tocher meshes.
- (v4) Strong imposition: Scott's `h^{1/2}` is reproduced only on a mesh whose alternate wall vertices lie in two triangles (mesh (a): rate 0.56 at N = 128, vertex-gradient error → 2); on the production meshes (three triangles per wall vertex) strong imposition shows no locking (rates 1.93 → 1.61, error ≈ 2.5× GS). GS/CNS* are unaffected by the mesh change. Observations only.
- (v4) Pressure: CNS* L² rate → 3/2 (1.80 at N = 96, decreasing); GS with G > 0 → 1/2, carried by the first layer.
- (v4) Growing penalty and flux correction behave as predicted (§14.10).

### Corrections in v4
- **Credit for the mean-free device.** Earlier drafts presented CNS* as our construction. Its velocity coincides with that of the zero-net-flux formulation of Frachon, Nilsson and Zahedi, BIT 64 (2024) 39 (§5.1); the paper now calls CNS* "the mean-free correction, cf. Frachon et al." (Remark 7.3) and claims as new only the analysis. Liu–Neilan–Otus (JNM 31 (2023), Remark 3.2), Burman–Hansbo–Larson, Freund–Stenberg, Becker, John et al. (SIAM Rev. 2017) and others are now cited.
- **The "4.5·10⁻³" divergence.** v3 said a bordered solve of the naive form returned discrete divergence ≈ 4.5·10⁻³, "consistent with incompatibility". That value is `|c||Ω_h|^{1/2} = 4.54·10⁻³` of the L²₀ variant (test A, N = 16, μ = 100): a mean-fixing bordered solve computes that variant, whose divergence is the constant c. Prop. 7.1 now states the L²₀ variant; the original runs were not archived, so the identification is by value.
- **Bibliography.** The h^{3/2} Nitsche observation is credited to Scott's prize statement (not Gjerde–Scott 2024); the §13 (3D) Worsey–Farin caveat wrongly attributed to Farrell–Mitchell–Scott Table 2.1 was removed; the Neilan 2015 description now follows FMS (a conditional k ≥ 6 SV result); the Galdi locator for the star-shaped Bogovskiĭ constant is now Lemma III.3.1 (per the audit; unverified against the book); Cavalcante's author field is "Cavalcante, Claudemir" (Zenodo). Unresolved items carry `%% VERIFY` comments in the sources.
- Corollary 6.2 now also requires `γ^{1/2}h_Γ ≪ G`; the text states that the transposed-gradient term is O(h_Γ) pointwise and O(h_Γ^{3/2}) only through odd symmetry (Lemma 3.5) or the penalty.

### Corrections in v3
The v2 numerics section and its penalty table said the meshes have "N equal boundary chords" and "ρ ≈ 3.3". Both are wrong. The chords are the radial projections of equispaced points on the square, so max|e|/min|e| is 1.43 (N = 16) to 1.97 (N = 256). The quantity reported as ρ is `h/max|e| = h/h_Γ`, not `ρ = h/min|e|` (which is 4.70, 5.73, 6.35 at N = 16, 32, 64). The reported γ* = μ*h_Γ/h was correct as defined. v3 relabels it (Section 13; v2 Table 4 and Figure 2 are Table 7 and Figure 3 in v3) and replaces the guess "a limit near 22" by the patch-computation prediction ≈ 21.4. v3 also removes the v2 claim that the code reproduces Cavalcante's benchmark table to 1e-5 on four meshes, because no archived run backs it; only the threshold μ* = 60.343 vs the published 60.338 (N = 20) is kept.

## Layout

```
paper/            main.tex, refs.bib, sections/*.tex (one file per section, 01–13 plus appendices), figures/*.pdf, main.pdf
code/             solver and every script behind a number in the paper
  svn.py            P4 Scott–Vogelius + Nitsche solver (mesh, assembly, GS / CNS* coupling, SuperLU solve, errors, tests A/B/P/C)
  job.py            one study or threshold job per process
  run_all.py        memory/CPU-aware scheduler for the full campaign (the RunPod run)
  report.py         tables for results/pod        report_local.py   tables for results/local (tests P, C)
  plots.py          paper figures                 penalty.py        coercivity threshold μ* by pivot inertia + bisection
  lemma_tests.py    exact discrete dual norms (Lemmas 3.5, 6.3)
  lamstar_exact.py, lamstar_exact_v2.py   exact rational certificates for Prop. 8.1 (stars S3 and S4′)
  lamstar.py        floating-point λ* on other fans     meshcheck.py   mesh sanity check
  -- added in v3 --
  svn_k.py          degree-general solver (P_k, k ≥ 2; standard or Clough–Tocher mesh; SuperLU or PARDISO)
  job_k.py, penalty_k.py, report_v3.py, validate_k.py   jobs, thresholds, tables, validation for results/v3
  lamstar_exact_ct.py   exact certificates on Clough–Tocher stars (§9)
  lamstar_family.py     single-triangle witness, shape-family certificates, patch thresholds (§8.1)
  h1_limit.py           H¹ leak constant: annulus bracket, series solution, FEM comparison (§6.1)
  lower_bound_tests.py  (M3) star constants, assembled test field, dual norms (§7.1)
  plot_leak.py          leak-field figure
  -- added in v4 --
  rescue.py, rescue_run.sh, rescue_report.py   growing penalty, μ-sweeps, flux-corrected data, conditioning (§7.2, §14.10)
  pressure_err.py       pressure errors of GS / CNS* (§14.9)
  strong_bc.py          strongly imposed no-slip, and GS/CNS* on the alternating mesh (a) (§14.8)
  report_v5.py          tables in results/v5/letter_tables.txt (and thresh_rho.jsonl)
  regen_logs.sh     regenerates logs/
  calibration.json  memory model used by run_all.py on the pod
results/
  pod/              the full campaign (32 vCPU / 125 GB RunPod, 36.4 min): study/*.jsonl (74 cases × GS, CNS*, D),
                    thresh/*.json (16 threshold runs), report.txt, run.log, err/ (per-job stdout)
  local/            tests P (pure pressure, 3 μ) and C (non-polynomial pressure), N = 16…96: study/*.jsonl,
                    report_local.txt, localrun.log
  v3/               k = 5, 6 and Clough–Tocher k = 2, 3 runs: study/*.jsonl, thresh/*.json, report_v3.txt,
                    validation.txt, jobs.txt, sched.py (scheduler), run.log, README.md
  v4/rescue/        growing-penalty / μ-sweep / flux-correction / conditioning records and summary.md (§7.2, §14.10)
  v5/               pressure/ and strong/ records, letter_tables.txt, thresh_rho.jsonl, job lists, README.md (§14.8–14.9)
logs/             certificate, lemma-test, v3 and v4 logs, all regenerated by code/regen_logs.sh (v3 adds lamstar_exact_ct,
                    lamstar_family_{bubble,certify,mesh,limit,limitstar,levels}, h1_limit_*, lower_bound_*; v4 adds rescue_*)
verification/     REFEREE_ROUND2.md (independent hostile review + point-by-point response), VALIDATION.md
notes/            proof.md (the Markdown working version of v2), archive/ (earlier drafts and scratch data),
                    v3/ (the v3 section drafts, with a report per draft stating exactly what is proved and what is open),
                    v4/ (rescue.tex, pressure_drag.tex, bib audit), v5/ (literature report, blind verification)
Makefile          pdf, report, figures, smoke, certificates, campaign, local
requirements.txt
LICENSE           (placeholder; no licence chosen yet)
SHA256SUMS
```

## Building the PDF

Building requires `pdflatex`, `bibtex` and `latexmk`.

    make pdf          # = cd paper && latexmk -pdf main.tex

## Rerunning

You need Python ≥ 3.11 and `pip install -r requirements.txt`. Run scripts from `code/`; the Makefile does this for you.

    make smoke          # < 1 min: N = 16 GS/CNS* solves for tests B1 and P (CNS* must return u_h = 0 on P)
    make report         # regenerates results/pod/report.txt, results/local/report_local.txt, results/v3/report_v3.txt,
                        # results/v5/letter_tables.txt and results/v4/rescue/summary.md from the raw records, and diffs
    make figures        # paper/figures/*.pdf
    make certificates   # logs/: all certificate, lemma-test and v3 logs (about an hour)
    make local          # tests P and C, N ≤ 96 (fits in 7 GB)
    make campaign       # full A/B campaign N ≤ 256 + thresholds: ~36 min on 32 cores / 125 GB

New runs write to `results/rerun/` (set `SVN_OUT` to change this), so the shipped results are never overwritten.

A single job:

    cd code && python3 job.py study 64 B1 100      # GS, CNS*, and D records for N = 64, test B1, μ = 100

v3 runs (degree-general solver; needs `pypardiso` for the larger sizes):

    cd code && python3 job_k.py study 5 std 64 B1 300   # k = 5, standard mesh, N = 64, test B1, μ = 300
    cd code && python3 job_k.py thresh 2 ct 32          # coercivity threshold, Clough–Tocher k = 2, N = 32
    cd code && python3 report_v3.py                     # rebuilds results/v3/report_v3.txt

## Runtimes and memory

| Step | Time | Memory |
|---|---|---|
| One study job, N ≤ 96 | seconds to ~2 min | ≤ ~3 GB |
| One study job, N = 192 | ~10 min | tens of GB (one was OOM-killed on the first try on the pod; the retry succeeded) |
| One study job, N = 256 | ~23 min | scheduler reserved up to 55 GB |
| Full campaign (74 study + 16 threshold jobs) | 36.4 min on 32 cores | 125 GB machine |
| Exact certificates (SymPy) | a few minutes each | small |

## Trust base

1. The written arguments of the paper, plus the cited results: Scott–Vogelius (1985) divergence range, the Guzmán–Scott (2019) inf-sup constant, Galdi's Bogovskii results (the locators Lemma III.3.4 / III.3.1 / Thm III.3.3 / Lemma II.1.3 are not verified against the 2nd edition), Argyris interpolation, and the first Korn inequality.
2. Exact rational arithmetic in SymPy for the Proposition 8.1 certificate. Floating point enters only in *finding* the witness; the check `2B − A − 9C > 0` is exact. The same holds for the v3 certificates (Clough–Tocher stars, the shape-family certificates of §8.1, which use exact rational interval arithmetic, and the k = 3 local-surjectivity rank).
3. v3 additionally cites: Clough–Tocher/HCT interpolation (Ciarlet), Scott–Zhang interpolation, Guzmán–Neilan and Arnold–Qin (not load-bearing), and for 3D the published Alfeld-split inf-sup results (Zhang 2005, Guzmán–Neilan 2018), whose constants are taken to depend only on shape regularity without re-audit.
4. The computations of Section 14 are evidence, not part of any proof. The (M3) condition of Theorem D′ for 4 ≤ k ≤ 6 is verified numerically only.

## Validation

Every v2 table in the paper is regenerated byte-for-byte from the shipped raw JSON by `make report`. The v3 tables of Section 14.11 come from `results/v3/report_v3.txt` (`cd code && python3 report_v3.py`); the v3 runs used PARDISO with a 1e-12 pressure-mass regularisation, validated against SuperLU to ≤ 1e-9 relative (`results/v3/validation.txt`). The solver reproduces the pod records from this layout. An independent second-round referee review and the response to each of its points are in `verification/`. See `verification/VALIDATION.md`.

## Use of generative AI

Anthropic's Claude Opus 5.5 was used as a research-support tool. It helped develop and check the proofs, implement, run and validate the code and certificates, review the arguments as a referee, and prepare the manuscript and this repository. The author directed the project and is responsible for all statements, computations, references and conclusions. The declaration is also at the end of the paper.

## Licence

TODO: no licence has been chosen. See `LICENSE`.

## Open items before submission

- Choose a licence and a venue.
- `refs.bib` (merged in v4 from `notes/v4/refs_verified.bib` and `notes/v5/refs_letter.bib`): still unverified are the Ciarlet §6.1 locator, the content of Neilan 2015, the Galdi theorem numbers, and several v5 entries marked `% VERIFY P/S` in the .bib (Freund–Stenberg, Becker, Babuška 1973 pages, Berger–Scott–Strang, Scott 1975).
- v3 open mathematics: (M3) for 4 ≤ k ≤ 6 on general meshes; convergence to the H¹ constant for γh_Γ → ∞; a uniform necessity constant for every nonsingular boundary star; Navier–Stokes beyond small data; penalty necessity and sharpness in 3D; general domains with Σ = ∅. See the reports in `notes/v3/`.
- v3 numerics caveats: μ = 100 is below the coercivity threshold for k = 5 (N ≥ 32), k = 6 and Clough–Tocher k = 3 (μ = 300/1000 used instead); Clough–Tocher k = 2, 3 are pre-asymptotic on test B1.
- Done in v4: the flux-corrected `g̃_I` (Definition 2.2) removes the flux term and the `h_Γ ≥ h_Ω⁴` proviso (Theorem 7.22).
- v4 open items: resolve the `%% VERIFY` comments (Galdi locators, Neilan 2015 content, Freund–Stenberg/Becker content, Gjerde–Scott 2022 slip form, Gjerde–Scott equation/table numbers used in §8, Cavalcante's μ* and scope, the original 4.5e-3 runs); decide whether to acknowledge the informal overlap in GitHub `woahwhattheheck/commons` PR #15236 (see notes/v5/lit_report.md); a lower bound for the CNS* pressure; the GS drag leading term for Navier–Stokes.
