# Scott–Vogelius–Nitsche on a polygonally approximated boundary

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23007081.svg)](https://doi.org/10.5281/zenodo.23007081)

**Preprint:** J. S. Padhi, *Scott–Vogelius–Nitsche on a polygonally approximated boundary: convergence, the missing pressure traction, and the penalty threshold*, Zenodo (2026). Published version **v2** (all results for every fixed k ≥ 4): [doi:10.5281/zenodo.23007081](https://doi.org/10.5281/zenodo.23007081). v1 (k = 4): [doi:10.5281/zenodo.23006372](https://doi.org/10.5281/zenodo.23006372).

> **This branch (`v3-draft`) is an unreleased draft of v3.** `paper/main.tex` and `paper/main.pdf` here are the v3 draft, not the published v2. v3 adds: every k ≥ 2 on Clough–Tocher refinements, general smooth curved domains, the matching H¹ lower bound (sharpness), the identification of the H¹ leak constant, steady Navier–Stokes, three dimensions, broader penalty-necessity certificates, and new numerics for k = 5, 6 and Clough–Tocher k = 2, 3. It also corrects a labelling error in the v2 numerics (see *Corrections* below). The draft source notes and their reports are in `notes/v3/`.

This repository holds a paper, the code behind it and all the raw results, answering Ridgway Scott's **zero-gradient prize question (PPL 115)**. The method in question uses exactly divergence-free P4 velocities, with the no-slip condition imposed weakly by Nitsche's method on an inscribed polygon, as in Gjerde–Scott (2024). The prize asks whether the error behaves like `h_Γ^{3/2} + h_Ω^k` for general Stokes data, *and if not, why not*.

The core theorems hold for every fixed polynomial degree k ≥ 4 on general meshes, and (v3) for every k ≥ 2 on Clough–Tocher refinements. Most computations use k = 4, as Gjerde–Scott do; v3 adds k = 5, 6 and Clough–Tocher k = 2, 3.

**Short answer:** for the method as printed, no. The printed Nitsche form has no `−pn` traction term, so whenever the pressure varies along the wall, the penalty balances the pressure instead of a traction. The discrete velocity then leaks through the wall with normal velocity `(h/μ)(p − p̄)`. A mean-corrected pressure-consistent variant restores Scott's rate.

## Main results

Notation: `h = h_Ω` is the bulk mesh size, `h_Γ` the longest boundary edge, `ρ = h/min|e|`, `γ = μh_Γ/h` and `G = ‖p − p̄‖_{L²(Γ)}`. Section numbers refer to the v3 draft.

| Result | Statement |
|---|---|
| **Theorem A** | The printed method converges, with energy-norm error `≍ (h/μ)^{1/2} G + O(h_Γ^{3/2})`. Upper and lower bounds are proved. |
| **Theorem C** | The H¹ seminorm error is `≤ C_p h/μ + C(1+γ^{1/2})h_Γ^{3/2} + …`. The rate is 1 (not 1/2) because div v = 0 turns the normal load into a tangential derivative. |
| **Theorem B** | Slip `≍ (h/μ)G` and H¹ error `≍ h/μ`, so Theorem C is sharp. |
| **Corollary B′** | `ũ − u_h = (h/μ)·v_φ + o(·)` for an explicit field `v_φ` with boundary values `≈ −(p − p̄)n`. The slip and energy errors divided by `(h/μ)G` and `(h/μ)^{1/2}G` tend to **exactly 1**. |
| **Theorem B″** (v3, §6.1) | `(μ/h)(ũ − u_h) → U` in H¹, where U is the Stokes flow driven by the leak as Dirichlet data. So the H¹ constant is `K = ‖∇U‖/G`. It depends on the domain: `K = 2.7565157` for test P in the box L = 2.5 (series solution, inside a rigorous bracket [2.0746, 3.0526]), and `K → 3/√7` as L → ∞. Proved for fixed γ; the large-γh_Γ regime is observed only. |
| **Prop. 7.1** | The naive consistent form `+⟨p_h, v·n⟩` is singular, with kernel `(0, 1)`. |
| **Theorem D** | The corrected method `CNS*`, which uses `+⟨p_h − p̄_Γ(p_h), v·n⟩`, is well posed and has `‖ũ − u_h‖_{H¹} ≤ C(h_Γ^{3/2} + h_Ω^k)` for general data and **every fixed k ≥ 4** (with boundary edges not too short: h_Γ ≥ h_Ω⁴ for even k). This is Scott's rate. |
| **Theorem D′** (v3, §7.1) | Matching lower bound `‖∇(ũ − u_h)‖ ≥ c‖W‖ h_Γ^{3/2} − Ch_Γ²` (W the wall shear) for CNS* and for GS, so the rate 3/2 is exact. **Unconditional for k ≥ 7. For 4 ≤ k ≤ 6 it is conditional on a local mesh hypothesis (M3)**, which is not proved; it was verified numerically at every boundary vertex of the production meshes (N ≤ 1024). |
| **Prop. E** (§8) | The penalty threshold scales like `ρ = h/min|e|`. `μ ≥ 4κC_I²ρ` suffices. Necessity (v3): on every mesh, `μ ≥ λ_T h/|e|` for each boundary triangle, and in exact rational arithmetic `λ_T > 1` (apex angle ≥ 35°) and `λ_T > 9.3` (apex ≥ 60°), base angles ≥ 5°. Also certified on stars near two reference stars. |
| **Theorem F** (v3, §9) | Clough–Tocher refinements: Theorems A–D and Corollary B′ for **every k ≥ 2**, with no angle or singular-vertex hypothesis; penalty necessity recertified exactly on split stars. |
| **Theorem G** (v3, §10) | General smooth curved boundaries: finitely many closed C³ curves (C^{k+3} for results using the leak field), nonconvex and multiply connected, with a nonempty polygonal remainder carrying strong data. The pressure constant must be the global boundary mean. |
| **Theorem H** (v3, §11) | Steady Navier–Stokes, assuming discrete stability of the linearisation (proved for small data): same leak and rates for GS, Scott's rate for CNS*. |
| **Theorem I** (v3, §12) | Three dimensions, assuming a uniform discrete inf-sup constant (known for Alfeld splits, k ≥ 3): Theorems A–D and Corollary B′ with the same rates. Penalty necessity, sharpness and unsplit meshes are open in 3D. |

The computations run up to N = 256 boundary chords.
- The `CNS*` H¹ rate falls monotonically toward 3/2 (1.57 at N = 256 for k = 4; 1.58–1.61 for k = 5, 6), independently of the pressure.
- GS has rate 1 whenever G > 0.
- A pure-pressure test isolates the traction response exactly, and its slip and energy constants reach 0.99–0.999 for every element tested (k = 4, 5, 6 and Clough–Tocher k = 2, 3).
- The penalty threshold is `γ* = μ*h_Γ/h ≈ 18–21` for k = 4, linear in `h/h_Γ` over more than a decade; it is larger for higher k (≈ 33 for k = 5, 48 for k = 6) and on Clough–Tocher meshes.

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
  regen_logs.sh     regenerates logs/
  calibration.json  memory model used by run_all.py on the pod
results/
  pod/              the full campaign (32 vCPU / 125 GB RunPod, 36.4 min): study/*.jsonl (74 cases × GS, CNS*, D),
                    thresh/*.json (16 threshold runs), report.txt, run.log, err/ (per-job stdout)
  local/            tests P (pure pressure, 3 μ) and C (non-polynomial pressure), N = 16…96: study/*.jsonl,
                    report_local.txt, localrun.log
  v3/               k = 5, 6 and Clough–Tocher k = 2, 3 runs: study/*.jsonl, thresh/*.json, report_v3.txt,
                    validation.txt, jobs.txt, sched.py (scheduler), run.log, README.md
logs/             certificate, lemma-test and v3 logs, all regenerated by code/regen_logs.sh (v3 adds lamstar_exact_ct,
                    lamstar_family_{bubble,certify,mesh,limit,limitstar,levels}, h1_limit_*, lower_bound_*)
verification/     REFEREE_ROUND2.md (independent hostile review + point-by-point response), VALIDATION.md
notes/            proof.md (the Markdown working version of v2), archive/ (earlier drafts and scratch data),
                    v3/ (the v3 section drafts, with a report per draft stating exactly what is proved and what is open)
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
    make report         # regenerates results/pod/report.txt, results/local/report_local.txt and results/v3/report_v3.txt from the raw JSON, and diffs
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

1. The written arguments of the paper, plus the cited results: Scott–Vogelius (1985) divergence range, the Guzmán–Scott (2019) inf-sup constant, Galdi's Lemma III.3.4 (Bogovskii on unions of star-shaped sets), Argyris interpolation, and the first Korn inequality.
2. Exact rational arithmetic in SymPy for the Proposition 8.1 certificate. Floating point enters only in *finding* the witness; the check `2B − A − 9C > 0` is exact. The same holds for the v3 certificates (Clough–Tocher stars, the shape-family certificates of §8.1, which use exact rational interval arithmetic, and the k = 3 local-surjectivity rank).
3. v3 additionally cites: Clough–Tocher/HCT interpolation (Ciarlet), Scott–Zhang interpolation, Guzmán–Neilan and Arnold–Qin (not load-bearing), and for 3D the published Alfeld-split inf-sup results (Zhang 2005, Guzmán–Neilan 2018), whose constants are taken to depend only on shape regularity without re-audit.
4. The computations of Section 13 are evidence, not part of any proof. The (M3) condition of Theorem D′ for 4 ≤ k ≤ 6 is verified numerically only.

## Validation

Every v2 table in the paper is regenerated byte-for-byte from the shipped raw JSON by `make report`. The v3 tables of Section 13.8 come from `results/v3/report_v3.txt` (`cd code && python3 report_v3.py`); the v3 runs used PARDISO with a 1e-12 pressure-mass regularisation, validated against SuperLU to ≤ 1e-9 relative (`results/v3/validation.txt`). The solver reproduces the pod records from this layout. An independent second-round referee review and the response to each of its points are in `verification/`. See `verification/VALIDATION.md`.

## Use of generative AI

Anthropic's Claude Opus 5.5 was used as a research-support tool. It helped develop and check the proofs, implement, run and validate the code and certificates, review the arguments as a referee, and prepare the manuscript and this repository. The author directed the project and is responsible for all statements, computations, references and conclusions. The declaration is also at the end of the paper.

## Licence

TODO: no licence has been chosen. See `LICENSE`.

## Open items before submission

- Choose a licence and a venue.
- `refs.bib`: Cavalcante's preprint has no arXiv number yet, and the edition and page numbers of Galdi's lemma should be double-checked against the book. Several v3 entries carry `note = {TODO verify ...}` (Arnold–Qin, Clough–Tocher, the Ciarlet section number, Zhang 2005/2011 issue numbers, Neilan 2015, Farrell–Mitchell–Scott volume/pages, Guzmán–Lischke–Neilan).
- v3 open mathematics: (M3) for 4 ≤ k ≤ 6 on general meshes; convergence to the H¹ constant for γh_Γ → ∞; a uniform necessity constant for every nonsingular boundary star; Navier–Stokes beyond small data; penalty necessity and sharpness in 3D; general domains with Σ = ∅. See the reports in `notes/v3/`.
- v3 numerics caveats: μ = 100 is below the coercivity threshold for k = 5 (N ≥ 32), k = 6 and Clough–Tocher k = 3 (μ = 300/1000 used instead); Clough–Tocher k = 2, 3 are pre-asymptotic on test B1.
- Optional: a flux-corrected `g_I` removes the flux term in Lemma 4.6 and the `h_Γ ≥ h_Ω⁴` proviso in Theorem D.
