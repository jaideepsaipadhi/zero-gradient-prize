# Literature report for the AML letter (v5), 2026-09-28

Files in this directory:

- `refs_letter.bib`: 36 entries, each with a `% VERIFY` line giving its status and source. It compiles cleanly with `amsplain` alongside `lit_framing.tex`; see the test at the end.
- `lit_framing.tex`: the history, positioning, contribution and scope paragraphs.

No existing file was edited, and nothing was committed.

**Access limits.**
- The bash sandbox cannot reach external hosts; Crossref and arXiv give a 403 at the proxy. Every check below was done with the web search and fetch tools.
- These sites refused access: AMS article PDFs (403), JSTOR (JS wall), ResearchGate (429), Wiley (403), epubs.siam.org (robots.txt) and ADS (robots.txt).
- **One page could not be read.** Springer's page for Dione–Urquiza (Numer. Math., DOI 10.1007/s00211-014-0646-9) was rate-limited (429), and the proxy said not to retry. That paper is therefore **not** in the bib.

## 1. Verification status of the handoff list

| Handoff item | Key | Status | Notes / corrections |
|---|---|---|---|
| Freund & Stenberg 1995 | FreundStenberg1995 | **P** | Entry quoted verbatim from Burman's 2012 SINUM reference list: *Proc. 9th Int. Conf. Finite Elements in Fluids*, M. (Morandi) Cecchi et al., eds., Univ. Padova, 1995, pp. 327–336. The full text was not seen. That it treats **Stokes** rests only on arXiv:2502.09550, which lists [FS95] among Nitsche works for Stokes/NS. |
| Becker 2002, Commun. Numer. Methods Eng. 18 | Becker2002 | **P** | DOI 10.1002/cnm.529 (Wiley URL). **18, 669–680** comes from the Numer. Algorithms 2023 reference list. The issue number is unconfirmed and was omitted. Content not read. |
| Burman & Hansbo (Nitsche/fictitious-domain Stokes) | BurmanHansbo2014 | **V** | M2AN **48(3) 859–874 (2014)**, DOI 10.1051/m2an/2013123 (numdam). |
| Juntunen & Stenberg, Math. Comp. 2009 | JuntunenStenberg2009 | **V** | **78(267) 1353–1374**, DOI 10.1090/S0025-5718-08-02183-2. It treats Poisson with Robin data, **not Stokes**. |
| Bramble–Dupont–Thomée 1972 | BrambleDupontThomee1972 | **V** | Math. Comp. **26(120) 869–879** (DGS 2025 reference list). No DOI added. |
| Burman–Hansbo–Larson, Math. Comp. 87 (2018) 633–657 | BurmanHansboLarson2018 | **V** | Issue **310**, DOI 10.1090/mcom/3240. Scalar problem. |
| cut-FEM boundary value correction for Stokes, arXiv:1801.07463 | BurmanHansboLarson2019Stokes | **V** | Published as an **ENUMATH 2017 chapter**: LNCSE 126, Springer 2019, pp. 183–192, DOI 10.1007/978-3-319-96415-7_15. Authors Burman, Hansbo, Larson. |
| Shifted boundary method (Main & Scovazzi) | MainScovazzi2018 | **V** | JCP **372 (2018) 972–995**, DOI 10.1016/j.jcp.2017.10.026. |
| arXiv:2105.10409 (closest SV work) | LiuNeilanOtus2023 | **V** | **Liu, Neilan, Otus**, JNM **31(2) 105–123 (2023)**, DOI 10.1515/jnma-2021-0125. Full text read (NSF-PAR); see §2. |
| Divergence-free cut-FEM, BIT 2024 (arXiv:2304.14230) | FrachonNilssonZahedi2024 | **V** | Authors are **Frachon, Nilsson, Zahedi**; Hansbo is not an author. BIT **64**, Paper 39 (2024), DOI 10.1007/s10543-024-01040-x. **Key prior art**; see §2. |
| Divergence-free cut-FEM, "M2AN 2022" | LiuNeilanOlshanskii2023 | **V (corrected)** | Liu, Neilan, Olshanskii, M2AN **57(1) 143–165, 2023** (not 2022), DOI 10.1051/m2an/2022072. |
| John–Linke–Merdon–Neilan–Rebholz, SIAM Rev. 59 (2017) | JohnEtAl2017 | **V** | 59(3) 492–544, DOI 10.1137/15M1047696. |
| Neilan's Stokes-complex review (2020) | Neilan2020 | **P** | Contemp. Math. **754**, AMS 2020, **pp. 141–158** (Neilan CV). The editors are Brenner, Shparlinski, Shu and Szyld (AMS page). The chapter DOI was not found. |
| Gjerde–Scott, Math. Comp. 91 (2022), slip | GjerdeScott2022 | **V (bibliographic)**; **content NOT verified** | 91(334) 597–622, DOI **10.1090/mcom/3682** (AMS page and GitHub README), online 5 Nov 2021. The AMS PDF returned 403 and there is no arXiv version, so **it is not known whether it includes the pressure/traction term or discusses exact divergence-freeness**. The abstract only covers normals/tangents and the Babuška–Sapondzhyan paradox. **Read it before submission** if the letter says anything about it. The current draft does not. |
| Dupont–Guzmán–Scott 2025 | DupontGuzmanScott2025 | **V** | JNM 33(1) 55–86. **Poisson only**; Stokes is never mentioned. |
| Eickmann–Scott–Tscherpel, arXiv:2509.17899 | EickmannScottTscherpel2025 | **V** | Submitted 22 Sep 2025, revised 24 Jun 2026. Covers polygonal domains, a flux-corrected Dirichlet interpolant (eq. 3.9) and L²₀ pressures. It does **not** treat Nitsche or curved boundaries. It shows that uncorrected data give ‖div u_h‖ ≈ 4e-5 vs 7e-11 when corrected. This is relevant to our "flux term". |

**Other entries added** (all with VERIFY lines):
- Stenberg 1995 (JCAM 63, DOI verified).
- Burman 2012 penalty-free Nitsche (P: pages from standard citation).
- Boiveau–Burman 2016 (V).
- Babuška 1973 (P: pages from snippet).
- Neilan–Otus 2021 (V).
- Durst–Neilan 2024 (V).
- Burman–Hansbo–Larson SINUM 62(2) 893–918 (2024), the Lagrange-multiplier divergence-free cut-FEM (V).
- Berger–Scott–Strang 1972 and Scott 1975 (S: as printed in GS2024's reference list; these are the scalar h^{3/2} results GS cite).
- `Cavalcante2026` (V).
- `PadhiV2`.

**Correction to `paper/refs.bib`.** The Zenodo record 22262466 lists the author as **"Cavalcante, Claudemir"** (published 3 Sep 2026, v1), not "de Souza Cavalcante". The record describes only the constant-pressure benchmark, with the bound ‖u−u_h‖_{H¹} ≤ C(h_Ω^k + h_Γ^{3/2}).

**Carried over from `notes/v4/bib_audit.md` (still applies).**
- The h^{3/2} observation must be credited to **ScottPrize**, not GS2024. The draft does this.
- GS2024 (arXiv) was checked again: (27) contains no pressure term. GS2024 cites **only Stenberg 1995** for Nitsche; it does not cite Freund–Stenberg, Juntunen–Stenberg, Becker or Burman–Hansbo. It never states that the method is inconsistent for Stokes.

## 2. Novelty findings, claim by claim

### Prior-art hits (honest summary)

**1. Frachon–Nilsson–Zahedi, BIT 64:39 (2024), §3 and §5.1. This is the most important hit.**
- **Their setting.** Divergence-free cut-FEM (BDM, H(div)-conforming) with a nonsymmetric Nitsche Dirichlet condition. The **boundary pressure term (u_h·n, p_h) appears in the momentum form only**, and the continuity form is the plain B₀. Their reason: *"We use B₀ since we do not want to perturb the divergence condition."*
- **§5.1, eqs (5.7)–(5.9).** When the pressure mean is fixed by a Lagrange multiplier α, **|Ω|α = −∫_{∂Ω} u_h·n and div u_h = −α**. With weak imposition the net flux is not zero, so exact incompressibility is lost.
- **Their fix.** Test the momentum equation with **V⁰ = {v : ∫_{∂Ω} v·n = 0}**, use the full pressure trial space, and fix the pressure mean separately.
- **Consequence for (iii) (proof sketch).**
  - Our setting uses test space V_R, with the outer boundary strong, so ∫_{∂Ω_h} v·n = ∫_{Γ_h} v·n.
  - For v with zero Γ_h-flux, ⟨p̄_Γ(p_h), v·n⟩ = 0. So on that subspace CNS\* coincides with the naive form, which is invariant under p ↦ p + const.
  - CNS\* adds exactly one equation, in a direction w with ∫_{Γ_h} w·n ≠ 0. There, B\*(p + c, w) = B\*(p, w) − c∫w·n, so that equation fixes the constant and nothing else.
  - Hence **the CNS\* velocity equals the FNZ-type (zero-flux test space) velocity in our setting**, and the pressures differ by a constant. Well-posedness of the two is equivalent.
- **Verdict on (iii).** The device (pressure term in the momentum row only, constant mode handled so that div u_h = 0 exactly) is **not new**. The mean-subtracted *square* form is a different packaging of it, whose pressure normalisation is automatically consistent (p_h ≈ p − p̄_{Γ_h}). The letter states this explicitly (paragraph C).
- **Handling of the "standard alternative" claim.** FNZ eq. (5.9) is the published version of that claim. The letter cites it.

**2. Liu–Neilan–Otus, JNM 31(2) (2023) (arXiv:2105.10409). Closest SV work.**
- **Method.** Scott–Vogelius on Clough–Tocher splits, k ≥ 2. The computational domain Ω_h is made of mesh elements inside Ω; its vertices are generally not on ∂Ω. The form a_h is Nitsche-type with a Taylor boundary correction S_h.
- **Normal component.** It is enforced by a **separate multiplier λ_h ∈ P_k(boundary edges), with mean zero**. The authors state that λ_h *"is an approximation to the pressure (modulo an additive constant) restricted to the computational boundary"*.
- **Zero-flux constraint.** The velocity space carries **∫_{∂Ω_h} v·n = 0**, and pressures are mean-zero.
- **Results.** Lemma 3.1 gives div u_h ≡ 0. Remark 3.2: *"If this constraint is not imposed in the Lagrange multiplier space, then in general (3.2) is ill-posed."* The velocity error is O(h^k) (optimal, above h^{3/2}), and the method is not pressure-robust.
- **Relevance to (ii).** This is an analogue: a constant-mode ill-posedness is noted and cured by a zero-flux or mean-zero pairing.
- **Relevance to (iii).** An exactly divergence-free, well-posed, consistent weak-BC Scott–Vogelius method on approximated curved domains **already exists**, with a higher rate. It uses a separate multiplier and boundary correction; CNS\* uses p_h itself and no correction. The letter must cite it, and says the Gjerde–Scott method, not the best possible method, is the object of study.

**3. Burman–Hansbo–Larson, ENUMATH 2017 (arXiv:1801.07463), Remark 1.**
- The Stokes cut-FEM uses b₁ = (q, div v) − (q, v·ν) in the momentum row and b₀ = (q, div v) in the mass row, i.e. the **same asymmetric placement** as CNS\*.
- The remark says the boundary term "is essential for consistency".
- So "the traction-free form is inconsistent" is known in general. The element there is Taylor–Hood, and the pressure is controlled by a stabilization.

**4. Liu–Neilan–Olshanskii, M2AN 57 (2023).**
- b(p, v) = −(p, div v) + ∫_Γ (v·n)p is used in **both** rows. Then div u_h = 0 only away from an O(h) boundary strip (Lemma 2.5), with grad-div scaling like h^{-1}.
- This is the published confirmation that the symmetric placement destroys exact incompressibility.

**5. Frachon–Hansbo–Nilsson–Zahedi, arXiv:2408.10089 (Darcy).**
- They note that a nonsymmetric penalty method with (p, v·n) needs the pressure mean prescribed, and that the multiplier "will perturb the divergence condition in the unfitted setting" (citing FNZ).
- Their Theorem 3.4 requires matching net flux for divergence preservation. The problem is Darcy, not Stokes.

**6. Informal overlap with (iv) and the missing-traction observation (not a publication).**
- GitHub `woahwhattheheck/commons` PR #15236 was merged on **2026-09-17**, before our v1 of 2026-09-28. It is a "partial research carrier" for this prize. It notes:
  - GS2024's V_h notation inconsistency (the same point as our footnote);
  - an O(h_Γ^{3/2}) consistency residual;
  - via exact rational arithmetic, that *"fixed global μ is not automatically uniform under independent boundary over-refinement"*. This is qualitatively part of (iv).
- It lists "missing traction/pressure consistency versus stress-consistent forms" as an **open** blocker. So it does **not** contain (i)–(iii).
- AML will not require citing it, but you may want to acknowledge it: priority on the qualitative part of (iv) is not ours alone.

### Claim-by-claim verdict

| Claim | Prior art found? | Verdict |
|---|---|---|
| (i) GS printed method: leak ≍ (h/μ)‖p−p̄‖_Γ, H¹ rate 1, energy 1/2, constant exactly 1; expected rate iff p is constant on the wall | No analysis of the Gjerde–Scott form (or any Nitsche–Stokes form without the pressure term) was found. The *necessity* of the term for consistency is stated in BHL 2019 Remark 1. The *mechanism* is the classical boundary-penalty consistency error (Babuška 1973), here acting on the pressure only, with penalty parameter h/μ → 0. | **Appears new** as a sharp, two-sided result. Do not present "the term is needed for consistency" as new. |
| (ii) naive ⟨p_h, v·n⟩ singular, kernel (0,1) | Analogues: LNO 2023 Remark 3.2 (ill-posed without the constraint); FNZ 2024 §5.1 (the constant mode must be fixed). Standard: a fully weak Nitsche–Stokes form always needs a pressure normalisation. | **Elementary / known in spirit.** Keep it as a one-line observation. |
| (iii) CNS\*: well posed, exactly divergence-free, h_Γ^{3/2}+h_Ω^k | The device is FNZ 2024 (the velocity is identical, see §2.1). LNO 2023 is an exactly divergence-free SV method with a separate multiplier, zero flux, boundary correction and O(h^k). | **The formulation is not new in substance; the analysis is.** New: the proof of Scott's conjectured rate for SV on an inscribed polygon, and the sharpness of h_Γ^{3/2} (Theorem D′). |
| (iv) penalty threshold ∝ ρ = h_Ω/min\|e\| | Sufficiency is a direct consequence of the standard trace-inverse inequality: Nitsche penalties must scale with the local boundary-element size (standard since Stenberg 1995 and Juntunen–Stenberg 2009). The non-uniformity of a fixed global μ under boundary over-refinement was noted informally in the GitHub PR above. | Sufficiency is **standard**. The explicit necessity (λ_T, exact arithmetic, stars) **appears new** among published work. |
| "Standard alternative": L²₀ pressure + ⟨p_h, v·n⟩ gives div u_h ≡ const | **FNZ 2024 eqs (5.7)–(5.9)**, in the fully weak setting. | Known. Cite FNZ. The precise form in our setting is below. |

### The "standard alternative": precise statement (checked algebraically and numerically)

Take pressures in Π_h ∩ L²₀ and the momentum term −(p_h, div v) + ⟨p_h, v·n⟩ for v ∈ V_R. The continuity row is tested by q ∈ Π_h ∩ L²₀; this is square and equivalent to a mean-zero multiplier α. Since div V_h = Π_h:

- div u_h ≡ c, a constant, with **c = |Ω_h|⁻¹ ∫_{∂Ω_h} u_h·n = |Ω_h|⁻¹ (∫_{∂R} g_I·n + ∫_{Γ_h} u_h·n)**, where n is the outward normal of Ω_h.
- For q ∈ L²₀ the pressure form equals B\*(q − p̄_{Γ_h}(q), v). Hence, if both problems are uniquely solvable, **c = 0 ⇔ the discrete net flux through ∂Ω_h vanishes ⇔ the solution coincides with that of CNS\* ⇔ p̄_{Γ_h}(p_h^{CNS\*}) = 0**.
- Flux-correcting g_I removes only the first term of c. Nothing forces ∫_{Γ_h} u_h·n = 0 under weak imposition, so c ≠ 0 generically.
- **c = 0 holds exactly** when all of ∂Ω_h carries strongly imposed normal data with zero net flux, which is FNZ's remark.

**Numerical check.** Scratch script `…/scratchpad/l20check.py` uses `code/svn.py` and SuperLU, with k = 4, μ = 100 and relative residual ≈ 6e-16.

| Test | N | CNS\* ‖div u_h‖ | L²₀ variant α | L²₀ variant ‖div u_h‖ | \|α\|·\|Ω_h\|^{1/2} | p̄_Γ(p_h^{CNS\*}) |
|---|---|---|---|---|---|---|
| B | 16 | 5e-12 | −5.98e-5 | 2.80e-4 | 2.80e-4 | −4.6e-3 |
| B | 32 | 9e-12 | −2.94e-5 | 1.37e-4 | — | −3.2e-3 |
| A (p = 0) | 16 | 4e-13 | 9.69e-4 | **4.54e-3** | 4.54e-3 | 7.5e-2 |

- The divergence of the L²₀ variant is exactly the constant α, confirming the claim.
- The H¹ errors of the two methods agree to 4 digits for test B.

**Side finding for the paper (07_consistent_method.tex, l. 15).** The paper reports that "a bordered direct solve" of the naive form gave discrete divergence ≈ 4.5·10⁻³, "consistent with incompatibility". That value coincides with the L²₀ / mean-bordered solution for test A at N = 16, μ = 100 (4.54e-3). The likelier explanation is that the bordered solve computed the mean-constrained solution, whose divergence is the constant c, rather than an incompatible system. Consider rewording that sentence to cite FNZ (5.9) and the formula for c. I have not confirmed which test and parameters produced the paper's number.

## 3. What the letter text does (lit_framing.tex)

| Paragraph | Content |
|---|---|
| A1 | Scott–Vogelius, Guzmán–Scott, pressure robustness (John et al.), split meshes, Neilan's review. The draft is careful to say that pressure robustness holds for fitted meshes with strong BCs; LNO and LNOl note that weakly imposed variants are not pressure-robust. |
| A2 | The gradient constraint, h^{1/2} vs the scalar h^{3/2}, isoparametric remedies, GS Nitsche, Scott's prize and "why not?". The prize wording is quoted from the PDF. |
| B | Nitsche–Stokes with traction (FS95, Stenberg, Becker, BH14, BHL19 "essential for consistency"); BVC / cut-FEM / SBM / DGS; divergence-free weak-BC works (LNOl, FNZ, LNO, BHL24, EST). |
| C | The contribution. It explicitly disclaims novelty of the (iii) device, states the velocity equivalence with FNZ, and puts the novelty in (i), the rate proof with sharpness, and (iv) necessity. |
| D | Scope: symmetric form only; the penalty-method mechanism (Babuška); Robin (Juntunen–Stenberg); multiplier and BVC methods reach higher order; the precise L²₀ statement with FNZ (5.9); element families not covered. |

Test compile: `pdflatex` + `bibtex` (amsplain) on a wrapper document produce no warnings, about 1,900 words including references. Paragraphs B–D may need trimming for AML's length limit.

## 4. To do before submission

1. Read GjerdeScott2022 (AMS PDF) and check whether its slip-Nitsche form has a pressure term.
2. Read Freund–Stenberg 1995 to confirm that it treats Stokes; otherwise drop it from "Stokes" citations and keep Stenberg 1995, Becker 2002 and BH14.
3. Confirm Becker 2002's issue number and content, and the pages of Babuška 1973 and Burman 2012.
4. Decide whether to acknowledge the commons PR #15236.
5. Reword the 4.5·10⁻³ sentence in §7, as above.
6. Update `paper/refs.bib` Cavalcante's author field.
