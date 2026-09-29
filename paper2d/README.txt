paper2d/ -- journal version of the long 2D paper (Paper B)
============================================================

Target: SIAM J. Numer. Anal. / IMA J. Numer. Anal. / Math. Comp. (plain amsart, 11pt).
Author: Jaideep Sai Padhi (Purdue).  %% TBD in main.tex: C. de S. Cavalcante may join (joint AML letter); not added.

Build
-----
  cd paper2d && pdflatex main && bibtex main && pdflatex main && pdflatex main
Result: main.pdf, 50 pages (about 46 pages of text, then declarations and references).
No undefined references or citations, no overfull boxes. Uses tikz (one schematic figure, Fig. 1) and
the four PDF figures in figures/ (copied from paper/figures/).
refs.bib = paper/refs.bib plus two entries: PadhiExt (the extended version, same record as PadhiExt in
letter/refs.bib; DOI placeholder, marked %% FILL v3 DOI) and ZeroGradRepo (the GitHub repository).
Nothing in paper/, paper3d/ or letter/ was modified. Nothing was committed.

Editorial feedback (Ramanujan J. editor) and where it is addressed
------------------------------------------------------------------
(1) AI-use statement: "Declarations" at the end (what Claude was used for: proofs, code, certificates,
    literature, drafting; how results were verified: independent referee passes, exact/interval
    certificates with negative controls, reproducible logs; author responsible; AI not an author).
(2) Readable exposition with history: Section 1 is prose (Scott-Vogelius 1985, Guzman-Scott, pressure
    robustness, Berger-Scott-Strang / Scott 1975, Bramble-Dupont-Thomee, isoparametric elements, the
    zero-gradient phenomenon, Nitsche 1971 / Babuska / Stenberg, Gjerde-Scott, the prize question quoted
    verbatim, the leak law and the repair explained before any notation). Each section opens with the
    idea; Section 3.7 is a roadmap of the error analysis; Fig. 1 illustrates the chordal bump, the leak
    and the two-triangle lock. Internal hypothesis labels are reduced to (M0)-(M2) (defined in words);
    (M3) became the named "star condition" (Definition 5.10); (D) became Assumption 3.5.
(3) Broader related work: Section 2 (curved-domain SV: Neilan-Otus, Durst-Neilan, Chen-Liu,
    Liu-Neilan-Otus; cut/Nitsche: Liu-Neilan-Olshanskii, Frachon-Nilsson-Zahedi, Burman-Hansbo-Larson,
    Eickmann-Scott-Tscherpel; Nitsche for Stokes; nearly singular vertices: Grassle-Bohne-Sauter, Park,
    Bohne-Grassle-Sauter, Kean-Neilan-Schneier; broader context incl. 3D families, local estimates,
    nonsingular branches).
(4) Explicit novelty statement: Section 2.6 ("Known" vs "New in this paper"); Cavalcante 2026 treated
    in Section 2.1 (constant-pressure case, same prize, our proofs independent).

What is in the core paper (with proofs or proof sketches)
---------------------------------------------------------
 Sec. 3  setting; chord geometry; transposed-gradient lemma (odd symmetry); coercivity; approximation
         lemma (full proof); error equation.                          [ext. Secs. 2-5]
 Sec. 4  printed method GS: energy rate 1/2 (two-sided), normal-load cancellation, H1 rate 1, slip and
         H1 lower bounds, exact leading term (constants 1), leak flow and K (full proof), annulus
         bracket, series value K=2.7565157, K -> 3/sqrt7.          [ext. Sec. 6]
 Sec. 5  singular naive form; CNS* well posed with Scott's rate (full proof); FNZ remark; penalty-uniform
         H1 bound with flux-corrected data (full proof via lift + discrete Dirichlet lemmas); lower
         bound c h_G^{3/2} (full proof), k>=7 single-triangle field, explicit P4 witness for k>=4
         (certificates referenced), why the transposed-gradient route fails; growing penalty
         (H1 rescue; energy lower bound with proof; rate table).       [ext. Secs. 7, 7a, 7e]
 Sec. 6  penalty: sufficiency; single-triangle necessity (explicit forms, proof); exact shape-family
         certificates (method); rho-scaling corollary; exact star counterexamples; resolved-region and
         every-(M0)-mesh necessity (statements + idea); stacked columns (proof); measured threshold and
         its gap.                                                     [ext. Sec. 8]
 Sec. 7  strong imposition: local structure (proof), two-sided h_G^{3/2} under (M2) (proofs), counting
         bound (proof), locked-star lemma (proof), constrained/sharp inf-sup and two-sided rate with
         locks (sketch), needle-angle law (proof of the local law; velocity bound cited).
                                                                      [ext. Sec. 8b]
 Sec. 9  numerics: selected tables (CNS* rates, B100, rate crossover, pressure, test P constants, K_h/d_h,
         other elements, dual norms, thresholds, growing penalty, flux, higher degree, strong vs Nitsche,
         drift, locks); mu=100 margin table and caveat; corner caveat.   [ext. Sec. 13 and tables]
 Sec. 10 conclusions, practical advice, open problems (prose).        [ext. Sec. 14]

What went to the supplement (cited as \cite{PadhiExt}; statements only in Sec. 8 or omitted)
-------------------------------------------------------------------------------------------
 - Penalty-uniform leak limit (ext. 7d), pressure/stress/drag incl. cell problem (7b), L2 theory (7c),
   Clough-Tocher k>=2 (9), general curved domains (10), Navier-Stokes (11, 11c), pressure robustness
   (11b), nearly singular stars / needle pairs / slit-domain inverse (8b), 3D pointer (PadhiThreeD).
 - Proof details: uniform Korn/inf-sup gluing (Lemma gd:lem:bog), full M3 witness algebra and interval
   certificates (07e), Argyris explicit constants and the (M0) necessity proof (08), the gamma*_infty
   cell limit (modulo LOC/CUT), conditioning proposition, rigorous bracket computations, the vertex
   lemma for strong imposition, needle-pair lemmas A1/A2.
 - Conditional results (boundary pressure-moment hypothesis, cell decay, transplant lemma, LOC/CUT) are
   mentioned only plainly in Sec. 8 (with the status table, Table 6) and Sec. 10; main-text theorems
   are all fully proved (computer-assisted where stated).
 - Appendix A (reproducibility file list) replaced by the data-availability statement and repository
   citation.

Open items for the author (marked in the sources)
-------------------------------------------------
 - %% TBD main.tex: co-author (Cavalcante) decision.
 - %% FILL Z_declarations.tex / refs.bib: Zenodo DOI of the extended version (v3).
 - %% TBD Z_declarations.tex: disclose the joint AML letter if submitted; further acknowledgements.
 - %% VERIFY Z_declarations.tex: funding statement.
 - %% VERIFY 03_setting.tex: BCDG2016 theorem number (checked only in arXiv/author copies).
 - VERIFY comments carried in refs.bib (FreundStenberg1995 / Becker2002 content not read, etc.).
