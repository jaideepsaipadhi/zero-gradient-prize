jsc/ -- Journal of Scientific Computing (Springer) submission of the long 2D paper (Paper B)
==========================================================================================

Copied from paper2d/ on 2026-09-29 and converted to the Springer Nature LaTeX template (sn-jnl.cls,
reference style sn-mathphys-num: numbered, as the journal requires). paper2d/, paper/, paper3d/ and
letter/ were not modified. Nothing was committed.

Build
  cd jsc && ./build.sh          (= python3 flatten.py; pdflatex; bibtex; pdflatex; pdflatex)
  Result: main.pdf, 66 pages; no LaTeX errors, no undefined references or citations, no overfull boxes.
  (Remaining warnings: hyperref "Token not allowed in a PDF string" x6, identical to paper2d.)

Files
  main.tex         GENERATED single-file manuscript (the template asks for one .tex file, no \input).
                   Upload this one. Do not edit it; edit main.src.tex or sections/ and rerun build.sh.
  main.src.tex     frame: class, packages, theorem styles, title page, abstract, keywords, MSC.
  flatten.py       inlines sections/01..10 verbatim, puts the AI statement in Section 3.8, builds the
                   Declarations from sections/Z_declarations.tex.
  paper2d_main.tex paper2d/main.tex (source of the macro block, copied verbatim by flatten.py).
  sections/        copies of paper2d/sections (edits listed below).
  sn-jnl.cls, sn-mathphys-num.bst   Springer Nature template files (source: GUIDELINES.txt item 1).
  GUIDELINES.txt   the journal's requirements, with URLs and quotations.
  cover_letter.txt, SUBMISSION_CHECKLIST.txt.

What changed relative to paper2d (no mathematical change)
  - Class amsart -> sn-jnl; geometry dropped (class layout); amsthm loaded before the class so the class
    defines its theorem styles; theorem environments keep their names and the shared per-section counter.
  - (Before the 2026-09-29 reframing, below:) all 96 labels resolved to the same numbers as in paper2d (checked from the .aux files); the same 64
    references are cited. One new label: sec:ai.
  - Abstract shortened from ~410 to 244 words (limit 150-250); no new claims. Keywords reduced from 7 to 6
    (dropped "polygonal approximation"). MSC via \pacs[MSC Classification].
  - Generative-AI statement (wording unchanged) moved from Declarations to a new final subsection of
    Section 3, "3.8 Use of generative AI", because the journal's LLM policy asks for it in the Methods
    section; Declarations and Author contributions point to it.
  - Declarations in the template's order: Funding, Competing interests (as paper2d), Ethics/Consent (not
    applicable), Data availability (as paper2d), Materials (not applicable), Code availability (GitHub +
    PadhiExt, doi:10.5281/zenodo.23044296), Author contributions (sole author), Use of generative AI
    (pointer). Acknowledgements (as paper2d) as \bmhead before the Declarations.
  - Layout-only edits in sections/ because the template's text block is narrower: seven wide displays
    split over two lines (\[..\] -> gather*/align*, same formulas, none labelled): 03_setting (two),
    04_printed, 05_corrected (two), 06_penalty, 07_strong; \tabcolsep reduced in three tables (01_intro
    Table 1, 09_numerics mu=100 margin table and Table "flux"); one "\end{tabular}\\[6pt]" ->
    "\end{tabular}\par\vspace{6pt}" (09_numerics, threshold table; the \\ was an error in sn-jnl).

Companion-paper reframing (2026-09-29; author's plan: the joint AML letter, letter/letter.tex, with
C. de S. Cavalcante is submitted FIRST; this sole-author paper is its companion). No mathematical change,
no proof deleted.
  - refs.bib: new entry CavalcantePadhiLetter (title from letter/letter.tex; %% FILL status/DOI).
  - New title (main.src.tex, cover letter, checklist); abstract rewritten (232 words) to lead with what the
    letter does not contain: strong imposition / zero-gradient lock, penalty theory beyond one triangle,
    H^1 leak constant; then "complete proofs of the letter's results".
  - Section order changed in flatten.py: Setting -> Strong imposition (07_strong, now Sec. 4) -> Penalty
    (06_penalty, Sec. 5) -> Printed (04_printed, Sec. 6) -> Corrected (05_corrected, Sec. 7) -> Further,
    Numerics, Conclusions. Labels unchanged, so all cross-references resolve; theorem NUMBERS changed
    (e.g. Theorem 4.x -> 6.x) relative to paper2d and to the letter's pointers into PadhiExt (those point
    to the Zenodo version and are unaffected).
  - Intro: new subsection 1.6 "This paper and the companion letter" (label sec:companion); main-results list
    reordered (lock, penalty, leak constant, complete proofs of the letter), each item says what is/is not
    in the letter. Related work: Sec. 2.6 (label sec:knownnew) now has three lists: known / announced in the
    joint letter / new here; Cavalcante's constant-pressure analysis still credited to his preprint.
  - Statements announced in the letter carry "; announced in [CavalcantePadhiLetter]" in their theorem
    headers (thm:A, thm:C, thm:B, cor:Bp, prop:singular (header shortened), thm:D, thm:lb, thm:pb1,
    prop:cert), and Sections 5, 6, 7 open with a "Relation to [letter]" sentence; Section 4 says none of
    it is in the letter. Transitional prose in 07_strong/06_penalty/03_setting adjusted for the new order.
  - Conclusions first paragraph reframed (two halves: letter's result proved in full; new material).
  - Author contributions (flatten.py) mention the joint letter; Z_declarations TBD comment replaced.
  - cover_letter.txt: related-manuscripts disclosure rewritten (letter supplied as supplementary file);
    SUBMISSION_CHECKLIST.txt: submit only after the letter; upload the letter PDF as related manuscript.

Open items (marked in sources): ORCID (%% FILL), funding (%% FILL), AML-letter submission date/number (%% FILL in cover letter) and
status/DOI (%% FILL in refs.bib), joint/sole attribution wording (%% CHECK in sections/02_related.tex), Neilan e-mail
(%% VERIFY), and the VERIFY items carried over from paper2d (BCDG2016 theorem number; refs.bib).
