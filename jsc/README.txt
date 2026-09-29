jsc/ -- Journal of Scientific Computing (Springer) submission of the long 2D paper (Paper B)
==========================================================================================

Copied from paper2d/ on 2026-09-29 and converted to the Springer Nature LaTeX template (sn-jnl.cls,
reference style sn-mathphys-num: numbered, as the journal requires). paper2d/, paper/, paper3d/ and
letter/ were not modified. Nothing was committed.

Build
  cd jsc && ./build.sh          (= python3 flatten.py; pdflatex; bibtex; pdflatex; pdflatex)
  Result: main.pdf, 64 pages; no LaTeX errors, no undefined references or citations, no overfull boxes.
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
  - All 96 labels resolve to the same numbers as in paper2d (checked from the .aux files); the same 64
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

Open items (marked in sources): ORCID (%% FILL), funding (%% FILL), co-author decision (%% TBD),
AML-letter disclosure (%% TBD in Declarations; %% AUTHOR CHECK in cover letter), Neilan e-mail
(%% VERIFY), and the VERIFY items carried over from paper2d (BCDG2016 theorem number; refs.bib).
