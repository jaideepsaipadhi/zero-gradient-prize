arXiv BUNDLES  (prepared 2026-09-30, branch v3-draft, not committed)
====================================================================

!! POSTING ORDER !!
The LETTER (letter/, letter_src.tar.gz) is joint work with Claudemir de Souza Cavalcante. As of 2026-10-08 it is
rebuilt from letter/letter.tex with all of his corrections (e-mail of Oct 2026) applied, and full-line %-comments
stripped. His details, CRediT roles and the corresponding author (Padhi) are filled in. Post it only once he has
also agreed to the arXiv posting itself, and confirmed the purpose clause for his AI use (see letter/CHANGES_v9.txt).
The EXTENDED paper is single-author and can be posted independently. Before posting it, fix the incorrect
statement about Cavalcante's preprint in paper/sections/09_clough_tocher.tex (see ../letter/CAVALCANTE_PREPRINT_NOTES.txt,
section 3), and then rebuild the bundle.

NOTE ON COMMENTS: arXiv serves the source to anyone who asks for it. All %-comments in the .tex files are therefore
public, including FILL/VERIFY/CHECK notes, internal file paths (notes/v6/...) and referee remarks. They were
deliberately left in, as instructed. Review them before uploading, or strip them (e.g. with arxiv_latex_cleaner).

CONTENTS
--------
letter/                  letter.tex (current joint version, unchanged), letter.bbl (fresh bibtex output),
                         elsarticle.cls, elsarticle-num.bst (v3.4c / v2.1; included so the build does not depend on the
                         version installed at arXiv). There are no figures.
letter_src.tar.gz        tarball of letter/ (files at top level).
extended/                main.tex, main.bbl, sections/*.tex (the 23 files actually \input; the unused
                         12_three_d.tex is left out), figures/*.pdf (4).
extended_src.tar.gz      tarball of extended/ (main.tex at top level).

VERIFICATION (done)
-------------------
Each tarball was extracted into an empty directory and compiled with pdflatex twice, without bibtex:
  letter:   6 pages, 0 errors, 0 undefined references/citations.
  extended: 156 pages, 0 errors, 0 undefined references/citations.
There are no absolute paths in any .tex or .bbl file (checked for /home, /root, /tmp, /mnt, /Users, C:\), and
\input, \includegraphics and \bibliography are all relative. The tarballs use relative member names (./...).
The letter compiles in 3p by default; the 5p variant is only switched on by an external \def\FIVEP.

REBUILDING
----------
letter:   copy letter.tex, refs.bib, elsarticle.cls and elsarticle-num.bst into a scratch directory, then run
          pdflatex, bibtex, pdflatex, pdflatex. Copy letter.tex and letter.bbl (plus the .cls and .bst) into
          arxiv/letter/, then  tar -czf letter_src.tar.gz -C letter .
          Strip full-line %-comments from letter.tex first (they are public on arXiv).
extended: the same, from paper/ (main.tex, refs.bib, sections/, figures/). Use style amsplain.

======================================================================
SUGGESTED arXiv METADATA -- EXTENDED PAPER
======================================================================
Title:    Scott--Vogelius--Nitsche on a polygonally approximated boundary: convergence, the missing pressure
          traction, and the penalty threshold
Authors:  Jaideep Sai Padhi
Primary:  math.NA  (cs.NA is cross-listed automatically, being the alias of math.NA)
Cross-list: physics.flu-dyn is NOT recommended. The paper is finite element error analysis, and the Stokes and
          Navier-Stokes content is a model problem rather than a fluid-dynamics result, so the moderators may
          remove the cross-list. Add it only if you want the fluids audience.
MSC:      65N30 65N12 76D07  (the paper itself also lists 65N15; arXiv accepts "65N30; 65N12, 65N15, 76D07")
Comments: 156 pages, 4 figures. Code and data: https://github.com/jaideepsaipadhi/zero-gradient-prize ;
          Zenodo DOI 10.5281/zenodo.23044296
Journal-ref / DOI: leave blank (the Zenodo DOI goes in Comments, not in the DOI field, which is for the
          published version).
Abstract (plain text; arXiv renders $...$ with MathJax):
We analyse the Scott-Vogelius-Nitsche method of Gjerde and Scott for two-dimensional Stokes flow whose no-slip curve is replaced by an inscribed polygon, answering Scott's zero-gradient prize question: is the error $O(h_\Gamma^{3/2}+h_\Omega^k)$, and if not, why not? For the method as printed the answer is no as soon as the wall pressure varies: the Nitsche form omits the traction $-pn$, and the discrete velocity leaks through the wall with normal velocity $(h/\mu)(p-\bar p)$. We prove the sharp $H^1$ rate $h/\mu$ and the exact energy rate $(h/\mu)^{1/2}$, with leading constant 1, and identify the $H^1$ constant with a Stokes leak flow. Subtracting the wall mean of the discrete pressure gives a well-posed, exactly divergence-free method (its velocity is that of the zero-net-flux formulation of Frachon, Nilsson and Zahedi); we prove it attains $h_\Gamma^{3/2}+h_\Omega^k$ for general data, for every $k\ge4$ under an angle condition and every $k\ge2$ on Clough-Tocher refinements, with a matching lower bound $c\,h_\Gamma^{3/2}$ whenever the wall shear is nonzero. We also give pressure, drag and $L^2$ rates. The penalty must scale like $h_\Omega/\min_e|e|$: sufficiency on every mesh, necessity certified in exact arithmetic for explicit shape classes. The theory extends to general smooth boundaries and to steady Navier-Stokes. Strongly imposed no-slip on the same polygon already attains the rate under our mesh condition; the zero-gradient lock comes from wall vertices in only two triangles. Computations up to 256 chords, for k=4,5,6 and Clough-Tocher k=2,3, confirm the rates and constants.
  [1610 characters; limit 1920]

======================================================================
SUGGESTED arXiv METADATA -- LETTER  (post ONLY after Cavalcante agrees)
======================================================================
Title:    The missing pressure traction in Scott--Vogelius--Nitsche methods on polygonal approximations of a
          curved wall
Authors:  Claudemir de Souza Cavalcante, Jaideep Sai Padhi   (alphabetical, as in the manuscript)
Primary:  math.NA   (cross-list: none; see the note above on physics.flu-dyn)
MSC:      65N30 65N12 76D07
Comments: 6 pages. Submitted to Applied Mathematics Letters. Companion to arXiv:XXXX.XXXXX (the extended version,
          once it is posted). Code and data: https://github.com/jaideepsaipadhi/zero-gradient-prize
          (Zenodo DOI left out at Cavalcante's request until the archive/version correspondence is confirmed.)
          (Journal policy: Elsevier allows preprints on arXiv. Mention the submission only if both authors agree.)
Abstract (the manuscript abstract, v9, in plain text):
We study smooth two-dimensional Stokes flow with a curved no-slip wall approximated by an inscribed polygon, using Scott-Vogelius elements of degree $k\ge4$ and Nitsche's method. For the printed Gjerde-Scott formulation, the pressure traction is absent. Under the approximation assumptions of Theorem 1(c), with fixed effective penalty, bounded mesh ratio and nonconstant wall pressure, the leading normal leak is $(h_\Omega/\mu)(p-\bar p)$, with asymptotic constant 1. The energy and $H^1$ errors have orders $(h_\Omega/\mu)^{1/2}$ and $h_\Omega/\mu$, respectively, with matching lower bounds. A mean-free pressure-traction correction has the same velocity as a zero-net-flux formulation of Frachon, Nilsson and Zahedi and attains $O(h_\Gamma^{3/2}+h_\Omega^k)$ under the stated assumptions and a sufficiently large, bounded effective penalty. Under condition (S) and nonzero wall shear, the exponent 3/2 is sharp for every $k\ge4$. A penalty proportional to $h_\Omega/\min_e|e|$ is sufficient and, under a triangle-shape condition, necessary. Quartic computations support the predicted rates, leak constant and penalty thresholds.
  [1132 characters; limit 1920]

SHA-256 of the tarballs:
c6b2cf0cdc3d1b856a79003b864a4350f1e29dd511006c89f6a99d57babbadfb  letter_src.tar.gz
32f9e8c092e3d7e71431e50cc35222e1041868060845f0f1b9c0f4cc2819f73e  extended_src.tar.gz
