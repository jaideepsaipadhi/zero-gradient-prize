# Applied Mathematics Letters submission

## Files

| file | contents |
|---|---|
| `letter.tex` | Manuscript source (elsarticle). |
| `refs.bib` | The 18 cited entries only, built from `notes/v5/refs_letter.bib` and `paper/refs.bib` (MorganScott1975, NeilanOtus2021), plus ChenLiu2026, with the VERIFY comments kept. Cavalcante's author field is "de Souza Cavalcante, Claudemir" (record: "Claudemir de Souza Cavalcante"). `PadhiExt` is the extended version, which will be posted as Zenodo **version 3**. Its DOI is the placeholder `10.5281/zenodo.XXXXXXXX`, marked `%% FILL v3 DOI`. |
| `highlights.txt` | 5 highlights, 64–84 characters each. |
| `cover_letter.txt` | Cover letter to the editor, joint (first person plural), with a FILL line for Cavalcante. |
| `letter.pdf` | The `[3p,times]` build. This is the page-count reference. |
| `letter_5p.pdf` | The `[5p,times]` two-column build, for checking only. |
| `elsarticle.cls`, `elsarticle-num.bst` | Elsevier class v3.4c (2025/01/11) and bst v2.1. They were not in this TeX Live install and CTAN was blocked, so they were taken from `github.com/quarto-journals/elsevier`. |
| `build.sh` | Builds both PDFs. |

## Build

```sh
cd letter && ./build.sh
```

This runs pdflatex, bibtex, then pdflatex twice, for `letter`. It then does the same for jobname `letter_5p` with `\def\FIVEP{}\input{letter}`. In the 5p build the tables become `table*` and one display gets a line break.

## Page counts and length

- **3p (`letter.pdf`): 6 pages**, including references and declarations. The limit is 6, and the last page is full (v6): anything Cavalcante adds must be offset by cuts.
- **5p (`letter_5p.pdf`): 5 pages.**
- **Abstract:** 230 words (limit 250).
- **Keywords:** 5.
- **Build:**
  - No undefined references or citations, and no bibtex warnings.
  - The 3p build has no overfull boxes and no Type 3 (bitmap) fonts; typewriter text uses `txtt`.
  - The 5p build has two overfull boxes of at most 8.3 pt, both in displays that are fine in 3p.

## Checklist against the AML guide for authors

- [x] At most 6 journal pages (3p estimate: 6).
- [x] Abstract of at most 250 words (230). It stands alone: it defines the setting, k, h_Γ and h_Ω in words.
- [x] 1–7 keywords (5).
- [x] Highlights file: 3–5 bullets of at most 85 characters each, including spaces (5 bullets, 64–84).
- [x] Numbered references in elsarticle-num style; only cited works are listed (18).
- [x] Declaration of competing interest. It states that the question is the subject of a prize offered by L. R. Scott and that the authors intend to submit this work for it; there are no other competing interests.
- [x] Data availability cites the GitHub repository and the Zenodo v3 record, which will also archive the code and data. The DOI is a placeholder.
- [x] Generative-AI declaration:
  - it uses Elsevier's current heading ("...in the manuscript preparation process") and wording ("...of the published article");
  - it sits before the references;
  - one sentence at the end of Section 1 also refers to it.
- [ ] Authors: Cavalcante (affiliation, e-mail, CRediT roles: FILL) and Padhi, alphabetical. See `CHANGES_v6.txt`.
- [x] Pointers into the extended version (macro `\ext`) follow the current numbering in `paper/main.aux` (penalty §9, pressure §8, rescue §7.2). Every one was checked against that file.
- [ ] **Author to do: insert the v3 DOI.** Replace `10.5281/zenodo.XXXXXXXX` in `refs.bib` (entry `PadhiExt`) and in the Data availability section of `letter.tex`. Both places are marked `%% FILL v3 DOI`.
- [ ] **Author to do: re-check pointers.** If the long paper is renumbered again, re-check the `\ext` pointers.
- [ ] **Author to do: declarations `.docx`.** AML requires the Elsevier declarations-tool `.docx` (competing interests) at submission. Its content must match the manuscript statement about Scott's prize.
- [ ] **Author to do: open literature items.** See `notes/v5/lit_report.md` §4.
  - Freund–Stenberg 1995 is not cited, because its bib entry is only partly verified; the text uses "see, e.g., [Burman–Hansbo 2014]".
  - GjerdeScott2022 is not cited.

## v6 changes

See `CHANGES_v6.txt` (joint authorship, unconditional sharpness, strong-imposition sentence, Chen–Liu and Neilan–Otus citations, Table 2 folded into text).

## Structure after the referee round

**What the theorems contain.** Theorems 1 and 3 and Propositions 2 and 4 now contain only what the letter proves or sketches:
- Theorem 3(c) is the H¹ bound for bounded γ, which follows from (b).

**Remark 1** ("Further results, proved in [14]") collects the results proved only in the extended version:
- the GS pressure;
- the γ-uniform CNS* H¹ bound and its pressure estimate;
- sharpness (the H¹ lower bound);
- the growing-penalty rescue, with its hypotheses and the conditioning cost.

## What was cut to fit 6 pages

**Material:**
- The explicit forms A_T, B_T and C_T, the star certificates, and the account of the gap between the certified 9.3 and the measured γ* ≈ 21. These now have pointers to §10.
- The H¹ leak constant K.
- The proofs of the lower bound, the γ-uniform bound and the pressure estimates. These are pointers only.

**References:** Berger–Scott–Strang, Scott 1975, Neilan's review, Durst–Neilan, Neilan–Otus 2021 (with its sentence on fitted elements), the shifted boundary method, Burman–Hansbo–Larson 2018, Bramble–Dupont–Thomée, Dupont–Guzmán–Scott, Freund–Stenberg and Becker.

**Tables and figures:**
- The first-layer pressure table and the N = 16 rows.
- The column for strong imposition on the standard mesh, now one sentence.
- The refinement half of the threshold table, now a caption line.
- The optional theorem-header pointers for Lemmas 2–3 and Proposition 2.
- There is no figure.
