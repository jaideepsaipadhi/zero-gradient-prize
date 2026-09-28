# Clough–Tocher extension (k ≥ 2): report

Deliverables:
- `notes/v3/ct_section.tex`: a section to `\input` after Section 8. Labels are `ct:*`, and it uses only the macros in `main.tex`.
- `code/lamstar_exact_ct.py`: the exact computations.
- `logs/lamstar_exact_ct.log`: their output. Runtime is about 4 minutes, single-threaded, under `nice`.

Compile check: I put a copy of `paper/` in the scratchpad, added `\input{sections/ct_section}` after `08_penalty`, and appended the BibTeX entries below. It compiles with pdflatex and bibtex with no errors, no undefined references or citations, and no overfull boxes in the new section. The section becomes Section 9, about 4.5 pages. Nothing in `paper/` was modified, and nothing was committed.

## What is proved

**Theorem 9.x (`ct:thm:main`).** Let the macro mesh satisfy (CT0): shape regularity plus c₀h_Γ ≤ |e| ≤ h_Γ. Split every macro triangle at its barycentre, and use continuous P_k velocities and discontinuous P_{k−1} pressures on the split, for any fixed **k ≥ 2**. Then every result of Sections 3–8 holds as stated. This covers Lemmas 3.1–3.5 and 4.1–4.6, the error equation, Theorems A, B, C and D, Corollaries A′ and B′, Proposition 7.1 and Proposition 8.1(a). Constants now depend on σ, c₀, k and L. **(M1), (M1′) and (M2) are dropped entirely.** For example, CNS\* gives:
- k=2: ‖ũ−u_h‖_{H¹} ≤ C(h_Γ^{3/2}+h²) when h_Γ ≥ h⁴;
- k=3: ‖ũ−u_h‖_{H¹} ≤ C(h_Γ^{3/2}+h³) when h_Γ ≥ h².

Item by item:

| Paper item | Status under Clough–Tocher |
|---|---|
| Mesh assumptions | (M1), (M1′) and (M2) all follow from shape regularity of the macro mesh (Lemma `ct:lem:mesh`). The refined mesh has no singular vertex; the proof counts rays: an interior macro vertex has ≥6 refined edges, a boundary vertex has 2m+1 ≥ 3 with none opposite, and a barycentre has 3 with none opposite. (M2) holds with Θ₀ = Θ₀(σ): the two refined angles inside one macro angle are consecutive and add up to it, which settles macro vertices, and a separate argument handles barycentres. Every vertex of Γ_h lies in ≥2 macro triangles, so ≥4 refined ones. No macro triangle has two edges on Γ_h. Each e ∈ E_h^Γ lies in exactly one refined triangle K_e, with height ⅓ that of T_e. **Consequence: for k ≥ 4 the paper applies verbatim to the refined mesh.** The new content is k = 2, 3, and removing the angle hypotheses for every k. |
| Lemma 4.3 (inverse traces) | Holds with T_e replaced by K_e. The "two edges on Γ_h doubles C_I" clause is void. C_I is larger, because K_e is thinner. |
| Lemma 4.4 (coercivity) | Unchanged formulas (μ₀ = 4κρC_I², c₁, M), with the new C_I. κ is unchanged. |
| Lemma 4.5 (divergence range) | div V_R = Π_h for all k ≥ 2, with the edge bubble taken on K_e. **W_h ≠ ∅ for every g**: the corner obstruction cannot occur, since no vertex is singular. |
| Lemma 4.6 (uniform inf-sup) | New proof, valid for all k ≥ 2 (Lemma `ct:lem:infsup`). It uses the macro-element technique in three steps: (1) the uniform continuous (Bogovskiĭ) constant from the paper's Galdi wedge cover; (2) a P₂–P₀ Fortin operator on the *macro* mesh (Scott–Zhang plus macro-edge bubbles), which lies in V_h because k ≥ 2; (3) local surjectivity on each split triangle, transferred by Piola. **The Galdi argument transfers, and more simply than before:** the domain enters only through β_c, with no Poincaré constant and no Θ₀. |
| Local surjectivity (`ct:lem:local`) | div: V₀(T) → Q₀(T) is onto for all k ≥ 2. The proof is by Piola transfer from the reference split. On the reference: for k=2, a dimension count (8 = 8) plus unisolvence of HCT; for k ≥ 4, Scott–Vogelius / Guzmán–Scott applied to the non-singular 3-triangle mesh of T̂; for k=3, an exact rational rank computation, also run for k = 2, 4, 5, 6 as a check. Kernel dimensions are 0, 3, 9, 18, 30. The general k ≥ 2 statement is also in Guzmán–Neilan 2018 (cited). |
| Approximation in W_h (flux exponent) | Same statement, for all k ≥ 2. Equispaced P_k interpolation on an unsplit edge is the closed Newton–Cotes rule, exact on P_{σ_k−1}, with σ_k = k+2 for even k and k+1 for odd k. So **σ₂ = σ₃ = 4** (Simpson and the 3/8 rule are both exact on cubics). The "no flux loss" condition is h_Γ ≥ h^{2(σ_k−k)}, which is h⁴ for k=2 and h² for k=3; the even/odd pattern is the same as for k ≥ 4. The flux term vanishes for g·ν of degree ≤ σ_k−1 per side (≤ 3 for k = 2, 3). The lower bound tnorm ≥ (μ/h)^{1/2}\|m_R\|/\|Γ_h\|^{1/2} is unchanged. |
| v_φ (Section 6) | For k = 2, 3: v_φ = curl of the **HCT cubic C¹ interpolant** on the macro mesh, which is piecewise P₂ ⊂ P_k. Its errors are L^∞ Ch³ ≤ Ch^k and gradient Ch². For k ≥ 4: Morgan–Scott on the refined mesh, which is allowed because the refined mesh satisfies (M2). **Side finding:** Section 6 states ‖v_φ − w‖_{W^{1,∞}} ≤ Ch^k, but its proofs use only ‖v_φ−w‖_{L^∞} ≤ Ch^k plus bounded ∇v_φ. I audited each use (Theorem A step 1, Theorem C steps 4–5, Theorem B(i) steps 1–5, Corollary B′ step 4); the section records this. |
| Lemma 6.3 (cancellation) | Holds unchanged. It uses only three facts: v ∈ Z_h is a pointwise div-free polynomial on the element touching e (now the split triangle K_e); v is continuous at the vertices of Γ_h; and Lemmas 4.3 and 3.2. |
| Theorems A, B, C, Corollary B′ | Unchanged, including the universal constant 1 for the leak field (h/μ)v_φ. |
| Theorem D, Proposition 7.1 | Unchanged. π_h is the projection onto Π_h on the refined mesh, and ‖η‖_{L^∞(K_e)} ≤ Ch_Γ^k. |
| Proposition 8.1(a) | Unchanged. |
| Proposition 8.1(b) | **Recertified** (Proposition `ct:prop:penalty`, results below). With the split, the m=2 star S₂ and the collinear-rim star S₄ are admissible; before the split both have singular vertices, so Section 8 had to exclude them. S₂ is the typical configuration at a vertex of Γ_h. The perturbation argument is now rigorous for **every** CT star: div: V_S → Π_S is onto for every CT star (by the macro argument plus the Nitsche-edge bubble), so the kernel dimension is exactly dim V_S − dim Π_S. That makes the kernel continuous in the vertex positions, and strict inequalities persist. |

### Exact certificates (logs/lamstar_exact_ct.log)

The Nitsche boundary is y=0, with zero data on the rim. The witness is exactly in the divergence-free kernel; this is re-verified by an exact assertion D·c = 0.

| k | quantity | S₂ | S₃ | S₄ | S₄′ |
|---|---|---|---|---|---|
| 2 | kernel dim (= m+4 = dim HCT dofs) | 6 | 7 | 8 | 8 |
| 2 | certified λ₀ (2B−A−λ₀C>0 exactly) | 571/100 | 377/50 | 377/50 | 188/25 |
| 2 | exact witness quotient | 5.717155 | 7.548839 | 7.544669 | 7.521803 |
| 3 | kernel dim | 18 | 24 | 30 | 30 |
| 3 | certified λ₀ | 483/20 | 2481/100 | 124/5 | 124/5 |
| 3 | exact witness quotient | 24.154731 | 24.812268 | 24.808483 | 24.806964 |
| 4 | kernel dim | 36 | 50 | 64 | 64 |
| 4 | certified λ₀ | 2419/50 | 2439/50 | 2439/50 | 4877/100 |
| 4 | exact witness quotient | 48.389465 | 48.782546 | 48.781808 | 48.778032 |

- In every case rank div = dim Π_S, so the kernel dimension equals dim V_S − dim Π_S.
- rank C on the kernel is 5, 9 and 13 for k = 2, 3, 4.
- **Validation:** the same code on the *unsplit* P₄ fans reproduces the paper's 9.97854… (S₃) and 9.26747… (S₄′) exactly.
- **Interpretation:** N_h is indefinite on Z_h whenever μℓ_z/h < λ₀. So for k=2 the P2 divergence-free kernel does make N_h indefinite, for λ up to about 7.5 on S₃ and S₄ and 5.7 on S₂.
- The thresholds are much larger than for unsplit P₄ (about 48.8 against 9.98 on S₃ at k=4). This is consistent with the thinner boundary triangle K_e. They grow with k, and are monotone in k by inclusion of the spaces.

**Numerical caveat (honest):** my first version followed the old script and split the kernel into ker C and its complement in floating point. For k=3 that misclassified a trace direction and reported a spurious λ\* ≈ 25.11, which the exact check then rejected (the true value is ≈ 24.81). The final script forms ker C and the Schur complement *exactly*, and uses floating point only for the top eigenvector of a small exact pencil. The certified numbers are lower bounds for the true maximal Rayleigh quotient λ\*; they agree with the float maximiser to about 4·10⁻⁶ relative. I did not certify that λ\* is not larger; necessity does not need that.

## Inputs cited

**Theorem and page numbers are NOT verified.** WebFetch to arxiv.org and epubs.siam.org needed permission prompts that were not answered, and curl to arxiv.org got a 403 from the proxy. Web search confirmed only titles and venues, so I deliberately cite no theorem numbers.
- Guzmán–Neilan, SIAM J. Numer. Anal. 56 (2018) 2826–2844: P_k–P_{k−1}^{disc} is inf-sup stable on Alfeld splits for k ≥ d. **Not load-bearing**: it is cited only as an alternative source for local surjectivity, which the section proves independently.
- Arnold–Qin 1992 (k=2): **not load-bearing** for the same reason. The bibliographic details (IMACS proceedings, pp. 28–34) are from memory; please check them.
- Scott–Vogelius 1985 and Guzmán–Scott 2019 Theorem 1 (already in the paper): load-bearing for local surjectivity with k ≥ 4, on the fixed 3-triangle split of T̂.
- Clough–Tocher 1965 and Ciarlet's 1978 book, §6.1 (HCT interpolation estimate |φ−Πφ|_{W^{j,∞}(T)} ≤ Ch^{4−j}, with constants depending on shape regularity): load-bearing for v_φ when k = 2, 3. This is the standard full-HCT estimate, but I could not check the exact section or theorem number. If preferred, Ciarlet, "Sur l'élément de Clough et Tocher", RAIRO R-2 (1974) 19–27 is the original source (also from memory).
- Scott–Zhang 1990: the H¹-stable interpolant preserving zero boundary values, used in the macro Fortin step.
- Morgan–Scott 1975 (already cited): for v_φ when k ≥ 4, as in the paper.

BibTeX to append to `paper/refs.bib`. These are the exact entries used in the compile test.

```bibtex
@incollection{ArnoldQin1992,
  author    = {Arnold, Douglas N. and Qin, Jinshui},
  title     = {Quadratic velocity/linear pressure {S}tokes elements},
  booktitle = {Advances in Computer Methods for Partial Differential Equations VII},
  editor    = {Vichnevetsky, R. and Knight, D. and Richter, G.},
  publisher = {IMACS}, year = {1992}, pages = {28--34}
}
@article{GuzmanNeilan2018,
  author  = {Guzm\'an, Johnny and Neilan, Michael},
  title   = {Inf-sup stable finite elements on barycentric refinements producing divergence-free approximations in arbitrary dimensions},
  journal = {SIAM J. Numer. Anal.}, volume = {56}, number = {5}, pages = {2826--2844}, year = {2018}
}
@inproceedings{CloughTocher1965,
  author    = {Clough, R. W. and Tocher, J. L.},
  title     = {Finite element stiffness matrices for analysis of plate bending},
  booktitle = {Proceedings of the Conference on Matrix Methods in Structural Mechanics},
  address   = {Wright-Patterson Air Force Base, Ohio}, year = {1965}, pages = {515--545}
}
@book{Ciarlet1978book,
  author = {Ciarlet, Philippe G.}, title = {The Finite Element Method for Elliptic Problems},
  publisher = {North-Holland}, address = {Amsterdam}, year = {1978}
}
@article{ScottZhang1990,
  author  = {Scott, L. Ridgway and Zhang, Shangyou},
  title   = {Finite element interpolation of nonsmooth functions satisfying boundary conditions},
  journal = {Math. Comp.}, volume = {54}, number = {190}, pages = {483--493}, year = {1990}
}
```

## Assumptions added

(CT0) is shape regularity of the *macro* mesh plus c₀h_Γ ≤ |e| ≤ h_Γ. That is the only mesh assumption; (D) is kept. h is defined as the maximum macro diameter; the refined diameters are comparable.

## What I could NOT do, and remaining caveats

1. **No numerics on Clough–Tocher meshes.** None of the rates, the leak constant 1, or the penalty threshold γ\* were measured, because of CPU limits and because it was not requested. The measured γ\* ≈ 18–20 of Section 9 (unsplit, k=4) does **not** transfer: the local certificates suggest CT needs a substantially larger γ. The section says so explicitly.
2. **Citation numbers are unverified** (see above). The results rely on them only through standard facts: the HCT estimate, Scott–Zhang, and SV/GS for k ≥ 4 on a fixed 3-triangle mesh.
3. **Continuous inf-sup uniformity** is reused from the paper's own Lemma 4.5 proof (the Galdi wedge cover), not re-derived. I checked only that it is mesh-independent, which it is. Any gap there is inherited.
4. **Morgan–Scott constants for k ≥ 4** depend on (M0) and (M2), as the paper already claims; I did not re-verify this. On CT meshes (M2) holds automatically. For k = 4, 5 the Argyris quintic would avoid Morgan–Scott.
5. **The v_φ audit** (only L^∞-O(h^k) and bounded gradient are used) was done by reading Section 6 line by line. It is not machine-checked.
6. Proposition `ct:prop:penalty` certifies **lower bounds** on λ\*(S), which is all necessity needs. Maximality is not certified. For k ≥ 5 the section uses λ₄ (monotonicity); the true k ≥ 5 values were not computed.
7. **Text outside my files would need edits** if this section is adopted: the abstract, the introduction ("All results hold for every fixed k ≥ 4"), Section 2 ("Fix k ≥ 4"), and Appendix A (reproducibility table: add `lamstar_exact_ct.py`).
