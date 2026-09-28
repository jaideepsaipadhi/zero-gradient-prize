# referee8 scratch
- ref8_vertex_fe.py/.log: FE (patchlib) non-flat sphere stars, sphere weights -> free-rank 3 in 2/3 cases => "mod D" test vacuous there (same for td_k3 generic-w cases with D-rank 3). Quad (faces 0,2; rank 2) non-vacuous, passes.
- ref8_vertex_fe_eqw.log: non-flat stars, distinct per-face equal weights (D=0): FE moment == independent formula to 1e-15 (3 stars incl. tau~0.5 hexagon). Non-vacuous confirmation of Thm vertex.
- ref8_lemmaK.py: sympy re-derivation of Lemma K (K = -wab/720 I mod D) and jet lemma (3 random exact geometries) -- OK.
- Criterion: 1 - tau(tau+rho) bound is for t in n_z^perp; for t || F extra error ~ sin^2 tau (1+rho) not covered by <=1/2 (e.g. tau=.7, rho~0). Asymptotically harmless.
- td_lift Ex.1/2: 4 points on unit sphere whose tangential projections form a square centred on the axis are equidistant from the axis => same height => coplanar. P' degenerate. Fix: off-centre / perturb one vertex.
