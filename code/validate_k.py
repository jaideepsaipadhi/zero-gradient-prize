"""Checks that svn_k reproduces svn at k = 4, and that the PARDISO path agrees with SuperLU.
  python3 validate_k.py  -> prints max relative differences (writes nothing)"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import svn, svn_k

ROOT = os.path.dirname(os.path.abspath(__file__))
Q = ("H1", "L2", "bL2", "bdn", "energy")


def rel(a, b):
    return max(abs(a[q] - b[q]) / max(abs(a[q]), 1e-300) for q in Q)


def absd(a, b):
    return max(abs(a[q] - b[q]) for q in Q)


print("1) svn_k (k=4, std, SuperLU) vs svn.py (in-process) and vs stored JSONL (results/pod, results/local)")
for N in (16, 24):
    for kind, store in (("B1", "pod"), ("P1", "local")):
        mu = 100.0
        ex = svn.make_exact(kind)
        v, t, c, o = svn.make_mesh(N)
        S0 = svn.Space(v, t, c, o); s0 = svn.assemble(S0, mu, ex["f"], None)
        S1 = svn_k.build(N, 4, "std"); s1 = svn_k.assemble(S1, mu, ex["f"])
        stored = [json.loads(l) for l in open(os.path.join(ROOT, "..", "results", store, "study",
                                                          f"N{N}_{kind}_mu{mu:g}.jsonl"))]
        for th in (0, 1):
            u0, *_ = svn.solve_lu(S0, s0, ex["u"], th); E0 = svn.errors(S0, u0, ex, mu)
            u1, *_ = svn_k.solve_lu(S1, s1, ex["u"], th); E1 = svn_k.errors(S1, u1, ex, mu)
            u2, *_ = svn_k.solve_pardiso(S1, s1, ex["u"], th); E2 = svn_k.errors(S1, u2, ex, mu)
            Es = stored[th]
            f = absd if (kind == "P1" and th == 1) else rel     # CNS* on P is round-off: compare absolutely
            lab = "abs" if f is absd else "rel"
            print(f"   N={N} {kind} theta={th}: svn_k-vs-svn {f(E0, E1):.1e}  svn_k-vs-stored {f(Es, E1):.1e}"
                  f"  pardiso-vs-lu {f(E1, E2):.1e}  ({lab})")
print("2) PARDISO vs SuperLU at k=6, N=32, B1, mu=1000")
ex = svn_k.make_exact("B1"); S = svn_k.build(32, 6, "std"); s = svn_k.assemble(S, 1000.0, ex["f"])
for th in (0, 1):
    u1, *_ = svn_k.solve_lu(S, s, ex["u"], th); u2, *_ = svn_k.solve_pardiso(S, s, ex["u"], th)
    print(f"   theta={th}: rel diff {rel(svn_k.errors(S, u1, ex, 1000.0), svn_k.errors(S, u2, ex, 1000.0)):.1e}")
