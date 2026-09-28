"""Exact certificates for notes/v6/penalty/REPORT.md.  Usage: python certs.py [neg|pos|ct|high]"""
import sys
from fractions import Fraction as Fr
from exact import certify_negative, certify_positive, fan_patch, ct_patch


def rect_star(H):
    return fan_patch([(1, 0), (1, H), (0, H), (-1, 0)])


def stars_union(H, W):
    """union of the stars of the boundary vertices (0,0),...,(W-1,0) of the column mesh
    (columns [i,i+1]x[0,H] split by the forward diagonal)."""
    P = []; ix = {}
    def v(p):
        if p not in ix:
            ix[p] = len(P); P.append(p)
        return ix[p]
    tris = []
    for i in range(W):
        for t in (((i - 1, 0), (i, 0), (i, H)), ((i, 0), (i + 1, 0), (i + 1, H)), ((i, 0), (i + 1, H), (i, H))):
            tt = tuple(v(p) for p in t)
            if tt not in tris:
                tris.append(tt)
    nit = {(v((i, 0)), v((i + 1, 0))) for i in range(-1, W)}
    return P, tris, nit


mode = sys.argv[1]
if mode == "neg":
    certify_negative(*rect_star(5), 4, lam1=0, label="rect star H=5")
    certify_negative(*rect_star(5), 4, lam1=-1, label="rect star H=5")
    certify_negative(*rect_star(Fr(9, 2)), 4, lam1=Fr(-1, 2), label="rect star H=9/2")
    certify_negative(*stars_union(8, 2), 4, lam1=-1, label="edge patch (2 stars) H=8")
if mode == "pos":
    certify_positive(*stars_union(5, 3), 4, Fr(1), label="3 stars H=5")
    certify_positive(*stars_union(8, 4), 4, Fr(1, 4), label="4 stars H=8")
if mode == "high":
    certify_negative(*rect_star(8), 5, lam1=-1, label="rect star H=8")
    certify_negative(*rect_star(12), 6, lam1=-1, label="rect star H=12")
if mode == "ct":
    certify_negative(*ct_patch(*rect_star(2)), 2, lam1=-1, label="CT rect star H=2")
    certify_negative(*ct_patch(*rect_star(4)), 3, lam1=Fr(-1, 2), label="CT rect star H=4")
    certify_negative(*ct_patch(*rect_star(8)), 4, lam1=-1, label="CT rect star H=8")
if mode == "controls":
    # negative controls: must print False
    certify_negative(*rect_star(5), 4, lam1=Fr(-6, 5), label="CONTROL rect star H=5 (lambda*=-1.1804 > -1.2)")
    certify_negative(*fan_patch([(1, 0), (1, 1), (-1, 1), (-1, 0)]), 4, lam1=0, label="CONTROL S3 (lambda*=9.98)")
    certify_negative(*rect_star(4), 4, lam1=0, label="CONTROL rect star H=4 (lambda*=+0.038)")
if mode == "edge6":
    certify_negative(*stars_union(6, 2), 4, lam1=Fr(-1, 10), label="edge patch (2 stars) H=6")
    certify_negative(*stars_union(6, 2), 4, lam1=Fr(-3, 20), label="CONTROL edge patch H=6 (lambda*=-0.1404 > -0.15)")
