"""referee11: INDEPENDENT exact check of robust6's corrected (P^T) fields v_1, v_2 on Q = T1 u T2 (k=4).
Only the BB coefficients are taken from the object under test (robust5/pt_fields.star_fields with UREF := robust6/uref_fix.get()).
Everything else is my own: Fraction polynomial arithmetic in (x,y), barycentrics by Cramer, triangle integrals by affine
pull-back + int_{T^} s^i t^j = i! j!/(i+j+2)!, edge integrals by parametrisation. No sympy integration, no BB integration formulas.
usage: python verify_PT.py x1 y1 x2 y2   (rationals)"""
import sys
from fractions import Fraction as Fr
from math import factorial as fa
sys.path.insert(0, '../robust5'); sys.path.insert(0, '../robust6')
import pt_fields as PF
import uref_fix

# ---------- polynomials: dict (i,j) -> Fr  in x^i y^j
def padd(*ps):
    o = {}
    for p in ps:
        for k, v in p.items(): o[k] = o.get(k, 0) + v
    return {k: v for k, v in o.items() if v != 0}
def pscale(p, c): return {k: v * c for k, v in p.items() if v * c != 0}
def pmul(p, q):
    o = {}
    for (a, b), u in p.items():
        for (c, d), w in q.items(): o[(a + c, b + d)] = o.get((a + c, b + d), 0) + u * w
    return {k: v for k, v in o.items() if v != 0}
def ppow(p, n):
    o = {(0, 0): Fr(1)}
    for _ in range(n): o = pmul(o, p)
    return o
def dx(p): return {(i - 1, j): v * i for (i, j), v in p.items() if i > 0}
def dy(p): return {(i, j - 1): v * j for (i, j), v in p.items() if j > 0}
def subst(p, X, Y):
    """p(X(s,t), Y(s,t)), X,Y polys in (s,t) stored as (i,j) dicts"""
    o = {}
    for (i, j), v in p.items(): o = padd(o, pscale(pmul(ppow(X, i), ppow(Y, j)), v))
    return o
def aff(c0, cs, ct): return {k: v for k, v in {(0, 0): Fr(c0), (1, 0): Fr(cs), (0, 1): Fr(ct)}.items() if v != 0}

def int_tri(p, P0, P1, P2):
    X = aff(P0[0], P1[0] - P0[0], P2[0] - P0[0]); Y = aff(P0[1], P1[1] - P0[1], P2[1] - P0[1])
    det = (P1[0] - P0[0]) * (P2[1] - P0[1]) - (P2[0] - P0[0]) * (P1[1] - P0[1])
    q = subst(p, X, Y)
    return abs(det) * sum(v * Fr(fa(i) * fa(j), fa(i + j + 2)) for (i, j), v in q.items())
def restrict(p, P, Q):
    """p on segment P + t(Q-P): univariate poly as dict (i,0)"""
    return subst(p, aff(P[0], Q[0] - P[0], 0), aff(P[1], Q[1] - P[1], 0))
def int_seg_param(p, P, Q):
    """int_0^1 p(P+t(Q-P)) dt (multiply by |Q-P| for ds)"""
    return sum(v * Fr(1, i + 1) for (i, j), v in restrict(p, P, Q).items())

def bary(P0, P1, P2):
    """barycentric coordinates as affine polys in x,y (Cramer)"""
    det = (P1[0] - P0[0]) * (P2[1] - P0[1]) - (P2[0] - P0[0]) * (P1[1] - P0[1])
    Ps = [P0, P1, P2]; out = []
    for i in range(3):
        A, B = Ps[(i + 1) % 3], Ps[(i + 2) % 3]
        # lambda_i = area(x, A, B)/area(P_i, A, B) ; area(x,A,B) = (A-x)x(B-x)
        c0 = A[0] * B[1] - B[0] * A[1]; cx = A[1] - B[1]; cy = B[0] - A[0]
        out.append(pscale({(0, 0): c0, (1, 0): cx, (0, 1): cy}, Fr(1) / det))
    return out

def bb_to_poly(coef, L, n=4):
    o = {}
    for al, v in coef.items():
        if v == 0: continue
        m = Fr(fa(n), fa(al[0]) * fa(al[1]) * fa(al[2]))
        o = padd(o, pscale(pmul(pmul(ppow(L[0], al[0]), ppow(L[1], al[1])), ppow(L[2], al[2])), m * v))
    return o

def main(x1, y1, x2, y2):
    # sanity of my integrators
    e0, e1 = (Fr(0), Fr(0)), (Fr(1), Fr(0))
    assert int_seg_param({(2, 1): Fr(1)}, e0, e1) == 0, "edge integral of y*x^2 on y=0 must vanish"
    assert int_seg_param({(3, 0): Fr(1)}, e0, e1) == Fr(1, 4)
    assert int_tri({(0, 0): Fr(1)}, (0, 0), (2, 0), (0, 3)) == 3
    assert int_tri({(1, 1): Fr(1)}, (0, 0), (1, 0), (0, 1)) == Fr(1, 24)

    PF.UREF = uref_fix.get()
    G, V, info = PF.star_fields(x1, y1, x2, y2, PF.NumExact)
    s1, s2 = x1 + y1, x2 + y2
    z = (Fr(0), Fr(0)); a = (Fr(1), Fr(0)); c = (x1 / s1, 1 / s1)
    b = ((x2 * c[0] - c[1]) / s2, (x2 * c[1] + c[0]) / s2)
    # sanity: cotangents of recomputed geometry
    def cot(P, A, B):
        u = (A[0] - P[0], A[1] - P[1]); w = (B[0] - P[0], B[1] - P[1])
        return (u[0] * w[0] + u[1] * w[1]) / (u[0] * w[1] - u[1] * w[0])
    assert cot(z, a, c) == x1 and cot(a, c, z) == y1 and cot(z, c, b) == x2 and cot(c, b, z) == y2, "geometry"
    L1 = bary(z, a, c); L2 = bary(z, c, b)
    B = [pscale(pmul(ppow({(1, 0): Fr(1)}, j), ppow({(0, 0): Fr(1), (1, 0): Fr(-1)}, 4 - j)), fa(4) // (fa(j) * fa(4 - j))) for j in range(5)]
    targets = [padd(B[1], pscale(B[2], -1)), padd(B[1], pscale(B[3], -1))]
    allok = True
    for f, F in enumerate(V):
        v1 = [bb_to_poly(F[1][d], L1) for d in range(2)]
        v2 = [bb_to_poly(F[2][d], L2) for d in range(2)]
        R = {}
        R["nonzero"] = any(v1) or any(v2)
        R["cont [z,c]"] = all(padd(restrict(p, z, c), pscale(restrict(q, z, c), -1)) == {} for p, q in zip(v1, v2))
        R["0 on [a,c]"] = all(restrict(p, a, c) == {} for p in v1)
        R["0 on [c,b]"] = all(restrict(p, c, b) == {} for p in v2)
        R["0 on [z,b]"] = all(restrict(p, z, b) == {} for p in v2)
        R["div=0 on T2"] = padd(dx(v2[0]), dy(v2[1])) == {}
        div1 = padd(dx(v1[0]), dy(v1[1]))
        vn = pscale(v1[1], -1)                      # n = (0,-1) outward on e
        rs = [{(p, q): Fr(1)} for p in range(4) for q in range(4 - p)]
        pb = []; pbg = []; green = []
        for r in rs:
            lhs = int_tri(pmul(div1, r), z, a, c)
            rhs = int_seg_param(pmul(vn, r), z, a)  # |e| = 1
            pb.append(lhs == rhs)
            vg = int_tri(padd(pmul(v1[0], dx(r)), pmul(v1[1], dy(r))), z, a, c)
            pbg.append(vg == 0)
            # Green on T1 with ALL three edges (outward normals, unnormalised: n ds = (dy,-dx) dt along ccw boundary)
            bd = 0
            for P, Q in ((z, a), (a, c), (c, z)):
                nx, ny = Q[1] - P[1], -(Q[0] - P[0])
                bd += int_seg_param(pmul(padd(pscale(v1[0], nx), pscale(v1[1], ny)), r), P, Q)
            green.append(lhs == -vg + bd)
        R["PB T1: (div v,r)=<v.n,r>_e, 10 r"] = all(pb)
        R["PB T1: (v,grad r)=0"] = all(pbg)
        R["Green consistency"] = all(green)
        # T2 moments directly too
        div2 = padd(dx(v2[0]), dy(v2[1]))
        R["PB T2 moments"] = all(int_tri(pmul(div2, r), z, c, b) == 0 for r in rs)
        R["(M) int_e v = 0"] = all(int_seg_param(p, z, a) == 0 for p in v1)
        R["(D) int_e d_n v = 0"] = all(int_seg_param(pscale(dy(p), -1), z, a) == 0 for p in v1)
        R["trace"] = padd(restrict(vn, z, a), pscale(targets[f], -1)) == {}
        # zero flux on all of dQ
        R["flux e = 0"] = int_seg_param(vn, z, a) == 0
        print(f"shape {x1},{y1},{x2},{y2} v_{f+1}: ", R)
        allok &= all(R.values())
    print("mA,mC,dA,dC =", info["mA"], info["mC"], info["dA"], info["dC"], " s1 =", s1)
    print("ALL OK" if allok else "FAILURE", flush=True)
    return allok

if __name__ == "__main__":
    args = [Fr(q) for q in sys.argv[1:5]]
    main(*args)
