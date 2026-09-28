"""referee11: independent check of nearsing Thm N2: d_z^2 = min_t sum_i |T_i| |H_i(t) - a P_nu|_F^2 for m=3 wall stars,
families I (needle at the wall) and II (split edge). Own formulation: kernel x of [vecP1..vecP4] (report's route), wall lines at
angles 0 and phi (line angles), bisector normal nu, a=1, ray lengths r_j (random in [c_r,1]). mpmath for tiny eps."""
import mpmath as mp, random
mp.mp.dps = 50
def Pv(psi):  # vec of P = n n^T for line at angle psi (normal n = (-sin, cos)), isometric vec
    s, c = mp.sin(psi), mp.cos(psi)
    return [s*s, c*c, -mp.sqrt(2)*s*c]
def det3(a, b, c):
    return a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0])
def dz(lines, thetas, r, phi):
    V = [Pv(p) for p in lines]
    x = [det3(V[1], V[2], V[3]), -det3(V[0], V[2], V[3]), det3(V[0], V[1], V[3]), -det3(V[0], V[1], V[2])]
    H = [[x[0]*v for v in V[0]], [x[0]*V[0][k] + x[1]*V[1][k] for k in range(3)], [-x[3]*v for v in V[3]]]
    # check closure
    H3b = [x[0]*V[0][k] + x[1]*V[1][k] + x[2]*V[2][k] for k in range(3)]
    assert max(abs(H3b[k] - H[2][k]) for k in range(3)) < mp.mpf(10)**-40 * (1 + max(abs(q) for q in x))
    g = Pv(phi/2)  # P_nu, nu = bisector normal
    area = [r[i]*r[i+1]*mp.sin(thetas[i])/2 for i in range(3)]
    # min_t sum A_i |t H_i - g|^2
    hh = sum(A*sum(h*h for h in Hi) for A, Hi in zip(area, H)); hg = sum(A*sum(h*q for h, q in zip(Hi, g)) for A, Hi in zip(area, H))
    gg = sum(area) * sum(q*q for q in g)
    return mp.sqrt(gg - hg*hg/hh), mp.sqrt(sum(area))
random.seed(1)
lo, hi = mp.inf, 0; worst = None
for fam in ("I", "II"):
    flo, fhi = mp.inf, 0
    for phi in [mp.mpf('0.08'), mp.mpf('0.01'), mp.mpf('1e-3'), mp.mpf('1e-5')]:
        for e in [mp.mpf('0.5'), mp.mpf('0.1'), phi, phi**1.5, phi**2, phi**2/10, phi**3, phi**4]:
            for gamma in [mp.pi/3, mp.pi/2, 2*mp.pi/3 - e]:
                for trial in range(3):
                    r = [mp.mpf(random.uniform(0.5, 1)) for _ in range(4)] if trial else [mp.mpf(1)]*4
                    if fam == "II":
                        lines = [0, gamma, gamma + e, phi]; th = [gamma, e, mp.pi + phi - gamma - e]
                    else:
                        lines = [0, gamma, phi - e, phi]; th = [gamma, mp.pi + phi - e - gamma, e]
                    if min(th) < e/2: continue
                    d, _ = dz(lines, th, r, phi)
                    w = min(1, phi/mp.sqrt(e))
                    q = d/w
                    if q < flo: flo = q
                    if q > fhi: fhi = q
    print(f"family {fam}: d_z / min(1, phi/sqrt(eps)) in [{mp.nstr(flo,4)}, {mp.nstr(fhi,4)}]  (a=1, h~1, rays in [0.5,1], phi in [1e-5,0.08], eps in [phi^4, 0.5])")
# compare with max-norm style phi/eps weight to show it is the wrong one in H^1
phi = mp.mpf('1e-3')
for e in [phi, phi**1.5, phi**2]:
    d, _ = dz([0, mp.pi/2, mp.pi/2+e, phi], [mp.pi/2, e, mp.pi/2+phi-e], [1]*4, phi)
    print(f"II phi=1e-3 eps={mp.nstr(e,3)}: d={mp.nstr(d,4)}  phi/sqrt(eps)={mp.nstr(phi/mp.sqrt(e),4)}  min(1,phi/eps)={mp.nstr(min(1,phi/e),4)}")
# exact kernel identities of the N2 proof
for (g, e, p) in [(mp.mpf('1.1'), mp.mpf('0.013'), mp.mpf('0.004')), (mp.mpf('0.9'), mp.mpf('1e-4'), mp.mpf('0.02'))]:
    V = [Pv(q) for q in (0, g, g+e, p)]
    x = [det3(V[1], V[2], V[3]), -det3(V[0], V[2], V[3]), det3(V[0], V[1], V[3]), -det3(V[0], V[1], V[2])]
    print("II |x2/x1| =", mp.nstr(abs(x[1]/x[0]), 12), " formula", mp.nstr(mp.sin(g+e)*mp.sin(p)/(mp.sin(e)*mp.sin(g-p)), 12), " |x4/x1+1| =", mp.nstr(abs(x[3]/x[0]+1), 6))
    V = [Pv(q) for q in (0, g, p-e, p)]
    x = [det3(V[1], V[2], V[3]), -det3(V[0], V[2], V[3]), det3(V[0], V[1], V[3]), -det3(V[0], V[1], V[2])]
    print("I |-x4/x1-1| =", mp.nstr(abs(-x[3]/x[0]-1), 12), " formula", mp.nstr(mp.sin(p)*abs(mp.sin(g-e))/(mp.sin(g-p)*mp.sin(e)), 12), " |x2/x1| =", mp.nstr(abs(x[1]/x[0]), 6))
