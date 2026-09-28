"""Assemble the explicit lower bound for Y(2B-A)/C of v = curl(I_Argyris F) and solve for kappa_0(theta).
Reads consts.json (argyris_consts.py). All formulas are those of REPORT.tex, Sections 2-3; every quantity is an
upper bound computed in floating point from rigorous rational/sup-norm inputs (float rounding is ~1e-15 relative,
and every final kappa_0 is rounded DOWN to 3 significant digits, which dominates it)."""
import json, math, sys, os

C = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "consts.json")))
m = {int(k): v for k, v in C["m"].items()}
Ah, Bh, Ch = C["Ahat"], C["Bhat"], C["Chat"]
Pv, Ne = C["Pv"], C["Ne"]
cHf, cHh = C["cH_full"], C["cH_half"]


def E(theta):
    """E_1, E_2 with |D^j(F - I F)| <= E_j M6 h_T^{6-j} on a triangle with min angle >= theta."""
    st = math.sin(theta)
    sa = min(math.sqrt(3) / 2, math.sin(2 * theta))    # sin of the largest angle
    Jinv_h = math.sqrt(2) / (st * sa)                   # h_T * ||J^{-1}||
    g_hyp = 2 / 120 + 2 * cHf                           # |gamma_hyp| / (M6 h^6)
    g_leg = 0.5**5 / 120 + (math.sqrt(2) / st) * cHh    # |gamma_leg| / (M6 h^6)
    out = []
    for j in (1, 2):
        q = Pv[j - 1] + g_hyp * Ne[j - 1][0] + g_leg * (Ne[j - 1][1] + Ne[j - 1][2])
        out.append(1 / math.factorial(6 - j) + Jinv_h**j * q)
    return out


def Qlb(kappa, theta, Y):
    E1, E2 = E(theta)
    m6 = m[6]
    Xp = (2 + kappa) * Y                                 # X' = X + kappa Y, X = 2Y
    ell = kappa * Y                                      # local boundary edges have length <= kappa Y
    if ell > 1:
        return -1e9
    r = 1 - ell**2 / 4
    if Xp / r >= 1:
        return -1e9
    wX = math.asin(Xp / r)
    delta = ell**2 / 4 + wX**2 / 2
    sig = math.tan(wX + math.asin(ell / 2))
    if delta / Y > 1 / 8:
        return -1e9
    S = 2 * Xp * ((1 + kappa) * Y + delta)              # area bound of the support patch
    L = 2 * Xp * math.sqrt(1 + sig**2)                   # length bound of Gamma_h in the patch
    m1, m2, m3 = m[1], m[2], m[3]
    # geometry (continuum field w on the true Omega_h, Gamma_h), in units where A,B dimensionless, C ~ Y
    A_sl = 16 * m2**2 * delta / Y
    eB = 4 * (sig * m1 * m2 + (delta / Y) * (m1 * m3 + m2**2))
    eC = 4 * Y * (m1**2 * (math.sqrt(1 + sig**2) - 1) + 2 * math.sqrt(1 + sig**2) * m1 * m2 * delta / Y)
    Chmax = Y * Ch + eC
    # interpolation
    etaA = 2 * E2 * m6 * kappa**4 * math.sqrt(S / Y**2)
    etaC = E1 * m6 * kappa**5 * math.sqrt(L)
    eBi = E2 * m6 * kappa**4 / Y * math.sqrt(L) * (math.sqrt(Chmax) + etaC) + m2 / Y * math.sqrt(L) * etaC
    num = 2 * (Bh - eB - eBi) - (math.sqrt(Ah + A_sl) + etaA)**2
    den = (math.sqrt(Chmax) + etaC)**2
    return Y * num / den


def kappa0(theta, Y0, target=15.0):
    # Qlb is decreasing in kappa and (for fixed kappa) in Y on the range used; check the monotonicity on a grid.
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        ok = Qlb(mid, theta, Y0) >= target
        lo, hi = (mid, hi) if ok else (lo, mid)
    return lo


if __name__ == "__main__":
    print("E_j(theta):")
    for deg in (10, 15, 20, 30, 45, 60):
        print(f"  theta={deg:2d}  E1={E(math.radians(deg))[0]:.4g}  E2={E(math.radians(deg))[1]:.4g}")
    print("\nQlb(kappa -> 0, Y -> 0) =", Qlb(1e-9, math.radians(30), 1e-12), "(continuum 15.4873)")
    for Y0 in (1e-12, 1e-6, 3e-6):
        print(f"\nY0={Y0:g}:  kappa_0(theta) with Qlb >= 15")
        for deg in (10, 15, 20, 30, 45, 60):
            th = math.radians(deg)
            k0 = kappa0(th, Y0)
            # monotonicity check in Y below Y0 and kappa below k0
            chk = min(Qlb(k0 * a, th, Y0 * b) for a in (0.1, 0.5, 0.9, 1.0) for b in (1e-6, 0.01, 0.5, 1.0))
            print(f"  theta={deg:2d}  kappa_0={k0:.4g}   min Qlb on grid below = {chk:.4f}   layers ~3/kappa0 = {3/max(k0,1e-300):.0f}")
