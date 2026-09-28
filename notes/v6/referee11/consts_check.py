"""referee11: independent spot checks of consts/REPORT.tex: Hermite midpoint formula and c_H, half-plane values,
E_j(theta) from the printed formula, and Q_lb at (kappa_0(theta), Y0=1e-6) recomputed from the report's printed formulas."""
import sympy as sp, math
s = sp.symbols('s'); f = sp.Function('f')
c = sp.symbols('c0:6'); p = sum(ci*s**i for i, ci in enumerate(c))
f0, f1, d0, d1, e0, e1 = sp.symbols('f0 f1 d0 d1 e0 e1')
sol = sp.solve([p.subs(s,0)-f0, p.subs(s,1)-f1, sp.diff(p,s).subs(s,0)-d0, sp.diff(p,s).subs(s,1)-d1, sp.diff(p,s,2).subs(s,0)-e0, sp.diff(p,s,2).subs(s,1)-e1], c)
Hd = sp.expand(sp.diff(p, s).subs(sol).subs(s, sp.Rational(1, 2)))
print("Hermite (Hf)'(1/2) =", Hd, "  report:", sp.expand(sp.Rational(15,8)*(f1-f0)-sp.Rational(7,16)*(d0+d1)+sp.Rational(1,32)*(e1-e0)) == Hd)
cH = sp.Rational(15,8)*2/720 + sp.Rational(7,16)*2/120 + sp.Rational(1,32)*2/24
print("c_H =", cH, "(report 29/1920)", cH == sp.Rational(29,1920))
# half-plane values
x, y = sp.symbols('x y', real=True)
g = y*(1-y)**6; ph = (1-(x/2)**2)**6; F = ph*g
w = (sp.diff(F, y), -sp.diff(F, x))
Gr = [[sp.diff(w[i], v) for v in (x, y)] for i in range(2)]
D = [[Gr[i][j]+Gr[j][i] for j in range(2)] for i in range(2)]
A = sp.integrate(sp.integrate(sp.expand(sum(D[i][j]**2 for i in range(2) for j in range(2))/2), (y, 0, 1)), (x, -2, 2))
wn = [q.subs(y, 0) for q in w]; dnw = [-sp.diff(q, y).subs(y, 0) for q in w]
B = sp.integrate(sp.expand(sum(a*b for a, b in zip(dnw, wn))), (x, -2, 2)); C = sp.integrate(sp.expand(sum(a*a for a in wn)), (x, -2, 2))
print("Ahat,Bhat,Chat =", A, B, C, float(A), float(B), float(C), " Qhat =", float((2*B-A)/C))
Ah, Bh, Ch = float(A), float(B), float(C)
m1, m2, m3, m6 = 3.6635, 32.187, 214.88, 18918.2
P1, P2 = 0.19597, 1.40308; N1 = (1.8857, 3.7713, 3.7713); N2 = (16.866, 33.731, 33.731); cHf = 29/1920
def E(j, th):
    sth = math.sin(th); st = min(math.sqrt(3)/2, math.sin(2*th)); P = P1 if j == 1 else P2; N = N1 if j == 1 else N2
    return 1/math.factorial(6-j) + (math.sqrt(2)/(sth*st))**j * (P + (1/60+2*cHf)*N[0] + (1/3840 + math.sqrt(2)*cHf/(2*sth))*(N[1]+N[2]))
for d in (10, 15, 20, 30, 45, 60): print(f"theta={d}: E1={E(1, math.radians(d)):.4g} E2={E(2, math.radians(d)):.4g}")
def Qlb(k, th, Y):
    Xp = (2+k)*Y; om = math.asin(Xp/(1-k*k*Y*Y/4)); dl = k*k*Y*Y/4 + om*om/2; sg = math.tan(om + math.asin(k*Y/2))
    S = 2*Xp*((1+k)*Y + dl); Lam = 2*Xp*math.sqrt(1+sg*sg)
    Asl = 16*m2**2*dl/Y
    eB = 4*(sg*m1*m2 + dl/Y*(m1*m3 + m2**2)); eC = 4*Y*(m1**2*(math.sqrt(1+sg*sg)-1) + 2*math.sqrt(1+sg*sg)*m1*m2*dl/Y)
    e1, e2 = E(1, th), E(2, th)
    etaA = 2*e2*m6*k**4*math.sqrt(S)/Y; etaC = e1*m6*k**5*math.sqrt(Lam)
    eBp = e2*m6*k**4*math.sqrt(Lam)/Y*(math.sqrt(Y*Ch+eC)+etaC) + m2*math.sqrt(Lam)/Y*etaC
    return Y*(2*(Bh-eB-eBp) - (math.sqrt(Ah+Asl)+etaA)**2)/(math.sqrt(Y*Ch+eC)+etaC)**2
for d, k0 in ((10, .00379), (15, .00594), (20, .00803), (30, .0117), (45, .0143), (60, .0161)):
    print(f"theta={d}: Qlb(kappa0={k0}, Y0=1e-6) = {Qlb(k0, math.radians(d), 1e-6):.5f};  Qlb(1.05*kappa0) = {Qlb(1.05*k0, math.radians(d), 1e-6):.4f}; Qlb(kappa0, Y=1e-9)={Qlb(k0, math.radians(d), 1e-9):.4f}")
