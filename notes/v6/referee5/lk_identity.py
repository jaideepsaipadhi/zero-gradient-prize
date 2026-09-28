"""Referee5: exact check of <eta, L_k^e>_e = |e|^{k+1} k!/(2k+1)! * Beta(k+1,k+1)-avg of d_tau^k phi,
eta = phi - pi_T phi (L2(T) projection onto P_{k-1}); several third vertices C for the same edge e."""
import sympy as sp, random
random.seed(7)
import sys; KS=[int(a) for a in sys.argv[1].split(",")]; NT=int(sys.argv[2])
x, y, s, t, u = sp.symbols('x y s t u')
R = lambda: sp.Rational(random.randint(-9, 9), random.randint(1, 7))
def proj(phi, V, k):
    A, B, C = [sp.Matrix(v) for v in V]
    X = A + s*(B-A) + u*(C-A); J = abs((B-A).row_join(C-A).det())
    mons = [x**i*y**j for i in range(k) for j in range(k-i)]
    sub = lambda f: sp.expand(f.subs({x: X[0], y: X[1]}))
    I = lambda f: J*sp.integrate(sp.integrate(sub(f), (u, 0, 1-s)), (s, 0, 1))
    G = sp.Matrix(len(mons), len(mons), lambda i, j: I(mons[i]*mons[j]))
    b = sp.Matrix([I(phi*m) for m in mons])
    c = G.LUsolve(b)
    return sum(ci*m for ci, m in zip(c, mons))
ok = True
for k in KS:
    phi = sum(R()*x**i*y**j for i in range(k+2) for j in range(k+2-i))
    A, B = (R(), R()), (R()+3, R()+1)
    Lk = sp.expand(sp.diff((t**2-t)**k, t, k)/sp.factorial(k))      # shifted Legendre, L_k(1)=1
    xt = {x: A[0]+t*(B[0]-A[0]), y: A[1]+t*(B[1]-A[1])}
    rhs = sp.integrate(sp.diff(sp.expand(phi.subs(xt)), t, k)*t**k*(1-t)**k, (t, 0, 1))/sp.factorial(k)
    # |e|^{k+1} k!/(2k+1)! * (2k+1)!/(k!)^2 * int d_tau^k phi t^k(1-t)^k  == |e| * rhs  (d/dt = |e| d_tau)
    vals = []
    for trial in range(NT):
        C = (sp.Rational(random.randint(-20, 20), 7), sp.Rational(random.randint(5, 30), 4))
        pi = proj(phi, [A, B, C], k)
        eta = sp.expand((phi-pi).subs(xt))
        lhs = sp.integrate(eta*Lk, (t, 0, 1))                          # /|e|
        etaL2 = sp.integrate(eta**2, (t, 0, 1))
        vals.append((lhs, etaL2))
        ok &= (sp.simplify(lhs-rhs) == 0)
    print(f"k={k}: lhs={[str(v[0]) for v in vals]} rhs={rhs}  |eta|^2 differ: {len(set(v[1] for v in vals))==NT}  Lk(0)={Lk.subs(t,0)}")
print("ALL EXACT:", ok)
