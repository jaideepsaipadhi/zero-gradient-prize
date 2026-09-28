"""referee9: independent exact check of robust4 Lemma dk / Lemma KT for k=4,5 on generic rational triangles.
T = L(0,0), R(1,0), c(a,r). eta_H = (H - Pi_{P_{k-1}(T)} H)|_e, e param s: point (s,0), lambda_L=1-s, lambda_R=s."""
import sympy as sp
x,y,s,u,w=sp.symbols('x y s u w')
def J(k):
    c=sp.symbols('c0:%d'%k); P=s**k+sum(c[i]*s**i for i in range(k))
    sol=sp.solve([sp.integrate(P*(1-s)*s**j,(s,0,1)) for j in range(k)],c,dict=True)[0]
    return sp.expand(P.subs(sol))
for k in (4,5):
    Jk=J(k); d=sp.expand(Jk.subs(s,1-s)-Jk); beta=sp.integrate(Jk,(s,0,1))
    g=sum(sp.Symbol('g%d'%i)*s**i for i in range(k+1))
    lhs=sp.expand(sp.integrate(d*g,(s,0,1))); rhs=sp.expand(beta*(g.subs(s,0)-g.subs(s,1)))
    print('k',k,'J(0)',Jk.subs(s,0),'beta',beta,'beta==|J(0)|/(k+1)',beta==abs(Jk.subs(s,0))/(k+1),'identity',sp.simplify(lhs-rhs)==0)
    for (a,r) in [(sp.Rational(1),sp.Rational(1)),(sp.Rational(2,7),sp.Rational(3,5)),(sp.Rational(-1,3),sp.Rational(5,4))]:
        def intT(f):
            return sp.integrate(sp.integrate(sp.expand(f.subs({x:u+w*a,y:w*r},simultaneous=True)*r),(u,0,1-w)),(w,0,1))
        mons=[x**i*y**j for i in range(k) for j in range(k-i)]
        G=sp.Matrix(len(mons),len(mons),lambda i,j:intT(mons[i]*mons[j]))
        def eta(H):
            cc=G.LUsolve(sp.Matrix([intT(H*m) for m in mons]))
            return sp.expand((H-sum(cc[i]*mons[i] for i in range(len(mons)))).subs(y,0).subs(x,s))
        # linear parts of barycentrics
        lL=-(x-a*y/r)*1 - y/r*0  # lambda_L = 1 - x + (a-1)... compute properly
        M=sp.Matrix([[0,1,a],[0,0,r],[1,1,1]])  # columns L,R,c ; rows x,y,1
        Minv=M.inv()  # lambda = Minv*[x,y,1]
        lam=[sp.expand(Minv[i,0]*x+Minv[i,1]*y) for i in range(3)]  # linear parts
        eL,eR,ec=[eta(l**k) for l in lam]
        print('  (a,r)=',(a,r),' eta(lc^k)==J(0):',sp.simplify(ec-Jk.subs(s,0))==0,
              ' eta(lL^k-lR^k)==d_k:',sp.simplify(eL-eR-d)==0)
