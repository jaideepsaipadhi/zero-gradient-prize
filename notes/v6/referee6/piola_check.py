# Check Lemma A:piola(b): int_F mu_i mu_j d_nu v dA = (hhat/h_T)^2 A int_Fhat mu_i mu_j ghat dAhat
import sympy as sp, random
x,y,z=sp.symbols('x y z'); X=sp.Matrix([x,y,z])
random.seed(3)
R=lambda: sp.Rational(random.randint(-9,9),random.randint(1,5))
A=sp.Matrix(3,3,lambda i,j:R())
while A.det()<=0: A=sp.Matrix(3,3,lambda i,j:R())
a=sp.Matrix([R(),R(),R()])
# vhat vanishing on zhat=0
vh=sp.Matrix([z*(1+2*x-y+z**2), z*(3*x*y-1+z), z*(x**2+y)])
Ainv=A.inv()
xh=Ainv*(X-a)
v=(A*vh.subs({x:sp.Symbol('X1'),y:sp.Symbol('X2'),z:sp.Symbol('X3')}))/A.det()
v=v.subs({sp.Symbol('X1'):xh[0],sp.Symbol('X2'):xh[1],sp.Symbol('X3'):xh[2]})
# physical face: image of zhat=0 triangle; inward normal nu
P=[A*sp.Matrix(p)+a for p in [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]]
nrm=(P[2]-P[1]).cross(P[3]-P[1]) if False else (P[1]-P[0]).cross(P[2]-P[0])
if nrm.dot(P[3]-P[0])<0: nrm=-nrm
nu=nrm/sp.sqrt(nrm.dot(nrm)); area=sp.sqrt(nrm.dot(nrm))/2
hT=abs((P[3]-P[0]).dot(nu))
s,t=sp.symbols('s t')
mu=[1-s-t,s,t]
pt=P[0]+s*(P[1]-P[0])+t*(P[2]-P[0])
dnu=(v.jacobian(X)*nu).subs({x:pt[0],y:pt[1],z:pt[2]})
gh=vh.diff(z).subs(z,0).subs({x:s,y:t})
for (i,j) in [(0,1),(0,2),(1,2)]:
    lhs=sp.Matrix([sp.integrate(sp.integrate(mu[i]*mu[j]*dnu[c],(t,0,1-s)),(s,0,1)) for c in range(3)])*2*area
    rhs=(1/hT)**2*A*sp.Matrix([sp.integrate(sp.integrate(mu[i]*mu[j]*gh[c],(t,0,1-s)),(s,0,1)) for c in range(3)])
    print((i,j),sp.simplify(lhs-rhs).T)
