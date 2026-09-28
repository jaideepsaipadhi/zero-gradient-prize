# Referee 8: independent sympy check of Lemma K (reference triangle), and the jet Lemma (G_i) symbolically
import sympy as sp
x,y=sp.symbols('x y'); wza,wzb,wab,e1,e2=sp.symbols('wza wzb wab e1 e2')
mons=[1,x,y,x*x,x*y,y*y]; c=sp.symbols('c0:12')
g=[sum(c[k]*mons[k] for k in range(6)), sum(c[6+k]*mons[k] for k in range(6))]
ev=lambda p,X,Y: [gi.subs({x:X,y:Y}) for gi in p]
H=sp.Rational(1,2)
eqs=ev(g,1,0)+ev(g,0,1)
gm=ev(g,H,H); eqs.append(gm[0]+gm[1])          # parallel to b-a=(-1,1): g.(1,1)=0
g0=ev(g,0,0); eqs+= [g0[0]-e1,g0[1]-e2]
eqs.append(ev(g,H,0)[1]); eqs.append(ev(g,0,H)[0])  # parallel to za, zb
I=lambda f: sp.integrate(sp.integrate(f,(y,0,1-x)),(x,0,1))
eqs+=[I(g[0]),I(g[1])]
sol=sp.solve(eqs,c,dict=True)[0]; free=[s for s in c if s not in sol]
print('free params:',free)
mz,ma,mb=1-x-y,x,y
Q=sp.Rational(1,2)*(wza*mz*ma+wzb*mz*mb+wab*ma*mb)
mom=[sp.expand(I(Q*gi.subs(sol))) for gi in g]
base=[sp.simplify(m.subs({f:0 for f in free})) for m in mom]
Dir=[sp.simplify(sp.diff(m,free[0])) for m in mom]
print('moment at free=0:',base); print('free direction (D):',Dir)
tgt=[-wab/720*e1,-wab/720*e2]
diff=[sp.simplify(base[i]-tgt[i]) for i in range(2)]
print('base - (-wab/720 e):',diff,' parallel to D:',sp.simplify(diff[0]*Dir[1]-diff[1]*Dir[0])==0)
# jet lemma: solve delta_i from jet system and compare with G/4
import random
Rv=lambda: sp.Matrix([sp.Rational(random.randint(-9,9),7) for _ in range(3)])
def bg(V):
    E=sp.Matrix.hstack(V[1]-V[0],V[2]-V[0],V[3]-V[0]); Ei=E.inv().T
    g1,g2,g3=Ei[:,0],Ei[:,1],Ei[:,2]; return [-(g1+g2+g3),g1,g2,g3]
random.seed(3)
for t in range(3):
    z=Rv();a=Rv();b=Rv();q=Rv(); d=Rv(); B=(q+z+a+b)/4
    G1=bg([B,z,a,b]);G2=bg([B,q,z,b]);G3=bg([B,q,z,a])
    di=sp.Matrix(sp.symbols('u0:3'))
    s=sp.solve([di.dot(G1[0]), d.dot(G2[1])+di.dot(G2[0]), d.dot(G3[1])+di.dot(G3[0])],list(di))
    di=di.subs(s)
    N=(a-z).cross(b-z); n=N/sp.sqrt(N.dot(N)); P=sp.eye(3)-n*n.T
    r=(a-z)+(b-z)+P*(q-z); Gd=P*d-(n.dot(d)/(q-z).dot(n))*r
    print('jet lemma trial',t,':',sp.simplify(di-Gd/4)==sp.zeros(3,1))
