"""Reference triangle: z=(0,0),a=(1,0),b=(0,1).  V0 = U cap {nu_a=nu_b=0} cap mean0 (dim 3):
moment K(e) of r in V0 with r(z)=e, modulo the single-macro line D (moment of S1 cap mean0)."""
import sympy as sp
s,t=sp.symbols('s t'); cs=sp.symbols('c0:12'); mons=[1,s,t,s*s,s*t,t*t]
g=[sum(cs[i]*mons[i] for i in range(6)), sum(cs[6+i]*mons[i] for i in range(6))]
def I(f):
    p=sp.Poly(sp.expand(f),s,t); return sum(c*sp.factorial(i)*sp.factorial(j)/sp.factorial(i+j+2) for (i,j),c in p.terms())
at=lambda f,S,T: f.subs({s:S,t:T})
h=sp.Rational(1,2)
condU=[at(g[0],1,0),at(g[1],1,0),at(g[0],0,1),at(g[1],0,1),at(g[0]+g[1],h,h)]
nu=[at(g[1],h,0), at(g[0],0,h)]   # normal comps at M_za (edge along x -> normal y) and M_zb
S1=[at(g[0],0,0),at(g[1],0,0)]+nu  # g(z)=0 and nu=0 -> S1
mean=[I(g[0]),I(g[1])]
w1,w2,w3=sp.symbols('w_za w_zb w_ab',positive=True)
mz,ma,mb=1-s-t,s,t
Qw=sp.Rational(1,2)*(w1*mz*ma+w2*mz*mb+w3*ma*mb)
mom=[I(Qw*g[0]),I(Qw*g[1])]
def lin(eqs): return sp.Matrix([[sp.diff(e,c) for c in cs] for e in eqs])
# S1 cap mean0 line
K1=lin(condU+S1+mean).nullspace(); assert len(K1)==1
D=sp.simplify(lin(mom)*K1[0]); print('single-macro direction D =',sp.factor(D.T))
e1,e2=sp.symbols('e1 e2')
sol={}
A=lin(condU+nu+mean+[at(g[0],0,0),at(g[1],0,0)]); rhs=sp.Matrix([0]*(5+2+2)+[e1,e2])
x=A.gauss_jordan_solve(rhs)[0]
x=x.subs({p:0 for p in x.free_symbols if str(p).startswith('tau')})
Kv=sp.simplify(lin(mom)*x)
print('K(e) (a particular representative) =',sp.factor(Kv.T))
KM=Kv.jacobian([e1,e2]); print('K matrix:',sp.factor(KM))
