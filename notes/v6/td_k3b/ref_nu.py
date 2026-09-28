"""Reference triangle z=(0,0),a=(1,0),b=(0,1). U = {g in P2(F;R^2): g(a)=g(b)=0, g(m_ab).(1,1)=0}.
Moment of mean-zero representatives with prescribed pi-data (e, alpha, beta), modulo the S-line D."""
import sympy as sp
s,t=sp.symbols('s t'); cs=sp.symbols('c0:12'); mons=[1,s,t,s*s,s*t,t*t]
g=[sum(cs[i]*mons[i] for i in range(6)), sum(cs[6+i]*mons[i] for i in range(6))]
def I(f):
    p=sp.Poly(sp.expand(f),s,t); return sum(c*sp.factorial(i)*sp.factorial(j)/sp.factorial(i+j+2) for (i,j),c in p.terms())
at=lambda f,S,T: f.subs({s:S,t:T})
h=sp.Rational(1,2)
condU=[at(g[0],1,0),at(g[1],1,0),at(g[0],0,1),at(g[1],0,1),at(g[0]+g[1],h,h)]
nu=[at(g[1],h,0), at(g[0],0,h)]
gz=[at(g[0],0,0),at(g[1],0,0)]
mean=[I(g[0]),I(g[1])]
w1,w2,w3=sp.symbols('w_za w_zb w_ab')
mz,ma,mb=1-s-t,s,t
Qw=sp.Rational(1,2)*(w1*mz*ma+w2*mz*mb+w3*ma*mb)
mom=[I(Qw*g[0]),I(Qw*g[1])]
def lin(eqs): return sp.Matrix([[sp.diff(e,c) for c in cs] for e in eqs])
e1,e2,al,be=sp.symbols('e1 e2 alpha beta')
A=lin(condU+nu+mean+gz); rhs=sp.Matrix([0]*5+[al,be]+[0,0]+[e1,e2])
x,params=A.gauss_jordan_solve(rhs)
x=x.subs({p:0 for p in params})
M=sp.expand(lin(mom)*x)
D=sp.Matrix([w1-w3, w3-w2])/720
print('D =',D.T)
for v in [al,be,e1,e2]:
    print(v, sp.factor(M.diff(v).T))
# reduce modulo D: write each column c as c - lambda D with lambda chosen to kill 2nd comp (when D2!=0) -- just print det with D
for v in [al,be]:
    c=M.diff(v); print('det[',v,',D] =',sp.factor(sp.Matrix.hstack(c,D).det()))
