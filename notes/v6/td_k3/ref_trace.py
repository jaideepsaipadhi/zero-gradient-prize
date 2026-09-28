"""Reference-triangle computation (exact, sympy).  F = conv(z,a,b), z=(0,0), a=(1,0), b=(0,1),
mu_z = 1-s-t, mu_a = s, mu_b = t.  U = {g in P2(F;R^2): g(a)=g(b)=0, g.(1,1)=0 at mid(ab)} (dim 7)
U0 = U cap {g(z)=0} (dim 5).  E(g) = (int mu_z mu_a g, int mu_z mu_b g, int mu_a mu_b g) in R^6.
Weighted moment for weights w=(w_za,w_zb,w_ab): sum w_ij E_ij(g) in R^2."""
import sympy as sp
s,t=sp.symbols('s t')
mons=[1,s,t,s*s,s*t,t*t]
cs=sp.symbols('c0:12')
g=[sum(cs[i]*mons[i] for i in range(6)), sum(cs[6+i]*mons[i] for i in range(6))]
def I(f):
    p=sp.Poly(sp.expand(f),s,t); r=0
    for (i,j),c in p.terms(): r+=c*sp.factorial(i)*sp.factorial(j)/sp.factorial(i+j+2)
    return r
mz,ma,mb=1-s-t,s,t
condU=[g[0].subs({s:1,t:0}),g[1].subs({s:1,t:0}),g[0].subs({s:0,t:1}),g[1].subs({s:0,t:1}),
       (g[0]+g[1]).subs({s:sp.Rational(1,2),t:sp.Rational(1,2)})]
condz=[g[0].subs({s:0,t:0}),g[1].subs({s:0,t:0})]
mean=[I(g[0]),I(g[1])]
E=[I(mz*ma*g[0]),I(mz*ma*g[1]),I(mz*mb*g[0]),I(mz*mb*g[1]),I(ma*mb*g[0]),I(ma*mb*g[1])]
def lin(eqs): return sp.Matrix([[sp.diff(e,c) for c in cs] for e in eqs])
def rank_E_on(conds):
    K=sp.Matrix.hstack(*lin(conds).nullspace())
    ME=lin(E)*K
    return K.shape[1], ME.rank(), ME
for nm,conds in [('U',condU),('U cap mean0',condU+mean),('U0',condU+condz),('U0 cap mean0',condU+condz+mean)]:
    d,r,ME=rank_E_on(conds)
    ann=ME.T.nullspace()
    print(f'{nm:14s} dim {d}  rank E {r}  annihilator of E-image:',[list(v.T) for v in ann])
# weighted rank on U cap mean0 and U0 cap mean0 for symbolic weights
w1,w2,w3=sp.symbols('w_za w_zb w_ab')
for nm,conds in [('U cap mean0',condU+mean),('U0 cap mean0',condU+condz+mean)]:
    d,r,ME=rank_E_on(conds)
    W=sp.Matrix([[w1,0,w2,0,w3,0],[0,w1,0,w2,0,w3]])*ME
    # 2x2 minors
    minors=set()
    for i in range(W.shape[1]):
        for j in range(i+1,W.shape[1]):
            m=sp.factor(W[:,[i,j]].det())
            if m!=0: minors.add(m)
    print(nm,'weighted 2x2 minors (factored):',minors)
