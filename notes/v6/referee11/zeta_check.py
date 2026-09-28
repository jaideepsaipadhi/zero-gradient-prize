# Riesz representer zeta in P3(T^) of r -> int_e r (e = {y=0}); check its trace on e is constant; also for P_{k-1}, k=4..7
import sympy as sp
x,y=sp.symbols('x y')
def iT(p): return sp.integrate(sp.integrate(sp.expand(p),(y,0,1-x)),(x,0,1))
for k in range(4,8):
    mons=[x**p*y**q for p in range(k) for q in range(k-p)]
    cs=sp.symbols('c0:%d'%len(mons)); z=sum(c*m for c,m in zip(cs,mons))
    sol=sp.solve([iT(z*m)-sp.integrate(m.subs(y,0),(x,0,1)) for m in mons],cs,dict=True)[0]
    tr=sp.expand(z.subs(sol).subs(y,0)); print(k, 'trace of zeta on e:', tr, ' expected k(k+1)/H*? H=1 ->', k*(k+1))
