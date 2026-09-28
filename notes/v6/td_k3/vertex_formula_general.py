"""EXACT check (reduced model, rational) of the general (non-flat) vertex-move formula:
  witness with jet d, nu=0 (p_j=-d, dp=-d/2, Phi=0), mean zero, S1 parts chosen canonically:
  sum of face moments  T(d) = - sum_i fac_i w_ab^(i)/360 * ghat_i(z),
  ghat_i(z) = 1/4 [ P_i d - (n_i.d / (q-z).n_i) ((a_i-z)+(b_i-z)+P_i(q-z)) ],   fac_i = 4|F_i|/H_i,
checked modulo span{D_i} with generic weights, and exactly with w_za=w_zb=w_ab (D_i=0)."""
from reduced import Reduced
from fractions import Fraction as Fr
import sympy as sp, random, math
from flat_formula2 import solve_total
random.seed(11)
def R3(v): return sp.Matrix([sp.Rational(c) for c in v])
def predict(link,q,z,gam,W,d):
    m=len(link); q=R3(q); z=R3(z); d=R3(d); tot=sp.zeros(3,1)
    for i in gam:
        a=R3(link[i]); b=R3(link[(i+1)%m])
        N=(a-z).cross(b-z); area=sp.sqrt(N.dot(N))/2; n=N/sp.sqrt(N.dot(N))
        H=abs((q-z).dot(n)); Pd=d-d.dot(n)*n; Pq=(q-z)-(q-z).dot(n)*n
        gz=(Pd-(n.dot(d)/(q-z).dot(n))*((a-z)+(b-z)+Pq))/4
        tot+= -(4*area/H)*W[i][2]/360*gz
    return sp.simplify(tot)
def full_total(link,q,z,gam,W,d):
    R=Reduced(link,q,z,gam,weights=W); m=R.m; names=R.names
    rows=[(r,0) for r in R.eqs]+[({('d',c):1},d[c]) for c in range(3)]
    for j in range(m): rows+=[({('p',j,c):1},-d[c]) for c in range(3)]
    rows+=[({('dp',c):1},-Fr(d[c])/2) for c in range(3)]+[({('Phi',):1},0)]
    for i in gam: rows+=[(r,0) for r in R.mean_rows(i)]
    A=sp.Matrix([[sp.Rational(e.get(n,0)) if n in e else 0 for n in names] for e,_ in rows])
    b=sp.Matrix([sp.Rational(v) for _,v in rows])
    sol,params=A.gauss_jordan_solve(b)
    tot=sp.zeros(3,1)
    for i in gam:
        rr=R.mom_rows_face(i)
        for c in range(3): tot[c]+=sum(sp.Rational(v)*sol[R.ix[k]] for k,v in rr[c].items())
    return tot,params
cases=[]
# octahedral fan (degenerate) and off-axis, curved random rings (all link points on a sphere through z)
ring4=[(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)]
cases.append(('octa fan',ring4,(0,0,2),(0,0,1),[0,1,2,3]))
cases.append(('octa off-axis',ring4,(Fr(1,5),Fr(1,7),2),(0,0,1),[0,1,2,3]))
for t in range(4):
    m=random.randint(3,6)
    angs=sorted(random.uniform(0,2*math.pi) for _ in range(m))
    if max((angs[(i+1)%m]-angs[i])%(2*math.pi) for i in range(m))>0.9*math.pi: continue
    gam=sorted(random.sample(range(m),random.randint(1,m)))
    link=[]
    for j,th in enumerate(angs):
        x=Fr(math.cos(th)).limit_denominator(40); y=Fr(math.sin(th)).limit_denominator(40)
        on=(j in gam) or ((j-1)%m in gam)
        link.append((x,y,-Fr(random.randint(1,6),10)*(x*x+y*y) if on else Fr(random.randint(2,8),10)))
    cases.append((f'curved random m={m}',link,(Fr(random.randint(-2,2),10),Fr(random.randint(-2,2),10),Fr(1)),(0,0,0),gam))
if __name__=="__main__":
 for nm,link,q,z,gam in cases:
     for eq in (True,False):
         W={i:(Fr(k:=random.randint(1,5)),)*3 for i in gam} if eq else {i:(Fr(random.randint(1,9)),Fr(random.randint(1,9)),Fr(random.randint(1,9))) for i in gam}
         for d in [(1,0,0),(0,1,0),(0,0,1)]:
             tot,params=full_total(link,q,z,gam,W,d)
             pr=predict(link,q,z,gam,W,d)
             base=tot.subs({p:0 for p in params}); Ds=sp.Matrix.hstack(*[tot.diff(p) for p in params]) if params else sp.zeros(3,0)
             diff=sp.simplify(base-pr)
             ok=(Ds.rank()==sp.Matrix.hstack(Ds,diff).rank())
             print(f'{nm:20s} {"equal w" if eq else "generic w"} d={d}: formula holds (mod D-span, D-rank {Ds.rank()}): {ok}',flush=True)
