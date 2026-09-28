"""Sharper check of the flat vertex formula: random weights with w_za=w_zb=w_ab per face (so D_i=0 and the
formula is an EXACT identity):  total moment = -(1/1440) sum_i fac_i w_ab^(i) d,  fac_i = 4|F_i|/H,
for d = e_x and d = e_y; also generic weights (checked modulo D-span).  Exact rational."""
from reduced import Reduced
from fractions import Fraction as Fr
import sympy as sp, random, math
random.seed(7)
def solve_total(R,link,gam,d):
    m=R.m; names=R.names
    rows=[(r,0) for r in R.eqs]+[({('d',c):1},d[c]) for c in range(3)]
    for j in range(m): rows+=[({('p',j,c):1},-d[c]) for c in range(3)]
    rows+=[({('dp',c):1},-d[c]/2) for c in range(3)]+[({('Phi',):1},0)]
    for i in gam: rows+=[(r,0) for r in R.mean_rows(i)]
    A=sp.Matrix([[sp.Rational(e.get(n,0)) if n in e else 0 for n in names] for e,_ in rows])
    b=sp.Matrix([sp.Rational(v) for _,v in rows])
    sol,params=A.gauss_jordan_solve(b)
    tot=sp.zeros(2,1)
    for i in gam:
        rr=R.mom_rows_face(i)
        for c in range(2): tot[c]+=sum(sp.Rational(v)*sol[R.ix[k]] for k,v in rr[c].items())
    return tot,params
if __name__=="__main__":
    nt=0; bad=0
    while nt<10:
        m=random.randint(3,7)
        angs=sorted(random.uniform(0,2*math.pi) for _ in range(m))
        if max((angs[(i+1)%m]-angs[i])%(2*math.pi) for i in range(m))>0.9*math.pi: continue
        nt+=1
        gam=sorted(random.sample(range(m),random.randint(1,m)))
        link=[]
        for j,t in enumerate(angs):
            r=random.uniform(0.7,1.3)
            x=Fr(r*math.cos(t)).limit_denominator(50); y=Fr(r*math.sin(t)).limit_denominator(50)
            on=(j in gam) or ((j-1)%m in gam)
            link.append((x,y,Fr(0) if on else Fr(random.randint(2,8),10)))
        H=Fr(random.randint(6,15),10); q=(Fr(random.randint(-3,3),10),Fr(random.randint(-3,3),10),H)
        W={i:(Fr(k:=random.randint(1,9)),Fr(k),Fr(k)) for i in gam}
        R=Reduced(link,q,(0,0,0),gam,weights=W)
        area=lambda i: abs(link[i][0]*link[(i+1)%m][1]-link[i][1]*link[(i+1)%m][0])/2
        kap=sum(4*area(i)/H*W[i][2] for i in gam)/1440
        for d in [(Fr(1),Fr(0),Fr(0)),(Fr(0),Fr(1),Fr(0))]:
            tot,params=solve_total(R,link,gam,d)
            pred=sp.Matrix([-kap*d[0],-kap*d[1]])
            ok=sp.simplify(tot-pred)==sp.zeros(2,1)
            bad+=not ok
            print(f'm={m} gam={gam} d={d[:2]} exact identity: {ok}  (free params: {len(params)})')
    print('bad',bad)
