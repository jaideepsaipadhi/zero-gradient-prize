"""Which 'moves' give rank 2 in flat multi-face rings?  (reduced model, exact)
  S1      : only single-macro fields (d=0, p_j=0 ... i.e. nu=0 and g(z)=0)
  S1+d    : plus the vertex move (p_j = -d so nu unchanged)
  S1+nu   : d=0, p_j free
  all"""
from reduced import Reduced
from fractions import Fraction as Fr
import sympy as sp
def variant(R,kind):
    m=R.m
    extra=[]
    if kind in ('S1','S1+nu'):
        extra+= [{('d',c):1} for c in range(3)]
    if kind in ('S1','S1+d'):
        for j in range(m):
            for c in range(3): extra.append({('p',j,c):1,('d',c):1})
    R.eqs=R.eqs+extra
    return R.analyse([(1,0,0),(0,1,0)])[1]
u=(1,-1,0); w=(0,1,-1); n=(1,1,1)
def Lt(a,b,h=0): return tuple(a*u[c]+b*w[c]+h*n[c] for c in range(3))
hexl=[Lt(1,0),Lt(1,1),Lt(0,1),Lt(-1,0),Lt(-1,-1),Lt(0,-1)]
cases=[('flat square fan',[(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)],(0,0,1),(0,0,0),range(4),[(1,0,0),(0,1,0)]),
 ('hex equilateral fan q axis',hexl,Lt(0,0,1),(0,0,0),range(6),[(1,-1,0),(1,1,-2)]),
 ('hex equilateral fan q off',hexl,Lt(Fr(1,3),Fr(1,5),1),(0,0,0),range(6),[(1,-1,0),(1,1,-2)]),
 ('hex eq, 2 adjacent Gamma',[Lt(1,0),Lt(1,1),Lt(0,1),Lt(-1,0,Fr(1,2)),Lt(-1,-1,Fr(1,2)),Lt(0,-1,Fr(1,3))],Lt(0,0,1),(0,0,0),[0,1],[(1,-1,0),(1,1,-2)]),
 ('hex eq, 3 adjacent Gamma',[Lt(1,0),Lt(1,1),Lt(0,1),Lt(-1,0),Lt(-1,-1,Fr(1,2)),Lt(0,-1,Fr(1,3))],Lt(0,0,1),(0,0,0),[0,1,2],[(1,-1,0),(1,1,-2)]),
]
for nm,link,q,z,gam,fr in cases:
    out=[]
    for kind in ('S1','S1+d','S1+nu','all'):
        R=Reduced(link,q,z,gam); m=R.m
        extra=[]
        if kind in ('S1','S1+nu'): extra+=[{('d',c):1} for c in range(3)]
        if kind in ('S1','S1+d'):
            for j in range(m):
                for c in range(3): extra.append({('p',j,c):1,('d',c):1})
        R.eqs=R.eqs+extra
        out.append(f'{kind}:{R.analyse(fr)[1]}')
    print(f'{nm:30s}',' '.join(out))
