"""EXACT check (reduced model of td_k3, rational arithmetic) of the PER-FACE moment formula, any geometry, any weights:
for v in Y^1(R) with interior-face data (delta, p_q, p_j, Phi) [fluxes equal], on every Gamma-face F_i=(z,a,b):
  M_i == (8|F_i|/H_i) * [ -(w_ab/2880) G_i delta  -  ((w_za-w_zb)/2880) * ( (P_a.u_a) (b-z) - (P_b.u_b) (a-z) ) ]
  modulo R*D_i,  D_i = (w_za-w_ab)(a-z) + (w_ab-w_zb)(b-z),
with P_x = delta + p_x, u_a = grad lam_q on s3=[B,q,z,a], u_b = grad lam_q on s2=[B,q,z,b].
The check is per face (1-dim ambiguity), hence never vacuous.  Also checks the flat closed forms of u_a,u_b."""
import sys, random, math
sys.path.insert(0,'../td_k3')
from reduced import Reduced
from patchlib import bary_grads, cross, sub, add, scal, dot
from fractions import Fraction as Fr
import sympy as sp
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 1)
def V(x): return sp.Matrix([sp.Rational(c) for c in x])
def predict_face(R,i,W,d,P):
    m=R.m; z=R.z; q=R.q; a=R.link[i]; b=R.link[(i+1)%m]
    B=scal(Fr(1,4),add(add(q,z),add(a,b)))
    G2,_=bary_grads([B,q,z,b]); G3,_=bary_grads([B,q,z,a])
    ua=V(G3[1]); ub=V(G2[1])
    N=V(cross(sub(a,z),sub(b,z))); area=sp.sqrt(N.dot(N))/2; n=N/sp.sqrt(N.dot(N))
    qz=V(sub(q,z)); H=abs(qz.dot(n)); ea=V(sub(a,z)); eb=V(sub(b,z)); dd=V(d)
    Gd=(dd-dd.dot(n)*n)-(n.dot(dd)/qz.dot(n))*(ea+eb+(qz-qz.dot(n)*n))
    wza,wzb,wab=[sp.Rational(x) for x in W[i]]
    Pa=V(P[i]); Pb=V(P[(i+1)%m])
    M=(8*area/H)*(-(wab/2880)*Gd-((wza-wzb)/2880)*((Pa.dot(ua))*eb-(Pb.dot(ub))*ea))
    D=(wza-wab)*ea+(wab-wzb)*eb
    return M,D,ua,ub,ea,eb,qz,H,area
def face_moment(link,q,z,gam,W,d,dp,p,Phi):
    R=Reduced(link,q,z,gam,weights=W); m=R.m; names=R.names
    rows=[(r,0) for r in R.eqs]+[({('d',c):1},d[c]) for c in range(3)]+[({('dp',c):1},dp[c]) for c in range(3)]
    for j in range(m): rows+=[({('p',j,c):1},p[j][c]) for c in range(3)]
    rows+=[({('Phi',):1},Phi)]
    for i in gam: rows+=[(r,0) for r in R.mean_rows(i)]
    A=sp.Matrix([[sp.Rational(e.get(n,0)) if n in e else 0 for n in names] for e,_ in rows])
    b=sp.Matrix([sp.Rational(v) for _,v in rows])
    sol,params=A.gauss_jordan_solve(b)
    out={}
    for i in gam:
        rr=R.mom_rows_face(i); tot=sp.zeros(3,1)
        for c in range(3): tot[c]=sum(sp.Rational(v)*sol[R.ix[k]] for k,v in rr[c].items())
        out[i]=tot
    return R,out,params
def rnd(): return Fr(random.randint(-9,9),random.randint(1,5))
nok=0;nall=0
for trial in range(int(sys.argv[2]) if len(sys.argv)>2 else 6):
    m=random.randint(3,6)
    angs=sorted(random.uniform(0,2*math.pi) for _ in range(m))
    if max((angs[(i+1)%m]-angs[i])%(2*math.pi) for i in range(m))>0.85*math.pi: continue
    flat=(trial%2==0)
    gam=sorted(random.sample(range(m),random.randint(1,m)))
    link=[]
    for j,th in enumerate(angs):
        x=Fr(math.cos(th)).limit_denominator(9); y=Fr(math.sin(th)).limit_denominator(9)
        onG=(j in gam) or ((j-1)%m in gam)
        zc=0 if (flat and onG) else Fr(random.randint(-4,4),10)
        link.append((x,y,zc))
    q=(Fr(random.randint(-3,3),10),Fr(random.randint(-3,3),10),Fr(random.randint(6,14),10)); z=(0,0,0)
    W={i:(rnd(),rnd(),rnd()) for i in gam}
    d=(rnd(),rnd(),rnd()); dp=(rnd(),rnd(),rnd()); Phi=rnd()
    p=[]
    for j in range(m):
        A=scal(Fr(1,2),cross(sub(q,z),sub(link[j],z)))
        pj=(rnd(),rnd(),rnd())
        # enforce A.((d+dp)/30 + pj/60) = Phi by moving along A
        cur=dot(A,scal(Fr(1,30),add(d,dp)))+dot(A,pj)/60
        pj=add(pj,scal((Phi-cur)*60/dot(A,A),A)); p.append(pj)
    P=[add(d,p[j]) for j in range(m)]
    try:
        R,out,params=face_moment(link,q,z,gam,W,d,dp,p,Phi)
    except Exception as e:
        print('skip',e); continue
    for i in gam:
        M,D,ua,ub,ea,eb,qz,H,area=predict_face(R,i,W,d,P)
        diff=sp.simplify(out[i]-M)
        ok=sp.Matrix.hstack(D,diff).rank()<=1 and (D.norm()!=0 or diff.norm()==0)
        # flat closed forms of u_a,u_b (only meaningful when face in plane x3=0 and z=0)
        fl=''
        if ea[2]==0 and eb[2]==0:
            Jsa=sp.Matrix([-ea[1],ea[0],0]); Jsb=sp.Matrix([-eb[1],eb[0],0])
            gmb=Jsa/(2*area); gma=-Jsb/(2*area); qp=sp.Matrix([qz[0],qz[1],0])
            ua_p=-gmb+sp.Matrix([0,0,(1+gmb.dot(qp))/H]); ub_p=-gma+sp.Matrix([0,0,(1+gma.dot(qp))/H])
            fl=f' flat u-forms ok: {sp.simplify(ua-ua_p)==sp.zeros(3,1) and sp.simplify(ub-ub_p)==sp.zeros(3,1)}'
        nall+=1; nok+=ok
        print(f'trial {trial} {"Gamma-flat" if flat else "curved"} m={m} I={gam} face {i}: per-face formula holds mod R D_i: {ok}{fl}',flush=True)
print(f'TOTAL {nok}/{nall} faces agree')
