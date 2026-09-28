"""EXACT test of the hypothesis (H): a flat edge star is degenerate for h (tangential total moment rank<2)
only if EVERY Gamma-face of the star is individually degenerate at z (h(e_ab)=0 and h(e_za)=h(e_zb))."""
import random, math, time
from fractions import Fraction as Fr
from degen import moment_mats, degenerate_h, hq
random.seed(int(__import__('sys').argv[1]) if len(__import__('sys').argv)>1 else 3)
def facedeg(link,i,h):
    m=len(link); a,b=link[i],link[(i+1)%m]
    H=((h[0],h[1]),(h[1],h[2]))
    return hq(H,(b[0]-a[0],b[1]-a[1]))==0 and hq(H,a[:2])==hq(H,b[:2])
t0=time.time(); n=0
while time.time()-t0<240:
    m=random.randint(3,6)
    angs=sorted(random.uniform(0,2*math.pi) for _ in range(m))
    if max((angs[(i+1)%m]-angs[i])%(2*math.pi) for i in range(m))>0.85*math.pi: continue
    link=[]
    for th in angs:
        r=random.uniform(0.6,1.4)
        link.append((Fr(round(3*r*math.cos(th)),3),Fr(round(3*r*math.sin(th)),3),0))
    if len(set(link))<m or any(l[0]==0 and l[1]==0 for l in link): continue
    from degen import hq as _h
    bad=False
    for j in range(m):
        a,b=link[j],link[(j+1)%m]
        if a[0]*b[1]-a[1]*b[0]<=0: bad=True
    if bad: continue
    k=random.randint(1,m); gam=sorted(random.sample(range(m),k))
    q=(Fr(random.randint(-3,3),10),Fr(random.randint(-3,3),10),Fr(random.randint(5,15),10))
    Ms=moment_mats(link,q,(0,0,0),gam)
    dg=degenerate_h(Ms); n+=1
    ok=True
    if isinstance(dg,str): ok=False; print('!!',dg)
    else:
        for lam,hs in dg:
            for h in hs:
                allf=all(facedeg(link,i,h) for i in gam)
                if not allf: ok=False
                print(f'  m={m} I={gam} degenerate h={h} lam={lam}: all faces individually degenerate: {allf}')
    print(f'case {n}: m={m} |I|={k}  degeneracies={0 if isinstance(dg,str) else sum(len(x[1]) for x in dg)}  H consistent: {ok}',flush=True)
