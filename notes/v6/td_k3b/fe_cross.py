"""Full-FE cross-check (patchlib, mod-p ranks = lower bounds for Q-ranks) of the degeneracy theorem:
sagitta weights w = h(e,e) with an indefinite / semidefinite h for which sum_i c_i = 0."""
import sys
sys.path.insert(0,'../td_k3')
from patchlib import *
from patchlib import _mom_face_rows
from fractions import Fraction as Fr
def hq(H,e): return H[0]*e[0]*e[0]+2*H[1]*e[0]*e[1]+H[2]*e[1]*e[1]
def run(nm,link,q,z,gam,H):
    m=len(link); P=Patch([[q,z,link[j],link[(j+1)%m]] for j in range(m)],[(j,(1,2,3)) for j in gam],3)
    tot=[dict(),dict()]
    for fi,i in enumerate(gam):
        a,b=link[i],link[(i+1)%m]
        w=(hq(H,sub(a,z)),hq(H,sub(b,z)),hq(H,sub(b,a)))
        rows=_mom_face_rows(P,fi,w)
        for c in range(2):
            for k,v in rows[c].items(): tot[c][k]=tot[c].get(k,0)+v
    rdc,r=rank_on_Y1_mod(P,tot)
    sc=sum(hq(H,sub(link[(i+1)%m],link[i])) for i in gam)
    print(f'{nm:46s} h={H}  sum_i w_ab={sc}  FE tangential total-moment rank (mod p) = {r}',flush=True)
sq=[(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)]; q=(Fr(1,5),Fr(1,7),1); z=(0,0,0)
run('square fan, all 4 faces (sum c_i = 0)',sq,q,z,[0,1,2,3],(0,1,0))
run('square fan, all 4 faces, h=(x+y)^2',sq,q,z,[0,1,2,3],(1,1,1))
run('square fan, faces 0,2 (both degenerate at z)',sq,q,z,[0,2],(1,1,1))
run('square fan, faces 0,1 (face 1 non-degenerate)',sq,q,z,[0,1],(1,1,1))
run('square fan, face 0 only, h=(x+y)^2 (degenerate)',sq,q,z,[0],(1,1,1))
run('square fan, face 0 only, h=x^2-y^2',sq,q,z,[0],(1,0,-1))
# Gamma-flat star with off-plane non-Gamma link vertex
lk=[(1,0,0),(0,1,0),(-1,Fr(1,3),Fr(2,5)),(-1,-1,Fr(-1,3)),(Fr(1,2),-1,0)]
run('Gamma-flat, off-plane non-Gamma vertices, I={4,0}',lk,q,z,[4,0],(0,1,0))
