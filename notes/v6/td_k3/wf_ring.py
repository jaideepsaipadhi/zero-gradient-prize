"""(Optional) WF-split rings around [z,q], k=3 and k=4: total moment rank (sphere weights l^2).
Interior macro faces: face point = intersection of the segment joining the two interior points with the face
(the WF collinearity condition); boundary macro faces: face barycentre."""
from wf import GPatch, wf_subs
from patchlib import F3, sub, add, scal, dot, cross
from fractions import Fraction as Fr
import sys
def seg_plane(B0,B1,P):
    N=cross(sub(P[1],P[0]),sub(P[2],P[0])); t=dot(N,sub(P[0],B0))/dot(N,sub(B1,B0))
    x=add(B0,scal(t,sub(B1,B0)))
    # inside check via barycentrics
    nn=dot(N,N); l1=dot(cross(sub(P[1],x),sub(P[2],x)),N)/nn; l2=dot(cross(sub(P[2],x),sub(P[0],x)),N)/nn
    assert 0<t<1 and l1>0 and l2>0 and 1-l1-l2>0, 'WF face point outside'
    return x
def wf_ring(link,q,gam,k,z=(0,0,0)):
    link=[F3(x) for x in link]; q=F3(q); z=F3(z); m=len(link)
    T=[[q,z,link[j],link[(j+1)%m]] for j in range(m)]
    B=[scal(Fr(1,4),add(add(t[0],t[1]),add(t[2],t[3]))) for t in T]
    subs=[]; zero=[]; gamlist=[]
    for j,t in enumerate(T):
        fp=[]
        for i in range(4):
            fv=[t[l] for l in range(4) if l!=i]
            if i==3:   # face [q,z,x_j] = f_j shared with T_{j-1}
                fp.append(seg_plane(B[j-1],B[j],fv))
            elif i==2: # face [q,z,x_{j+1}] = f_{j+1} shared with T_{j+1}
                fp.append(seg_plane(B[j],B[(j+1)%m],fv))
            else:
                fp.append(scal(Fr(1,3),add(add(fv[0],fv[1]),fv[2])))
        for (s,i) in wf_subs(t,B[j],fp):
            subs.append(s)
            if i in (0,1): zero.append(s[1:])   # faces opposite q (i=0) and z (i=1) are ring boundary
            if i==0 and j in gam: gamlist.append((len(subs)-1,[t[1],t[2],t[3]]))
    return GPatch(subs,zero,gamlist,k)
if __name__=='__main__':
    cases=[('m=3 single Gamma',[(1,0,0),(0,1,0),(-1,-1,Fr(1,2))],(Fr(1,5),Fr(-1,7),1),[0]),
           ('m=4 single Gamma',[(1,0,0),(0,1,0),(-1,Fr(1,3),Fr(1,2)),(Fr(1,4),-1,Fr(1,2))],(Fr(1,9),Fr(1,11),1),[0]),
           ('flat square fan',[(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)],(Fr(1,9),Fr(1,11),1),[0,1,2,3]),
           ('flat pentagon fan',[(1,0,0),(Fr(3,10),1,0),(-Fr(4,5),Fr(3,5),0),(-Fr(4,5),-Fr(3,5),0),(Fr(3,10),-1,0)],(Fr(1,7),Fr(1,9),1),[0,1,2,3,4])]
    for k in (3,4):
        for nm,link,q,gam in cases:
            P=wf_ring(link,q,gam,k); print(f'k={k} {nm:20s}',P.ranks([(1,0,0),(0,1,0)]),flush=True)
