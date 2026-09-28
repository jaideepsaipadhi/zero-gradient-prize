"""Icosphere macro meshes (level L), one prism layer of thickness delta*h extruded radially, prisms split by
global vertex order.  For every boundary face F: ring R(z0, q_F) (z0 lowest vertex of F, q_F = lift of z0).
Report: max tilt tau, max rho, max tau(tau+rho) (explicit criterion <= 1/2), and min over rings of
sigma_min(V_R)/sum_i c_i (flat value 1)."""
import numpy as np, itertools
from vertex_matrix import vertex_matrix
def icosphere(L):
    t=(1+5**.5)/2
    V=[(-1,t,0),(1,t,0),(-1,-t,0),(1,-t,0),(0,-1,t),(0,1,t),(0,-1,-t),(0,1,-t),(t,0,-1),(t,0,1),(-t,0,-1),(-t,0,1)]
    V=[np.array(v)/np.linalg.norm(v) for v in V]
    F=[(0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),(1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),(3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),(4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1)]
    for _ in range(L):
        mid={}; NF=[]
        def mp(i,j):
            k=(min(i,j),max(i,j))
            if k not in mid:
                v=V[i]+V[j]; V.append(v/np.linalg.norm(v)); mid[k]=len(V)-1
            return mid[k]
        for (a,b,c) in F:
            ab,bc,ca=mp(a,b),mp(b,c),mp(c,a); NF+=[(a,ab,ca),(b,bc,ab),(c,ca,bc),(ab,bc,ca)]
        F=NF
    return np.array(V),F
for L in range(0,5):
    V,F=icosphere(L)
    h=max(np.linalg.norm(V[a]-V[b]) for f in F for a,b in itertools.combinations(f,2))
    for delta in (0.5,1.0):
        lift=lambda i: V[i]*(1+delta*h)
        # faces at each vertex, with apex = lift of the lowest vertex
        byz={}
        for f in F:
            z0=min(f); byz.setdefault(z0,[]).append(f)
        worst=0; minrat=9; tmax=0; rmax=0
        for z0,fs in byz.items():
            # ring around [z0, lift z0]: link = cyclic neighbours of z0; Gamma faces = faces with lowest vertex z0
            # vertex_matrix only needs the Gamma faces (a,b) pairs in the ring:
            link=[]; gam=[]
            for k,f in enumerate(fs):
                others=[v for v in f if v!=z0]
                # orient (a,b) so that normal points outward (n.z>0)
                a,b=V[others[0]],V[others[1]]
                if np.cross(a-V[z0],b-V[z0])@V[z0]<0: a,b=b,a
                link+= [a,b]; gam.append(2*k)
            Vm,cs,t,r=vertex_matrix(link,lift(z0),V[z0],gam,V[z0])
            s=np.linalg.svd(Vm,compute_uv=False)[-1]/cs
            minrat=min(minrat,s); tmax=max(tmax,t); rmax=max(rmax,r); worst=max(worst,t*(t+r))
        print(f'L={L} h={h:.3f} delta={delta}: max tau={tmax:.3f}  max rho={rmax:.2f}  max tau(tau+rho)={worst:.3f}  min sigma_min(V)/sum c={minrat:.3f}')
