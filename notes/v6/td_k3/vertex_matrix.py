"""The vertex-move matrix of a ring (closed form, exact):
   V_R := sum_i c_i P_T G_i |_T ,  c_i = |F_i| w_ab^(i) / H_i,
   G_i d = P_i d - (n_i.d / (q-z).n_i) r_i,  r_i = (a_i-z)+(b_i-z)+P_i(q-z),
T = tangent plane of Gamma at z (for the sphere: z-perp).  The witness of the report has tangential total
moment  P_T T(d) = -(1/360) V_R d.  Also the tilt/shape quantities tau, rho of the explicit criterion
tau (tau + rho) <= 1/2."""
import numpy as np
def vertex_matrix(link,q,z,gam,nz,weights=None):
    link=[np.array(x,float) for x in link]; q=np.array(q,float); z=np.array(z,float); nz=np.array(nz,float); nz/=np.linalg.norm(nz)
    m=len(link)
    # tangent basis
    e1=np.cross(nz,[1,0,0] if abs(nz[0])<0.9 else [0,1,0]); e1/=np.linalg.norm(e1); e2=np.cross(nz,e1)
    Tb=np.array([e1,e2])
    V=np.zeros((2,2)); taus=[]; rhos=[]; csum=0
    for i in gam:
        a=link[i]; b=link[(i+1)%m]
        N=np.cross(a-z,b-z); area=np.linalg.norm(N)/2; n=N/np.linalg.norm(N)
        if n@nz<0: n=-n
        Hs=(q-z)@n; H=abs(Hs)
        w=weights[i][2] if weights else (a-b)@(a-b)
        P=np.eye(3)-np.outer(n,n)
        r=(a-z)+(b-z)+P@(q-z)
        G=P-np.outer(r,n)/Hs
        c=area*w/H; csum+=c
        V+=c*Tb@G@Tb.T
        taus.append(np.arccos(min(1,n@nz))); rhos.append(np.linalg.norm(r)/H)
    return V,csum,max(taus),max(rhos)
if __name__=='__main__':
    ring4=[(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)]
    V,cs,t,r=vertex_matrix(ring4,(0,0,2),(0,0,1),[0,1,2,3],(0,0,1))
    print('octa fan: V =',V.round(12).tolist(),' tau=%.3f rho=%.3f tau(tau+rho)=%.3f'%(t,r,t*(t+r)))
    V,cs,t,r=vertex_matrix(ring4,(0.2,1/7,2),(0,0,1),[0,1,2,3],(0,0,1))
    print('octa off-axis: V sv =',np.linalg.svd(V,compute_uv=False).round(5),' (sum c = %.3f)'%cs)
