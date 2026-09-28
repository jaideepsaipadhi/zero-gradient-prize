# Checks Lemma td:faces(d): on a sphere-inscribed triangle, the face mean of n(x*)-n_F
# (normals oriented away from the centre) equals (x_c - o_F) + O(h^2): zero only for equilateral faces.
# Run: python3 notes/v3/face_mean_check.py
import numpy as np
def tri_on_sphere(angles_deg, h):
    # triangle in plane z=d with circumradius r=h, vertices at given polar angles
    r=h; d=np.sqrt(1-r*r)
    return np.array([[r*np.cos(np.radians(a)), r*np.sin(np.radians(a)), d] for a in angles_deg])
def face_mean(V, n=60):
    a,b,c=V; nF=np.cross(b-a,c-a); nF/=np.linalg.norm(nF)
    if nF@a<0: nF=-nF
    acc=np.zeros(3); area=0
    for i in range(n):
        for j in range(n-i):
            # centroid-rule on sub-triangles (two orientations)
            for (s,t,w) in [((i+1/3)/n,(j+1/3)/n,1)]+([((i+2/3)/n,(j+2/3)/n,1)] if i+j<n-1 else []):
                x=a+s*(b-a)+t*(c-a); nx=x/np.linalg.norm(x)
                acc+=w*(nx-nF); area+=w
    mean=acc/area
    xc=V.mean(0); o=np.array([0,0,V[0,2]])
    return mean, xc-o
for h in [0.2,0.1,0.05]:
    for ang in [(90,210,330),(90,200,340),(60,180,330)]:
        m,dx=face_mean(tri_on_sphere(ang,h))
        print(h,ang,"mean/h=",np.round(m/h,4)," (xc-o)/h=",np.round(dx/h,4))
