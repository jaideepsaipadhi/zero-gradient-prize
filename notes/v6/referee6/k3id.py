import numpy as np
rng=np.random.default_rng(1)
for _ in range(3):
    z=rng.normal(size=(3,2))
    e=[z[0]-z[1],z[1]-z[2],z[2]-z[0]]  # cyclic
    S=sum(np.dot(v,v)*v for v in e)
    xc=z.mean(0)
    # circumcentre
    A=2*np.array([z[1]-z[0],z[2]-z[0]]); b=np.array([z[1]@z[1]-z[0]@z[0],z[2]@z[2]-z[0]@z[0]])
    o=np.linalg.solve(A,b)
    area=abs(np.cross(z[1]-z[0],z[2]-z[0]))/2
    d=xc-o; Rd=np.array([-d[1],d[0]])
    print(S, 12*area*Rd, S/(12*area*Rd))
