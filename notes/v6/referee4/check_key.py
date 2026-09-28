# Independent check of Lemma K(d): K = sigma (w-z) lam_w lam_z^2 on a random locked star
import sympy as sp, random
x,y=sp.symbols('x y')
def bary(P,Q,R):
    M=sp.Matrix([[P[0],Q[0],R[0]],[P[1],Q[1],R[1]],[1,1,1]]); Mi=M.inv()
    v=Mi*sp.Matrix([x,y,1]); return [sp.expand(v[i]) for i in range(3)]
random.seed(3)
for t in [sp.Rational(1,10),sp.Rational(1,100),sp.Rational(1,1000)]:
    # wall: chords of unit circle approx: z=(0,0), z-=(-t, -a t^2), z+=(t, -b t^2) (turning angle ~ t)
    a=sp.Rational(random.randint(3,7),10); b=sp.Rational(random.randint(3,7),10)
    z=(0,0); zm=(-t,-a*t**2); zp=(t,-b*t**2)
    w=(sp.Rational(random.randint(-3,3),10)*t, sp.Rational(random.randint(7,12),10)*t)
    T1=(zm,z,w); T2=(z,zp,w)
    l1=bary(*T1); l2=bary(*T2)   # order: zm,z,w / z,zp,w
    res=[]
    for (lz,lw) in [(l1[1],l1[2]),(l2[0],l2[2])]:
        K=[ (w[0]-z[0])*lw*lz**2, (w[1]-z[1])*lw*lz**2 ]
        d=sp.diff(K[0],x)+sp.diff(K[1],y)
        res.append((d.subs({x:0,y:0}), [ [sp.diff(K[i],v).subs({x:0,y:0}) for v in (x,y)] for i in range(2)]))
    # flux of K through s: K parallel to w-z on s => K.n_s=0
    # divergence integral over T1 via boundary: K=0 on other edges
    print('t=',t,' div at z T1,T2:',res[0][0],res[1][0])
