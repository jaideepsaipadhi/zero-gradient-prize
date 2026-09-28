"""Flat case, one face F=[z,a,b] (z=0, a,b in the plane x3=0) in macro [q,z,a,b], q=(q1,q2,h):
jet map d -> ghat(z)=d_1 (tangential), symbolic."""
import sympy as sp
a1,a2,b1,b2,q1,q2,h=sp.symbols('a1 a2 b1 b2 q1 q2 h',real=True)
z=sp.Matrix([0,0,0]); a=sp.Matrix([a1,a2,0]); b=sp.Matrix([b1,b2,0]); q=sp.Matrix([q1,q2,h])
B=(q+z+a+b)/4
def grads(V):
    M=sp.Matrix.hstack(*[sp.Matrix([*v,1]) for v in V]).T  # rows (v,1)
    Minv=M.inv()   # columns? lam = coeffs: lam_i(x)= row i of (M^T)^{-1} [x;1]
    L=(M.T).inv()
    return [sp.Matrix(L[i,0:3]).T for i in range(4)]
s2=[B,q,z,b]; s3=[B,q,z,a]
G2=grads(s2); G3=grads(s3)
d=sp.Matrix(sp.symbols('d1:4')); x,y=sp.symbols('x y')
d1=sp.Matrix([x,y,0])
sol=sp.solve([ (G2[0].T.dot(d1)+G2[1].T.dot(d)), (G3[0].T.dot(d1)+G3[1].T.dot(d))],[x,y],dict=True)[0]
gz=sp.Matrix([sp.simplify(sol[x]),sp.simplify(sol[y])])
print('ghat(z) =');sp.pprint(sp.simplify(gz))
J=gz.jacobian(d); print('L matrix (2x3):'); sp.pprint(sp.simplify(J))
