"""SYMBOLIC check of the flat reduction (Gamma-face in the plane x3=0, z=0, q=(q1,q2,H)):
 (i)  P.u_a = -sigma_a/(2|F|) + nu_a/H,  P.u_b = sigma_b/(2|F|) + nu_b/H  with sigma_x = J s_x.(pi_x - nu_x q'/H)
      (u_a = grad lam_q on [B,q,z,a], u_b = grad lam_q on [B,q,z,b]; computed from scratch here)
 (ii) flux: A_x.P = (H/2) sigma_x,  A_x = (q - z) x (x - z)/2
 (iii) sigma_x over the flux-admissible family = gamma - s_x ^ eta (affine in s_x)."""
import sympy as sp
a1,a2,b1,b2,q1,q2,H=sp.symbols('a1 a2 b1 b2 q1 q2 H',real=True)
p1,p2,nu=sp.symbols('p1 p2 nu')
z=sp.Matrix([0,0,0]); a=sp.Matrix([a1,a2,0]); b=sp.Matrix([b1,b2,0]); q=sp.Matrix([q1,q2,H])
B=(q+a+b)/4
def gradlam(verts,k):
    # gradient of the barycentric of verts[k] in simplex verts (4 points)
    M=sp.Matrix([[1]+list(v) for v in verts]); Mi=M.inv()
    return sp.Matrix(Mi[1:,k])
ua=gradlam([B,q,z,a],1); ub=gradlam([B,q,z,b],1)
area2=a1*b2-a2*b1   # = 2|F| (ccw)
J=lambda v: sp.Matrix([-v[1],v[0],0])
qp=sp.Matrix([q1,q2,0])
P=sp.Matrix([p1,p2,nu])
sig_a=J(a).dot(P-nu*qp/H); sig_b=J(b).dot(P-nu*qp/H)
print('(i) a:',sp.simplify(P.dot(ua)-(-sig_a/area2+nu/H))==0,' b:',sp.simplify(P.dot(ub)-(sig_b/area2+nu/H))==0)
A=(q-z).cross(a-z)/2
print('(ii) flux:',sp.simplify(A.dot(P)-H/2*sig_a)==0)
g,e1,e2,psi1,psi2,psin,Phi=sp.symbols('gamma e1 e2 psi1 psi2 psin Phi')
psi=sp.Matrix([psi1,psi2,psin])
sig=(2/H)*(60*Phi-A.dot(psi))      # from A.P = 60 Phi - A.psi, psi = delta + 2 p_q
eta=sp.Matrix([psi1-psin*q1/H,psi2-psin*q2/H])
print('(iii) sigma = 120Phi/H - (a ^ eta):',sp.simplify(sig-(120*Phi/H-(a1*eta[1]-a2*eta[0])))==0)
