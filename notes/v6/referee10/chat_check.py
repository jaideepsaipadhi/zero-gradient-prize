# Independent check: representer of r -> int_{e^} r on P_{k-1}(T^), T^={(0,0),(1,0),(0,1)}, e^=[0,1]x{0}
import sympy as sp, sys
x,y,t=sp.symbols('x y t')
def check(k):
    basis=[x**a*y**b for a in range(k) for b in range(k-a)]
    n=len(basis)
    def ip(f,g): return sp.integrate(sp.integrate(f*g,(x,0,1-y)),(y,0,1))
    G=sp.Matrix(n,n,lambda i,j: ip(basis[i],basis[j]))
    rhs=sp.Matrix([sp.integrate(b.subs(y,0),(x,0,1)) for b in basis])
    c=G.LUsolve(rhs)
    r=sp.expand(sum(c[i]*basis[i] for i in range(n)))
    g=sp.expand(2*sp.diff(sp.legendre(k,t),t).subs(t,1-2*y))
    ck=sp.integrate(r.subs(y,0),(x,0,1))
    dk=sp.integrate(r.subs(y,0)**2,(x,0,1))
    return sp.simplify(r-g)==0, ck, dk, ip(r,r)
for k in range(int(sys.argv[1]),int(sys.argv[2])+1):
    ok,ck,dk,nn=check(k)
    print(k, ok, ck, ck==k*(k+1), dk==(k*(k+1))**2, nn==ck, flush=True)
