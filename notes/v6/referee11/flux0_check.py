# independent check of robust6 Lemma flux0 (k=4,5): for g in P_k(e^), int g = 0, the system
# p in lam_z lam_a P_{k-2}^2, (p, grad r)_T^ = -<g, r>_e^ for r in P_{k-1} is solvable and int_e p.n = 0 for every solution.
import sympy as sp
x,y=sp.symbols('x y')
def iT(p): return sp.integrate(sp.integrate(sp.expand(p),(y,0,1-x)),(x,0,1))
for k in (4,5):
    lz,la=1-x-y,x
    mons=[x**i*y**j for i in range(k-1) for j in range(k-1-i)]
    amb=[(lz*la*m,sp.Integer(0)) for m in mons]+[(sp.Integer(0),lz*la*m) for m in mons]
    rs=[x**i*y**j for i in range(k) for j in range(k-i)]
    A=sp.Matrix([[iT(u0*sp.diff(r,x)+u1*sp.diff(r,y)) for (u0,u1) in amb] for r in rs])
    flux_row=sp.Matrix([[sp.integrate((-u1).subs(y,0),(x,0,1)) for (u0,u1) in amb]])   # n=(0,-1)
    # flux functional in row space of A  <=> flux constant on solution sets; express flux = lam^T A
    lam=sp.symbols('l0:%d'%len(rs)); sol=sp.solve(list(sp.Matrix([lam]).T.T*A - flux_row) if False else list((sp.Matrix([lam])*A-flux_row)),lam,dict=True)
    print('k=%d rank A=%d rows=%d; flux in rowspace(A): %s'%(k,A.rank(),len(rs),bool(sol)))
    if sol:
        s=sol[0]; lam_v=[s.get(l,0) for l in lam]; lam_v=[sp.sympify(v).subs({l:0 for l in lam}) for v in lam_v]
        # flux(p) = lam^T A p = -sum lam_r <g,r>_e ; check this vanishes for all mean-zero g in P_k
        gs=[x**j - sp.Rational(1,j+1) for j in range(1,k+1)]
        vals=[sp.simplify(-sum(l*sp.integrate((g*r).subs(y,0),(x,0,1)) for l,r in zip(lam_v,rs))) for g in gs]
        print('   flux of the solution for basis of mean-zero g:',vals)
        # solvability: rhs in column space
        for g in gs:
            b=sp.Matrix([-sp.integrate((g*r).subs(y,0),(x,0,1)) for r in rs])
            assert A.rank()==A.row_join(b).rank()
        print('   solvable for all mean-zero g: True')
