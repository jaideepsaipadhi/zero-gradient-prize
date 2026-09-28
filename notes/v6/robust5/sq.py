import sympy as sp
bz,ba=sp.symbols('beta_z beta_a'); bc=1-bz-ba
n=5
I5=[(i,j,n-i-j) for i in range(n,-1,-1) for j in range(n-i,-1,-1)]
d1={a:sp.Symbol('d1_%d%d%d'%a) for a in I5}   # T1: (z,a,c)
d2={a:sp.Symbol('d2_%d%d%d'%a) for a in I5}   # T2: (z,c,b)
eqs=[]
for a in I5:
    iz,ia,ic=a
    if iz<=1 or ic==0: eqs.append(d1[a])
    jz,jc,jb=a
    if jz<=1 or jc<=1: eqs.append(d2[a])
for i in range(n+1):
    j=n-i; eqs.append(d1[(i,0,j)]-d2[(i,j,0)])
for i in range(n):
    j=n-1-i
    eqs.append(d2[(i,j,1)]-(bz*d1[(i+1,0,j)]+ba*d1[(i,1,j)]+bc*d1[(i,0,j+1)]))
unk=[d1[a] for a in I5]+[d2[a] for a in I5]
M=sp.Matrix([[sp.diff(e,u) for u in unk] for e in eqs])
ns=M.nullspace(simplify=True)
print('dim S(Q) =',len(ns))
for v in ns:
    print({str(u):sp.factor(q) for u,q in zip(unk,v) if q!=0})
