"""(Optional) Worsey-Farin split, k=3.  Generic sub-tet patch code: fields C^0 piecewise P_k on a given list of
sub-tetrahedra, zero on given boundary triangles, div-free.  Gamma macro face F (split into 3 sub-faces
[m_F,v_i,v_j]); constraints: zero mean of d_nu v on EACH sub-face (they are the Gamma_h faces);
moment: int_F Q_F d_nu v with Q_F = 1/2 sum l_ij^2 mu_i mu_j (mu = barycentrics of the macro face F)."""
from fractions import Fraction as Fr
import numpy as np, itertools, sys
from patchlib import F3, sub, add, scal, dot, cross, bary_grads, mindex, multinom, tri_mono_int, rank_mod, P_MOD

class GPatch:
    def __init__(self, subs, zero_tris, gam, k=3):
        """subs: list of 4-point lists (vertex 0 must be the one opposite a Gamma sub-face if used);
        zero_tris: list of 3-point triangles on which v=0 (patch boundary); gam: list of
        (sub index, macro face (P1,P2,P3)) -- the sub-face is the face of subs[i] opposite subs[i][0]."""
        self.k=k; self.subs=[[F3(p) for p in S] for S in subs]
        zero=set()
        for tri in zero_tris:
            P=[F3(p) for p in tri]
            for b in mindex(k,3): zero.add(tuple(sum(Fr(b[j],k)*P[j][c] for j in range(3)) for c in range(3)))
        self.key={}; self.sidx=[]
        for V in self.subs:
            d={}
            for a in mindex(k,4):
                pt=tuple(sum(Fr(a[j],k)*V[j][c] for j in range(4)) for c in range(3))
                if pt in zero: d[a]=None
                else:
                    if pt not in self.key: self.key[pt]=len(self.key)
                    d[a]=self.key[pt]
            self.sidx.append(d)
        self.nunk=3*len(self.key)
        self.grads=[bary_grads(V) for V in self.subs]
        self.div_rows=[]
        for V,d,(G,D) in zip(self.subs,self.sidx,self.grads):
            for b in mindex(k-1,4):
                r={}
                for j in range(4):
                    a=tuple(b[m]+(1 if m==j else 0) for m in range(4)); n=d[a]
                    if n is None: continue
                    for c in range(3):
                        if G[j][c]!=0: r[3*n+c]=r.get(3*n+c,0)+k*G[j][c]
                self.div_rows.append(r)
        self.cons=[]; self.mom=[]
        for (si,Fm) in gam:
            V=self.subs[si]; G,D=self.grads[si]; P=V[1:]
            Fm=[F3(p) for p in Fm]
            N=cross(sub(P[1],P[0]),sub(P[2],P[0])); fac=(dot(N,N)/4)/(abs(D)/2)
            # macro-face barycentrics of the sub-face vertices
            NF=cross(sub(Fm[1],Fm[0]),sub(Fm[2],Fm[0])); nn=dot(NF,NF)
            def mbary(x):
                l1=dot(cross(sub(Fm[1],x),sub(Fm[2],x)),NF)/nn; l2=dot(cross(sub(Fm[2],x),sub(Fm[0],x)),NF)/nn
                return (l1,l2,1-l1-l2)
            MB=[mbary(p) for p in P]   # MB[k][i] = mu_i(p_k)
            l2={(0,1):dot(sub(Fm[0],Fm[1]),sub(Fm[0],Fm[1])),(0,2):dot(sub(Fm[0],Fm[2]),sub(Fm[0],Fm[2])),(1,2):dot(sub(Fm[1],Fm[2]),sub(Fm[1],Fm[2]))}
            # Q in sub-face barycentrics nu: coefficient of nu_k nu_l (k<=l)
            Qc={}
            for (i,j),lv in l2.items():
                for kk in range(3):
                    for ll in range(3):
                        key=tuple(sorted((kk,ll))); Qc[key]=Qc.get(key,0)+Fr(1,2)*lv*MB[kk][i]*MB[ll][j]
            d=self.sidx[si]; crow=[dict() for _ in range(3)]; mrow=[dict() for _ in range(3)]
            for b in mindex(k-1,3):
                n=d[(1,)+b]
                if n is None: continue
                cb=multinom(b); i0=cb*tri_mono_int(b); i1=0
                for (kk,ll),qc in Qc.items():
                    g=list(b); g[kk]+=1; g[ll]+=1; i1+=qc*cb*tri_mono_int(tuple(g))
                for c in range(3):
                    crow[c][3*n+c]=fac*2*k*i0; mrow[c][3*n+c]=fac*2*k*i1
            self.cons.append(crow); self.mom.append(mrow)
    def dense(self,rows):
        M=np.zeros((len(rows),self.nunk),dtype=np.int64)
        for i,r in enumerate(rows):
            for c,v in r.items(): M[i,c]=(v.numerator%P_MOD)*pow(v.denominator%P_MOD,P_MOD-2,P_MOD)%P_MOD
        return M
    def ranks(self,frame):
        D=self.dense(self.div_rows); C=self.dense(sum(self.cons,[]))
        Mr=[]
        for t in frame:
            r={}
            for mrow in self.mom:
                for c in range(3):
                    for col,v in mrow[c].items(): r[col]=r.get(col,0)+Fr(t[c])*v
            Mr.append(r)
        M=self.dense(Mr)
        rd=rank_mod(D); rdc=rank_mod(np.vstack([D,C])); ra=rank_mod(np.vstack([D,C,M]))
        return dict(nunk=self.nunk,rank_div=rd,div_rows=len(self.div_rows),rank_cons=rdc-rd,moment_rank=ra-rdc)

def wf_subs(T, Bp, fpts):
    """WF split of macro T (4 pts) with interior point Bp and face points fpts[i] (face opposite T[i]).
    returns list of (sub-tet, info) with sub-tet = [Bp, m_i, v_a, v_b]."""
    out=[]
    for i in range(4):
        fv=[T[j] for j in range(4) if j!=i]
        for a,b in [(0,1),(1,2),(0,2)]:
            out.append(([Bp,fpts[i],fv[a],fv[b]],i))
    return out

if __name__=='__main__':
    k=int(sys.argv[1]) if len(sys.argv)>1 else 3
    # single macro, general shape, Gamma face opposite vertex 0
    T=[F3(p) for p in [(0,0,3),(0,0,0),(5,1,0),(1,4,0)]]
    Bp=scal(Fr(1,4),add(add(T[0],T[1]),add(T[2],T[3])))
    fp=[scal(Fr(1,3),add(add(*[T[j] for j in range(4) if j!=i][:2]),[T[j] for j in range(4) if j!=i][2])) for i in range(4)]
    S=wf_subs(T,Bp,fp)
    subs=[s for s,_ in S]
    zero=[s[1:] for s in subs]   # all sub-faces opposite the interior point = boundary of T
    gam=[(n,[T[1],T[2],T[3]]) for n,(s,i) in enumerate(S) if i==0]
    P=GPatch(subs,zero,gam,k)
    print('WF single macro k=%d:'%k,P.ranks([(1,0,0),(0,1,0)]))
