"""First-moment image M^ = int_F^ g (x) x on Y^1 (Alfeld reference), k=3 and k=4: prints sample image matrices
[[g1.x1, g1.x2],[g2.x1, g2.x2]] (exact). At k=3 the image is spanned by the skew matrix. usage: python3 leak_k3_image.py"""
import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import refalfeld as ra, ratlin as rl
for k in (3,4):
    A = ra.Alfeld(k,'all'); rk,N = ra.nullspace(A.div_rows(), A.nunk)
    W=[('m',(0,0,0)),('m2',(0,1,0)),('m3',(0,0,1))]
    FR=ra.face_moment_rows(A,W)
    L={key: ra.apply_rows([r],N)[0] for key,r in FR.items()}
    K1=rl.nullspace([L[('m',0)],L[('m',1)]])
    restr=lambda key:[sum(a*b for a,b in zip(L[key],v)) for v in K1]
    # M^ = int g (x) x : rows c (component of g), cols (mu2->x1, mu3->x2)
    rows=[restr(('m2',0)),restr(('m3',0)),restr(('m2',1)),restr(('m3',1))]
    print(k,'rank',rl.rank(rows))
    # find nonzero column
    for j in range(len(K1)):
        col=[r[j] for r in rows]
        if any(col):
            print(' sample image [[g1x1,g1x2],[g2x1,g2x2]] =',[str(c) for c in col]); 
            if k==4: break
