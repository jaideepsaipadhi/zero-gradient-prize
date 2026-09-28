"""Flat fans with EQUAL weights (1,1,1) on every face: single-macro fields contribute nothing,
so this isolates the moment carried by fields that cross interior macro faces."""
from patchlib import *
from fractions import Fraction as Fr
def fan(ring,q,z=(0,0,0)):
    m=len(ring); return [[q,z,ring[j],ring[(j+1)%m]] for j in range(m)],[(j,(1,2,3)) for j in range(m)]
def run(name,ring,q,frame=[(1,0,0),(0,1,0)]):
    mac,gf=fan(ring,q); m=len(ring)
    P=Patch(mac,gf,3,weights={j:(1,1,1) for j in range(m)})
    print(f'{name:45s}', moment_rank_mod(P,frame), kappa_float(P,frame).round(6))
run('square, apex on axis (singular edge)',[(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)],(0,0,1))
run('square, apex off axis',[(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)],(Fr(1,5),Fr(1,7),1))
run('kite 4-star straight diagonals, apex (1/5,1/7,1)',[(1,0,0),(0,2,0),(-3,0,0),(0,-1,0)],(Fr(1,5),Fr(1,7),1))
run('4-star non-straight, apex on axis',[(1,0,0),(Fr(1,5),1,0),(-1,Fr(1,3),0),(0,-1,0)],(0,0,1))
run('pentagon-ish, apex on axis',[(1,0,0),(Fr(3,10),1,0),(-Fr(4,5),Fr(3,5),0),(-Fr(4,5),-Fr(3,5),0),(Fr(3,10),-1,0)],(0,0,1))
run('hexagon-ish, apex on axis',[(1,0,0),(Fr(1,2),Fr(7,8),0),(-Fr(1,2),Fr(7,8),0),(-1,0,0),(-Fr(1,2),-Fr(7,8),0),(Fr(1,2),-Fr(7,8),0)],(0,0,1))
run('3-star, apex on axis',[(1,0,0),(-Fr(1,2),Fr(7,8),0),(-Fr(1,2),-Fr(7,8),0)],(0,0,1))
