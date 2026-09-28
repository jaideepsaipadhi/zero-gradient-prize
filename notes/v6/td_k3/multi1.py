"""Rings with 2 or more Gamma faces: total tangential moment rank (weights l^2), flat and curved."""
from patchlib import *
from fractions import Fraction as Fr
import numpy as np, random
def ringP(link,gam,q,z=(0,0,0)):
    m=len(link)
    mac=[[q,z,link[j],link[(j+1)%m]] for j in range(m)]
    return Patch(mac,[(j,(1,2,3)) for j in gam],3)
fr=[(1,0,0),(0,1,0)]
# 2 adjacent Gamma faces [z,x0,x1],[z,x1,x2] flat, others lifted
L=[(1,0,0),(0,1,0),(-1,0,0),(-Fr(1,2),-1,Fr(1,2)),(Fr(1,2),-1,Fr(1,3))]
for q in [(0,0,1),(Fr(1,5),Fr(1,7),1)]:
    P=ringP(L,[0,1],q); print('2 adjacent flat, q',q,moment_rank_mod(P,fr),kappa_float(P,fr).round(5))
# 2 adjacent, equilateral-ish? use lattice plane x+y+z=0 : faces [0,u,u+w],[0,u+w,w]
u=(1,-1,0); w=(0,1,-1); n=(1,1,1)
def Lt(a,b,h=0): return tuple(a*u[c]+b*w[c]+h*n[c] for c in range(3))
fr2=[(1,-1,0),(1,1,-2)]
link=[Lt(1,0),Lt(1,1),Lt(0,1),Lt(-1,0,Fr(1,2)),Lt(-1,-1,Fr(1,2)),Lt(0,-1,Fr(1,3))]
for q in [Lt(0,0,1),Lt(Fr(1,3),Fr(1,3),1),Lt(Fr(1,5),-Fr(1,7),1)]:
    P=ringP(link,[0,1],q); print('2 adjacent equilateral flat, q',q,moment_rank_mod(P,fr2),kappa_float(P,fr2).round(5))
    P=ringP(link,[0],q); print('   1 face (same ring)          ',moment_rank_mod(P,fr2),kappa_float(P,fr2).round(5))
# full hexagonal equilateral fan (6 Gamma faces), apex on axis / off
hexl=[Lt(1,0),Lt(1,1),Lt(0,1),Lt(-1,0),Lt(-1,-1),Lt(0,-1)]
for q in [Lt(0,0,1),Lt(Fr(1,3),Fr(1,3),1)]:
    P=ringP(hexl,range(6),q); print('hex fan equilateral, q',q,moment_rank_mod(P,fr2),kappa_float(P,fr2).round(5))
