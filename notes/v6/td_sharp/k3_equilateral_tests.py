"""k=3: total circumsphere-weighted moment on Y^1 for fans of EQUILATERAL faces inscribed in a sphere."""
import sympy as sp
from alfeld_patch import analyse
R = sp.Rational
print('octahedral fan, off-axis apex (1/5,1/7,2):')
ring = [(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)]
q = (R(1,5), R(1,7), 2)
analyse([[q,(0,0,1),ring[j],ring[(j+1)%4]] for j in range(4)], [(j,(1,2,3)) for j in range(4)], 3, [(1,0,0),(0,1,0)])
print('regular-tetrahedron fan (3 equilateral faces), apex 2*(1,1,1):')
z = (1,1,1); rr = [(1,-1,-1),(-1,1,-1),(-1,-1,1)]
analyse([[(2,2,2),z,rr[j],rr[(j+1)%3]] for j in range(3)], [(j,(1,2,3)) for j in range(3)], 3, [(1,-1,0),(1,1,-2)])
print('octahedral fan, apex (0,0,2), k=4 (control):')
analyse([[(0,0,2),(0,0,1),ring[j],ring[(j+1)%4]] for j in range(4)], [(j,(1,2,3)) for j in range(4)], 4, [(1,0,0),(0,1,0)])
