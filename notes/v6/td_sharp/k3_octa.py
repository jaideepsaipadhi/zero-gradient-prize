"""k=3, EQUILATERAL boundary faces (octahedron inscribed in the unit sphere): single macro-tet
fields have zero tangential moment there (their moment is prop. to sum l_ij^2 e_ij).
Does the vertex-fan patch (4 macro-tets sharing the edge [e3, 2 e3]) supply both directions?"""
import sys
from alfeld_patch import analyse
k = int(sys.argv[1]) if len(sys.argv) > 1 else 3
z = (0, 0, 1); qa = (0, 0, 2)
ring = [(1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0)]
macro = [[qa, z, ring[j], ring[(j + 1) % 4]] for j in range(4)]
gf = [(j, (1, 2, 3)) for j in range(4)]
print('single macro-tet on an equilateral face:')
analyse([macro[0]], [(0, (1, 2, 3))], k, [(1, 0, 0), (0, 1, 0)])
print('octahedral vertex fan (4 equilateral faces):')
analyse(macro, gf, k, [(1, 0, 0), (0, 1, 0)])
