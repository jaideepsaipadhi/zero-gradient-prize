"""structured strip: W columns of width 1, J layers of height H each, cells split by forward diagonal
(or alternating). Nitsche on y=0, Dirichlet elsewhere."""
import numpy as np, sys
from patch import lamstar
def strip(W, J, H, alt=False, k=4, shear=0.0):
    idx = lambda i, j: j*(W+1)+i
    P = [(i + shear*j*H, j*H) for j in range(J+1) for i in range(W+1)]
    tris = []
    for j in range(J):
        for i in range(W):
            a,b,c,d = idx(i,j), idx(i+1,j), idx(i+1,j+1), idx(i,j+1)
            if alt and (i+j)%2: tris += [(a,b,d),(b,c,d)]
            else: tris += [(a,b,c),(a,c,d)]
    nit = {(idx(i,0), idx(i+1,0)): (0,-1) for i in range(W)}
    return lamstar(P, tris, nit, k)
if __name__ == "__main__":
    for H in [0.5, 0.75, 1, 1.5, 2, 3, 4, 6]:
        row = []
        for (W,J) in [(2,1),(4,1),(8,1),(4,2),(8,2),(8,3)]:
            row.append(strip(W,J,H)[0])
        print(H, np.round(row,3), flush=True)
