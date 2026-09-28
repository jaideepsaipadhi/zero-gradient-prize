"""NUMERICAL: as graded_check.py but with mesh size TRIPLING per layer (row j>=1: squares of side 3^j over three
row-(j-1) squares; two hanging bottom nodes n1,n2; triangles (bl,n1,tl),(n1,n2,tl),(n2,tr,tl),(n2,br,tr);
min angle about 15.3 deg)."""
import numpy as np, math
from graded_check import quotient


def triadic_mesh(J, W):
    T = []
    for i in range(-W, W):
        a, b, c, d = (i, 0), (i + 1, 0), (i + 1, 1), (i, 1)
        T += [(a, b, c), (a, c, d)]
    y0 = 1
    for j in range(1, J + 1):
        s = 3**j
        for i in range(-W // s, W // s):
            x0 = i * s
            bl, br, tr, tl = (x0, y0), (x0 + s, y0), (x0 + s, y0 + s), (x0, y0 + s)
            n1, n2 = (x0 + s / 3, y0), (x0 + 2 * s / 3, y0)
            T += [(bl, n1, tl), (n1, n2, tl), (n2, tr, tl), (n2, br, tr)]
        y0 += s
    return [np.array(t, float) for t in T]


if __name__ == "__main__":
    T = triadic_mesh(5, 243)
    ang = 180
    for P in T:
        for i in range(3):
            u = P[(i + 1) % 3] - P[i]; w = P[(i + 2) % 3] - P[i]
            ang = min(ang, math.degrees(math.acos(u @ w / np.linalg.norm(u) / np.linalg.norm(w))))
    print("min angle", ang)
    for Y in (1, 2, 3, 4, 6, 9, 27):
        q, k = quotient(T, Y)
        print(f"Y={Y:4.1f}  Y(2B-A)/C = {q:9.4f}   (max diam/Y {k:5.2f})", flush=True)
