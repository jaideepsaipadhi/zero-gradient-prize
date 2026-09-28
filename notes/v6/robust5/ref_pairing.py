"""Reference pairing matrix P[f][j] = < (I - pi_T^) lam_z^j lam_a^(4-j) , g_f >_e  (exact), T^ = (0,0),(1,0),(0,1),
g_1 = B_1 - B_2, g_2 = B_1 - B_3 (quartic Bernstein in s = lam_a on e).  By affine invariance of (I - pi_T) and the
Piola invariance of the normal trace this is the pairing of m_j = l_z^j l_a^(4-j) with the fields v_1, v_2 on EVERY T1."""
import sympy as sp
x, y, s = sp.symbols('x y s')
lz, la = 1 - x - y, x
mons = [x**i * y**j for i in range(4) for j in range(4 - i)]
intT = lambda f: sp.integrate(sp.integrate(sp.expand(f), (y, 0, 1 - x)), (x, 0, 1))
Gm = sp.Matrix(len(mons), len(mons), lambda i, j: intT(mons[i] * mons[j]))
def eta(H):
    c = Gm.LUsolve(sp.Matrix([intT(H * m) for m in mons]))
    return sp.expand((H - sum(c[i] * mons[i] for i in range(len(mons)))).subs({y: 0}).subs(x, s))
B = [sp.binomial(4, i) * s**i * (1 - s)**(4 - i) for i in range(5)]
g = [B[1] - B[2], B[1] - B[3]]
P = sp.Matrix(2, 5, lambda f, j: sp.integrate(eta(lz**j * la**(4 - j)) * g[f], (s, 0, 1)))
print("P =", P)
# kernel: Ker_T = span{m_0 = l_a^4, m_4 = l_z^4, l_c^4 = (l_z + l_a)^4}
kc = sp.Matrix([sp.binomial(4, j) for j in range(5)])
print("P m_0, P m_4, P l_c^4:", list(P[:, 0]), list(P[:, 4]), list(P * kc))
print("P_W (cols m_1, m_2):", P[:, 1:3], " det =", P[:, 1:3].det())
