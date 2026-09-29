import mpmath as mp
# Stokes Dirichlet corner, velocity exponent lambda: sin^2(lambda*w) = lambda^2 sin^2 w
for w,guess in [(mp.pi/2,2.7+1.1j),(mp.pi/2,4.8+1.46j),(mp.pi/2,6.8+1.7j),(3*mp.pi/2,0.54)]:
    f=lambda l: mp.sin(l*w)**2-l**2*mp.sin(w)**2
    r=mp.findroot(f,mp.mpc(guess)); print(w, r, abs(f(r)))
