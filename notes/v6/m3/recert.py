"""Re-run the v4 interval certificates (code/m3_certify.py, unchanged) for given targets, as a reproducibility check.
Run: python3 recert.py 30 10        (about 3-5 min on 1 core)
Pod (all v4 targets, ~20 min): cd code && python3 m3_certify.py certify"""
import sys, time
from fractions import Fraction as Fr
sys.path.insert(0, '/home/claude/zero-gradient-prize/code')
import m3_certify as mc
td, dd = int(sys.argv[1]), int(sys.argv[2])
for fam, k0 in [(f, k) for (a, b, f, k) in mc.TARGETS if (a, b) == (td, dd)]:
    t0 = time.time(); ok, info = mc.branch_and_bound(td, dd, fam, k0, progress=False)
    print(f"{'F3' if fam=='m3' else 'Fs'}({td},{dd}) kappa_w >= {float(k0)}: CERTIFIED={ok} {info} [{time.time()-t0:.0f}s]", flush=True)
