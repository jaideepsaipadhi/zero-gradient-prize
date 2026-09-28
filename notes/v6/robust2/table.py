"""Tabulate normalised quantities from leak_*.jsonl (REPORT.tex Sec. 5)."""
import json, glob, sys
rows = []
for f in sorted(glob.glob('leak_*.jsonl')):
    for l in open(f):
        if l.startswith('{'):
            r = json.loads(l)
            if 'Wgam' in r: rows.append(r)
seen = set(); uniq = []
for r in rows:
    key = (r['phi'], r['mu'], r['N'])
    if key not in seen: seen.add(key); uniq.append(r)
uniq.sort(key=lambda r: (r['phi'], r['mu'], r['N']))
print("phi    mu      N  gamma   |eta|/hG^k  slip/(e|eta|) osc/(e|eta|) H1/(e hG^-1/2|eta|) En/(e^1/2|eta|)  leak    H1")
for r in uniq:
    eps = r['h'] / r['mu']; e = r['eta']; k = r['k']
    print(f"{r['phi']:5s} {r['mu']:7.0f} {r['N']:3d} {r['gamma']:7.1f}  {e/r['hG']**k:.4f}     {r['Wgam']/(eps*e):7.3f}     {r['Wosc']/(eps*e):7.3f}      {r['H1']/(eps*r['hG']**-0.5*e):7.3f}         {r['energy']/(eps**0.5*e):7.3f}     {r['leak']:.3f}  {r['H1']:.3e}")
