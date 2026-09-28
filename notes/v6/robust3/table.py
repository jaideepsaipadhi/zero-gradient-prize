"""Print the certificate tables of REPORT.tex from the logs."""
import json, glob
rows = []
for f in ["cert_sin_8.jsonl", "cert_sin_16_32.jsonl", "cert_sin_64.jsonl", "cert_sin_128.jsonl", "cert_sin_256.jsonl", "cert_radial.jsonl", "cert_sin_32_gamma.jsonl"]:
    for l in open(f):
        r = json.loads(l); r["file"] = f; rows.append(r)
seen = set()
for r in rows:
    key = (r["phi"], r["N"], r["mu"])
    if key in seen: continue
    seen.add(key)
    s = r["hG"]**4.5 / r["gamma"]
    H1 = r.get("H1"); 
    print(f'{r["phi"]:6s} N={r["N"]:4d} mu={r["mu"]:7.0f} gam={r["gamma"]:8.2f} LB={r["LB"]:.3e} LB/(hG^4.5/gam)={r["LB"]/s:.4f} LB_n={r["LB_n"]:.3f}'
          + (f' H1={H1:.3e} H1/(hG^4.5/gam)={H1/s:.4f} H1_n={r["H1_n"]:.3f} H1/LB={r["ratio"]:.2f} id={r["identity_lhs"]/r["identity_rhs"]-1:.1e}' if H1 else ''))
