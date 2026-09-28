"""Table of the uniform-mesh runs (robust4 REPORT Sec. 6): H1n = ||grad W|| gamma / hG^4.5 and the gamma-normalised
observed order p = 4.5 + log(H1n_i/H1n_j)/log(hG_i/hG_j).   python table4.py [files...]"""
import sys, json, math, glob
files = sys.argv[1:] or ["u_8_64.jsonl", "u_mid.jsonl"] + glob.glob("u_pod.jsonl")
rows = []
for f in files:
    try:
        rows += [json.loads(l) for l in open(f)]
    except FileNotFoundError:
        pass
rows = [r for r in rows if r["mu"] == 1000]
for p in ["sin", "Q", "B", "A", "R"]:
    q = sorted({r["N"]: r for r in rows if r["phi"] == p}.values(), key=lambda r: r["N"])
    print(p.ljust(3), "N     :", " ".join("%7d" % r["N"] for r in q))
    print("    H1n   :", " ".join("%7.4f" % r["H1n"] for r in q))
    print("    order :", "        " + " ".join("%7.2f" % (4.5 + math.log(q[i]["H1n"] / q[i + 1]["H1n"]) / math.log(q[i]["hG"] / q[i + 1]["hG"])) for i in range(len(q) - 1)))
    print("    osc(eta):", " ".join("%7.3f" % r["eta_osc_frac"] for r in q))
    print("    eta/h^4 :", " ".join("%7.4f" % r["eta_h4"] for r in q))
