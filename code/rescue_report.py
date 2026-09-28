"""Tables for notes/v4/rescue.tex / rescue_report.md from results/v4/rescue/.   python3 rescue_report.py"""
import os, json, glob, math
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "v4", "rescue")
NS = [16, 24, 32, 48, 64, 96]


def load(path):
    return [json.loads(l) for l in open(path)]


def rate(e0, e1, h0, h1):
    return math.log(e0 / e1) / math.log(h0 / h1)


def study_table(kind, alpha, c, theta=0):
    rows = []
    for N in NS:
        p = os.path.join(ROOT, "study", f"{kind}_N{N}_a{alpha:g}_c{c:g}.jsonl")
        if not os.path.exists(p):
            continue
        r = [d for d in load(p) if d["theta"] == theta][0]
        rows.append(r)
    return rows


def fmt_study(kind, theta=0):
    out = [f"\n### {kind}, {'GS' if theta == 0 else 'CNS*'}: H1 error / energy error (rate in h_Gamma)\n",
           "| N | " + " | ".join(f"a={a:g},c={c:g}: mu, H1 (rate), energy (rate)" for a, c in CASES) + " |",
           "|---|" + "---|" * len(CASES)]
    tabs = [study_table(kind, a, c, theta) for a, c in CASES]
    for i, N in enumerate(NS):
        cells = []
        for t in tabs:
            if i >= len(t):
                cells.append("--"); continue
            r = t[i]
            s = f"{r['mu']:.0f}, {r['H1']:.3e}"
            if i > 0:
                s += f" ({rate(t[i-1]['H1'], r['H1'], t[i-1]['hG'], r['hG']):.2f})"
            s += f", {r['energy']:.3e}"
            if i > 0:
                s += f" ({rate(t[i-1]['energy'], r['energy'], t[i-1]['hG'], r['hG']):.2f})"
            cells.append(s)
        out.append(f"| {N} | " + " | ".join(cells) + " |")
    return "\n".join(out)


CASES = [(0, 100), (0.5, 100), (0.5, 1000), (1, 100)]


def fmt_P():
    """Normalised leak constants on test P (G = 1.658...)."""
    G = math.sqrt(7 * math.pi / 8)
    out = ["\n### P1 (GS): H1*mu/(hG) and energy*(mu/h)^{1/2}/G\n", "| N | " + " | ".join(f"a={a:g},c={c:g}" for a, c in CASES) + " |",
           "|---|" + "---|" * len(CASES)]
    tabs = [study_table("P1", a, c, 0) for a, c in CASES]
    for i, N in enumerate(NS):
        cells = []
        for t in tabs:
            if i >= len(t):
                cells.append("--"); continue
            r = t[i]
            cells.append(f"{r['H1'] * r['mu'] / (r['h'] * G):.4f} / {r['energy'] * math.sqrt(r['mu'] / r['h']) / G:.4f}")
        out.append(f"| {N} | " + " | ".join(cells) + " |")
    return "\n".join(out)


def fmt_musweep():
    out = ["\n### mu sweep (fixed N): H1, energy, energy/(gamma^{1/2} hG^{3/2}), slip*(mu/h)^{1/2}\n"]
    for p in sorted(glob.glob(os.path.join(ROOT, "musweep", "*.jsonl"))):
        R = load(p)
        out.append(f"\n{os.path.basename(p)}\n\n| mu | gamma | method | H1 | energy | energy/(g^.5 hG^1.5) | relres |\n|---|---|---|---|---|---|---|")
        for r in R:
            out.append(f"| {r['mu']:.0e} | {r['gamma']:.3g} | {'GS' if r['theta'] == 0 else 'CNS'} | {r['H1']:.4e} | "
                       f"{r['energy']:.4e} | {r['energy'] / (math.sqrt(r['gamma']) * r['hG'] ** 1.5):.4f} | {r['relres']:.1e} |")
    return "\n".join(out)


def fmt_flux():
    out = ["\n### flux correction (CNS*): slip norm ||u~-u_h||_{Gamma_h} vs floor |m_R|/|Gamma_h|^{1/2} and ||u~||_{Gamma_h}\n"]
    for kind in ("F", "F2"):
        out.append(f"\n{kind}\n\n| N | mu | m_R | slip floor | ||u~||_Gh | slip uncorr | slip corr | slipE/floor uncorr | "
                   f"energy uncorr | energy corr | H1 uncorr | H1 corr | flux_G uncorr | flux_G corr |\n" + "|---" * 14 + "|")
        for N in NS:
            p = os.path.join(ROOT, "flux", f"{kind}_N{N}.jsonl")
            if not os.path.exists(p):
                continue
            R = load(p)
            for mu in sorted(set(r["mu"] for r in R)):
                u = [r for r in R if r["mu"] == mu and r["corrected"] == 0][0]
                c = [r for r in R if r["mu"] == mu and r["corrected"] == 1][0]
                slipE = math.sqrt(mu / u["h"]) * u["bL2"]
                out.append(f"| {N} | {mu:.0e} | {u['mR']:.3e} | {u['slip_floor']:.3e} | {u['utrace']:.3e} | {u['bL2']:.4e} | "
                           f"{c['bL2']:.4e} | {slipE / u['floor']:.4f} | {u['energy']:.4e} | {c['energy']:.4e} | "
                           f"{u['H1']:.6e} | {c['H1']:.6e} | {u['fluxG']:.2e} | {c['fluxG']:.1e} |")
    return "\n".join(out)


def fmt_cond():
    out = ["\n### conditioning (k = 4, Lagrange basis, L2-orthonormal pressure basis)\n",
           "| N | h | mu | gamma | lmax(A) | lmin(A) | kappa(A) | lmin(S) | (|G_h|/|O_h|) h/mu | 2nd eig S | lmax(S) | min|eig K| | kappa(K) |",
           "|---" * 13 + "|"]
    for p in sorted(glob.glob(os.path.join(ROOT, "cond", "N*.jsonl")), key=lambda s: int(os.path.basename(s)[1:-6])):
        for r in load(p):
            ratio = 2 * math.pi / (25 - math.pi)          # |Gamma_h|/|Omega_h| ~ 2 pi / (4 L^2 - pi), L = 2.5
            out.append(f"| {r['N']} | {r['h']:.3f} | {r['mu']:.0e} | {r['gamma']:.3g} | {r['lmaxA']:.4g} | {r['lminA']:.4g} | "
                       f"{r['condA']:.3e} | {r['schur_min']:.3e} | {ratio * r['h'] / r['mu']:.3e} | {r['schur_min2']:.3e} | "
                       f"{r['schur_max']:.3f} | {r['minabsK']:.3e} | {r['condK']:.3e} |")
    return "\n".join(out)


if __name__ == "__main__":
    parts = [fmt_study(k) for k in ("A", "P1", "B1", "B100")]
    parts += [fmt_study(k, 1) for k in ("A", "B1")]
    parts += [fmt_P(), fmt_musweep(), fmt_flux(), fmt_cond()]
    txt = "\n".join(parts)
    with open(os.path.join(ROOT, "summary.md"), "w") as f:
        f.write(txt + "\n")
    print(txt)
