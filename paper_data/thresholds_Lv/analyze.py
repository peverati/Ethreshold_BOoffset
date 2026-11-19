"""Analysis of thresholds_Lv exact runs -> RESULT.md (and fits_<target>.txt).

Written 2026-11-15.  Run from the project root:  python3 paper_data/thresholds_Lv/analyze.py
Works on partial data (missing stages are reported as missing).

Per target (J1v1, L3A, L4S, L1S):
  * proven upper bound: smallest fixed-basis crossing (tol 1e-12, margin 0; rigorous up to rounding);
    rounding-safe bound: margin 1e-11 Eh on the largest basis.
  * extrapolated M_c: best crossing + geometric extrapolation of the compact-size sequence
    N_D - 2, N_D, N_D + 2 (all diffuse sets), error = |extrapolation| (as paper_data/rotational).
    Diffuse convergence: last two diffuse increments.
  * threshold law (L > 0: binding linear in M - M_c): quadratic fit of E_b(M) over x <= 0.05,
    root = independent estimate (biased by basis incompleteness at the smallest x).
  * comparison with the frozen blind prediction (sha256 checked) and with BO, BO+ad, literature.
"""
import csv
import hashlib
import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TARGETS = ["J1v1", "L3A", "L4S", "L1S"]
PRED_KEY = {"J1v1": "J1v1", "L3A": "J3v0", "L4S": "J4v0"}
NAMES = {"J1v1": "(J=1, v=1): L=1 'A', root 1", "L3A": "(J=3, v=0): L=3 'A', root 0",
         "L4S": "(J=4, v=0): L=4 'S', root 0", "L1S": "2p sigma_u J=1: L=1 'S', root 0"}
BOK = {"J1v1": ("g L=1", 1), "L3A": ("g L=3", 0), "L4S": ("g L=4", 0), "L1S": ("u L=1", 0)}
LIT = {"J1v1": "Frolov & Wardlaw arXiv:1008.3010: bound up to m_mu ~ 214.72 m_e, M_c ~ 17.09 (basis-limited upper estimate)",
       "L3A": "ddmu (M = 17.752) has no bound J=3 level, ttmu (26.585) has one: 17.752 < M_c < 26.585",
       "L1S": "none checked"}


def rows(path):
    return list(csv.DictReader(open(path))) if os.path.exists(path) else []


def geom_extrap(v):
    """Three values of a nested sequence -> (limit, |correction|, ratio) or None."""
    if len(v) < 3:
        return None
    d1, d2 = v[1] - v[0], v[2] - v[1]
    if d1 == 0 or d2 / d1 <= 0 or d2 / d1 >= 1:
        return (v[2], abs(d2), float("nan"))
    q = d2 / d1
    corr = d2 * q / (1 - q)
    return (v[2] + corr, abs(corr), q)


def analyze(t, out):
    P = out.append
    cr = [r for r in rows(os.path.join(HERE, f"{t}_crossings.csv"))]
    bd = rows(os.path.join(HERE, f"{t}_binding.csv"))
    P(f"\n## {NAMES[t]}\n")
    if not cr:
        P("No exact data yet.")
        return None
    std = [r for r in cr if float(r["tol"]) == 1e-12 and float(r["margin"]) == 0.0]
    fin = [(r["label"], int(r["size"]), int(r["kept"]), float(r["M_star"])) for r in std]
    P("| basis | size | kept | M* (upper bound) |\n|---|---|---|---|")
    seen = set()
    for lab, n, k, m in fin:
        if lab in seen:
            continue
        seen.add(lab)
        P(f"| {lab} | {n} | {k} | {m:.10f} |")
    good = [m for *_, m in fin if np.isfinite(m)]
    res = dict(target=t)
    if not good:
        P("\nNo finite crossing.")
        return res
    best = min(good)
    res["best"] = best
    P(f"\n* **best proven upper bound (up to rounding): {best:.9f}**")
    safe = [float(r["M_star"]) for r in cr if float(r["margin"]) == 1e-11 and np.isfinite(float(r["M_star"]))]
    if safe:
        res["safe"] = min(safe)
        P(f"* rounding-safe bound (margin 1e-11 Eh): {min(safe):.9f}")
    tl = [(r["label"], float(r["M_star"])) for r in cr if float(r["tol"]) == 1e-10 and float(r["margin"]) == 0.0]
    for lab, m in tl:
        P(f"* tol 1e-10 on '{lab}': {m:.10f} (difference from tol 1e-12 bound: {m - best:+.1e})")
    # compact-size check with all diffuse sets
    full = {}
    for lab, n, k, m in fin:
        if "diffuse" in lab:
            N = int(lab.split("N=")[1].split()[0])
            nd = int(lab.split("+ ")[1].split()[0])
            full.setdefault(nd, {})[N] = m
    ext = None
    if full:
        ndmax = max(full)
        seq = [full[ndmax][N] for N in sorted(full[ndmax])]
        if len(seq) >= 3:
            ext = geom_extrap(seq[-3:])
            res["extrap"] = ext[0]
            res["extrap_err"] = ext[1]
            P(f"* compact-size sequence with {ndmax} diffuse sets, N = {sorted(full[ndmax])}: "
              + ", ".join(f"{x:.9f}" for x in seq) + f"; ratio {ext[2]:.3f}")
            P(f"* **extrapolated M_c = {ext[0]:.8f} +/- {max(ext[1], 1e-9):.1e}** (geometric in N; assumes geometric "
              f"convergence, the bound above is rigorous)")
        # diffuse convergence at N_D
        ND = None
        for lab, n, k, m in fin:
            if "diffuse" in lab:
                ND = int(lab.split("N=")[1].split()[0])
                break
        dseq = [(nd, full[nd].get(ND)) for nd in sorted(full) if full[nd].get(ND) is not None]
        if len(dseq) >= 2:
            P("* diffuse increments at N_D: " + ", ".join(f"{a}:{b:.10f}" for a, b in dseq)
              + f"; last increment {dseq[-1][1] - dseq[-2][1]:+.1e}")
            res["diffuse_last"] = dseq[-1][1] - dseq[-2][1]
    if bd:
        M = np.array([float(r["M"]) for r in bd])
        Eb = np.array([float(r["binding"]) for r in bd])
        o = np.argsort(M)
        M, Eb = M[o], Eb[o]
        sel = M <= best * 1.051
        P("\nBinding energies (largest basis):\n\n| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |\n|---|---|---|---|")
        for m, e in zip(M, Eb):
            P(f"| {m:.8f} | {m / best - 1:.4f} | {e:.6e} | {e / (m - best):.4e} |")
        if sel.sum() >= 3:
            p = np.polyfit(M[sel], Eb[sel], 2)
            r = [x.real for x in np.roots(p) if abs(x.imag) < 1e-12 and abs(x.real - best) < 0.05 * best]
            if r:
                rr = min(r, key=lambda x: abs(x - best))
                res["fit_root"] = rr
                P(f"\n* quadratic fit (x <= 0.05) root {rr:.8f}; slope at root {np.polyval(np.polyder(p), rr):.4e} Eh per mass unit "
                  f"(linear threshold law for L > 0)")
    return res


def main():
    bo = json.load(open(os.path.join(HERE, "bo_brackets.json")))
    out = ["# Result: exact critical masses of L > 0 levels of {X+, X+, e-} (thresholds_Lv)", "",
           f"*Last updated: {time.strftime('%Y-%m-%d %H:%M')}*", "",
           "Written by `analyze.py` from the run_exact.py CSV files. Crossings M* of fixed bases are rigorous upper "
           "bounds to M_c (variational root below E_th proves binding), up to rounding."]
    pj = os.path.join(HERE, "prediction_J1v1.json")
    pred, okhash = None, None
    if os.path.exists(pj):
        pred = json.load(open(pj))["predictions"]
        want = open(os.path.join(HERE, "prediction_J1v1.sha256")).read().split()[0]
        have = hashlib.sha256(open(os.path.join(HERE, "PREDICTION_J1v1.md"), "rb").read()).hexdigest()
        okhash = want == have
        out += ["", "## Blind-test integrity", "",
                f"* sha256(PREDICTION_J1v1.md) = {have}: {'matches' if okhash else 'DOES NOT MATCH'} prediction_J1v1.sha256.",
                "* Frozen 2026-11-15 20:17. Caveat (FREEZE_NOTES.md): the J1v1 scale optimization (stage 1, N = 12) ran 20:11-20:13, "
                "before the freeze; it computed no crossing."]
    summary = []
    for t in TARGETS:
        r = analyze(t, out)
        key, v = BOK[t]
        Mbo = bo.get(f"{key} BO", [None] * 2)[v]
        Mad = bo.get(f"{key} BO+ad", [None] * 2)[v] if f"{key} BO+ad" in bo else None
        out.append(f"\n* BO {Mbo:.6f}" + (f", BO+ad {Mad:.6f}" if Mad else " (BO+ad not available for 2p sigma_u)"))
        if t in LIT:
            out.append(f"* literature / consistency: {LIT[t]}")
        if r and "best" in r:
            Mx = r.get("extrap", r["best"])
            line = f"| {t} | {Mbo:.5f} | {Mad if Mad is None else f'{Mad:.5f}'} | {r['best']:.7f} | {Mx:.7f} |"
            if pred and t in PRED_KEY:
                q = pred[PRED_KEY[t]]
                d = Mx - q["pf_nonpert"]
                ok = abs(d) <= q["unc"] + r.get("extrap_err", 0)
                line += f" {q['pf_nonpert']:.4f} +/- {q['unc']:.4f} | {d:+.4f} | {'PASS' if ok else 'FAIL'} |"
                if Mad:
                    out.append(f"* non-adiabatic offset M_c(BO+ad) - M_c(exact) = {Mad - Mx:.5f} "
                               f"(predicted {Mad - q['pf_nonpert']:.5f}); offset vs BO {Mbo - Mx:.5f}")
            else:
                line += " | | |"
            summary.append(line)
    out += ["", "## Summary", "",
            "| target | BO | BO+ad | best upper bound | extrapolated M_c | prediction (frozen) | exact - pred | pass |",
            "|---|---|---|---|---|---|---|---|"] + summary
    open(os.path.join(HERE, "RESULT.md"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
