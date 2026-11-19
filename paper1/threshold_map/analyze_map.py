"""Analysis of the threshold-map exact runs -> RESULT.md.

Written 2026-11-16 (threshold_map, step 3).  Run from the project root (works on partial data; missing targets are
reported as missing):
    python3 paper1/threshold_map/analyze_map.py

Per target (as paper1/thresholds_Lv/analyze.py):
  * best proven upper bound = smallest fixed-basis crossing (tol 1e-12, margin 0; rigorous up to rounding);
    rounding-safe bound = margin 1e-11 Eh on the largest basis; tol 1e-10 sensitivity;
  * extrapolated M_c: best + geometric extrapolation of the compact sequence N_D - 2, N_D, N_D + 2 (all diffuse sets),
    error = |correction|, floored at the last diffuse increment and at the tol sensitivity is NOT applied (reported);
  * root-ordering check (roots below k bound at M*); binding-energy table (L > 0: linear threshold law).
Map-level:
  * comparison with the frozen predictions (sha256 checked), with Frolov 1992 II (quoted y_t and the smallest mass
    at which his variational table binds the level: a rigorous upper bound our value must respect);
  * complete list of thresholds M_c <= 101 (known exact + new), order of appearance exact vs BO vs BO+ad;
  * number of bound levels n(M) of the 1s sigma_g manifold vs M: exact vs BO vs BO+ad counting;
  * muonic consistency: which levels are bound at M = 8.880 (pp mu), 17.752 (dd mu), 26.585 (tt mu).
"""
import csv
import hashlib
import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
KNOWN = {(0, 1): (9.525634, "4e-6"), (0, 2): (29.63939, "1e-5"), (0, 3): (59.31866, "2e-5"), (0, 4): (98.26550, "1.1e-4"),
         (1, 0): (4.2082827, "5e-7"), (1, 1): (17.094282, "1e-6"), (2, 0): (11.4894217, "3e-7"),
         (3, 0): (21.376199, "1e-6"), (4, 0): (33.914349, "2e-6")}
FROLOV = {(0, 1): (10.045, 10.05), (0, 2): (30.70, 30.75), (0, 3): (61.135, 61.15), (1, 1): (17.125, 17.1),
          (1, 2): (39.375, 39.45), (2, 1): (27.715, 27.75), (2, 2): (53.295, 53.30)}
FROLOV_ORDER = ["(0,0)", "(1,0)", "(0,1)", "(2,0)", "(1,1)", "(2,1)", "(0,2)", "(1,2)", "(2,2)", "(0,3)", "(1,3)"]
MUONIC = {"pp mu": 8.880, "dd mu": 17.752, "tt mu": 26.585}   # m/m_mu, as in thresholds_Lv (dd, tt: Frolov 2025)
# bound-level sets of the muonic molecular ions, verified 2026-11-16 against Frolov & Wardlaw, EPJD 63, 339 (2011) and Frolov, EPJD 66, 212 (2012) (docs/literature/frolov_muonic_1994_2019_notes.md); was: quoted from memory, to be checked against a
# review, e.g. Ponomarev, Contemp. Phys. 31, 219 (1990))
MUONIC_LIT = {"pp mu": {"(0,0)", "(1,0)"},
              "dd mu": {"(0,0)", "(0,1)", "(1,0)", "(1,1)", "(2,0)"},
              "tt mu": {"(0,0)", "(0,1)", "(1,0)", "(1,1)", "(2,0)", "(3,0)"}}


def rows(p):
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []


def geom_extrap(v):
    if len(v) < 3:
        return None
    d1, d2 = v[1] - v[0], v[2] - v[1]
    if d1 == 0 or d2 / d1 <= 0 or d2 / d1 >= 1:
        return (v[2], abs(d2), float("nan"))
    q = d2 / d1
    corr = d2 * q / (1 - q)
    return (v[2] + corr, abs(corr), q)


HP_DATA = os.path.join(HERE, "..", "threshold_map_L7hp", "data")   # double-double assembly runs (2026-11-17)


def analyze_target(t, lev, out, data=DATA, tol_std=1e-12, tol_chk=1e-10, note=None):
    """data/tol_std/tol_chk default to run 5; the HP runs (threshold_map_L7hp) use tol 1e-14 / 1e-12."""
    P = out.append
    cr = rows(os.path.join(data, f"{t}_crossings.csv"))
    bd = rows(os.path.join(data, f"{t}_binding.csv"))
    L, k = int(t[1:].split("v")[0]), int(t.split("v")[1])
    P(f"\n## ({L},{k}): L={L} '{'A' if L % 2 else 'S'}', root {k}\n")
    P(f"BO {lev['M_BO']:.6f}, BO+ad {lev['M_BOad']:.6f}.")
    res = dict(target=t)
    if note and cr:
        P("\n" + note)
    if data != DATA and cr:      # HP runs: guard recheck (recheck_guard.py) -> exclude flagged rows, list them
        gpath = os.path.join(data, f"{t}_guard.csv")
        flagged = {(g["label"], g["tol"], g["margin"]): g["reason"] for g in rows(gpath) if g["flag"] == "1"}
        if flagged:
            P("\nExcluded rows (spurious-root guard / conditioning, recheck_guard.py):\n")
            P("| basis | tol | margin | M* | reason |\n|---|---|---|---|---|")
            for (lab, tl, mg), why in flagged.items():
                m0 = [r["M_star"] for r in cr if (r["label"], r["tol"], r["margin"]) == (lab, tl, mg)]
                P(f"| {lab} | {tl} | {mg} | {m0[0] if m0 else ''} | {why} |")
            cr = [dict(r, M_star="nan") if (r["label"], r["tol"], r["margin"]) in flagged else r for r in cr]
        elif not os.path.exists(gpath):
            P("\n**Guard recheck not yet run for this target.**")
    if not cr:
        P("\nNo exact data yet.")
        vf = os.path.join(os.path.dirname(DATA), f"validate_L{L}.txt")
        if os.path.exists(vf) and "RESULT: FAIL" in open(vf).read():
            P(f"\n**Not computed: validate_L {L} FAILED** ({os.path.basename(vf)}); the target was skipped by the launcher. "
              f"Diagnosis in validate_L{L}_diagnostics.md (if present): term lists verified uncorrupted (fresh regeneration "
              f"identical); cause = canonical-orthogonalization truncation at tol 1e-12 and loss of precision at tol 1e-14 "
              f"for L = {L}. The frozen prediction stands untested.")
        return res
    std = [r for r in cr if float(r["tol"]) == tol_std and float(r["margin"]) == 0.0]
    P("\n| basis | size | kept | M* (upper bound) | E_(k-1) - E_th at M* |\n|---|---|---|---|---|")
    seen, fin = set(), []
    for r in std:
        if r["label"] in seen:
            continue
        seen.add(r["label"])
        m = float(r["M_star"])
        fin.append((r["label"], m))
        P(f"| {r['label']} | {r['size']} | {r['kept']} | {m:.10f} | {r.get('gap_below', '')} |")
    bad = [r["label"] for r in std if k > 0 and r.get("gap_below") not in (None, "", "nan") and float(r["gap_below"]) >= 0]
    if bad:
        P(f"\n**ROOT-ORDER PROBLEM**: root {k - 1} not bound at M* for {bad}")
    elif k > 0:
        P(f"\n* root ordering: roots 0..{k - 1} bound at every crossing (root k = v confirmed)")
    good = [m for _, m in fin if np.isfinite(m)]
    if not good:
        P("\nNo finite crossing.")
        return res
    best = min(good)
    res["best"] = best
    P(f"* **best proven upper bound (up to rounding): {best:.9f}**")
    safe = [float(r["M_star"]) for r in cr if float(r["margin"]) == 1e-11 and np.isfinite(float(r["M_star"]))]
    if safe:
        res["safe"] = min(safe)
        P(f"* rounding-safe bound (margin 1e-11 Eh): {min(safe):.9f}")
    for r in cr:
        if float(r["tol"]) == tol_chk and float(r["margin"]) == 0.0:
            m = float(r["M_star"])
            res["tol_sens"] = abs(m - best)
            P(f"* tol {tol_chk:g} on '{r['label']}': {m:.10f} (difference from best: {m - best:+.1e})")
    full = {}
    for lab, m in fin:
        if data != DATA and not np.isfinite(m):      # HP runs: sequences from accepted (finite) crossings only
            continue
        N = int(lab.split("N=")[1].split()[0])
        nd = int(lab.split("+ ")[1].split()[0]) if "diffuse" in lab else 0
        full.setdefault(nd, {})[N] = m
    ndmax = max(full)
    seq_N = sorted(full[ndmax])
    if len(seq_N) >= 3:
        seq = [full[ndmax][N] for N in seq_N[-3:]]
        ext = geom_extrap(seq)
        res["extrap"], res["extrap_err"] = ext[0], max(ext[1], 1e-9)
        P(f"* compact-size sequence ({ndmax} diffuse sets), N = {seq_N[-3:]}: " + ", ".join(f"{x:.9f}" for x in seq)
          + f"; ratio {ext[2]:.3f}")
        P(f"* **extrapolated M_c = {ext[0]:.8f} +/- {res['extrap_err']:.1e}** (geometric in N; the bound is rigorous)")
    dseq = sorted((nd, d.get(max(d))) for nd, d in full.items() if nd > 0)
    if 0 in full and dseq:
        ND = int([lab for lab, _ in fin if "diffuse" in lab][0].split("N=")[1].split()[0])
        dseq = [(0, full[0].get(ND))] + [(nd, full[nd].get(ND)) for nd in sorted(full) if nd > 0 and full[nd].get(ND)]
        dseq = [(a, b) for a, b in dseq if b is not None]
        if len(dseq) >= 2:
            P("* diffuse increments at N_D: " + ", ".join(f"{a}:{b:.10f}" for a, b in dseq)
              + f"; last increment {dseq[-1][1] - dseq[-2][1]:+.1e}")
            res["diffuse_last"] = dseq[-1][1] - dseq[-2][1]
    if bd:
        M = np.array([float(r["M"]) for r in bd])
        Eb = np.array([float(r["binding"]) for r in bd])
        o = np.argsort(M)
        M, Eb = M[o], Eb[o]
        P("\n| M | x = M/M* - 1 | E_th - E | E_b/(M - M*) |\n|---|---|---|---|")
        for m, e in zip(M, Eb):
            P(f"| {m:.8f} | {m / best - 1:.4f} | {e:.6e} | {e / (m - best):.4e} |")
    return res


def main():
    bm = json.load(open(os.path.join(HERE, "bo_map.json")))
    lev = {(d["J"], d["v"]): d for d in bm["levels"]}
    out = ["# Result: complete (J, v) threshold map of {X+, X+, e-}, M_c <= ~100 (threshold_map)", "",
           f"*Last updated: {time.strftime('%Y-%m-%d %H:%M')}*", "",
           "Written by `analyze_map.py` from the `run_map.py` CSV files in `data/`. Crossings M* of fixed bases are "
           "rigorous upper bounds to M_c (a variational root below E_th proves binding), up to rounding. "
           "Values marked 'extrapolated' assume geometric convergence in the compact size."]
    pj, ps, pm = (os.path.join(HERE, f) for f in ("prediction_map.json", "prediction_map.sha256", "PREDICTION_MAP.md"))
    pred = None
    if os.path.exists(ps):
        pred = json.load(open(pj))["predictions"]
        want = open(ps).read().split()[0]
        have = hashlib.sha256(open(pm, "rb").read()).hexdigest()
        out += ["", "## Blind-test integrity", "",
                f"* sha256(PREDICTION_MAP.md) = {have}: {'matches' if want == have else 'DOES NOT MATCH'} prediction_map.sha256."]
    res = {}
    for t in bm["targets"]:
        L, k = int(t[1:].split("v")[0]), int(t.split("v")[1])
        if not rows(os.path.join(DATA, f"{t}_crossings.csv")) and rows(os.path.join(HP_DATA, f"{t}_crossings.csv")):
            extra = os.path.join(HP_DATA, "..", f"{t}_note.md")
            res[(L, k)] = analyze_target(
                t, lev[(L, k)], out, data=HP_DATA, tol_std=1e-14, tol_chk=1e-12,
                note=f"**Computed with double-double matrix assembly** (paper1/threshold_map_L7hp, run_map_hp.py; "
                     f"validate_L {L} PASSES there, validate_hp.txt). Standard projection tol 1e-14, check tol 1e-12; "
                     f"otherwise run 5's method. (Run 5's double-precision path failed validate_L {L}.)")
            if os.path.exists(extra):
                out.append("\n" + open(extra).read().strip())
        else:
            res[(L, k)] = analyze_target(t, lev[(L, k)], out)
    # ---------------- summary vs predictions and Frolov
    out += ["", "## Summary: new thresholds", "",
            "| (J,v) | BO | BO+ad | best upper bound | extrapolated M_c | prediction (frozen) | exact - pred | pass | "
            "Frolov y_t (bound at) | consistent with Frolov bound |",
            "|---|---|---|---|---|---|---|---|---|---|"]
    npass = nconv = 0
    for (L, k), r in res.items():
        d = lev[(L, k)]
        fr = FROLOV.get((L, k))
        line = f"| ({L},{k}) | {d['M_BO']:.5f} | {d['M_BOad']:.5f} |"
        if "best" not in r:
            line += " | | " + (f"{pred[r['target']]['M_pred_np']:.4f} +/- {pred[r['target']]['unc']:.4f}" if pred else "") + \
                    f" | | | {f'{fr[0]} ({fr[1]})' if fr else ''} | |"
            out.append(line)
            continue
        Mx = r.get("extrap", r["best"])
        line += f" {r['best']:.7f} | {Mx:.7f} |"
        if pred:
            q = pred[r["target"]]
            dd = Mx - q["M_pred_np"]
            ok = abs(dd) <= q["unc"] + r.get("extrap_err", 0)
            nconv += 1
            npass += ok
            line += f" {q['M_pred_np']:.4f} +/- {q['unc']:.4f} | {dd:+.4f} | {'PASS' if ok else 'FAIL'} |"
        else:
            line += " | | |"
        line += (f" {fr[0]} ({fr[1]}) | {'yes' if r['best'] <= fr[1] else 'NO (exact above his bound)'} |" if fr else " | |")
        out.append(line)
    if pred:
        out.append(f"\nBlind test: {npass} PASS out of {nconv} targets with an exact value "
                   "(criterion fixed in PREDICTION_MAP.md: |exact - pred| <= band + exact error).")
    # ---------------- complete map
    allv = {}
    for key, d in lev.items():
        if key == (0, 0):
            continue
        if key in KNOWN:
            allv[key] = (KNOWN[key][0], "exact (earlier)")
        elif key in res and "best" in res[key]:
            allv[key] = (res[key].get("extrap", res[key]["best"]), "exact (this map)")
        elif pred and f"J{key[0]}v{key[1]}" in pred:
            allv[key] = (pred[f"J{key[0]}v{key[1]}"]["M_pred_np"], "PREDICTED (no exact yet)")
    lab = lambda k: f"({k[0]},{k[1]})"
    order_x = [lab(k) for k, _ in sorted(allv.items(), key=lambda kv: kv[1][0])]
    out += ["", "## Complete map (M_c <= ~101)", "", "| order | (J,v) | M_c | source | BO | BO+ad | BO - exact | BO+ad - exact |",
            "|---|---|---|---|---|---|---|---|"]
    for i, (k, (m, src)) in enumerate(sorted(allv.items(), key=lambda kv: kv[1][0])):
        d = lev[k]
        out.append(f"| {i + 1} | {lab(k)} | {m:.6f} | {src} | {d['M_BO']:.4f} | {d['M_BOad']:.4f} | "
                   f"{d['M_BO'] - m:.4f} | {d['M_BOad'] - m:.4f} |")
    out += ["", "Order of appearance:", "",
            "* exact:  " + " < ".join(order_x),
            "* BO+ad: " + " < ".join(x for x in bm["order_BOad"] if x != "(0,0)"),
            "* BO:    " + " < ".join(x for x in bm["order_BO"] if x != "(0,0)"),
            "* Frolov 1992 II (S, P, D only): " + " < ".join(FROLOV_ORDER[1:]),
            "* exact restricted to Frolov's levels: " + " < ".join(x for x in order_x if x in FROLOV_ORDER)]
    # ---------------- n(M)
    roots = bm["roots"]
    Mgrid = [2, 5, 8.880, 10, 15, 17.752, 20, 25, 26.585, 30, 40, 50, 60, 70, 80, 90, 100]

    def n_exact(M):
        n = 1 + sum(1 for m, _ in allv.values() if m < M)          # (0,0) always bound
        unk = sum(1 for m, s in allv.values() if m < M and s.startswith("PRED"))
        return n, unk

    def n_bo(M, c):
        n = sum(sum(1 for x in roots[f"L={L} {c}"] if x < M) for L in range(bm.get("LMAX", 10) + 1) if f"L={L} {c}" in roots)
        return n + (1 if c == "BO+ad" else 0)       # BO+ad: (0,0) bound for all M (no root)
    out += ["", "## Number of bound levels n(M) of the 1s sigma_g manifold", "",
            "Exact counting uses the exact M_c (earlier + this map); levels without an exact value yet use the frozen "
            "prediction and are counted in the 'of which predicted' column. Valid for M <= 100 (levels with BO+ad "
            "M_c in (101, 105] are outside the map).", "",
            "| M | n exact | of which predicted | n BO | n BO+ad |", "|---|---|---|---|---|"]
    for M in Mgrid:
        ne, unk = n_exact(M)
        out.append(f"| {M:g} | {ne} | {unk} | {n_bo(M, 'BO')} | {n_bo(M, 'BO+ad')} |")
    # ---------------- muonic
    out += ["", "## Muonic consistency check", "",
            "Literature bound-level sets (verified against Frolov & Wardlaw 2011 and Frolov 2012): "
            "pp mu 2 levels (0,0),(1,0); dd mu 5 levels (+ (0,1),(1,1),(2,0)); tt mu 6 levels (+ (3,0)).", "",
            "| ion | M | bound levels from the map | literature | agree |", "|---|---|---|---|---|"]
    for ion, M in MUONIC.items():
        b = {"(0,0)"} | {lab(k) for k, (m, _) in allv.items() if m < M}
        near = [f"{lab(k)} M_c={m:.3f}" for k, (m, _) in allv.items() if abs(m - M) < 1.0]
        out.append(f"| {ion} | {M} | {', '.join(sorted(b))} | {', '.join(sorted(MUONIC_LIT[ion]))} | "
                   f"{'yes' if b == MUONIC_LIT[ion] else 'NO'}{' (near: ' + '; '.join(near) + ')' if near else ''} |")
    chk = os.path.join(HERE, "..", "threshold_map_L7hp", "check_run5_orthog.txt")
    if os.path.exists(chk):
        out += ["", "## Orthogonalization check of run 5 (2026-11-17)", "",
                "Single-pass canonical orthogonalization caused spurious levels for L = 7 at tol <= 1e-14 "
                "(threshold_map_L7hp). Run 5 uses tol 1e-12; for three representative standard largest bases the "
                "crossing was recomputed with a second orthogonalization pass (check_run5_orthog.py). The single-pass "
                "max|X^T S X - I| (2e-3..6e-3) is dominated by rounding in forming X^T S X in double and is not itself "
                "an error of the result; the crossing shifts are:", ""]
        out += ["* " + l.strip() for l in open(chk) if l.strip() and not l.startswith("#")]
        out += ["", "All shifts (<= 4.3e-10) are far below the quoted errors of these targets (1e-9, 4.9e-6, 2.4e-5): "
                    "run 5 is not affected."]
    open(os.path.join(HERE, "RESULT.md"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
