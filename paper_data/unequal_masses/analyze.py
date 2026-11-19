"""Analysis of the unequal_masses exact runs -> RESULT.md; blind comparison with PREDICTION_UNEQUAL.md.

Written 2026-11-16.  Run from the project root (works on partial data; missing stages are reported):
    python3 paper_data/unequal_masses/analyze.py

Per target:
  * proven upper bound: smallest fixed-basis crossing (tol 1e-12, margin 0); rounding-safe bound (margin 1e-11);
  * extrapolated M_c: best crossing + geometric extrapolation of the compact sequence N_D-2, N_D, N_D+2 (all
    diffuse sets); error = |extrapolation| + |diffuse-enlarged change| + (L = 0) last ladder step + tol change;
  * threshold-law fit from the stage-5 binding table: L = 0 sqrt(E_b) = c1 dM + c2 dM^2 (+ c3 dM^3),
    L > 0 E_b = c1 dM + c2 dM^2 (independent estimate; biased by basis incompleteness at the smallest x);
  * blind comparison: sha256 of PREDICTION_UNEQUAL.md / prediction_unequal.json vs prediction_unequal.sha256,
    freeze time vs first non-dry run entry; theory (2ch+m) pass/fail (|exact - pred| <= unc); null models NS, NT
    and persistence NP (offset(v) = offset(v-1), offset = M_BOad - M_c) with the same tolerance.
"""
import csv
import hashlib
import json
import os
import re
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ORDER = ["J0v1", "J0v2", "J0v3", "J0v4", "J1v0", "J2v0", "S15J0v1", "R1J0v1"]
PKEY = {"J0v1": "r2_J0v1", "J0v2": "r2_J0v2", "J0v3": "r2_J0v3", "J0v4": "r2_J0v4", "J1v0": "r2_J1v0",
        "J2v0": "r2_J2v0", "S15J0v1": "r1.5_J0v1"}
SYM_EXACT_R1 = 9.525634          # paper_data/critical_masses (symmetric family, v = 1)


def rows(path):
    return list(csv.DictReader(open(path))) if os.path.exists(path) else []


def geom_extrap(v):
    """Geometric extrapolation of a monotone 3-term sequence; returns (limit, |correction|) or (last, nan)."""
    if len(v) < 3 or not all(np.isfinite(v)):
        return (v[-1] if v else np.nan), np.nan
    d1, d2 = v[1] - v[0], v[2] - v[1]
    if d1 == 0 or d2 / d1 <= 0 or d2 / d1 >= 1:
        return v[-1], abs(d2)
    q = d2 / d1
    corr = d2 * q / (1 - q)
    return v[-1] + corr, abs(corr)


def blind_check():
    md, js, sha = (os.path.join(HERE, f) for f in ("PREDICTION_UNEQUAL.md", "prediction_unequal.json",
                                                   "prediction_unequal.sha256"))
    if not os.path.exists(sha):
        return False, "NO FROZEN PREDICTION (prediction_unequal.sha256 missing)"
    want = dict(reversed(l.split()) for l in open(sha).read().splitlines() if l.strip())
    ok = True
    msg = []
    for f in (md, js):
        h = hashlib.sha256(open(f, "rb").read()).hexdigest()
        good = want.get(os.path.basename(f)) == h
        ok &= good
        msg.append(f"{os.path.basename(f)} sha256 {'OK' if good else 'MISMATCH'}")
    frozen = re.search(r"Frozen (\d{4}-\d\d-\d\d \d\d:\d\d)", open(md).read())
    first = None
    for t in ORDER:
        for line in (open(os.path.join(HERE, f"{t}_log.txt")).read().splitlines()
                     if os.path.exists(os.path.join(HERE, f"{t}_log.txt")) else []):
            if "--dry" in line or "calibration" in line or line[20:].startswith("  "):
                continue
            if "=== run_exact.py" in line or "stage" in line or "scales" in line:
                ts = line[:16]
                first = ts if first is None or ts < first else first
                break
    if frozen and first:
        good = frozen.group(1) <= first
        ok &= good
        msg.append(f"frozen {frozen.group(1)}, first exact run entry {first}: {'OK' if good else 'RUN BEFORE FREEZE'}")
    return ok, "; ".join(msg)


def fit_threshold(t, L):
    B = rows(os.path.join(HERE, f"{t}_binding.csv"))
    if len(B) < 4:
        return np.nan, "binding table missing/short"
    M = np.array([float(r["M"]) for r in B])
    Eb = np.array([float(r["binding"]) for r in B])
    o = np.argsort(M)
    M, Eb = M[o], Eb[o]
    sel = Eb > 0
    M, Eb = M[sel], Eb[sel]
    if len(M) < 4:
        return np.nan, "too few bound points"
    out = []
    for order in (2, 3):
        y = np.sqrt(Eb) if L == 0 else Eb
        n = min(len(M), 6 + order)
        c = np.polyfit(M[:n], y[:n], order)
        rts = [x.real for x in np.roots(c) if abs(x.imag) < 1e-9 and x.real < M[0] and x.real > 0.8 * M[0]]
        if rts:
            out.append(max(rts))
    if not out:
        return np.nan, "no root"
    return float(np.mean(out)), f"fit roots {['%.6f' % x for x in out]} ({'sqrt(E_b)' if L == 0 else 'E_b'} polynomial, lowest points)"


ENLARGE_NOTE = ["", "Conclusions (2026-11-16 17:30):", "",
    "* r = 1 (R1J0v1): the production basis gives 9.5873; the symmetric exact value is 9.525634. Cause confirmed: at "
    "r = 1 the 'closed' channel (B far from A e) is degenerate with the open one and lacked the diffuse ladder. "
    "Mirror basis M: 9.52565573 (N = 14/16/18: 9.52565574, 9.52565573, 9.52565550; n_d + 20/40 at margin 1e-11: "
    "unchanged to 2e-10); M + (n_c + 2, n_d + 20) (MB, 14000 functions): 9.52564891. Remaining difference to "
    "9.525634: +2.2e-5 (M), +1.5e-5 (MB). The last enlargement step was -6.8e-6, so the residual is about 2x that "
    "step: convergence toward the exact value, but not yet within the variant's own error, which would be ~7e-6 "
    "(compact 2e-7 + enlargement 6.8e-6). The remaining slow direction is the atom polarization in the "
    "diffuse blocks (n_c). Polarization alone without the mirror (P6, P8) does not help (9.5873, 9.5871).",
    "* r = 1, MP6 (mirror + n_c 6): both the margin-0 and margin-1e-11 searches ended on a flat near-threshold "
    "level (gap ~1e-11, below M ~ 9.3) and were rejected by guard (c). Unverified interpretation: the larger "
    "unsymmetrized basis binds a very weak u-type threshold state that the symmetric code excludes, so root k = 1 "
    "is ambiguous at r = 1 only. Not relevant for r > 1, where there is no degeneracy.",
    "* r = 1.5 (S15J0v1): shifts M -3.5e-8, P6 +1.2e-6, MP6 +9.6e-8; r = 2 (J0v1): M +6.1e-7, P6 +7.7e-7, MP6 "
    "+3.8e-7. All below 1e-5, and positive (tol 1e-12 drops more functions in the larger bases, so not strictly "
    "variational at the 1e-6 level). The production r = 2 and r = 1.5 results stand; no rerun was needed. These "
    "~1e-6 shifts lie within the quoted J0v1 and S15J0v1 errors (7.7e-6, 7.1e-6). The test is L = 0 only: J1v0 and "
    "J2v0 (quoted errors 2e-8, 9e-7) are not covered by it.",
]
SPURIOUS = []   # (target, label, tol, margin, M*, reason): every row excluded from the analysis
FLAGGED = []    # same, rows flagged by one soft rule but kept (guard-checked)
NOISE = 1.63e-13    # max build-to-build |dE| at genuine J0v1 crossings (noise_check.py, BLAS threads 2/5/10)
CERT_MARGIN = max(1e-11, 10 * NOISE)    # certified-bound margin (2026-11-16): = 1e-11
TOL_JUMP = 1e-3     # tol 1e-12 vs 1e-10 crossing of the same basis differing by more -> tol 1e-12 rows excluded
OUTLIER = 1e-2      # a diffuse-basis crossing below all other diffuse-basis crossings by more -> excluded


def screen(t, C):
    """Consistency screen (2026-11-16, after the J0v1 noise crossing).  Returns the rows kept; every excluded
    row is recorded in SPURIOUS with its reason.  Rules: (1) rows flagged spurious=1 by the run_exact.py guard;
    (2) tol 1e-12 vs tol 1e-10 crossing of the same basis (margin 0) differing by > TOL_JUMP: the tol 1e-12 rows
    of that basis are flagged; (3) a tol-1e-12 margin-0 crossing of a basis with diffuse functions lying below
    all other such crossings of the target by > OUTLIER (all rows of that basis and tol are flagged);
    (4) a margin-1e-11 row below the margin-0 row of the same basis and tol (a margin shift can only raise M*).
    Revision 2026-11-16 17:00: rules (1) and (4) exclude; rules (2) and (3) exclude a row on their own only if it
    was written before guard (c) existed (spurious column empty); a row checked by the guard is excluded only when
    (2) and (3) both flag it, otherwise it is kept and listed as "flagged, kept" (R1J0v1: the enlarged-diffuse
    row 9.5590 passed the guard and is a genuine, better upper bound; the r = 1 basis was not converged)."""
    out = {}
    soft = {}

    def drop(r, why):
        k = (r["label"], r["tol"], r["margin"])
        if k not in out:
            out[k] = why
            SPURIOUS.append((t, r["label"], r["tol"], r["margin"], r["M_star"], why))

    def flag(r, why):
        if r.get("spurious", "") == "":
            drop(r, why + " [pre-guard row]")
        else:
            soft.setdefault((r["label"], r["tol"], r["margin"]), []).append((r, why))

    fin = lambda r: np.isfinite(float(r["M_star"]))
    for r in C:
        if r.get("spurious", "") == "1":
            drop(r, "run_exact.py guard: " + r.get("reason", ""))
    m0 = {(r["label"], r["tol"]): float(r["M_star"]) for r in C if r["margin"] == repr(0.0) and fin(r)}
    for (lab, tol), v in m0.items():
        v10 = m0.get((lab, repr(1e-10)))
        if tol == repr(1e-12) and v10 is not None and abs(v - v10) > TOL_JUMP:
            for r in C:
                if r["label"] == lab and r["tol"] == tol:
                    flag(r, f"tol 1e-12 crossing {v:.8f} vs tol 1e-10 {v10:.8f} of the same basis (|diff| > {TOL_JUMP:g})")
    dif = {lab: v for (lab, tol), v in m0.items() if tol == repr(1e-12) and "diffuse" in lab}
    for lab, v in dif.items():
        others = [w for l2, w in dif.items() if l2 != lab]
        if others and v < min(others) - OUTLIER:
            for r in C:
                if r["label"] == lab and r["tol"] == repr(1e-12):
                    flag(r, f"crossing {v:.8f} below all other diffuse-basis crossings (min {min(others):.8f}) by > {OUTLIER:g}")
    for r in C:
        v0 = m0.get((r["label"], r["tol"]))
        if r["margin"] == repr(1e-11) and fin(r) and v0 is not None and float(r["M_star"]) < v0 - 1e-9:
            drop(r, f"margin-1e-11 crossing below the margin-0 crossing {v0:.10f} of the same basis")
    for k, lst in soft.items():
        if k in out:
            continue
        r = lst[0][0]
        why = "; ".join(w for _, w in lst)
        if len(lst) >= 2:
            drop(r, "two screen rules agree: " + why)
        else:
            FLAGGED.append((t, r["label"], r["tol"], r["margin"], r["M_star"], why))
    return [r for r in C if (r["label"], r["tol"], r["margin"]) not in out]


def analyze(t):
    C = [r for r in rows(os.path.join(HERE, f"{t}_crossings.csv")) if r["target"] == t]
    if not C:
        return None
    n_all = len(C)
    C = screen(t, C)
    ok = lambda r, tol=1e-12, m=0.0: r["tol"] == repr(tol) and r["margin"] == repr(m) and np.isfinite(float(r["M_star"]))
    std = [r for r in C if ok(r)]
    best = min(std, key=lambda r: float(r["M_star"])) if std else None
    safe = [float(r["M_star"]) for r in C if ok(r, 1e-12, 1e-11)]
    m0lab = {r["label"]: float(r["M_star"]) for r in std}
    cert = [(r["label"], m0lab[r["label"]], float(r["M_star"])) for r in C
            if ok(r, 1e-12, CERT_MARGIN) and r["label"] in m0lab]
    tol10 = [float(r["M_star"]) for r in C if ok(r, 1e-10, 0.0)]
    labs = {r["label"]: float(r["M_star"]) for r in std}
    # compact-size sequence with all diffuse sets
    base = [l for l in labs if l.startswith("compact N=") and ("diffuse" in l) and "enlarged" not in l]
    seqN = {}
    for l in base:
        m = re.match(r"compact N=(\d+) \+ (.*)", l)
        seqN.setdefault(m.group(2), {})[int(m.group(1))] = labs[l]
    seq_key = max(seqN, key=lambda k: len(seqN[k])) if seqN else None
    seq = [seqN[seq_key][n] for n in sorted(seqN[seq_key])] if seq_key else []
    extrap, corr = geom_extrap(seq) if len(seq) >= 3 else ((seq[-1] if seq else np.nan), np.nan)
    enl = [v for l, v in labs.items() if "enlarged" in l]
    d_enl = abs(min(enl) - best_val) if (enl and best) and (best_val := float(best["M_star"])) else np.nan
    ladder = sorted([(l, v) for l, v in labs.items() if "alpha_s=" in l and "enlarged" not in l], key=lambda x: x[1])
    d_lad = abs(ladder[1][1] - ladder[0][1]) if len(ladder) >= 2 else np.nan
    return dict(n=n_all, best=float(best["M_star"]) if best else np.nan, best_label=best["label"] if best else "",
                safe=min(safe) if safe else np.nan, cert=cert, tol10=min(tol10) if tol10 else np.nan, seq=seq,
                extrap=extrap, corr=corr, d_enl=d_enl, d_lad=d_lad)


def main():
    pred = json.load(open(os.path.join(HERE, "prediction_unequal.json"))) if os.path.exists(
        os.path.join(HERE, "prediction_unequal.json")) else None
    blind_ok, blind_msg = blind_check()
    now = time.strftime("%Y-%m-%d %H:%M")
    L = [f"# unequal_masses: exact critical masses and blind comparison", "", f"*Last updated: {now}*", "",
         f"Generated by analyze.py. Blind check: **{'OK' if blind_ok else 'NOT OK'}** ({blind_msg}).", "",
         "M_c = mass of the light nucleus A (heavy B = rM), threshold -rM/(2(rM+1)).", "",
         "| target | runs | best upper bound | rounding-safe | tol 1e-10 | compact seq (N_D-2, N_D, N_D+2) | extrap. M_c | error | threshold-law fit |",
         "|---|---|---|---|---|---|---|---|---|"]
    res = {}
    for t in ORDER:
        a = analyze(t)
        if a is None:
            L.append(f"| {t} | not run | | | | | | | |")
            continue
        Lq = {"J1v0": 1, "J2v0": 2}.get(t, 0)
        fit, fmsg = fit_threshold(t, Lq)
        err = np.nansum([a["corr"], a["d_enl"], a["d_lad"] if Lq == 0 else 0.0, abs(a["tol10"] - a["best"]) if np.isfinite(a["tol10"]) else 0.0])
        Mc = a["extrap"] if np.isfinite(a["extrap"]) else a["best"]
        res[t] = dict(Mc=Mc, err=err, best=a["best"], fit=fit, cert=a["cert"])
        L.append(f"| {t} | {a['n']} | {a['best']:.8f} | {a['safe']:.8f} | {a['tol10']:.8f} | "
                 f"{', '.join('%.7f' % x for x in a['seq'])} | {Mc:.6f} | {err:.1e} | {fit:.6f} ({fmsg}) |")
    L += ["", "## Excluded rows (spurious-root guard and consistency screen)", ""]
    if SPURIOUS:
        L += ["| target | basis | tol | margin | M* | reason |", "|---|---|---|---|---|---|"]
        L += [f"| {t} | {lab} | {tol} | {m} | {v} | {why} |" for t, lab, tol, m, v, why in SPURIOUS]
    else:
        L += ["None."]
    if FLAGGED:
        L += ["", "Flagged by one consistency rule but kept (rows checked by guard (c); revision 2026-11-16 17:00: the "
                  "tol-jump and outlier rules exclude a guard-checked row only when both agree):", "",
              "| target | basis | tol | margin | M* | flag |", "|---|---|---|---|---|---|"]
        L += [f"| {t} | {lab} | {tol} | {m} | {v} | {why} |" for t, lab, tol, m, v, why in FLAGGED]
    L += ["", "## Certified bounds (rounding noise)", "",
          f"Build-to-build noise of E_k - E_th at genuine crossings (J0v1, alpha_s = 1e-3 and 1e-4, rebuilt with BLAS "
          f"thread counts 2, 5, 10 and 5 again; noise_check.py / noise_check_output.jsonl): max |dE| = {NOISE:.2e} Eh "
          "(2.8e-14 at M*, 1.6e-13 at M*(1+1e-3)); identical thread counts give bit-identical results. The ~1e-8 "
          "effects seen during the J0v1 incident belong to the alpha_s=3e-4 artefact level only (rejected by the guard). "
          f"Certified margin = max(1e-11, 10 x noise) = {CERT_MARGIN:.0e} Eh, i.e. unchanged: the certified bound is the "
          "crossing of E_k = E_th - margin (tol 1e-12) of the stage-4 standard basis.", "",
          "dM = M*(margin) - M*(0) is the actual shift. For L = 0 it is NOT margin/slope: the gap starts quadratically "
          "at the crossing (J0v1, alpha_s=1e-4: gap(M*(1-1e-3)) = +5.8e-13, gap(M*(1+1e-3)) = -3.9e-8, secant slope "
          "1.7e-6 Eh per unit M, so margin/slope = 6e-6, but the actual shift is 1.7e-4). L > 0 gaps are linear.", "",
          "| target | basis | M* (margin 0) | certified bound M* (margin 1e-11) | dM |", "|---|---|---|---|---|"]
    for t in ORDER:
        for lab, v0, vc in (res.get(t, {}).get("cert") or []):
            L.append(f"| {t} | {lab} | {v0:.8f} | {vc:.8f} | {vc - v0:.2e} |")
        if t in res and not res[t].get("cert"):
            L.append(f"| {t} | (margin rows excluded by the screen) | | | |")
    L += ["", "## Diffuse-basis enlargement test (2026-11-16, enlarge_test.py)", "",
          "Variants of the stage-4 standard basis (compact N = 16, alpha_s = 1e-4): M = mirror (closed channel gets the "
          "open-channel sets 0.15, 0.03, 1e-4), P6/P8 = atom-coordinate degree n_c 4 -> 6/8 on all diffuse blocks, "
          "MP6 = M + P6, MB = M + (n_c + 2, n_d + 20), MD20/MD40 = M + n_d + 20/40, MN14/MN18 = M with compact N = 14/18; "
          "suffix @m = crossing at margin 1e-11 (noise-robust). shift = M* - production margin-0 value.", "",
          "| target | variant | size | kept | M* | shift | note |", "|---|---|---|---|---|---|---|"]
    for t in ORDER:
        for r in rows(os.path.join(HERE, f"enlarge_{t}.csv")):
            note = ("REJECTED by guard: " + r["reason"]) if r["spurious"] == "1" else ""
            L.append(f"| {t} | {r['variant']} | {r['size']} | {r['kept']} | {float(r['M_star']):.8f} | "
                     f"{float(r['shift']):+.2e} | {note} |")
    L += ENLARGE_NOTE
    L += ["", "Note on J0v1 (2026-11-16): the first J0v1 run (before the guard existed) had a noise crossing for "
              "\"compact N=16 + diffuse alpha_s=0.0003\" at tol 1e-12 (M* = 9.1976224, margin 1e-11: 9.2172467). It was not "
              "a low spurious E_0 but a near-threshold level whose E - E_th is flat in M and not reproducible between "
              "builds (ill-conditioned projection; check_guard.py). Because it was the lowest ladder crossing, the original "
              "stage 4 used alpha_s=0.0003 instead of the genuine best ladder member; stages 3-5 of J0v1 were discarded "
              "(moved to Trash) and redone with the guard (run_exact.py, check (c): crossing slope)."]
    if "R1J0v1" in res:
        L += ["", f"r = 1 check (unsymmetrized basis): M_c = {res['R1J0v1']['Mc']:.6f} vs symmetric exact {SYM_EXACT_R1} "
                  f"(diff {res['R1J0v1']['Mc'] - SYM_EXACT_R1:+.2e})."]
    if pred:
        L += ["", "## Blind comparison (frozen predictions)", "",
              "| target | exact M_c | theory 2ch+m | unc | exact - theory | PASS? (unc / 0.10) | null NS (diff) | null NT (diff) | "
              "offset M_BOad - M_c | persistence pred. | persistence diff | BO | BO+ad | BO+ad thr |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
        prev_off = None
        for t in ORDER:
            if t not in PKEY:
                continue
            q = pred["targets"][PKEY[t]]
            if t not in res:
                L.append(f"| {t} | not run | {q['M_pred']:.4f} | {q['unc']:.4f} | | | {q['null_NS']:.4f} | {q['null_NT']:.4f} | | | | "
                         f"{q['M']['BO']:.4f} | {q['M']['BOad']:.4f} | {q['M']['BOad_thr']:.4f} |")
                prev_off = None
                continue
            x = res[t]["Mc"]
            d = x - q["M_pred"]
            off = q["M"]["BOad"] - x
            if t in ("J0v2", "J0v3", "J0v4") and prev_off is not None:
                pers = q["M"]["BOad"] - prev_off
                pdiff = f"{x - pers:+.4f}"
                pstr = f"{pers:.4f}"
            else:
                pstr, pdiff = "", ""
            prev_off = off if t.startswith("J0v") else prev_off
            L.append(f"| {t} | {x:.6f} ± {res[t]['err']:.1e} | {q['M_pred']:.4f} | {q['unc']:.4f} | {d:+.4f} | "
                     f"{'PASS' if abs(d) <= q['unc'] else 'FAIL'} / {'PASS' if abs(d) <= 0.10 else 'FAIL'} | {q['null_NS']:.4f} ({x - q['null_NS']:+.3f}) | "
                     f"{q['null_NT']:.4f} ({x - q['null_NT']:+.3f}) | {off:.4f} | {pstr} | {pdiff} | {q['M']['BO']:.4f} | "
                     f"{q['M']['BOad']:.4f} | {q['M']['BOad_thr']:.4f} |")
        L += ["", "A null model 'passes' a level when |exact - null| <= the theory's frozen unc for that level.",
              "Persistence uses the exact offset of the previous J=0 level of the same family (needs both runs)."]
    open(os.path.join(HERE, "RESULT.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
