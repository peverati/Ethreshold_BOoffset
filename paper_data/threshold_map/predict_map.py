"""Blind, parameter-free predictions of the exact critical masses of every target (J, v) of the threshold map.

Written 2026-11-16 (threshold_map, step 2).  Same code path as paper_data/thresholds_Lv/predict.py (frozen J1v1 test,
passed) and paper_data/critical_masses_v4/predict.py (v = 4, passed): BO+ad zero-energy coupling problem
(paper_data/bo_offset/coupling.py, lam = M + 1/2), threshold weights w_rad = u'^2/(lam B), w_cen = L(L+1) u^2/(r^2 lam B),
folded with the first-principles mass functions m_v(R), m_r(R) (paper_data/bo_offset/massfunc/results_fine.csv),
first order and non-perturbatively:
    predicted M_c(exact) = M_c(BO+ad) - predicted non-adiabatic offset.
No parameter is fitted.  Uncertainty band FIXED IN ADVANCE: the largest |exact - predicted| over all nine levels
with a known exact M_c (it is L=1 v=0, +0.0295), applied uniformly to every target.

Run from the project root (~3-6 min):
    python3 paper_data/threshold_map/predict_map.py [--freeze]
--freeze writes PREDICTION_MAP.md + prediction_map.sha256 and makes the outputs read-only (chmod 444); refuses if
prediction_map.sha256 already exists.
"""
import hashlib
import json
import os
import stat
import sys
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "paper_data", "bo_offset"))
sys.path.insert(0, os.path.join(ROOT, "paper_data", "bo_offset", "massfunc"))
sys.path.insert(0, os.path.join(ROOT, "solvers"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import coupling as C  # noqa: E402
import analyze as A   # noqa: E402

# exact critical masses known before this prediction (calibration of the error bar only; nothing is fitted)
EXACT = {(0, 1): 9.525634, (0, 2): 29.63939, (0, 3): 59.31866, (0, 4): 98.26550,
         (1, 0): 4.2082827, (1, 1): 17.094282, (2, 0): 11.4894217, (3, 0): 21.376199, (4, 0): 33.914349}
# earlier frozen parameter-free values (thresholds_Lv/PREDICTION_J1v1.md): this code must reproduce them
OLD_FROZEN = {(1, 1): 17.1072, (3, 0): 21.3818, (4, 0): 33.9175}
FROLOV = {(1, 2): 39.375, (2, 1): 27.715, (2, 2): 53.295}
CUT = {}


def bo_ad_lam(L, v, M_ad):
    """BO+ad critical lam of root v, refined in a narrow bracket around the bo_map value."""
    lam0 = M_ad + 0.5
    r = C.critical(L, 1, lam0 - 0.05, lam0 + 0.05, n=5)
    assert len(r) == 1 and abs(r[0] - lam0) < 1e-6, (L, v, r, lam0)
    return float(r[0])


def grid_weights(lam, L):
    """As thresholds_Lv/predict.py (large-L guard: beyond the minimum of |u| for r > 30 the residual growing part
    of the numerical zero-energy solution is replaced by the decaying tail u(r_cut)(r_cut/r)^L)."""
    I, r, u, up = C.weights(lam, L, 1)
    if L > 0:
        sel = np.where(r > 30)[0]
        ic = sel[np.argmin(np.abs(u[sel]))]
        if ic < len(r) - 1:
            rc = r[ic]
            u = u.copy(); up = up.copy()
            u[ic:] = u[ic] * (rc / r[ic:]) ** L
            up[ic:] = -L * u[ic:] / r[ic:]
        CUT[(L, round(lam, 4))] = float(r[ic])
    B = np.trapz(-C.dU(r) * u ** 2, r)
    lB = lam * B
    return r, up ** 2 / lB, L * (L + 1) * u ** 2 / r ** 2 / lB


def predict_level(L, v, M_ad, mv, mr):
    lam = bo_ad_lam(L, v, M_ad)
    r, wr, wc = grid_weights(lam, L)
    s = dict(L=L, v=v, lam=lam, _grid=(r, wr, wc))
    p1 = A.fold_first_order(s, mv, mr)
    pnp = lam - A.lam_nonpert(s, mv, mr, p1)
    fcen = np.trapz(wc, r) / (np.trapz(wc, r) + np.trapz(wr, r))
    mvw = np.trapz(mv(r) * wr, r) / np.trapz(wr, r)
    mrw = np.trapz(mr(r) * wc, r) / np.trapz(wc, r) if L else float("nan")
    return dict(L=L, v=v, lam=lam, M_BOad=lam - 0.5, off_first=float(p1), off_nonpert=float(pnp), f_cen=float(fcen),
                mv_w=float(mvw), mr_w=float(mrw), M_pred_np=float(lam - 0.5 - pnp), M_pred_first=float(lam - 0.5 - p1))


def main():
    t0 = time.time()
    lines = []
    P = lambda s="": (print(s, flush=True), lines.append(s))
    P(f"# threshold_map/predict_map.py  {time.strftime('%Y-%m-%d %H:%M')}")
    df = pd.read_csv(os.path.join(ROOT, "paper_data", "bo_offset", "massfunc", "results_fine.csv")).drop_duplicates("R").sort_values("R")
    mv, mr = A.MassFunction(df.R.values, df.m_v.values), A.MassFunction(df.R.values, df.m_r.values)
    bm = json.load(open(os.path.join(HERE, "bo_map.json")))
    lev = {(d["J"], d["v"]): d for d in bm["levels"]}
    # ---- calibration (known exact levels)
    P("\n## calibration: levels with known exact M_c (not fitted; only the error band)")
    P("level   M_c(BO+ad)  off_1st  off_np   f_cen  <m_v>_w <m_r>_w  pred_np    exact      exact-pred")
    cal = {}
    for key in sorted(EXACT):
        q = predict_level(*key, lev[key]["M_BOad"], mv, mr)
        q["exact"] = EXACT[key]
        q["resid"] = EXACT[key] - q["M_pred_np"]
        cal[key] = q
        P(f"({key[0]},{key[1]})  {q['M_BOad']:10.6f}  {q['off_first']:7.4f}  {q['off_nonpert']:7.4f}  {q['f_cen']:5.3f}  "
          f"{q['mv_w']:6.4f}  {q['mr_w']:6.4f}  {q['M_pred_np']:9.4f}  {EXACT[key]:9.6f}  {q['resid']:+.4f}")
    for key, old in OLD_FROZEN.items():
        assert abs(cal[key]["M_pred_np"] - old) < 1.5e-4, (key, cal[key]["M_pred_np"], old)
    P("(reproduces the frozen thresholds_Lv values 17.1072, 21.3818, 33.9175 to < 1.5e-4: asserted)")
    unc = float(max(abs(q["resid"]) for q in cal.values()))
    kmax = max(cal, key=lambda k: abs(cal[k]["resid"]))
    rms = float(np.sqrt(np.mean([q["resid"] ** 2 for q in cal.values()])))
    P(f"uncertainty band (fixed now) = max |exact - pred| = {unc:.4f} (level {kmax}); rms of the 9 residuals {rms:.4f}")
    # ---- targets
    P("\n## PREDICTIONS (targets: BO+ad <= 101, no exact value)")
    P("target  sector         M_c(BO)     M_c(BO+ad)  off_np   f_cen  pred_np +/- unc       pred_1st   Frolov92")
    pred = {}
    for name in bm["targets"]:
        J, v = int(name[1:].split("v")[0]), int(name.split("v")[1])
        d = lev[(J, v)]
        q = predict_level(J, v, d["M_BOad"], mv, mr)
        q.update(name=name, sym=d["sym"], M_BO=d["M_BO"], unc=unc, frolov=FROLOV.get((J, v)))
        pred[name] = q
        fr = f"{FROLOV[(J, v)]:.3f}" if (J, v) in FROLOV else ""
        P(f"{name:6s}  L={J} '{d['sym']}' root {v}  {d['M_BO']:10.6f}  {q['M_BOad']:10.6f}  {q['off_nonpert']:7.4f}  "
          f"{q['f_cen']:5.3f}  {q['M_pred_np']:9.4f} +/- {unc:.4f}  {q['M_pred_first']:9.4f}  {fr}")
    P(f"\nweight cut radius r_cut (L, lam) -> r (RM = {C.RM:g} means no cut): {CUT}")
    P(f"run time {time.time() - t0:.0f} s")
    out = dict(date=time.strftime("%Y-%m-%d %H:%M"), unc=unc, unc_level=f"{kmax}", rms=rms,
               calibration={f"J{k[0]}v{k[1]}": v for k, v in cal.items()}, predictions=pred)
    json.dump(out, open(os.path.join(HERE, "prediction_map.json"), "w"), indent=1, default=float)
    open(os.path.join(HERE, "predict_map_output.txt"), "w").write("\n".join(lines) + "\n")
    if "--freeze" not in sys.argv:
        return
    md, sha = os.path.join(HERE, "PREDICTION_MAP.md"), os.path.join(HERE, "prediction_map.sha256")
    if os.path.exists(sha):
        print("prediction already frozen; not overwritten")
        return
    now = time.strftime("%Y-%m-%d %H:%M")
    cal_rows = "\n".join(f"| ({q['L']},{q['v']}) | {q['M_BOad']:.6f} | {q['off_first']:.4f} | {q['off_nonpert']:.4f} | "
                         f"{q['f_cen']:.3f} | {q['M_pred_np']:.4f} | {q['exact']} | {q['resid']:+.4f} |" for q in cal.values())
    pred_rows = "\n".join(f"| ({q['L']},{q['v']}) | L={q['L']} '{q['sym']}' root {q['v']} | {q['M_BO']:.6f} | {q['M_BOad']:.6f} | "
                          f"{q['off_nonpert']:.4f} | {q['f_cen']:.3f} | **{q['M_pred_np']:.4f} +/- {unc:.4f}** | "
                          f"{q['M_pred_first']:.4f} | {q['frolov'] if q['frolov'] else ''} |" for q in pred.values())
    txt = f"""# Blind predictions: exact critical masses of the (J, v) threshold map of {{X+, X+, e-}}

*Last updated: {now}*

**Frozen {now} (local time), before any exact three-body calculation of these targets was run.**
sha256 of this file: `prediction_map.sha256`; this file, `prediction_map.json`, `predict_map_output.txt` and the
sha256 file are chmod 444. Produced by `predict_map.py` (BO input `bo_map.json`). Must not be edited after the
exact runs start. Honesty notes on what was seen before the freeze: `PROGRESS.md`.

M_c = mass M (electron masses) of the two identical X+ at which the level becomes bound below
E_th = -M/(2(M+1)). Level (J, v) of the 1s sigma_g manifold = root k = v of the natural-parity L = J sector with
exchange symmetry (-1)^J (odd J 'A', even J 'S').

## Method (no fitted parameter)

BO+ad zero-energy coupling problem (lam = M + 1/2); threshold weights folded with the first-principles mass functions
m_v(R) (radial) and m_r(R) (centrifugal) from `paper_data/bo_offset/massfunc/results_fine.csv`, solved
non-perturbatively. Same code path as the two earlier blind tests (v = 4: passed; J1v1, J3v0, J4v0: passed).

## Calibration (nine levels with known exact M_c; used ONLY for the band)

| level | M_c(BO+ad) | offset 1st | offset non-pert | f_cen | predicted | exact | exact - pred |
|---|---|---|---|---|---|---|---|
{cal_rows}

Uncertainty band, fixed now and applied to every target: **+/- {unc:.4f}** = largest |exact - pred| of the nine
(level {kmax}); rms of the nine residuals {rms:.4f}.

## PREDICTIONS

| level | sector | BO | BO+ad | offset non-pert | f_cen | **parameter-free (non-pert.)** | first order | Frolov 1992 II |
|---|---|---|---|---|---|---|---|---|
{pred_rows}

**Pass/fail (set now).** Each target passes if |M_c(exact) - prediction| <= {unc:.4f} (+ the exact error bar).
Summary statistic, also set now: the number of passes out of all targets that reach a converged exact value.
"""
    open(md, "w").write(txt)
    h = hashlib.sha256(open(md, "rb").read()).hexdigest()
    open(sha, "w").write(f"{h}  PREDICTION_MAP.md\n")
    for f in (md, sha, os.path.join(HERE, "prediction_map.json"), os.path.join(HERE, "predict_map_output.txt")):
        os.chmod(f, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
    print("FROZEN", h)


if __name__ == "__main__":
    main()
