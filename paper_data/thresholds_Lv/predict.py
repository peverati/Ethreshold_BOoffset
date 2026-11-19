"""Blind, parameter-free prediction of the exact critical masses of the 1s sigma_g-manifold levels
(J=1, v=1) [primary], (J=3, v=0) and (J=4, v=0) [secondary] of {X+, X+, e-}.

Written 2026-11-15.  Same method as paper_data/critical_masses_v4/predict.py (frozen v = 4 test, passed), now
with L > 0: BO+ad zero-energy coupling problem (paper_data/bo_offset/coupling.py, lam = M + 1/2), threshold
weights w_rad = u'^2/(lam B), w_cen = L(L+1) u^2/(r^2 lam B), folded with the first-principles
vibrational and rotational mass functions m_v(R), m_r(R) (paper_data/bo_offset/massfunc/results_fine.csv),
first order and non-perturbatively (massfunc/analyze.py: fold_first_order, lam_nonpert):
    predicted M_c(exact) = M_c(BO+ad) - predicted non-adiabatic offset.
Calibration (not fitted, only used for the error bar): the exact L=1 v=0 (4.2082827) and L=2 v=0
(11.4894217) critical masses (paper_data/rotational) and L=0 v=1..4 (paper_data/critical_masses, critical_masses_v4).

Run from the project root (~2-4 min):
    python3 paper_data/thresholds_Lv/predict.py [--freeze]
--freeze writes PREDICTION_J1v1.md + prediction_J1v1.sha256 (chmod 444) unless the sha256 file exists.
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

# exact critical masses known before this prediction (calibration of the error bar only)
EXACT = {(0, 1): 9.525634, (0, 2): 29.63939, (0, 3): 59.31866, (0, 4): 98.26550,
         (1, 0): 4.2082827, (2, 0): 11.4894217}
# (L, v): bracket in lam for the BO+ad root search, index of the root inside it
STATES = {(1, 0): (3.0, 8.0, 0), (2, 0): (9.0, 15.0, 0), (0, 1): (5.0, 20.0, 0),
          (1, 1): (8.0, 22.0, 1), (3, 0): (12.0, 45.0, 0), (4, 0): (20.0, 70.0, 0)}
# note: for (1, 1) the bracket [8, 22] contains the v=0 root? no: v=0 is at lam ~5.07 < 8, so the
# first root in [8, 22] is v=1; index fixed below after the search (asserted).


def bo_ad_lam(L, v):
    lo, hi, _ = STATES[(L, v)]
    r = C.critical(L, 1, lo, hi, n=60)
    return float(r[0]), r


def grid_weights(lam, L):
    """As massfunc/analyze.py export_weights, plus a guard for large L: at the critical lam the residual
    growing part A r^(L+1) (A ~ root tolerance) can dominate u ~ r^-L at r ~ 1e4 for L >= 4 and spoil the
    max-normalization. Beyond the minimum of |u| for r > 30 (r_cut; = RM when no growth), u is replaced by
    its decaying tail u(r_cut)(r_cut/r)^L. For L <= 3 r_cut = RM (no change; reported)."""
    I, r, u, up = C.weights(lam, L, 1)
    if L > 0:
        sel = np.where(r > 30)[0]
        ic = sel[np.argmin(np.abs(u[sel]))]
        if ic < len(r) - 1:
            rc = r[ic]
            u = u.copy(); up = up.copy()
            u[ic:] = u[ic] * (rc / r[ic:]) ** L
            up[ic:] = -L * u[ic:] / r[ic:]
            # (overall normalization cancels: the weights are u^2 / int(-dU u^2))
        CUT[L] = float(r[ic])
    B = np.trapz(-C.dU(r) * u ** 2, r)
    lB = lam * B
    return r, up ** 2 / lB, L * (L + 1) * u ** 2 / r ** 2 / lB


CUT = {}


def main():
    t0 = time.time()
    lines = []
    P = lambda s="": (print(s, flush=True), lines.append(s))
    P(f"# thresholds_Lv/predict.py  {time.strftime('%Y-%m-%d %H:%M')}")
    df = pd.read_csv(os.path.join(ROOT, "paper_data", "bo_offset", "massfunc", "results_fine.csv")).drop_duplicates("R").sort_values("R")
    mv, mr = A.MassFunction(df.R.values, df.m_v.values), A.MassFunction(df.R.values, df.m_r.values)
    bo = json.load(open(os.path.join(HERE, "bo_brackets.json")))
    res = {}
    P("state   lam_c(BO+ad)  M_c(BO+ad)  pred_1st  pred_np   f_cen   <m_v>_w  <m_r>_w   exact     resid_np")
    for (L, v) in [(0, 1), (1, 0), (2, 0), (1, 1), (3, 0), (4, 0)]:
        lam, allr = bo_ad_lam(L, v)
        if (L, v) == (1, 1):
            assert lam > 8.0 and abs(lam - 0.5 - bo["g L=1 BO+ad"][1]) < 1e-6, (lam, bo["g L=1 BO+ad"])
        r, wr, wc = grid_weights(lam, L)
        s = dict(L=L, v=v, lam=lam, _grid=(r, wr, wc))
        p1 = A.fold_first_order(s, mv, mr)
        pnp = lam - A.lam_nonpert(s, mv, mr, p1)
        fcen = np.trapz(wc, r) / (np.trapz(wc, r) + np.trapz(wr, r))
        mvw = np.trapz(mv(r) * wr, r) / np.trapz(wr, r)
        mrw = np.trapz(mr(r) * wc, r) / np.trapz(wc, r) if L else float("nan")
        ex = EXACT.get((L, v))
        resid = (lam - 0.5 - ex) - pnp if ex else float("nan")
        res[(L, v)] = dict(lam=lam, M_BOad=lam - 0.5, pred_first=p1, pred_nonpert=pnp, f_cen=fcen, mv_w=mvw,
                           mr_w=mrw, exact=ex, resid_nonpert=resid, resid_first=((lam - 0.5 - ex) - p1) if ex else float("nan"),
                           M_pred_np=lam - 0.5 - pnp, M_pred_first=lam - 0.5 - p1)
        P(f"L={L} v={v}  {lam:11.7f}  {lam - 0.5:10.6f}  {p1:8.5f}  {pnp:8.5f}  {fcen:6.3f}  {mvw:7.4f}  {mrw:7.4f}  "
          f"{ex if ex else float('nan'):9.6f}  {resid:+.5f}")
    # consistency with massfunc/analysis.txt (L=1: 0.3328, L=2: 0.2778, L=0 v=1: 0.5281 non-pert)
    for key, ref in (((1, 0), 0.3328), ((2, 0), 0.2778), ((0, 1), 0.5281)):
        assert abs(res[key]["pred_nonpert"] - ref) < 2e-4, (key, res[key]["pred_nonpert"], ref)
    P(f"weight cut radius r_cut by L (RM = {C.RM:g} means no cut): {CUT}")
    P("(L=1 v=0, L=2 v=0, L=0 v=1 reproduce massfunc/analysis.txt to < 2e-4: asserted)")
    # uncertainty: largest |residual| among the L > 0 calibration levels (non-perturbative)
    rl = [res[k]["resid_nonpert"] for k in ((1, 0), (2, 0))]
    unc = float(max(abs(x) for x in rl))
    # empirical: residual * lam roughly constant for L > 0 (L=1: lam*resid, L=2: lam*resid); mean
    kap = float(np.mean([res[k]["resid_nonpert"] * res[k]["lam"] for k in ((1, 0), (2, 0))]))
    P(f"\nL>0 calibration residuals (exact - predicted offset): L=1 v=0 {rl[0]:+.4f}, L=2 v=0 {rl[1]:+.4f}; "
      f"lam*resid = {res[(1, 0)]['resid_nonpert'] * res[(1, 0)]['lam']:.4f}, {res[(2, 0)]['resid_nonpert'] * res[(2, 0)]['lam']:.4f} "
      f"(mean {kap:.4f})")
    pred = {}
    for key, name in (((1, 1), "J1v1"), ((3, 0), "J3v0"), ((4, 0), "J4v0")):
        q = res[key]
        emp = q["M_pred_np"] - kap / q["lam"]
        pred[name] = dict(M_BO=bo[f"g L={key[0]} BO"][key[1]], M_BOad=q["M_BOad"], pf_nonpert=q["M_pred_np"],
                          pf_first=q["M_pred_first"], unc=unc, empirical=emp, emp_shift=-kap / q["lam"])
        P(f"{name}: BO {pred[name]['M_BO']:.6f}  BO+ad {q['M_BOad']:.6f}  parameter-free (non-pert) {q['M_pred_np']:.5f} "
          f"+/- {unc:.4f}  (first order {q['M_pred_first']:.5f});  empirical (resid ~ {kap:.3f}/lam) {emp:.5f}")
    P(f"\nrun time {time.time() - t0:.0f} s")
    out = dict(date=time.strftime("%Y-%m-%d %H:%M"), states={f"L{k[0]}v{k[1]}": v for k, v in res.items()},
               predictions=pred, kappa_emp=kap, unc=unc)
    json.dump(out, open(os.path.join(HERE, "prediction_J1v1.json"), "w"), indent=1, default=float)
    open(os.path.join(HERE, "predict_output.txt"), "w").write("\n".join(lines) + "\n")
    md, sha = os.path.join(HERE, "PREDICTION_J1v1.md"), os.path.join(HERE, "prediction_J1v1.sha256")
    if "--freeze" not in sys.argv:
        return
    if os.path.exists(sha):
        print("prediction already frozen; not overwritten")
        return
    p = pred["J1v1"]
    now = time.strftime("%Y-%m-%d %H:%M")
    rows = "\n".join(f"| L={L} v={v} | {q['M_BOad']:.6f} | {q['pred_first']:.5f} | {q['pred_nonpert']:.5f} | {q['f_cen']:.3f} | "
                     f"{q['mr_w']:.4f} | {q['exact'] if q['exact'] else '?'} | "
                     f"{q['resid_nonpert']:+.4f} |" if q['exact'] else
                     f"| L={L} v={v} | {q['M_BOad']:.6f} | {q['pred_first']:.5f} | {q['pred_nonpert']:.5f} | {q['f_cen']:.3f} | "
                     f"{q['mr_w']:.4f} | ? | |" for (L, v), q in res.items())
    txt = f"""# Blind prediction: exact critical mass of the (J = 1, v = 1) level of {{X+, X+, e-}}

*Last updated: {now}*

**Frozen {now} (local time), before any exact J = 1, v = 1 (or L = 3, 4) three-body calculation was run.**
sha256 of this file: `prediction_J1v1.sha256`; files chmod 444. Produced by `predict.py` (output
`predict_output.txt`, numbers `prediction_J1v1.json`). Must not be edited after the exact runs start.

M_c = mass M (electron masses) of the two identical X+ at which the level becomes bound below
E_th = -M/(2(M+1)). (J=1, v=1) = second root of the L = 1 exchange-antisymmetric natural-parity sector.

## Method (no fitted parameter)
BO+ad zero-energy coupling problem (lam = M + 1/2) for L = 1, second root; threshold weights folded with
m_v(R) (radial) and m_r(R) (centrifugal) from `paper_data/bo_offset/massfunc/results_fine.csv`; first order
and non-perturbative. Same code path as the v = 4 blind test (passed: predicted 98.2627 +/- 0.0095, exact
98.26550 +/- 0.00011).

| state | M_c(BO+ad) | offset 1st | offset non-pert | f_cen | <m_r>_w | exact M_c | resid (exact - pred offset) |
|---|---|---|---|---|---|---|---|
{rows}

## PREDICTION

| level | BO | BO+ad | **parameter-free (non-pert.)** | first order | empirical |
|---|---|---|---|---|---|
""" + "\n".join(f"| {n} | {q['M_BO']:.6f} | {q['M_BOad']:.6f} | **{q['pf_nonpert']:.4f} +/- {q['unc']:.4f}** | "
                f"{q['pf_first']:.4f} | {q['empirical']:.4f} |" for n, q in pred.items()) + f"""

* **Primary: M_c(J=1, v=1) = {p['pf_nonpert']:.4f} +/- {p['unc']:.4f}** (parameter-free). J = 3 and J = 4 rows are
  secondary blind predictions made by the same code at the same time.
* Uncertainty = largest |residual| of the two L > 0 calibration levels (L = 1 v = 0: {rl[0]:+.4f};
  L = 2 v = 0: {rl[1]:+.4f}). Both residuals are positive (exact offset larger than predicted), i.e. the
  exact M_c tends to lie BELOW the parameter-free value; 'empirical' subtracts a residual
  {kap:.4f}/lam (mean of lam*resid of the two calibration levels; a fit, NOT parameter-free).
* Literature cross-check (not used): Frolov & Wardlaw (arXiv:1008.3010), (J=1, v=1) of dd-mu-like system bound up
  to m_mu ~ 214.72 m_e, i.e. M_c ~ 17.09 (basis-limited upper estimate; docs/literature/referee_check_2026-11-15.md).
* Consistency (not used): ttmu (M = 26.585) has a bound J = 3 level, ddmu (17.752) does not: 17.752 < M_c(J=3) < 26.585.

**Pass/fail (set now).** Parameter-free passes if |M_c(exact) - {p['pf_nonpert']:.4f}| <= {p['unc']:.4f}.
"""
    open(md, "w").write(txt)
    h = hashlib.sha256(open(md, "rb").read()).hexdigest()
    open(sha, "w").write(f"{h}  PREDICTION_J1v1.md\n")
    for f in (md, sha, os.path.join(HERE, "prediction_J1v1.json"), os.path.join(HERE, "predict_output.txt")):
        os.chmod(f, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
    print("FROZEN", h)


if __name__ == "__main__":
    main()
