"""Blind prediction of the exact critical mass of the 5th symmetric S level (v = 4) of {X+, X+, e-}.

Written 2026-11-14.  Run from the project root (about 1-2 min):

    python3 paper_data/critical_masses_v4/predict.py [--freeze]

Steps (all parameter-free except the clearly labelled empirical variants):
 1. BO and BO+ad critical masses of v = 0..4 from the zero-energy coupling problem of
    paper_data/bo_offset/coupling.py (lam = M + 1/2), cross-checked against an independent
    bo_h2plus count_levels bisection (Numerov, ppw 16 and 32, Richardson) for n = 5.
 2. Zero-energy BO+ad threshold weights w_rad(r) = u'^2/(lam B) of each level.
 3. Fold with the first-principles mass functions m_v(R) of paper_data/bo_offset/massfunc
    (results_fine.csv), first order and non-perturbative (same code as massfunc/analyze.py):
        predicted M_c(exact) = M_c(BO+ad) - predicted non-adiabatic offset.
    The v = 1..3 fold is recomputed here and must reproduce massfunc/fold_results.csv.
 4. Empirical variants: (a) parameter-free prediction + extrapolated residual trend in 1/lam;
    (b) direct extrapolation of the exact BO+ad offsets of v = 1..3 in 1/lam.
 5. Uncertainty from the v = 1..3 residuals.
Outputs: prediction.json, predict_output.txt, PREDICTION.md, prediction.sha256.
PREDICTION.md (+ sha256) is written only with --freeze, and never overwritten once prediction.sha256
exists (blind test), unless --force.
"""
import hashlib
import json
import math
import os
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
import analyze as A   # noqa: E402  (massfunc/analyze.py: MassFunction, fold_first_order, lam_nonpert)

EXACT = {1: (9.525634, 0.000004), 2: (29.63939, 0.00001), 3: (59.31866, 0.00002)}   # paper_data/critical_masses
BRACKETS = [(0.5, 5.0), (5.0, 20.0), (20.0, 45.0), (45.0, 80.0), (85.0, 120.0)]


def crit_lams(c):
    out = []
    for lo, hi in BRACKETS:
        r = C.critical(0, c, lo, hi, n=40)
        out.append(r)
    flat = [x for r in out for x in r]
    return flat


def bo_h2plus_check(n=5):
    """Independent route: bo_h2plus integer level count, bisection, ppw 16 and 32, Richardson."""
    from bo_h2plus import count_levels
    res = {}
    for ad in (False, True):
        def crit(ppw, lo, hi, rtol=1e-11):
            cnt = lambda M: count_levels(M, 0, adiabatic=ad, ppw=ppw)
            assert cnt(lo) < n <= cnt(hi), (ad, ppw, lo, hi, cnt(lo), cnt(hi))
            while hi / lo - 1 > rtol:
                mid = math.sqrt(lo * hi)
                if cnt(mid) >= n:
                    hi = mid
                else:
                    lo = mid
            return math.sqrt(lo * hi)
        M1 = crit(16.0, 85.0, 115.0)
        M2 = crit(32.0, M1 * 0.99, M1 * 1.01)
        res["BOad" if ad else "BO"] = (M2 + (M2 - M1) / 15, abs(M2 - M1) / 15, M1, M2)
    return res


def grid_weights(lam, c):
    I, r, u, up = C.weights(lam, 0, c)
    lB = lam * I["B"]
    return I, r, up ** 2 / lB, np.zeros_like(r)


def main():
    force = "--force" in sys.argv
    t0 = time.time()
    lines = []
    P = lambda s="": (print(s, flush=True), lines.append(s))
    P(f"# predict.py  {time.strftime('%Y-%m-%d %H:%M')}")

    # ---- 1. BO and BO+ad critical lam
    lam_bo, lam_ad = crit_lams(0), crit_lams(1)
    P(f"BO    lam_c: {lam_bo}")
    P(f"BO+ad lam_c: {lam_ad}  (BO+ad has v=0 bound already at M=0.05, so its list starts at v=1)")
    assert len(lam_bo) == 5 and len(lam_ad) == 4, (lam_bo, lam_ad)
    ref = pd.read_csv(os.path.join(ROOT, "paper_data", "critical_masses", "bo_critical_masses.csv"), comment="#",
                      header=None, names=["model", "n", "M_c", "err", "a", "b"])
    P("consistency with paper_data/critical_masses/bo_critical_masses.csv (n <= 4):")
    for _, q in ref.iterrows():
        n = int(q.n)
        if q.model == "BO":
            mine = lam_bo[n - 1] - 0.5
        elif q.model.startswith("BO+ad(own"):
            mine = lam_ad[n - 2] - 0.5
        else:
            continue
        P(f"   {q.model:22s} n={n}  csv {float(q.M_c):.10f}  here {mine:.10f}  diff {mine - float(q.M_c):+.1e}")
    chk = bo_h2plus_check(5)
    M_bo4, M_ad4 = lam_bo[4] - 0.5, lam_ad[3] - 0.5
    P(f"v=4 (n=5): BO    M_c = {M_bo4:.10f} (shooting);  bo_h2plus count_levels {chk['BO'][0]:.10f} +/- {chk['BO'][1]:.1e}"
      f"  diff {M_bo4 - chk['BO'][0]:+.1e}")
    P(f"v=4 (n=5): BO+ad M_c = {M_ad4:.10f} (shooting);  bo_h2plus count_levels {chk['BOad'][0]:.10f} +/- {chk['BOad'][1]:.1e}"
      f"  diff {M_ad4 - chk['BOad'][0]:+.1e}")

    # ---- 1b. first-order adiabatic piece (BO solution), as paper_data/bo_offset/offset_table.py
    P("\n## BO - BO+ad: exact (both models solved) vs first order <1/4 - D>/<-dU> on the BO zero-energy solution")
    ad_rows = {}
    for v in range(1, 5):
        I, *_ = C.weights(lam_bo[v], 0, 0)
        first = -I["Dad"] / I["B"]
        ex = lam_bo[v] - lam_ad[v - 1]
        ad_rows[v] = (ex, first)
        P(f"v={v}: exact {ex:.5f}   first order {first:.5f}   ratio {first / ex:.4f}")

    # ---- 2+3. weights and fold with m_v(R)
    df = pd.read_csv(os.path.join(ROOT, "paper_data", "bo_offset", "massfunc", "results_fine.csv")).drop_duplicates("R").sort_values("R")
    mv, mr = A.MassFunction(df.R.values, df.m_v.values), A.MassFunction(df.R.values, df.m_r.values)
    P(f"\n## Fold with m_v(R) (results_fine.csv, {len(df)} R points)")
    P("v   lam_c(BO+ad)  M_c(BO+ad)    exact_offset  pred_1st  pred_nonpert  <m_v>_w   frac(r>8)  resid_1st  resid_np")
    fold = {}
    old = pd.read_csv(os.path.join(ROOT, "paper_data", "bo_offset", "massfunc", "fold_results.csv"))
    for v in range(1, 5):
        lam = lam_ad[v - 1]
        I, r, wr, wc = grid_weights(lam, 1)
        s = dict(L=0, v=v, lam=lam, _grid=(r, wr, wc))
        p1 = A.fold_first_order(s, mv, mr)
        pnp = lam - A.lam_nonpert(s, mv, mr, p1)
        mvw = np.trapz(mv(r) * wr, r) / np.trapz(wr, r)
        g = mv(r) * wr
        far = np.trapz(g[r > 8], r[r > 8]) / np.trapz(g, r)
        ex = EXACT[v][0] if v in EXACT else None
        off = (lam - 0.5 - ex) if ex else float("nan")
        fold[v] = dict(lam=lam, M_BOad=lam - 0.5, exact_offset=off, pred_first=p1, pred_nonpert=pnp, mv_w=mvw,
                       frac_r_gt_8=far, resid_first=off - p1, resid_nonpert=off - pnp)
        P(f"{v}   {lam:12.7f}  {lam - 0.5:12.7f}  {off:12.5f}  {p1:8.5f}  {pnp:10.5f}    {mvw:.5f}   {far:.4f}"
          f"     {off - p1:+.5f}   {off - pnp:+.5f}")
        if v <= 3:
            o = old[old.state == f"L=0 v={v}"].iloc[0]
            assert abs(o.pred_first - p1) < 2e-4 and abs(o.pred_nonpert - pnp) < 2e-4, (v, o.pred_first, p1, o.pred_nonpert, pnp)
    P("(v = 1..3 reproduce massfunc/fold_results.csv to < 2e-4: checked by assertion)")

    # ---- 4. predictions
    f4 = fold[4]
    pred = {}
    pred["pf_first"] = f4["M_BOad"] - f4["pred_first"]
    pred["pf_nonpert"] = f4["M_BOad"] - f4["pred_nonpert"]
    x = np.array([1 / fold[v]["lam"] for v in (1, 2, 3)])
    x4 = 1 / f4["lam"]
    emp = {}
    for kind in ("first", "nonpert"):
        res = np.array([fold[v][f"resid_{kind}"] for v in (1, 2, 3)])
        c23 = np.polyfit(x[1:], res[1:], 1)           # smooth regime v = 2, 3 (line through two points)
        c123 = np.polyfit(x, res, 1)                   # least squares, all three
        emp[kind] = dict(resid_lin23=float(np.polyval(c23, x4)), resid_lin123=float(np.polyval(c123, x4)),
                         resid_mean23=float(res[1:].mean()))
    offs = np.array([fold[v]["exact_offset"] for v in (1, 2, 3)])
    o_lin23 = float(np.polyval(np.polyfit(x[1:], offs[1:], 1), x4))
    o_quad = float(np.polyval(np.polyfit(x, offs, 2), x4))
    o_lin123 = float(np.polyval(np.polyfit(x, offs, 1), x4))
    # residual = exact offset - predicted offset, and M_exact = M_BOad - offset, so M = pf - residual
    pred["emp_resid_nonpert_lin23"] = pred["pf_nonpert"] - emp["nonpert"]["resid_lin23"]
    pred["emp_resid_first_lin23"] = pred["pf_first"] - emp["first"]["resid_lin23"]
    pred["emp_resid_nonpert_lin123"] = pred["pf_nonpert"] - emp["nonpert"]["resid_lin123"]
    pred["emp_offset_lin23"] = f4["M_BOad"] - o_lin23
    pred["emp_offset_quad123"] = f4["M_BOad"] - o_quad
    pred["emp_offset_lin123"] = f4["M_BOad"] - o_lin123
    rf = np.array([fold[v]["resid_first"] for v in (1, 2, 3)])
    rn = np.array([fold[v]["resid_nonpert"] for v in (1, 2, 3)])
    unc = dict(first_max=float(np.abs(rf).max()), first_rms=float(np.sqrt((rf ** 2).mean())),
               nonpert_max=float(np.abs(rn).max()), nonpert_rms=float(np.sqrt((rn ** 2).mean())),
               nonpert_mean23=float(rn[1:].mean()))
    P("\n## Predictions for M_c(exact, v = 4)")
    P(f"BO     M_c = {M_bo4:.7f}")
    P(f"BO+ad  M_c = {f4['M_BOad']:.7f}")
    P(f"parameter-free, first order    : {pred['pf_first']:.5f}  (offset {f4['pred_first']:.5f})")
    P(f"parameter-free, non-perturbative: {pred['pf_nonpert']:.5f}  (offset {f4['pred_nonpert']:.5f})")
    P(f"  uncertainty from v=1..3 residuals: first max {unc['first_max']:.4f} rms {unc['first_rms']:.4f}; "
      f"nonpert max {unc['nonpert_max']:.4f} rms {unc['nonpert_rms']:.4f}")
    P(f"empirical (a) pf-nonpert - residual line in 1/lam through v=2,3 ({emp['nonpert']['resid_lin23']:+.5f}): {pred['emp_resid_nonpert_lin23']:.5f}")
    P(f"              pf-first   - residual line in 1/lam through v=2,3 ({emp['first']['resid_lin23']:+.5f}): {pred['emp_resid_first_lin23']:.5f}")
    P(f"              pf-nonpert - least-squares line over v=1..3 ({emp['nonpert']['resid_lin123']:+.5f}): {pred['emp_resid_nonpert_lin123']:.5f}")
    P(f"empirical (b) exact offset extrapolated in 1/lam: line v=2,3 -> offset {o_lin23:.5f}, M_c {pred['emp_offset_lin23']:.5f}; "
      f"quadratic v=1..3 -> {o_quad:.5f}, {pred['emp_offset_quad123']:.5f}; LSQ line v=1..3 -> {o_lin123:.5f}, {pred['emp_offset_lin123']:.5f}")
    P(f"\nrun time {time.time() - t0:.0f} s")

    out = dict(date=time.strftime("%Y-%m-%d %H:%M"), M_BO=M_bo4, M_BOad=f4["M_BOad"], lam_BO=lam_bo, lam_BOad=lam_ad,
               bo_h2plus_check=chk, adiabatic=ad_rows, fold=fold, predictions=pred, empirical=emp,
               offset_extrap=dict(lin23=o_lin23, quad123=o_quad, lin123=o_lin123), uncertainty=unc)
    json.dump(out, open(os.path.join(HERE, "prediction.json"), "w"), indent=1, default=float)
    open(os.path.join(HERE, "predict_output.txt"), "w").write("\n".join(lines) + "\n")
    md, sha = os.path.join(HERE, "PREDICTION.md"), os.path.join(HERE, "prediction.sha256")
    if "--freeze" not in sys.argv:
        print("dry prediction (no --freeze): PREDICTION.md not written")
        return
    if os.path.exists(sha) and not force:
        print("prediction.sha256 exists: PREDICTION.md left untouched (blind test). Use --force only before the exact run.")
        return
    write_md(md, out)
    h = hashlib.sha256(open(md, "rb").read()).hexdigest()
    open(sha, "w").write(f"{h}  PREDICTION.md\n")
    print("wrote PREDICTION.md, sha256", h)


def write_md(path, o):
    p, f4, u, e = o["predictions"], o["fold"][4], o["uncertainty"], o["empirical"]
    fd = o["fold"]
    stamp = o["date"]
    rows = "\n".join(
        f"| {v} | {fd[v]['M_BOad']:.7f} | {fd[v]['exact_offset']:.5f} | {fd[v]['pred_first']:.5f} | {fd[v]['pred_nonpert']:.5f} "
        f"| {fd[v]['resid_first']:+.5f} | {fd[v]['resid_nonpert']:+.5f} |" for v in (1, 2, 3))
    adr = "\n".join(f"| {v} | {o['adiabatic'][v][0]:.5f} | {o['adiabatic'][v][1]:.5f} |" for v in (1, 2, 3, 4))
    txt = f"""# Blind prediction: exact critical mass of the v = 4 symmetric S level of {{X+, X+, e-}}

*Last updated: {stamp}*

**Frozen {stamp} (local time), before any exact v = 4 calculation was run.** The sha256 of this file is in
`prediction.sha256`; `analyze.py` checks it. This file must not be edited after the exact run starts.
Produced by `predict.py` (full output: `predict_output.txt`, numbers: `prediction.json`).

M_c = mass M (electron masses) of the two identical X+ at which the 5th symmetric S level (v = 4, root k = 4)
becomes bound below E_th = -M/(2(M+1)).

## Inputs (no new fit)

| level | BO M_c | BO+ad M_c (own asymptote) |
|---|---|---|
| v = 4 | **{o['M_BO']:.7f}** | **{o['M_BOad']:.7f}** |

* Zero-energy coupling problem of `paper_data/bo_offset/coupling.py` (lam = M + 1/2). Cross-check with
  `bo_h2plus` count_levels bisection (ppw 16/32, Richardson): BO {o['bo_h2plus_check']['BO'][0]:.7f},
  BO+ad {o['bo_h2plus_check']['BOad'][0]:.7f}. The n <= 4 values reproduce `paper_data/critical_masses/bo_critical_masses.csv`
  (see `predict_output.txt`).
* Mass function m_v(R): `paper_data/bo_offset/massfunc/results_fine.csv` (first-principles, converged to 3e-11).

## Adiabatic piece (BO - BO+ad)

| v | exact (both models solved) | first order |
|---|---|---|
{adr}

## Non-adiabatic offset M_c(BO+ad) - M_c(exact): parameter-free fold, and the v = 1..3 test

| v | M_c(BO+ad) | exact offset | fold, 1st order | fold, non-pert. | residual 1st | residual non-pert. |
|---|---|---|---|---|---|---|
{rows}
| **4** | {f4['M_BOad']:.7f} | ? | **{f4['pred_first']:.5f}** | **{f4['pred_nonpert']:.5f}** | | |

Residual = exact offset - predicted offset. <m_v> weighted for v = 4: {f4['mv_w']:.4f}; fraction of the fold from r > 8 bohr: {f4['frac_r_gt_8']:.3f}.

## PREDICTION

| variant | predicted M_c(exact, v = 4) |
|---|---|
| **parameter-free, non-perturbative fold (primary)** | **{p['pf_nonpert']:.4f} ± {u['nonpert_max']:.4f}** |
| parameter-free, first-order fold | {p['pf_first']:.4f} ± {u['first_max']:.4f} |
| empirical (a): primary − extrapolated residual (residual line in 1/lam through v = 2, 3 = {e['nonpert']['resid_lin23']:+.4f}) | {p['emp_resid_nonpert_lin23']:.4f} ± 0.005 |
| empirical (a'): primary − least-squares residual line over v = 1..3 ({e['nonpert']['resid_lin123']:+.4f}) | {p['emp_resid_nonpert_lin123']:.4f} |
| empirical (b): exact BO+ad offsets of v = 1..3 extrapolated in 1/lam (line through v = 2, 3) | {p['emp_offset_lin23']:.4f} |
| empirical (b'): same, quadratic through v = 1..3 | {p['emp_offset_quad123']:.4f} |

**Uncertainties.**
* Parameter-free: the quoted ± is the largest |residual| of v = 1..3 (non-perturbative {u['nonpert_max']:.4f},
  rms {u['nonpert_rms']:.4f}; first order {u['first_max']:.4f}, rms {u['first_rms']:.4f}). The v = 2, 3 residuals are
  both negative (mean {u['nonpert_mean23']:+.4f}), so if that trend holds the exact M_c lies *above* the parameter-free
  value by about 0.005; the parameter-free number does not include it.
* Empirical (a): ± 0.005 is a judgement (the other empirical variants (a'), (b), (b') all lie within it); it is not a statistical error. Empirical variants are fits to the v = 1..3 data and are NOT parameter-free.
* The BO and BO+ad masses themselves are numerically exact to ~1e-8 and do not contribute.

**Pass/fail criterion (set now).** The parameter-free theory passes if |M_c(exact) - {p['pf_nonpert']:.4f}| <= {u['nonpert_max']:.4f}.
The empirical (a) extrapolation passes if |M_c(exact) - {p['emp_resid_nonpert_lin23']:.4f}| <= 0.005.
"""
    open(path, "w").write(txt)


if __name__ == "__main__":
    main()
