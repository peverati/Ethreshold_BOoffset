"""Blind, parameter-free predictions of the exact critical masses of {A+ (M), B+ (rM), e-} levels.

Written 2026-11-16.  Theory: THEORY.md.  All model critical couplings from twochannel.py (zero-energy problems,
lam = M_eff/mu_e, M = (1+r) lam/(2r) - 1/(1+r)).  No three-body result is used.

    python3 paper1/unequal_masses/predict.py            # compute -> predict_output.txt, bo_brackets.json
    python3 paper1/unequal_masses/predict.py --freeze   # + PREDICTION_UNEQUAL.md/.json, sha256, chmod 444
(~10-20 min with 6 worker processes, cached in model_roots.json; --recompute forces; --freeze refuses to overwrite an existing prediction_unequal.sha256.)

Targets: r = 2 (primary; HD+-like) and r = 3/2 (secondary), levels (J, v) = (0, 1..4), (1, 0), (2, 0).
Prediction = model '2ch+m' (two-channel g/u + 2x2 non-adiabatic mass matrices, non-perturbative).
Uncertainty = |2ch+m_asym - 2ch+m| + |2ch+m_rc3 - 2ch+m| + (symmetric-family residual of the same (L, v),
times (1+r)/(2r)), floor 0.03.
Null models (for comparison, also frozen):
    NS  symmetric transfer:            lam = lam_BOad(r) - dlam_sym(L, v)        (exact symmetric BO+ad offset)
    NT  symmetric transfer, true thr.: lam = lam_BOad_thr(r) - dlam_sym(L, v)
    NP  persistence (same family):     offset(v) = offset(v-1)  -> evaluated by analyze.py with the exact values;
        the theory's predicted level-to-level offset changes are listed here.
"""
import hashlib
import json
import os
import stat
import sys
import time
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import warnings  # noqa: E402

warnings.filterwarnings("ignore")
import numpy as np  # noqa: E402

import twochannel as T  # noqa: E402

R_LIST = (2.0, 1.5)
LEVELS = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (2, 0)]
KINDS = ["BO", "BOad", "BOad_thr", "BOad+m", "2ch", "2ch+m", "2ch+m_asym", "2ch+m_rc3"]
# lam scan per L: (lo, hi, step); the lowest root found is v = 0 (L = 0: BO v=0 lam_c = 1.465, so lo < 1)
SCAN = {0: (0.3, 165.0, 1.25), 1: (2.0, 45.0, 0.75), 2: (6.0, 55.0, 0.75)}
# exact symmetric-family BO+ad offsets (lam units = M units at r = 1) and residuals of the parameter-free fold
# (exact - predicted offset; paper1/critical_masses, critical_masses_v4, rotational, massfunc/analysis.txt,
#  thresholds_Lv/predict_output.txt).
SYM_OFFSET = {(0, 1): 0.5375, (0, 2): 0.5062, (0, 3): 0.4970, (0, 4): 0.4930, (1, 0): 0.3623, (2, 0): 0.2892}
SYM_RESID = {(0, 1): 0.0094, (0, 2): -0.0062, (0, 3): -0.0047, (0, 4): -0.0028, (1, 0): 0.0295, (2, 0): 0.0114}
WORKERS = int(os.environ.get("PREDICT_WORKERS", "6"))


def job(args):
    kind, r, L = args
    lo, hi, st = SCAN[L]
    if kind.startswith("2ch") and L == 0:
        lo = 1.0
    n = int(round((hi - lo) / st)) + 1
    t0 = time.time()
    try:
        roots = T.Model(kind, r, L).roots(lo, hi, n)[0]
        err = ""
    except Exception as exc:  # pragma: no cover
        roots, err = [], repr(exc)
    return kind, r, L, roots, time.time() - t0, err


def compute():
    """Model roots for every (kind, r, L); cached in model_roots.json (delete it or pass --recompute to redo)."""
    cache = os.path.join(HERE, "model_roots.json")
    if os.path.exists(cache) and "--recompute" not in sys.argv:
        d = json.load(open(cache))
        res = {(k, float(r), int(L)): v for key, v in d.items() for k, r, L in [key.split("|")]}
        if all((k, r, L) in res for r in R_LIST for L in SCAN for k in KINDS):
            print(f"  model roots read from {cache}")
            return res
    jobs = [(k, r, L) for r in R_LIST for L in SCAN for k in KINDS]
    res = {}
    with ProcessPoolExecutor(WORKERS) as ex:
        for kind, r, L, roots, dt, err in ex.map(job, jobs):
            res[(kind, r, L)] = roots
            print(f"  {kind:11s} r={r} L={L}: {len(roots)} roots {[round(x, 6) for x in roots]} ({dt:.0f}s) {err}", flush=True)
    json.dump({f"{k}|{r}|{L}": v for (k, r, L), v in res.items()}, open(cache, "w"), indent=1)
    return res


def level_lam(res, kind, r, L, v):
    """Root of level v: v = number of BOad_thr roots strictly below it (BOad_thr is the most shifted model;
    every model's v-th root lies in (lam_thr(v-1), lam_thr(v)], checked in predict_output.txt)."""
    thr = np.array(res[("BOad_thr", r, L)])
    cand = [x for x in res[(kind, r, L)] if int(np.sum(thr < x - 1e-9)) == v]
    if len(cand) != 1:
        return float("nan")
    return cand[0]


def main():
    t0 = time.time()
    lines = []
    P = lambda s="": (print(s, flush=True), lines.append(s))
    P(f"# unequal_masses/predict.py  {time.strftime('%Y-%m-%d %H:%M')}")
    res = compute()
    # sanity: BO lam_c independent of r and equal to the symmetric values (coupling.py)
    for r in R_LIST:
        assert abs(res[("BO", r, 0)][1] - 11.5814565935) < 1e-6, res[("BO", r, 0)]
        for k in KINDS:      # level assignment: every model has exactly one root per target level
            for (L, v) in LEVELS:
                assert np.isfinite(level_lam(res, k, r, L, v)), (k, r, L, v, res[(k, r, L)])
    out = dict(date=time.strftime("%Y-%m-%d %H:%M"), r_list=R_LIST, kinds=KINDS, levels=LEVELS,
               lam={f"{k}|{r}|{L}": v for (k, r, L), v in res.items()}, targets={})
    P("\nlam_c and M_c (light-nucleus mass) per model; v counted from the lowest root of each (L) sector")
    for r in R_LIST:
        conv = lambda lam: T.UQ.M_of_lam(lam, r)
        f = (1 + r) / (2 * r)
        P(f"\n## r = {r}  (beta = {T.UQ.beta_of(r):.6f}, M = {f:.4f} lam - {1 / (1 + r):.4f})")
        P("level  " + "  ".join(f"{k:>11s}" for k in KINDS))
        for (L, v) in LEVELS:
            P(f"J{L}v{v}   " + "  ".join(f"{conv(level_lam(res, k, r, L, v)):11.5f}" for k in KINDS))
        P("\nlevel   M_pred(2ch+m)   unc     M_BO      M_BOad    M_BOad_thr  null NS    null NT   "
          "off_BOad  off_thr  d(off_BOad) vs prev v")
        prev = None
        for (L, v) in LEVELS:
            lam = {k: level_lam(res, k, r, L, v) for k in KINDS}
            M = {k: conv(x) for k, x in lam.items()}
            spread = abs(M["2ch+m_asym"] - M["2ch+m"]) + abs(M["2ch+m_rc3"] - M["2ch+m"])
            unc = max(0.03, spread + abs(SYM_RESID[(L, v)]) * f)
            NS = conv(lam["BOad"] - SYM_OFFSET[(L, v)])
            NT = conv(lam["BOad_thr"] - SYM_OFFSET[(L, v)])
            off_ad = M["BOad"] - M["2ch+m"]
            off_thr = M["BOad_thr"] - M["2ch+m"]
            d = (off_ad - prev) if (prev is not None and L == 0 and v > 1) else float("nan")
            prev = off_ad if L == 0 else prev
            name = f"r{r:g}_J{L}v{v}"
            out["targets"][name] = dict(r=r, L=L, v=v, M_pred=M["2ch+m"], unc=unc, spread=spread,
                                        M=M, lam=lam, null_NS=NS, null_NT=NT, offset_BOad=off_ad,
                                        offset_BOad_thr=off_thr, d_offset_prev=d,
                                        mass_shift=M["2ch"] - M["2ch+m"], mixing_shift=M["BOad_thr"] - M["2ch"])
            P(f"J{L}v{v}   {M['2ch+m']:12.5f}  {unc:6.4f}  {M['BO']:9.5f} {M['BOad']:9.5f} {M['BOad_thr']:10.5f}  "
              f"{NS:9.5f} {NT:9.5f}  {off_ad:8.4f} {off_thr:8.4f}  {d:+.4f}")
    P(f"\nrun time {time.time() - t0:.0f} s")
    json.dump(out, open(os.path.join(HERE, "bo_brackets.json"), "w"), indent=1, default=float)
    open(os.path.join(HERE, "predict_output.txt"), "w").write("\n".join(lines) + "\n")
    if "--freeze" in sys.argv:
        freeze(out)


def freeze(out):
    md = os.path.join(HERE, "PREDICTION_UNEQUAL.md")
    sha = os.path.join(HERE, "prediction_unequal.sha256")
    js = os.path.join(HERE, "prediction_unequal.json")
    if os.path.exists(sha):
        print("prediction already frozen; not overwritten")
        return
    now = time.strftime("%Y-%m-%d %H:%M")
    json.dump(out, open(js, "w"), indent=1, default=float)
    rows = []
    for name, q in out["targets"].items():
        M = q["M"]
        rows.append(f"| {name} | {M['BO']:.5f} | {M['BOad']:.5f} | {M['BOad_thr']:.5f} | {M['2ch']:.5f} | "
                    f"**{q['M_pred']:.4f} ± {q['unc']:.4f}** | {q['null_NS']:.4f} | {q['null_NT']:.4f} | "
                    f"{q['offset_BOad']:+.4f} | {q['d_offset_prev']:+.4f} |")
    txt = f"""# Blind prediction: exact critical masses of {{A+ (M), B+ (rM), e-}}

*Last updated: {now}*

**Frozen {now} (local time), before any exact crossing search in this family.** sha256 of this file in
`prediction_unequal.sha256`; files chmod 444. Numbers: `prediction_unequal.json`; full model table:
`predict_output.txt`. Produced by `predict.py --freeze`. Exact three-body results seen before the freeze:
see PROGRESS.md (only HD+ validation energies at physical masses and an r = 1 energy comparison; no crossing).

M_c = mass M of the LIGHT nucleus (electron masses; heavy nucleus rM) at which the level becomes bound below
the exact threshold E_th = -mu_Be/2 = -rM/(2(rM+1)) (A+ + (B e)(1s)). Level (J, v) = root v of the
natural-parity L = J sector (no exchange symmetry).

## Method (no fitted parameter; THEORY.md)
Two-channel (1s sigma_g, 2p sigma_u) zero-energy problem with exact g/u couplings (a, K) and the exact
threshold built in, plus second-order non-adiabatic 2x2 mass matrices (Sternheimer, massfunc_gu.csv),
solved non-perturbatively ('2ch+m'). Uncertainty = model-variant spread (u-row masses atomic / frozen
below r = 3) + the symmetric-family residual of the same (L, v) scaled by (1+r)/(2r); floor 0.03.

## PREDICTION (M_c of the light nucleus)

| target | BO | BO+ad (own asymptote) | BO+ad (exact threshold) | 2ch (no masses) | **theory 2ch+m** | null NS | null NT | offset vs BO+ad | change vs previous v |
|---|---|---|---|---|---|---|---|---|---|
""" + "\n".join(rows) + f"""

* Null NS: BO+ad (own asymptote) minus the exact symmetric-family offset of the same (L, v) (lam units).
  Null NT: the same applied to BO+ad measured from the exact threshold. Null NP (persistence): offset(v) =
  offset(v-1) of this family, evaluated after the exact runs; the theory predicts the level-to-level change
  listed in the last column (offset = M_BOad - M_c).
* **Pass/fail (set now).** A target passes if |M_c(exact) - M_pred| <= unc. The test discriminates where
  |M_pred - null| > unc + (exact error bar): for NS and NT this holds for every level (columns above); for NP
  wherever the predicted change exceeds 0.03 + unc.
* **Secondary (lenient) criterion, set now:** |M_c(exact) - M_pred| <= 0.10. Reason: the floor 0.03 is
  calibrated on the symmetric family only; the O(1/lam) terms of the g/u part (2p pi_u, higher Sigma states as
  explicit channels, potential-like second-order terms near the localization radius) are not calibrated. Even
  0.10 is far smaller than the theory-null differences (>= 0.4 for NT, >= 2.4 for NS, level-to-level offset
  changes 1.7-2.3 for NP).
* Offsets are M_BOad - M_c (negative here: the exact threshold lies below the BO+ad asymptote).
* Primary targets: r = 2, J0v1..J0v4, J1v0, J2v0. Secondary: r = 1.5 (exact run only if time permits).
"""
    open(md, "w").write(txt)
    h = hashlib.sha256(open(md, "rb").read()).hexdigest()
    open(sha, "w").write(f"{h}  PREDICTION_UNEQUAL.md\n")
    hj = hashlib.sha256(open(js, "rb").read()).hexdigest()
    open(sha, "a").write(f"{hj}  prediction_unequal.json\n")
    for f in (md, sha, js, os.path.join(HERE, "predict_output.txt")):
        os.chmod(f, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
    print("FROZEN", h)


if __name__ == "__main__":
    main()
