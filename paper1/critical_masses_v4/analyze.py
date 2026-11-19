"""Analysis of run_exact.py output: exact M_c(v = 4) with error bar, blind comparison with PREDICTION.md.

Written 2026-11-15.  Run from the project root after run_exact.py:
    python3 paper1/critical_masses_v4/analyze.py [--dir DIR]

1. Verifies sha256(PREDICTION.md) against prediction.sha256 and that PREDICTION.md predates the first
   non-dry entry of run_log.txt.  A mismatch is reported loudly in RESULT.md (the comparison is then not blind).
2. observables_state4.csv = best (lowest) ladder member at each M; then the k = 1..3 fit script
   (solvers/threebody/examples/critical_mass_fit.py) writes fits_state4.txt and fig_state4.png.
3. M_c = order-3 fit of sqrt(E_b) = c1 x + c2 x^2 + c3 x^3 on x <= 0.1 (as for k = 1..3); error bar = full
   spread of the order-2 and order-3 fits on windows x <= 0.05 .. 0.3 and of the free power-law fits, plus a
   basis term from the convergence checks (energy change / slope dE_b/dM, at the check masses).
4. Upper bounds: best fixed-basis crossing (margin 0) and rounding-safe bound (margin 1e-11).
5. Comparison with every PREDICTION variant and the pass/fail criteria set in PREDICTION.md -> RESULT.md.
"""
import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "solvers", "threebody", "examples"))
import numpy as np  # noqa: E402

import critical_mass_fit as F  # noqa: E402

K = 4


def read(path):
    return F.read(path) if os.path.exists(path) else []


def blind_check(d):
    md, sha = os.path.join(HERE, "PREDICTION.md"), os.path.join(HERE, "prediction.sha256")
    out = []
    ok = True
    h = hashlib.sha256(open(md, "rb").read()).hexdigest()
    ref = open(sha).read().split()[0]
    if h == ref:
        out.append(f"sha256(PREDICTION.md) = {h} matches prediction.sha256: the prediction is unchanged.")
    else:
        ok = False
        out.append(f"**WARNING: sha256 MISMATCH** (file {h}, stored {ref}): PREDICTION.md was edited.")
    m = re.search(r"Frozen (\d{4}-\d\d-\d\d \d\d:\d\d)", open(md).read())
    frozen = m.group(1) if m else "?"
    first = None
    logp = os.path.join(d, "run_log.txt")
    if os.path.exists(logp):
        for line in open(logp):
            if "=== run_exact.py stages" in line:
                first = line[:16]
                break
    out.append(f"PREDICTION.md frozen {frozen}; first production entry of run_log.txt: {first}.")
    if first and frozen != "?" and first < frozen:
        ok = False
        out.append("**WARNING: the exact run started before the prediction was frozen.**")
    out.append("Short code tests of run_exact.py (2026-11-15, after the freeze) are described in FREEZE_NOTES.md.")
    return ok, out, md


def best_observables(d):
    rows = read(os.path.join(d, "observables_all.csv"))
    best = {}
    for r in rows:
        M = round(r["M"], 10)
        if M not in best or r["E"] < best[M]["E"]:
            best[M] = r
    out = []
    for M in sorted(best):
        r = best[M]
        mu = r["M"] * (r["M"] + 1) / (2 * r["M"] + 1)
        Eb = r["binding"]
        kap = np.sqrt(2 * mu * Eb) if Eb > 0 else float("nan")
        out.append(dict(state=K, M=r["M"], best_alpha_s=r["alpha_s"], E=r["E"], E_th=r["E_th"], binding=Eb, kappa=kap,
                        r01=r["r01"], r01_rms=np.sqrt(r["r01_2"]), rXe=r["rXe"],
                        one_over_2kappa=1 / (2 * kap) if kap > 0 else float("nan"), dE_dM=r["dE_dM"],
                        dEth_dM=-1 / (2 * (r["M"] + 1) ** 2)))
    with open(os.path.join(d, "observables_state4.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        for r in out:
            w.writerow({k: repr(v) if isinstance(v, float) else v for k, v in r.items()})
    return out


def fits(d):
    er = read(os.path.join(d, "energies_state4.csv"))
    cr = read(os.path.join(d, "crossings_state4.csv"))
    data = F.best_binding(er)
    ok = (data[:, 1] > 0) & (data[:, 2] < 1e-2 * data[:, 1])
    M, Eb = data[ok, 0], data[ok, 1]
    std = [r for r in cr if r.get("variant", "standard") == "standard" and np.isfinite(r["M_star"])]
    Mstar = min(r["M_star"] for r in std)
    res = []
    for xmax in (0.05, 0.1, 0.2, 0.3):
        sel = (M - Mstar) <= xmax + 1e-12
        for order in (2, 3):
            if sel.sum() < order + 2:
                continue
            try:
                p, e, r = F.fit_sqrt(M[sel], Eb[sel], order, Mstar)
                res.append(("sqrt", xmax, order, p[0], p[1:].tolist(), int(sel.sum())))
            except RuntimeError:
                pass
    for xmax in (0.05, 0.1, 0.3):
        sel = (M - Mstar) <= xmax + 1e-12
        if sel.sum() < 4:
            continue
        try:
            p, e = F.fit_power(M[sel], Eb[sel], Mstar)
            res.append(("power", xmax, p[2], p[0], [], int(sel.sum())))
        except (RuntimeError, ValueError):
            pass
    central = [r for r in res if r[0] == "sqrt" and r[2] == 3 and abs(r[1] - 0.1) < 1e-9]
    if not central:
        central = [r for r in res if r[0] == "sqrt" and r[2] == 2]
    Mc = central[0][3]
    allM = np.array([r[3] for r in res])
    c1 = central[0][4][0]
    return dict(Mc=Mc, lo=allM.min(), hi=allM.max(), fit_err=max(Mc - allM.min(), allM.max() - Mc), fits=res,
                Mstar=Mstar, c1=c1, cpar=central[0][4], n_points=int(ok.sum()), n_total=len(data),
                spread=data[:, 2].tolist(), crossings=cr)


def basis_error(d, fr):
    conv = read(os.path.join(d, "convergence_state4.csv"))
    out = []
    worst = 0.0
    for r in conv:
        if r["variant"].startswith("reference"):
            continue
        x = r["M"] - fr["Mc"]
        # dE_b/dM from the fit: E_b = (c1 x + c2 x^2 + ...)^2
        c = [fr["c1"]] + list(fr["cpar"][1:])
        s = sum(ci * x ** (i + 1) for i, ci in enumerate(c))
        ds = sum((i + 1) * ci * x ** i for i, ci in enumerate(c))
        slope = 2 * s * ds
        dM = abs(r["diff_from_ref"]) / slope if slope > 0 else float("nan")
        out.append((r["variant"], r["alpha_s"], r["M"], r["diff_from_ref"], dM))
        if not r["variant"].startswith("tol"):      # tol 1e-10 shrinks the space: not a convergence test
            worst = max(worst, dM if np.isfinite(dM) else 0.0)
    return worst, out



def crossing_line(d):
    """Margin-0 crossings, reported honestly. Near M_c the fixed-basis slope is only ~3e-11 Eh per mass
    unit, so the ~1e-14 Eh rounding floor moves a margin-0 crossing by ~3e-4; crossings whose E - E_th
    never exceeds ~1e-13 (the most diffuse ladder members) sit at that floor and are not bounds at all.
    (Fixed 2026-11-15: the shared fit script's 'best crossing' picked such a noise crossing, 98.2500.)"""
    cr = read(os.path.join(d, "crossings_state4.csv"))
    ok = [r for r in cr if np.isfinite(r["M_star"]) and r["alpha_s"] >= 1e-4]
    noise = [r for r in cr if r["alpha_s"] < 1e-4]
    lo, hi = min(r["M_star"] for r in ok), max(r["M_star"] for r in ok if r["alpha_s"] <= 1e-4)
    return (f"{lo:.5f} to {hi:.5f} for the alpha_s = 1e-4 bases (scatter ~3e-4 = rounding floor 1e-14 Eh / slope); "
            f"{len(noise)} more diffuse members (alpha_s < 1e-4) are at the rounding floor near M_c and give no bound")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=HERE)
    a = ap.parse_args()
    d = a.dir
    blind_ok, blind_lines, md = blind_check(d)
    best_observables(d)
    F.main(d, [K])                                     # fits_state4.txt, fig_state4.png (k = 1..3 script)
    fr = fits(d)
    berr, brows = basis_error(d, fr)
    err = fr["fit_err"] + berr
    ub = read(os.path.join(d, "upper_bounds.csv"))
    ub0 = [r for r in ub if r["margin"] == 0.0 and np.isfinite(r["M_bound"])]
    ub1 = [r for r in ub if abs(r["margin"] - 1e-11) < 1e-20 and np.isfinite(r["M_bound"])]
    pj = json.load(open(os.path.join(HERE, "prediction.json")))
    P = pj["predictions"]
    mdtxt = open(md).read()
    for k in ("pf_nonpert", "pf_first", "emp_resid_nonpert_lin23"):
        assert f"{P[k]:.4f}" in mdtxt, f"prediction.json value {k} not found in PREDICTION.md"
    unc = pj["uncertainty"]
    variants = [("parameter-free, non-perturbative (primary)", P["pf_nonpert"], unc["nonpert_max"]),
                ("parameter-free, first order", P["pf_first"], unc["first_max"]),
                ("empirical (a)", P["emp_resid_nonpert_lin23"], 0.005),
                ("empirical (a')", P["emp_resid_nonpert_lin123"], None),
                ("empirical (b)", P["emp_offset_lin23"], None),
                ("empirical (b')", P["emp_offset_quad123"], None),
                ("BO+ad (no non-adiabatic correction)", pj["M_BOad"], None),
                ("BO", pj["M_BO"], None)]
    stamp = time.strftime("%Y-%m-%d %H:%M")
    L = [f"# Result: exact critical mass of the v = 4 symmetric S level of {{X+, X+, e-}}", "",
         f"*Last updated: {stamp}*", "",
         "Written by `analyze.py` from the output of `run_exact.py`. Fits: `fits_state4.txt`, `fig_state4.png`.", "",
         "## Blind-test check", ""] + [f"* {x}" for x in blind_lines] + [""]
    if not blind_ok:
        L += ["**The comparison below is NOT blind.**", ""]
    L += ["## Exact M_c (v = 4)", "",
          f"**M_c = {fr['Mc']:.5f} ± {err:.5f}** (order-3 fit, x <= 0.1; {fr['n_points']} of {fr['n_total']} masses usable).", "",
          f"* fit spread (order 2/3, windows 0.05-0.3, power law): {fr['lo']:.6f} to {fr['hi']:.6f} (± {fr['fit_err']:.1e})",
          f"* basis term (largest |dE| of the convergence checks / slope dE_b/dM): {berr:.1e}",
          f"* threshold-law coefficients sqrt(E_b) = c1 x + c2 x^2 + c3 x^3: {np.array2string(np.array(fr['cpar']), precision=5)}",
          f"* fixed-basis crossings (margin 0, i.e. upper bounds only up to rounding): {crossing_line(d)}"]
    if ub1:
        L.append(f"* rounding-safe upper bound (margin 1e-11 Eh): {min(r['M_bound'] for r in ub1):.6f}")
    if ub0:
        sl = [r["slope_at_Mstar"] for r in ub0 if np.isfinite(r["slope_at_Mstar"])]
        if sl:
            L.append(f"* fixed-basis slope at M*: {np.min(np.abs(sl)):.1e} to {np.max(np.abs(sl)):.1e} Eh per mass unit")
    L += ["* The lower side is extrapolated through the threshold law (not rigorous), as for k = 1..3.", "",
          "Convergence checks (energy change vs the standard basis, and implied |dM|):", "",
          "| variant | alpha_s | M | dE (Eh) | implied dM |", "|---|---|---|---|---|"]
    for v, al, M, dE, dM in brows:
        L.append(f"| {v} | {al:g} | {M:.4f} | {dE:+.2e} | {dM:.1e} |")
    L += ["", "## Comparison with the blind prediction", "",
          "| prediction | value | stated ± | exact − predicted | within stated ±? |", "|---|---|---|---|---|"]
    for name, val, u in variants:
        dlt = fr["Mc"] - val
        ok = "" if u is None else ("**yes**" if abs(dlt) <= u + err else "**no**")
        L.append(f"| {name} | {val:.4f} | {'' if u is None else f'{u:.4f}'} | {dlt:+.4f} | {ok} |")
    off = pj["M_BOad"] - fr["Mc"]
    L += ["", f"Exact non-adiabatic offset M_c(BO+ad) − M_c(exact) = {off:.4f} (predicted {pj['fold']['4']['pred_nonpert']:.4f}); "
              f"residual exact − predicted offset = {off - pj['fold']['4']['pred_nonpert']:+.4f} "
              f"(v = 1..3: {', '.join(f'{pj['fold'][str(v)]['resid_nonpert']:+.4f}' for v in (1, 2, 3))}).",
          f"Offset vs BO: {pj['M_BO'] - fr['Mc']:.4f} (v = 1..3: 1.556, 1.528, 1.528).", "",
          "Pass criteria were fixed in PREDICTION.md; 'within' above also allows the exact error bar."]
    txt = "\n".join(L) + "\n"
    open(os.path.join(d, "RESULT.md"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
