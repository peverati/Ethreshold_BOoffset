"""Comparison with King, Longford & Cox 2013 and the location of the g(M) minimum.

Reads king2013_table3.csv, king2013_table2.csv, ours_table.csv, ours_lit.csv, ours_extra.csv,
fine_grid.csv, minimum.csv.  Writes comparison_table3.csv, comparison_table2.csv, gmin.txt,
fig_g_comparison.png.
"""
import csv
import json
import math
import os

import numpy as np
from scipy.optimize import least_squares, minimize_scalar

from common import HERE, read_table3, a3_of, M_of_a3, Eth_reduced, g_of

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def rcsv(name):
    p = os.path.join(HERE, name)
    if not os.path.exists(p):
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh))


def ours_by_key(name):
    return {r["key"]: r for r in rcsv(name)}


def fmt(x, n=12):
    return f"{x:.{n}f}"


# ------------------------------------------------------------------ Table III comparison
def table3(out_lines):
    ours = ours_by_key("ours_table.csv")
    rows = []
    for r in read_table3():
        key = f"{round(r['M'], 12):.12g}"
        if key not in ours:
            continue
        o = ours[key]
        mref = r["m3v"] if o["unit"] == "m3" else r["m1v"]
        Ec, err = float(o["E_ext"]), float(o["E_err"])
        E, Eerr = mref * Ec, mref * err
        Evar = mref * float(o["E_var"])
        if math.isinf(r["m3v"]):
            Eth = -0.5 * r["m1v"]
        else:
            Eth = -r["m1v"] * r["m3v"] / (2 * (r["m1v"] + r["m3v"]))
        g = (E - Eth) / Eth
        gerr = Eerr / abs(Eth)
        diff = E - r["E_their"]
        uc = 10.0 ** (-r["dec_certain"])                       # unit of the last certain digit
        uu = 10.0 ** (-(r["dec_certain"] + len(r["unc_digits"])))   # unit of the (last) uncertain digit
        # their certain digits are right if our value lies within one unit of the last certain digit
        within = abs(diff) <= uc + Eerr
        # stricter: truncation of our value reproduces their certain digits
        s_ours = f"{abs(E):.{r['dec_certain']+4}f}"[:-4]
        s_their = f"{abs(r['E_their']):.{r['dec_certain']+4}f}"[:-4]
        rows.append(dict(system=r["system"], M=r["M"], a3=a3_of(r["M"]) if r["M"] > 0 else 0.0,
                         E_their=r["E_printed"], E_ours=E, E_ours_err=Eerr, E_ours_var=Evar, diff=diff,
                         diff_in_last_certain=diff / uc, diff_in_uncertain_digit=diff / uu,
                         within_stated=within, certain_digits_match=(s_ours == s_their),
                         g_their=float(r["g_printed"]), g_ours=g, g_ours_err=gerr,
                         g_diff=g - float(r["g_printed"]), note=r["note"]))
    keys = list(rows[0].keys())
    with open(os.path.join(HERE, "comparison_table3.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for d in rows:
            w.writerow({k: (repr(v) if isinstance(v, float) else v) for k, v in d.items()})
    out_lines.append("Table III: ours - theirs (their energy units)")
    out_lines.append(f"{'system':15s} {'M':>12s} {'theirs':>22s} {'ours':>20s} {'err':>8s} {'diff':>9s} "
                     f"{'d/ulc':>7s} {'ok':>3s} {'trunc':>5s} {'g_ours':>13s} {'g_diff':>9s}")
    for d in rows:
        out_lines.append(f"{d['system']:15s} {d['M']:12.6g} {d['E_their']:>22s} {d['E_ours']:20.12f} "
                         f"{d['E_ours_err']:8.1e} {d['diff']:9.1e} {d['diff_in_last_certain']:7.2f} "
                         f"{'yes' if d['within_stated'] else 'NO':>3s} {'yes' if d['certain_digits_match'] else 'no':>5s} "
                         f"{d['g_ours']:13.10f} {d['g_diff']:9.1e}")
    n = len(rows)
    nw = sum(d["within_stated"] for d in rows)
    out_lines.append(f"within one unit of their last certain digit: {nw}/{n}")
    return rows


# ------------------------------------------------------------------ Table II comparison
T2MAP = {"H-": ("e", "p"), "Mu-": ("e", "mu"), "muH-": ("mu", "p"), "Ps-": ("e", "e"), "muD2+": ("d", "mu"),
         "Mu2+": ("mu", "e"), "H2+": ("p", "e"), "D2+": ("d", "e")}


def table2(out_lines):
    from common import mass_of
    ours = ours_by_key("ours_table.csv")
    lit = ours_by_key("ours_lit.csv")
    rows = []
    for r in rcsv("king2013_table2.csv"):
        m1, m3 = T2MAP[r["system"]]
        M = mass_of(m1) / mass_of(m3)
        o = ours.get(f"{round(M, 12):.12g}")
        if o is None:
            continue
        E, err = float(o["E_ext"]), float(o["E_err"])   # unit particle = e or mu, as in their table
        d = dict(system=r["system"], units=r["units"], M=M, E_ours=E, E_ours_err=err,
                 E_ours_var=float(o["E_var"]), E_their_2856=float(r["E_2856_2par"]),
                 diff_2856=E - float(r["E_2856_2par"]), E_their_1078=float(r["E_1078_2par"]),
                 diff_1078=E - float(r["E_1078_2par"]))
        litv = r["E_literature"]
        litE = float(litv) if not litv.startswith("0") else -float(litv)
        d["E_literature"] = litE
        lk = r["system"] + "@lit"
        if lk in lit:
            d["E_ours_at_lit_masses"] = float(lit[lk]["E_ext"])
            d["E_ours_at_lit_masses_err"] = float(lit[lk]["E_err"])
            d["diff_lit"] = d["E_ours_at_lit_masses"] - litE
        elif r["system"] == "Ps-":
            d["E_ours_at_lit_masses"], d["E_ours_at_lit_masses_err"] = E, err
            d["diff_lit"] = E - litE
        rows.append(d)
    keys = ["system", "units", "M", "E_ours", "E_ours_err", "E_ours_var", "E_their_2856", "diff_2856",
            "E_their_1078", "diff_1078", "E_literature", "E_ours_at_lit_masses", "E_ours_at_lit_masses_err",
            "diff_lit"]
    with open(os.path.join(HERE, "comparison_table2.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for d in rows:
            w.writerow({k: (repr(d[k]) if isinstance(d.get(k), float) else d.get(k, "")) for k in keys})
    out_lines.append("")
    out_lines.append("Table II (2856 terms) and their 'literature' column")
    out_lines.append(f"{'system':7s} {'ours (CODATA2012)':>20s} {'err':>8s} {'King 2856':>18s} {'diff':>9s} "
                     f"{'ours @lit masses':>20s} {'literature':>18s} {'diff':>9s}")
    for d in rows:
        ol = d.get("E_ours_at_lit_masses", float("nan"))
        out_lines.append(f"{d['system']:7s} {d['E_ours']:20.13f} {d['E_ours_err']:8.1e} {d['E_their_2856']:18.12f} "
                         f"{d['diff_2856']:9.1e} {ol:20.13f} {d['E_literature']:18.12f} {d.get('diff_lit', float('nan')):9.1e}")
    return rows


# ------------------------------------------------------------------ g(M) minimum
def eq8(p, a3):
    A, B, C = p
    return A - B * np.sqrt(1 - a3) + C * (1 - a3) ** 2


def eq9(c, a3):
    x = np.sqrt(1 - a3)
    return sum(ci * x ** i for i, ci in enumerate(c))


def fit_min(a3, g, form, nterms=11):
    if form == "eq8":
        p = least_squares(lambda p: eq8(p, a3) - g, [0.2, 0.24, 0.09]).x
        f = lambda a: eq8(p, a)  # noqa: E731
    else:
        X = np.vstack([np.sqrt(1 - a3) ** i for i in range(nterms)]).T
        p, *_ = np.linalg.lstsq(X, g, rcond=None)
        f = lambda a: eq9(p, a)  # noqa: E731
    res = minimize_scalar(f, bounds=(0.05, 0.4), method="bounded", options=dict(xatol=1e-10))
    ss = np.sum((f(a3) - g) ** 2)
    r2 = 1 - ss / np.sum((g - g.mean()) ** 2)
    return res.x, float(res.fun), p, r2, f


def gmin(out_lines, t3rows):
    L = out_lines
    L.append("")
    L.append("=" * 100)
    L.append("Minimum of g = (E0 - E_th)/E_th,  E_th = -M/(2(M+1)),  a3 = M/(M+2)")
    L.append("King et al.: a3 = 0.197 +- 0.001, M = 0.491 +- 0.003 (fit, Eq. 9 n = 5, plus alternative forms)")
    # (0) their printed Eq. 9 coefficients
    c_pr = [0.588641277, 0.326670402, 0.383300517, 2.04869070, 12.3665410, 46.8666865, 110.492237,
            162.933845, 146.384326, 73.3459246, 15.7455692]
    L.append(f"(0) their printed Eq. 9 coefficients: g(a3=1) = c0 = {eq9(c_pr, 1.0):.4f}, g(a3=0.2) = "
             f"{eq9(c_pr, 0.2):.1f}, g(a3=0) = {eq9(c_pr, 0.0):.1f} (data: 0.199, 0.045, 0.0555): the printed "
             f"coefficients do not reproduce their own data (signs and/or form misprinted).")
    # distinct finite-mass points
    uniq = {}
    for d in t3rows:
        uniq.setdefault(round(d["M"], 12), d)
    pts = sorted(uniq.values(), key=lambda d: d["M"])
    a3 = np.array([a3_of(d["M"]) for d in pts])
    g_o = np.array([d["g_ours"] for d in pts])
    g_t = np.array([d["g_their"] for d in pts])
    fin = np.array([d["M"] > 0 for d in pts])
    extra = rcsv("ours_extra.csv")
    a3x = np.array([a3_of(float(r["M"])) for r in extra])
    gx = np.array([g_of(float(r["E_ext"]), float(r["M"])) if float(r["M"]) >= 1 else
                   g_of(float(r["E_ext"]) * float(r["M"]), float(r["M"])) for r in extra]) if extra else np.array([])
    L.append(f"(1) fits in King's functional forms ({fin.sum()} distinct finite masses of Table III + infinite-mass H-"
             f"{', + their 7 extra a3 points 0.17..0.23 (ours)' if len(gx) else ''})")
    fits = {}
    for label, aa, gg in [("their g (printed), all Table III points", a3, g_t),
                          ("our g, Table III points", a3, g_o),
                          ("our g, Table III points + 7 extra", np.concatenate([a3, a3x]), np.concatenate([g_o, gx]))]:
        if len(aa) == len(a3) + 0 or len(gx):
            for form in ("eq8", "eq9"):
                am, gm, p, r2, f = fit_min(aa, gg, form)
                fits[(label, form)] = (am, gm, r2)
                L.append(f"    {form}  {label:45s}: a3_min = {am:.6f}  M_min = {M_of_a3(am):.6f}  g_min = {gm:.9f}  "
                         f"R^2 = {r2:.10f}")
                L.append("          params: " + ", ".join(f"{x:.6g}" for x in p))
            for nt in (7, 9, 13):   # sensitivity of Eq. 9 to its order (2n+1 terms)
                am, gm, p, r2, f = fit_min(aa, gg, "eq9", nt)
                L.append(f"    eq9 with {nt:2d} terms (n={(nt-1)//2}) {label:33s}: M_min = {M_of_a3(am):.6f}  R^2 = {r2:.10f}")
    # (2) fit-free, dense grid
    fg = rcsv("fine_grid.csv")
    res = {}
    if fg:
        fg = sorted(fg, key=lambda r: float(r["M"]))
        M = np.array([float(r["M"]) for r in fg])
        L.append(f"(2) fit-free: dense grid M = {M.min():.4f}..{M.max():.4f}, step {np.median(np.diff(M)):.4f}, "
                 f"{len(M)} points; scales re-optimized every 0.05 in M")
        for N in (28, 32):
            g = np.array([float(r[f"g_N{N}"]) for r in fg])
            dg = np.array([float(r[f"dgdM_N{N}"]) for r in fg])
            i = int(np.argmin(g))
            # parabola through the three lowest-neighbour points
            x, y = M[i - 1:i + 2], g[i - 1:i + 2]
            cpar = np.polyfit(x, y, 2)
            Mpar = -cpar[1] / (2 * cpar[0])
            # sign change of the exact HF derivative, linear interpolation
            j = np.where((dg[:-1] < 0) & (dg[1:] >= 0))[0]
            Mhf = [M[k] - dg[k] * (M[k + 1] - M[k]) / (dg[k + 1] - dg[k]) for k in j]
            res[N] = dict(Mgrid=M[i], gmin=g[i], Mpar=Mpar, Mhf=Mhf)
            L.append(f"    N={N}: grid minimum M = {M[i]:.4f} (g = {g[i]:.10f}; +- half step {0.5*np.median(np.diff(M)):.5f}); "
                     f"3-point parabola M = {Mpar:.6f}; dg/dM sign change(s) at M = "
                     + ", ".join(f"{m:.6f}" for m in Mhf))
        # panel-boundary continuity check
        jumps = []
        for k in range(1, len(fg)):
            if fg[k]["panel_center"] != fg[k - 1]["panel_center"]:
                jumps.append(abs(float(fg[k]["E_N32"]) - float(fg[k]["E_N28"])))
        L.append(f"    E(N=32) - E(N=28) on the grid: max {max(abs(float(r['E_N32']) - float(r['E_N28'])) for r in fg):.1e}")
        res["grid"] = (M, np.array([float(r["g_N32"]) for r in fg]))
    # (3) Hellmann-Feynman root
    mn = rcsv("minimum.csv")
    if mn:
        L.append("(3) fit-free, root of the exact (Hellmann-Feynman) dg/dM at fixed basis, Brent, xtol 1e-12")
        for r in mn:
            L.append(f"    alphas@M={float(r['M_alpha']):.2f} (N_opt {r['N_opt']}: {float(r['alpha_X']):.6f}, "
                     f"{float(r['alpha_e']):.6f})  N={r['N']:>2s} size={r['size']:>5s}: M_min = {float(r['M_min']):.10f}  "
                     f"g_min = {float(r['g_min']):.12f}  g'' = {float(r['d2g_dM2']):.5f}")
        best = [r for r in mn if r["N"] == "32"]
        Ms = np.array([float(r["M_min"]) for r in best])
        main = [r for r in mn if float(r["M_alpha"]) == 0.49 and r["N_opt"] == "12"]
        seq = [float(r["M_min"]) for r in sorted(main, key=lambda r: int(r["N"]))]
        L.append(f"    N=32 spread over the four scale choices: {Ms.min():.10f} .. {Ms.max():.10f} "
                 f"(range {Ms.max() - Ms.min():.1e});  sequence N=20..32 (alphas@0.49): "
                 + ", ".join(f"{m:.10f}" for m in seq))
        res["hf"] = dict(Ms=Ms, seq=seq, g=float(best[0]["g_min"]), d2=float(best[0]["d2g_dM2"]))
    return res, fits, (a3, g_o, g_t, fin)


# ------------------------------------------------------------------ figure
def figure(t3rows, res, fits_data):
    a3, g_o, g_t, fin = fits_data
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    Mo = np.array([d["M"] for d in t3rows])
    pos = Mo > 0
    ax[0].semilogx(Mo[pos], [d["g_ours"] for d, p in zip(t3rows, pos) if p], "o", mfc="none", ms=8, label="ours (exact)")
    ax[0].semilogx(Mo[pos], [d["g_their"] for d, p in zip(t3rows, pos) if p], "x", ms=6, label="King et al. 2013")
    g_inf = [d for d in t3rows if d["M"] == 0][0]
    ax[0].axhline(g_inf["g_ours"], color="0.6", ls=":", lw=1, label=r"$M\to0$ ($\infty$H$^-$)")
    ax[0].set_xlabel(r"$M = m_1/m_3$")
    ax[0].set_ylabel(r"$g = (E_0 - E_{th})/E_{th}$")
    ax[0].set_title("(a) relative excess binding")
    ax[0].legend(fontsize=8)
    if "grid" in res:
        M, g = res["grid"]
        ax[1].plot(M, g, "-", lw=1.2, label="ours, dense grid (N = 32)")
        sel = (Mo > 0.28) & (Mo < 0.82)
        ax[1].plot(Mo[sel], [d["g_their"] for d, s in zip(t3rows, sel) if s], "x", ms=8, label="King et al. Table III")
        if "hf" in res:
            Mm = float(np.median(res["hf"]["Ms"]))
            ax[1].plot([Mm], [res["hf"]["g"]], "v", color="C3", ms=9, label=f"our minimum M = {Mm:.5f}")
        ax[1].axvspan(0.488, 0.494, color="C1", alpha=0.25, label="King: 0.491 ± 0.003")
        ax[1].axvline(0.490402, color="0.4", ls="--", lw=1, label="Eq. 9 (n = 5) fit to our g: 0.4904")
        ax[1].set_xlabel(r"$M$")
        ax[1].set_ylabel(r"$g$")
        ax[1].set_title("(b) the minimum, fit-free")
        ax[1].legend(fontsize=8)
    dif = np.array([abs(d["diff"]) for d in t3rows])
    unc = np.array([10.0 ** -0 for d in t3rows])
    E = np.array([abs(d["E_ours"]) for d in t3rows])
    ulc = np.array([abs(d["diff"] / d["diff_in_last_certain"]) if d["diff_in_last_certain"] != 0 else np.nan
                    for d in t3rows])
    err = np.array([d["E_ours_err"] for d in t3rows])
    Mp = np.where(Mo > 0, Mo, 1e-5)
    ax[2].loglog(Mp, np.maximum(dif / E, 1e-17), "o", label=r"$|E_{ours}-E_{King}|/|E|$")
    ax[2].loglog(Mp, ulc / E, "_", ms=12, mew=2, label="their last certain digit / |E|")
    ax[2].loglog(Mp, np.maximum(err / E, 1e-17), ".", color="0.5", label="our error estimate / |E|")
    ax[2].set_xlabel(r"$M$ ($M=0$ plotted at $10^{-5}$)")
    ax[2].set_ylabel("relative energy difference")
    ax[2].set_title("(c) agreement with Table III")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_g_comparison.png"), dpi=150)


def main():
    L = []
    t3 = table3(L)
    table2(L)
    res, fits, fd = gmin(L, t3)
    figure(t3, res, fd)
    txt = "\n".join(L)
    with open(os.path.join(HERE, "results.txt"), "w") as fh:
        fh.write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
