"""Fits, tables and figures from the CSVs written by scan.py.

python analyze.py   ->  energies.csv, descriptors.csv, spacings.csv, fits.txt, density_crossover.txt,
                        fig1_energies_bo_error.png, fig2_lindemann_R.png, fig3_gR.png,
                        fig4_density_cm.png, fig5_mass_derivatives.png, fig6_spacings.png
"""
import csv
import json
import math
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
HARTREE_CM = 219474.6313705
RE_BO = 1.99719332            # bo_h2plus R_e for electron mass 1 (README)

plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "legend.fontsize": 8, "figure.dpi": 150,
                     "axes.grid": True, "grid.alpha": 0.25, "lines.linewidth": 1.4, "lines.markersize": 4})
ST = dict(exact=dict(color="k", marker="o", ls="-", label="exact (non-BO)"),
          BO=dict(color="tab:blue", marker="s", ls="--", label="BO"),
          BOad=dict(color="tab:orange", marker="^", ls="-.", label="BO + adiabatic"))


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def load(name):
    path = os.path.join(HERE, name)
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return [{k: (fnum(v) if k not in ("E0_flag", "E1_flag", "E1_opt_flag", "E_flag", "sym", "weights", "grid", "model")
                     and not k.endswith("_seq") else v) for k, v in r.items()} for r in csv.DictReader(fh)]


def col(rows, k):
    return np.array([r.get(k, float("nan")) for r in rows], dtype=float)


def write_csv(path, rows, fields):
    with open(os.path.join(HERE, path), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{r[k]:.15g}" if isinstance(r.get(k), float) else r.get(k, "")) for k in fields})


def powfit(x, y, sel):
    """ln|y| = a + p ln x on the selected points; returns p, sigma_p (stat.), A."""
    X, Y = np.log(x[sel]), np.log(np.abs(y[sel]))
    n = len(X)
    A = np.vstack([np.ones(n), X]).T
    c, res, *_ = np.linalg.lstsq(A, Y, rcond=None)
    if n > 2:
        s2 = np.sum((Y - A @ c) ** 2) / (n - 2)
        cov = s2 * np.linalg.inv(A.T @ A)
        sp = math.sqrt(cov[1, 1])
    else:
        sp = float("nan")
    return c[1], sp, math.exp(c[0])


def fit_report(name, kap, y, Ms, windows, out, expected=None):
    """Free power-law fits in kappa on several large-M windows; error = max(stat, window spread)."""
    lines = [f"--- {name}"]
    ps = []
    for Mmin in windows:
        sel = (Ms >= Mmin) & np.isfinite(y) & (y != 0)
        if sel.sum() < 2:
            continue
        p, sp, A = powfit(kap, y, sel)
        ps.append((Mmin, p, sp, A, int(sel.sum())))
        lines.append(f"  window M >= {Mmin:g} ({sel.sum()} pts): exponent in kappa p = {p:.4f} +/- {sp:.4f} (stat), "
                     f"prefactor {A:.5g}")
    # local slopes between consecutive masses
    sel = np.isfinite(y) & (y != 0)
    k2, y2, M2 = kap[sel], y[sel], Ms[sel]
    loc = [(M2[i], M2[i + 1], math.log(abs(y2[i + 1] / y2[i])) / math.log(k2[i + 1] / k2[i]))
           for i in range(len(M2) - 1) if M2[i] >= 20]
    lines.append("  local slopes d ln|y| / d ln kappa: " + ", ".join(f"[{a:g},{b:g}] {s:.3f}" for a, b, s in loc))
    if ps:
        main = ps[min(1, len(ps) - 1)]
        spread = max(abs(q[1] - main[1]) for q in ps)
        err = max(main[2] if np.isfinite(main[2]) else 0, spread)
        lines.append(f"  RESULT p = {main[1]:.3f} +/- {err:.3f} (window M >= {main[0]:g}; error = max(stat, window spread))"
                     + (f"; BO-theory expectation {expected}" if expected else ""))
        out[name] = (main[1], err, main[0])
    return "\n".join(lines)


LARGE_M = os.path.join(HERE, "..", "convergence_large_M")


def apply_large_M_override(ex, rot):
    """Replace the M = 1e4 extrapolated energies by the longer sequences of paper1/convergence_large_M
    (added 2026-11-14): old mass-scan sequence + new N, same fixed alphas, re-extrapolated."""
    import json
    import common  # noqa: F401  (puts solvers/threebody on sys.path)
    from threebody import extrapolate
    path = os.path.join(LARGE_M, "results.csv")
    if not os.path.exists(path):
        return
    new = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            new.setdefault(int(r["L"]), []).append((int(r["N"]), float(r["E"])))
    for r in ex:
        if round(r["M"]) == 10000 and 0 in new:
            seq = json.loads(r["E0_seq"]) + [e for _, e in sorted(new[0])]
            r["E0_ext"], r["E0_err"] = extrapolate(seq)[:2]
            r["E0"] = seq[-1]
    for r in rot:
        L = int(r["L"])
        if round(r["M"]) == 10000 and L in new:
            seq = json.loads(r["E_seq"]) + [e for _, e in sorted(new[L])]
            r["E_ext"], r["E_err"] = extrapolate(seq)[:2]
            r["E"] = seq[-1]


def main():
    bo = {round(r["M"], 6): r for r in load("bo.csv")}
    ex = sorted(load("exact_L0.csv"), key=lambda r: r["M"])
    rot = load("exact_rot.csv")
    apply_large_M_override(ex, rot)
    rotd = {(round(r["M"], 6), int(r["L"])): r for r in rot}
    if not ex:
        raise SystemExit("no exact data yet")
    Ms = col(ex, "M")
    kap = Ms ** -0.25
    B = lambda k: np.array([bo.get(round(m, 6), {}).get(k, float("nan")) for m in Ms], dtype=float)

    # ------------------------------------------------------------------ energies table
    E0 = col(ex, "E0_ext")
    E0err = col(ex, "E0_err")
    E0var = col(ex, "E0")
    Ebo, Ebad = B("BO_E00"), B("BOad_E00")
    dBO, dBOad = E0 - Ebo, E0 - Ebad
    Eth = col(ex, "E_th")
    eps = E0 / (Ms / (Ms + 1))
    rows = []
    for i, M in enumerate(Ms):
        rows.append(dict(M=M, kappa=kap[i], E_th=Eth[i], E0_variational=E0var[i], E0_extrap=E0[i], E0_err=E0err[i],
                         E0_flag=ex[i]["E0_flag"], E_BO=Ebo[i], E_BOad=Ebad[i], exact_minus_BO=dBO[i],
                         exact_minus_BOad=dBOad[i], E_over_muXe=eps[i], binding=Eth[i] - E0[i],
                         N_last=ex[i]["N_last"], size=ex[i]["size_last"], alpha_X=ex[i]["alpha_X"],
                         alpha_e=ex[i]["alpha_e"], weights=ex[i]["weights"], E0_seq=ex[i]["E0_seq"]))
    write_csv("energies.csv", rows, list(rows[0].keys()))

    fits = {}
    rep = ["Fits of the exact-minus-BO differences and descriptors at large M (kappa = M^-1/4).",
           "Exponents are free fits of ln|y| vs ln kappa; error = max(statistical, spread over windows M>=500/1000/1836).",
           ""]
    W = [500, 1000, 1836]
    rep.append(fit_report("E_exact - E_BO", kap, dBO, Ms, W, fits, "4 (adiabatic correction, O(1/M))"))
    rep.append(fit_report("E_exact - E_BOad", kap, dBOad, Ms, W, fits,
                          "6 (leading non-adiabatic term, O(M^-3/2): R-dependent vibrational mass)"))
    # fixed-exponent two-term fits
    sel = Ms >= 500
    A = np.vstack([kap[sel] ** 4, kap[sel] ** 6]).T
    c, *_ = np.linalg.lstsq(A, dBO[sel], rcond=None)
    rep.append(f"  two-term fit E-E_BO = a k^4 + b k^6 (M>=500): a = {c[0]:.6f}, b = {c[1]:.5f}  "
               f"(compare mu_e^2 D(R_e)/M -> D(2.0) = 0.1748 for the adiabatic correction at R ~ R_e)")
    A = np.vstack([kap[sel] ** 6, kap[sel] ** 8]).T
    c2, *_ = np.linalg.lstsq(A, dBOad[sel], rcond=None)
    rep.append(f"  two-term fit E-E_BOad = a k^6 + b k^8 (M>=500): a = {c2[0]:.6f}, b = {c2[1]:.5f}")
    rep.append("  M*(E-E_BO): " + ", ".join(f"{m:g}:{m * d:.5f}" for m, d in zip(Ms, dBO) if m >= 20))
    rep.append("  M^1.5*(E-E_BOad): " + ", ".join(f"{m:g}:{m ** 1.5 * d:.5f}" for m, d in zip(Ms, dBOad) if m >= 20))

    # ------------------------------------------------------------------ descriptors table
    keys = ["R", "R2", "lindemann", "invR", "invR2", "rXe", "T", "V", "virial", "T_X", "T_e", "dE_dlnM",
            "d2E_dlnM2", "d2E_dlnM2_fd_err", "g0_R3", "cm_curv_a2", "cm_s_max_over_ell", "cm_contrast"]
    drow = []
    for i, r in enumerate(ex):
        d = dict(M=r["M"], kappa=kap[i])
        for k in keys:
            d[k] = r[k]
            if k + "_prevN" in r and k not in ("d2E_dlnM2_fd_err",):
                d[k + "_basis_err"] = abs(r[k] - r[k + "_prevN"])
        d["cm_curv_a2_basis_err"] = abs(r["cm_curv_a2"] - r["cm_curv_a2_prevN"])
        d["d2E_dlnM2_err"] = d["d2E_dlnM2_fd_err"] + d["d2E_dlnM2_basis_err"]
        b = bo.get(round(r["M"], 6), {})
        for tag in ("BO", "BOad"):
            for k in ("R", "dR", "lindemann", "invR2", "Tnuc", "dE_dlnM", "d2E_dlnM2", "d2E_dlnM2_err"):
                d[f"{tag}_{k}"] = b.get(f"{tag}_{k}", float("nan"))
        drow.append(d)
    write_csv("descriptors.csv", drow, list(drow[0].keys()))
    dcol = lambda k: np.array([d[k] for d in drow], dtype=float)
    lind, R = dcol("lindemann"), dcol("R")
    rep.append(fit_report("Lindemann ratio dR/<R> (exact)", kap, lind, Ms, W, fits, "1 (harmonic zero-point amplitude)"))
    rep.append(fit_report("Lindemann ratio (BO)", kap, dcol("BO_lindemann"), Ms, W, fits, "1"))
    rep.append(fit_report("<R> - R_e (exact; R_e = 1.99719 bohr, BO limit)", kap, R - RE_BO, Ms, W, fits,
                          "2 (anharmonic shift ~ dR^2; plus O(1/M) from mu_e)"))
    rep.append(fit_report("-dE/dlnM = X kinetic energy (exact, Hellmann-Feynman)", kap, -dcol("dE_dlnM"), Ms, W, fits,
                          "2 (-dE/dlnM -> omega/4 + O(1/M))"))
    rep.append(fit_report("d2E/dlnM2 (exact, fixed basis)", kap, dcol("d2E_dlnM2"), Ms, W, fits,
                          "2 (-> omega/8)"))
    rep.append(fit_report("Lindemann exact - BO", kap, lind - dcol("BO_lindemann"), Ms, W, fits, None))
    rep.append(fit_report("<R> exact - BO", kap, R - dcol("BO_R"), Ms, W, fits, None))
    # other "crossover" markers: zero of d2E/dlnM2 (maximum of the X kinetic energy) and Lindemann maximum
    d2 = dcol("d2E_dlnM2")
    for i in range(len(Ms) - 1):
        if d2[i] < 0 <= d2[i + 1]:
            tt = -d2[i] / (d2[i + 1] - d2[i])
            Mz = math.exp(math.log(Ms[i]) + tt * math.log(Ms[i + 1] / Ms[i]))
            rep.append(f"--- d2E/dlnM2 changes sign (max of -dE/dlnM = X kinetic energy) between M={Ms[i]:g} "
                       f"and {Ms[i + 1]:g}: interpolated M = {Mz:.3f}  (BO+ad: see bo.csv; BO is unbound there)")
    il = int(np.nanargmax(lind))
    rep.append(f"--- Lindemann ratio maximum on the grid: {lind[il]:.5f} at M = {Ms[il]:g} (neighbours "
               f"{lind[max(il - 1, 0)]:.5f}, {lind[min(il + 1, len(lind) - 1)]:.5f}); 1/sqrt(3) = {1 / math.sqrt(3):.5f}")
    epsv = E0 / (Ms / (Ms + 1))
    ie = int(np.argmax(epsv))
    rep.append(f"--- H6: E/mu_Xe maximum on the grid {epsv[ie]:.6f} at M = {Ms[ie]:g} (non-monotonic: "
               f"{epsv[0]:.5f} at M={Ms[0]:g}, {epsv[-1]:.5f} at M={Ms[-1]:g})")
    # coalescence: ln g(0) vs kappa^-2 (only where resolvable)
    g0 = dcol("g0_R3")
    ok = (Ms >= 9) & (Ms <= 50)   # beyond M ~ 50 g(0) < 1e-8 <R>^-3 is below the basis resolution (non-monotonic)
    cc = np.polyfit(kap[ok] ** -2, np.log(g0[ok]), 1)
    loc = [(Ms[i], Ms[i + 1], (math.log(g0[i + 1]) - math.log(g0[i])) / (kap[i + 1] ** -2 - kap[i] ** -2))
           for i in range(len(Ms) - 1) if 1 <= Ms[i] and Ms[i + 1] <= 50]
    rep.append(f"--- coalescence: ln(<R>^3 g_XX(0)) vs kappa^-2 = sqrt(M), fit on M in [9.5, 50]: slope {cc[0]:.4f}, "
               f"intercept {cc[1]:.3f}; local slopes " + ", ".join(f"[{a:g},{b:g}] {s:.3f}" for a, b, s in loc) +
               ". Values (M >= 100 are basis noise, not physics): " +
               ", ".join(f"{m:g}:{g:.3e}" for m, g in zip(Ms, g0)))

    # ------------------------------------------------------------------ CM density crossover
    dens = load("density_cm.csv")
    crossover_txt = density_crossover(dens, ex, Ms)
    rep.append("")
    rep.append(crossover_txt)

    # ------------------------------------------------------------------ spacings
    srows = []
    for i, r in enumerate(ex):
        M = r["M"]
        b = bo.get(round(M, 6), {})
        e1, e1err, src = r["E1_ext"], r["E1_err"], "ground-state scales"
        if np.isfinite(r.get("E1_opt", float("nan"))) and r["E1_opt"] < r["E1"]:
            e1, e1err, src = r["E1_opt_ext"], r["E1_opt_err"], "root-1 scales"
        bound1 = e1 < r["E_th"] - 1e-10 and M > 9.5257
        L1 = rotd.get((round(M, 6), 1))
        L2 = rotd.get((round(M, 6), 2))
        d = dict(M=M, kappa=kap[i], E_th=r["E_th"], E00=E0[i], E00_err=E0err[i],
                 E10=e1 if bound1 else float("nan"), E10_err=e1err if bound1 else float("nan"), E10_source=src,
                 E01=L1["E_ext"] if L1 and L1["E_ext"] < r["E_th"] else float("nan"),
                 E01_err=L1["E_err"] if L1 else float("nan"),
                 E02=L2["E_ext"] if L2 and L2["E_ext"] < r["E_th"] else float("nan"),
                 E02_err=L2["E_err"] if L2 else float("nan"))
        d["vib"] = d["E10"] - d["E00"]
        d["vib_err"] = d["E10_err"] + d["E00_err"]
        d["d1"] = d["E01"] - d["E00"]
        d["d1_err"] = d["E01_err"] + d["E00_err"]
        d["d2"] = d["E02"] - d["E00"]
        d["d2_err"] = d["E02_err"] + d["E00_err"]
        d["d2_over_d1"] = d["d2"] / d["d1"]
        d["M_d1"] = M * d["d1"]
        d["vib_over_d1"] = d["vib"] / d["d1"]
        for tag in ("BO", "BOad"):
            e00 = b.get(f"{tag}_E00", float("nan"))
            d[f"{tag}_vib"] = b.get(f"{tag}_E10", float("nan")) - e00
            d[f"{tag}_d1"] = b.get(f"{tag}_E01", float("nan")) - e00
            d[f"{tag}_d2"] = b.get(f"{tag}_E02", float("nan")) - e00
            d[f"{tag}_d2_over_d1"] = d[f"{tag}_d2"] / d[f"{tag}_d1"]
            d[f"{tag}_vib_relerr"] = d[f"{tag}_vib"] / d["vib"] - 1
            d[f"{tag}_d1_relerr"] = d[f"{tag}_d1"] / d["d1"] - 1
            # drop relative errors not resolved by the exact calculation (error > 30% of the difference)
            for k, e in (("vib", "vib_err"), ("d1", "d1_err")):
                if not (d[e] < 0.3 * abs(d[f"{tag}_{k}"] - d[k])):
                    d[f"{tag}_{k}_relerr"] = float("nan")
        d["rigid_rotor_2invR2"] = 2 * r["invR2"] / M   # 2 B, B = <1/R^2>/(2 mu_N), mu_N = M/2
        srows.append(d)
    write_csv("spacings.csv", srows, list(srows[0].keys()))
    sc = lambda k: np.array([d[k] for d in srows], dtype=float)
    rep.append(fit_report("vibrational spacing E(v=1)-E(v=0) (exact)", kap, sc("vib"), Ms, W, fits, "2"))
    rep.append(fit_report("rotational spacing E(L=1)-E(L=0) (exact)", kap, sc("d1"), Ms, W, fits, "4"))
    rep.append(fit_report("BO+ad relative error of vib spacing", kap, sc("BOad_vib_relerr"), Ms, W, fits, None))
    rep.append(fit_report("BO+ad relative error of d1", kap, sc("BOad_d1_relerr"), Ms, W, fits, None))
    rep.append(fit_report("BO relative error of vib spacing", kap, sc("BO_vib_relerr"), Ms, W, fits, None))
    with open(os.path.join(HERE, "fits.txt"), "w") as fh:
        fh.write("\n".join(rep) + "\n")
    print("\n".join(rep))

    figures(ex, bo, drow, srows, dens, Ms, kap, E0, E0err, Ebo, Ebad, dBO, dBOad, eps)


def density_crossover(dens, ex, Ms):
    """Sign change of the curvature of rho_X(s) at s = 0 (about the CM), several fit windows."""
    prof = {}
    for d in dens:
        prof.setdefault(round(d["M"], 6), []).append((d["s_over_ell"], d["rho_over_rho0"]))
    a2 = {}
    for w in (0.05, 0.08, 0.1):
        a = []
        for M in Ms:
            p = np.array(sorted(prof.get(round(M, 6), [])))
            if len(p) == 0:
                a.append(np.nan)
                continue
            sel = p[:, 0] <= w + 1e-12
            A = np.vstack([np.ones(sel.sum()), p[sel, 0] ** 2, p[sel, 0] ** 4]).T
            c, *_ = np.linalg.lstsq(A, p[sel, 1], rcond=None)
            a.append(c[1] / c[0])
        a2[w] = np.array(a)
    lines = ["--- Density of X about the centre of mass: curvature a2 of rho(s)/rho(0) = 1 + a2 (s/l)^2 + ..., "
             "l = sqrt(<s^2>) (fits on s/l <= w)"]
    lines.append("  M: " + ", ".join(f"{m:g}:{x:+.4f}" for m, x in zip(Ms, a2[0.08]) if m <= 20))
    lines.append("  basis change of a2 (N vs N-4) and rho(0)/rho_max: " + ", ".join(
        f"{r['M']:g}:{abs(r['cm_curv_a2'] - r['cm_curv_a2_prevN']):.1e}/{1 / r['cm_contrast']:.1e}" for r in ex))
    est = []
    for w, a in a2.items():
        for i in range(len(Ms) - 1):
            if np.isfinite(a[i]) and np.isfinite(a[i + 1]) and a[i] < 0 <= a[i + 1]:
                # linear interpolation in ln M
                t = -a[i] / (a[i + 1] - a[i])
                Mc = math.exp(math.log(Ms[i]) + t * (math.log(Ms[i + 1]) - math.log(Ms[i])))
                est.append(Mc)
                lines.append(f"  window w={w}: sign change between M={Ms[i]:g} and {Ms[i + 1]:g}; "
                             f"interpolated (ln M, linear) M_c = {Mc:.4f}")
    basis = [abs(r["cm_curv_a2"] - r["cm_curv_a2_prevN"]) for r in ex if r["M"] <= 3]
    if est:
        lines.append(f"  RESULT (descriptor-dependent, reference point = CM): M_c = {np.mean(est):.3f}, spread over "
                     f"windows {max(est) - min(est):.3f}; basis change of a2 (N vs N-4) <= {max(basis):.1e} for M <= 3")
    return "\n".join(lines)


def figures(ex, bo, drow, srows, dens, Ms, kap, E0, E0err, Ebo, Ebad, dBO, dBOad, eps):
    big = Ms >= 1
    # ---------------- fig 1: energies and BO error
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.4))
    a = ax[0]
    a.plot(Ms, E0, **ST["exact"])
    a.plot(Ms, Ebo, **ST["BO"])
    a.plot(Ms, Ebad, **ST["BOad"])
    mm = np.geomspace(Ms.min(), Ms.max(), 200)
    a.plot(mm, -mm / (2 * (mm + 1)), color="gray", lw=0.8, ls=":", label="threshold $-M/2(M+1)$")
    a.set_xscale("log")
    a.set_xlabel("M (electron masses)")
    a.set_ylabel("E (hartree)")
    a.set_title("(a) ground state, L = 0")
    a.legend()
    a = ax[1]
    s = np.isfinite(dBO)
    a.loglog(kap[s], np.abs(dBO[s]), **dict(ST["BO"], label="|E − E$_{BO}$|"))
    s2 = np.isfinite(dBOad)
    a.loglog(kap[s2], np.abs(dBOad[s2]), **dict(ST["BOad"], label="|E − E$_{BO+ad}$|"))
    a.loglog(kap, E0err, color="gray", ls=":", marker=".", label="exact: extrap. error")
    kk = np.array([kap[Ms >= 500].min(), kap[Ms >= 500].max()])
    for p, ref, lab in ((4, np.abs(dBO[Ms >= 1836][0]) / kap[Ms >= 1836][0] ** 4, "κ$^4$"),
                        (6, np.abs(dBOad[Ms >= 1836][0]) / kap[Ms >= 1836][0] ** 6, "κ$^6$")):
        a.loglog(kk, ref * kk ** p, color="gray", lw=0.8)
        a.text(kk[1] * 1.05, ref * kk[1] ** p, lab, fontsize=8)
    for m, lab in ((1836.15267343, "H$_2^+$"), (206.768283, "Mu$_2^+$")):
        a.axvline(m ** -0.25, color="gray", lw=0.5, ls="--")
        a.text(m ** -0.25, 1.5e-11, lab, rotation=90, fontsize=7, va="bottom")
    a.set_xlabel("κ = M$^{-1/4}$")
    a.set_ylabel("|ΔE| (hartree)")
    a.set_title("(b) exact − BO (sign: E−E$_{BO}$>0, E−E$_{BO+ad}$<0)")
    a.legend(loc="upper left")
    a = ax[2]
    a.plot(Ms, eps, **ST["exact"])
    a.set_xscale("log")
    a.set_xlabel("M")
    a.set_ylabel("E / μ$_{Xe}$  (μ$_{Xe}$ = M/(M+1))")
    a.set_title("(c) normalized energy (H6)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig1_energies_bo_error.png"))
    plt.close(fig)

    # ---------------- fig 2: Lindemann ratio and <R>
    dc = lambda k: np.array([d[k] for d in drow], dtype=float)
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    a = ax[0]
    a.loglog(Ms, dc("lindemann"), **ST["exact"])
    a.loglog(Ms, dc("BO_lindemann"), **ST["BO"])
    a.loglog(Ms, dc("BOad_lindemann"), **ST["BOad"])
    a.axhline(1 / math.sqrt(3), color="gray", lw=0.6, ls=":")
    a.text(Ms.min(), 1 / math.sqrt(3) * 1.04, "1/√3 (two-body 1s)", fontsize=7)
    for m, lab in ((9.525634, "v=1 binds"), (4.2082827, "L=1 binds")):
        a.axvline(m, color="gray", lw=0.5, ls="--")
        a.text(m * 1.05, 0.65, lab, rotation=90, fontsize=7)
    a.set_xlabel("M")
    a.set_ylabel("Lindemann ratio ΔR/⟨R⟩")
    a.set_title("(a) X–X Lindemann ratio")
    a.legend(loc="lower left")
    a = ax[1]
    a.loglog(Ms, dc("R"), **ST["exact"])
    a.loglog(Ms, dc("BO_R"), **ST["BO"])
    a.loglog(Ms, dc("BOad_R"), **ST["BOad"])
    a.axhline(RE_BO, color="gray", lw=0.6, ls=":")
    a.set_xlabel("M")
    a.set_ylabel("⟨R⟩ (bohr)")
    a.set_title("(b) mean X–X distance (BO R$_e$ dotted)")
    a.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig2_lindemann_R.png"))
    plt.close(fig)

    # ---------------- fig 3: pair distribution
    gex = load("gR_exact.csv")
    gbo = load("gR_bo.csv")
    sel_M = [0.1, 0.5, 1.0, 3.0, 20.0, 206.768283, 1836.15267343, 1.0e4]
    cmap = plt.get_cmap("viridis")
    fig, ax = plt.subplots(1, 3, figsize=(14, 3.8))
    for j, M in enumerate(sel_M):
        c = cmap(j / (len(sel_M) - 1))
        pts = [(d["x"], d["P"]) for d in gex if abs(d["M"] - M) < 1e-6 * M and d["grid"] == "scaled"]
        if pts:
            p = np.array(pts)
            ax[0].plot(p[:, 0], p[:, 1], color=c, label=f"M = {M:g}")
        pts = [(d["x"], d["P"]) for d in gbo if abs(d["M"] - M) < 1e-6 * M and d["grid"] == "scaled"
               and d["model"] == "BOad"]
        if pts and M >= 20:
            p = np.array(pts)
            ax[0].plot(p[:, 0], p[:, 1], color=c, ls=":", lw=1.0)
        pts = [(d["x"], d["P"]) for d in gex if abs(d["M"] - M) < 1e-6 * M and d["grid"] == "fixed"]
        if pts and M >= 3:
            p = np.array(pts)
            ax[1].plot(p[:, 0], p[:, 1], color=c, label=f"M = {M:g}")
    ax[0].set_xlabel("R / ⟨R⟩")
    ax[0].set_ylabel("⟨R⟩ P(R)   (P = 4πR² g$_{XX}$)")
    ax[0].set_title("(a) X–X distance distribution, scaled\n(dotted: BO+ad, M ≥ 20; indistinguishable here)")
    ax[0].legend(fontsize=7)
    ax[1].set_xlabel("R (bohr)")
    ax[1].set_ylabel("P(R) (1/bohr)")
    ax[1].set_xlim(0, 8)
    ax[1].set_title("(b) unscaled, M ≥ 3")
    ax[1].legend(fontsize=7)
    for j, M in enumerate(sel_M):
        if M < 20 or M > 300:   # at M >~ 1000 BO DVR-grid interpolation noise exceeds the difference
            continue
        pe = np.array([(d["x"], d["P"]) for d in gex if abs(d["M"] - M) < 1e-6 * M and d["grid"] == "fixed"])
        for model, ls in (("BO", "--"), ("BOad", "-")):
            pb = np.array([(d["x"], d["P"]) for d in gbo if abs(d["M"] - M) < 1e-6 * M and d["grid"] == "fixed"
                           and d["model"] == model])
            if len(pe) and len(pb):
                ax[2].plot(pe[:, 0], (pe[:, 1] - pb[:, 1]) * M ** 0.75, color=cmap(j / (len(sel_M) - 1)), ls=ls,
                           label=f"M = {M:g}, " + ("BO" if model == "BO" else "BO+ad"))
    ax[2].set_xlim(0.5, 4.5)
    ax[2].set_xlabel("R (bohr)")
    ax[2].set_ylabel("[P$_{exact}$ − P$_{model}$] × M$^{3/4}$")
    ax[2].set_title("(c) exact minus BO pair distribution (M = 20, 206.8)")
    ax[2].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig3_gR.png"))
    plt.close(fig)

    # ---------------- fig 4: density about the CM
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    sel_M2 = [0.1, 0.3, 0.5, 0.6, 0.62, 0.66, 0.7, 0.8, 1.0, 1.5]
    for j, M in enumerate(sel_M2):
        pts = sorted((d["s_over_ell"], d["rho_over_rho0"]) for d in dens if abs(d["M"] - M) < 1e-6 * M)
        if pts:
            p = np.array(pts)
            ax[0].plot(p[:, 0], p[:, 1], color=cmap(j / (len(sel_M2) - 1)), label=f"M = {M:g}")
    ax[0].set_xlabel("s / ℓ   (s = |x$_X$ − x$_{CM}$|, ℓ = ⟨s²⟩$^{1/2}$)")
    ax[0].set_ylabel("ρ$_X$(s) / ρ$_X$(0)")
    ax[0].set_xlim(0, 2)
    ax[0].set_ylim(0, 1.3)
    ax[0].set_title("(a) one-particle density of X about the CM")
    ax[0].legend(fontsize=7, ncol=2)
    a2 = dc("cm_curv_a2")
    rel = np.abs(a2 - np.array([r["cm_curv_a2_prevN"] for r in ex])) / np.abs(a2)
    okc = (rel < 0.05) & (Ms <= 30)   # beyond M ~ 30 rho(0) is exponentially small: curvature is basis noise
    ax[1].plot(Ms[okc], a2[okc], **ST["exact"])
    ax[1].axhline(0, color="gray", lw=0.8)
    ax[1].set_xscale("log")
    ax[1].set_yscale("symlog", linthresh=1)
    ax[1].set_xlabel("M")
    ax[1].set_ylabel("curvature a$_2$ at s = 0")
    ax[1].axvspan(0.4, 0.8, color="tab:green", alpha=0.12, label="Mátyus et al. 2011 (0.4–0.8, secondary)")
    ax[1].axvline(0.491, color="tab:purple", lw=0.8, ls="--", label="King et al. 2013 (0.491, other criterion)")
    ax[1].axvline(0.98, color="tab:red", lw=0.8, ls=":", label="Salas et al. 2022 (≈0.98, other criterion)")
    ax[1].set_title("(b) curvature at the CM (a$_2$<0: maximum; >0: hole)")
    ax[1].legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig4_density_cm.png"))
    plt.close(fig)

    # ---------------- fig 5: mass derivatives
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    a = ax[0]
    a.loglog(Ms, -dc("dE_dlnM"), **dict(ST["exact"], label="exact −dE/dlnM (= X kinetic energy, HF)"))
    a.loglog(Ms, -dc("BO_dE_dlnM"), **ST["BO"])
    a.loglog(Ms, -dc("BOad_dE_dlnM"), **ST["BOad"])
    a.loglog(Ms, dc("BO_Tnuc"), color="tab:blue", ls=":", label="BO ⟨T$_R$⟩ (nuclear kinetic)")
    a.set_xlabel("M")
    a.set_ylabel("−dE/d ln M (hartree)")
    a.set_title("(a) first mass derivative")
    a.legend(fontsize=7)
    a = ax[1]
    y = dc("d2E_dlnM2")
    e = dc("d2E_dlnM2_err")
    a.errorbar(Ms, y, yerr=e, **{k: v for k, v in ST["exact"].items()})
    a.plot(Ms, dc("BO_d2E_dlnM2"), **ST["BO"])
    a.plot(Ms, dc("BOad_d2E_dlnM2"), **ST["BOad"])
    a.axhline(0, color="gray", lw=0.6)
    a.set_xscale("log")
    a.set_xlabel("M")
    a.set_ylabel("d²E/d(ln M)² (hartree)")
    a.set_title("(b) second derivative (exact: fixed basis, FD of HF)")
    a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig5_mass_derivatives.png"))
    plt.close(fig)

    # ---------------- fig 6: spacings
    sc = lambda k: np.array([d[k] for d in srows], dtype=float)
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
    a = ax[0]
    a.loglog(Ms, sc("vib"), **dict(ST["exact"], label="exact E(v=1)−E(v=0)"))
    a.loglog(Ms, sc("BO_vib"), **dict(ST["BO"], label="BO"))
    a.loglog(Ms, sc("BOad_vib"), **dict(ST["BOad"], label="BO+ad"))
    a.loglog(Ms, sc("d1"), color="k", marker="D", ls="-", mfc="w", label="exact E(L=1)−E(L=0)")
    a.loglog(Ms, sc("BO_d1"), color="tab:blue", marker="D", mfc="w", ls="--")
    a.loglog(Ms, sc("BOad_d1"), color="tab:orange", marker="D", mfc="w", ls="-.")
    a.loglog(Ms, sc("d2"), color="k", marker="v", ls="-", mfc="w", label="exact E(L=2)−E(L=0)")
    a.set_xlabel("M")
    a.set_ylabel("spacing (hartree)")
    a.set_title("(a) vibrational and rotational spacings")
    a.legend(fontsize=7)
    a = ax[1]
    a.semilogx(Ms, sc("d2_over_d1"), **dict(ST["exact"], label="exact d$_2$/d$_1$"))
    a.semilogx(Ms, sc("BO_d2_over_d1"), **dict(ST["BO"], label="BO"))
    a.semilogx(Ms, sc("BOad_d2_over_d1"), **dict(ST["BOad"], label="BO+ad"))
    a.axhline(3, color="gray", lw=0.6, ls=":")
    a2_ = a.twinx()
    a2_.loglog(Ms, sc("vib_over_d1"), color="tab:green", ls=":", marker="x", label="vib/d$_1$ (right axis)")
    a2_.set_ylabel("vib / d$_1$", color="tab:green")
    a2_.legend(fontsize=7, loc="lower right")
    a.set_xlabel("M")
    a.set_ylabel("ratio")
    a.set_title("(b) tower ratios (rigid rotor: 3)")
    a.legend(fontsize=7)
    a = ax[2]
    for k, st, lab in (("BO_vib_relerr", ST["BO"], "BO, vib"), ("BOad_vib_relerr", ST["BOad"], "BO+ad, vib"),
                       ("BO_d1_relerr", dict(ST["BO"], mfc="w"), "BO, rot d$_1$"),
                       ("BOad_d1_relerr", dict(ST["BOad"], mfc="w"), "BO+ad, rot d$_1$")):
        y = sc(k)
        s = np.isfinite(y)
        a.loglog(kap[s], np.abs(y[s]), **dict(st, label=lab + " (|rel. err.|)"))
    a.set_xlabel("κ = M$^{-1/4}$")
    a.set_ylabel("|spacing(model)/spacing(exact) − 1|")
    a.set_title("(c) relative error of BO spacings")
    a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig6_spacings.png"))
    plt.close(fig)


if __name__ == "__main__":
    main()
