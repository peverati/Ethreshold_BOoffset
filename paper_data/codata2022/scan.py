"""Exact (threebody) versus Born-Oppenheimer (bo_h2plus) mass scans of {X+, X+, e-}.

Usage (from anywhere):
    python scan.py bo            # BO and BO+adiabatic levels, J = 0, 1, 2 (fast, ~1 min)
    python scan.py exact [M...]  # exact L = 0 'S' roots 0, 1 + descriptors (restartable)
    python scan.py rot [M...]    # exact L = 1 'A' and L = 2 'S' lowest levels (restartable)
then  python analyze.py         # fits, tables and figures.

Every stage appends one row per mass to its CSV and skips masses already present, so an
interrupted run is resumed by running the same command again.
"""
import csv
import json
import math
import os
import sys
import time

from common import (DATA, HERE, MASSES, R_FIXED, X_SCALED, L1_CRIT, L2_CRIT, V1_CRIT,
                    density_about_cm, family_threshold, bo_imports)
import numpy as np
from scipy.interpolate import CubicSpline

from threebody import symmetric_family, PerimetricBasis, MatrixSet, solve, optimize_alphas, extrapolate
from threebody.rotational import optimize_scales, convergence_table_L, observables_L


def log(msg):
    line = time.strftime("%H:%M:%S ") + msg
    print(line, flush=True)
    with open(os.path.join(HERE, "run_log.txt"), "a") as fh:
        fh.write(line + "\n")


def done_masses(path, key="M"):
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return {round(float(r[key]), 9): r for r in csv.DictReader(fh)}


class Appender:
    def __init__(self, path, fields):
        self.path, self.fields = path, fields
        new = not os.path.exists(path)
        self.fh = open(path, "a", newline="")
        self.w = csv.DictWriter(self.fh, fieldnames=fields, extrasaction="ignore")
        if new:
            self.w.writeheader()

    def row(self, d):
        self.w.writerow({k: (repr(float(v)) if isinstance(v, (float, np.floating)) else v) for k, v in d.items()})
        self.fh.flush()


# ======================================================================= BO stage
def stage_bo(masses=MASSES):
    solve_levels, thresholds = bo_imports()
    path = os.path.join(HERE, "bo.csv")
    fields = ["M", "kappa", "E_thr_BO", "E_thr_BOad", "E_thr_exact"]
    for tag in ("BO", "BOad"):
        fields += [f"{tag}_n_J0", f"{tag}_E00", f"{tag}_err00", f"{tag}_E10", f"{tag}_E01", f"{tag}_E02",
                   f"{tag}_R", f"{tag}_dR", f"{tag}_lindemann", f"{tag}_invR2", f"{tag}_Tnuc",
                   f"{tag}_dE_dlnM", f"{tag}_dE_dlnM_err", f"{tag}_d2E_dlnM2", f"{tag}_d2E_dlnM2_err"]
    out = Appender(path, fields)
    have = done_masses(path)
    gpath = os.path.join(HERE, "gR_bo.csv")
    gout = Appender(gpath, ["M", "model", "grid", "x", "P"])

    def E0(M, ad):
        s = solve_levels(M, 0, adiabatic=ad, check=False, ppw=4.0)
        return s.levels[0].E if s.levels else float("nan")

    for M in masses:
        if round(M, 9) in have:
            continue
        t0 = time.time()
        thr_bo, thr_ex = thresholds(M)
        row = dict(M=M, kappa=M ** -0.25, E_thr_BO=thr_bo, E_thr_exact=thr_ex,
                   E_thr_BOad=thresholds(M, adiabatic=True)[0])
        for tag, ad in (("BO", False), ("BOad", True)):
            s0 = solve_levels(M, 0, adiabatic=ad, check=True)
            row[f"{tag}_n_J0"] = len(s0.levels)
            if not s0.levels:
                continue
            l0 = s0.levels[0]
            row.update({f"{tag}_E00": l0.E, f"{tag}_err00": l0.err_E, f"{tag}_R": l0.R_mean, f"{tag}_dR": l0.dR,
                        f"{tag}_lindemann": l0.lindemann, f"{tag}_invR2": l0.inv_R2, f"{tag}_Tnuc": l0.T_kin})
            if len(s0.levels) > 1:
                row[f"{tag}_E10"] = s0.levels[1].E
            for J in (1, 2):
                sJ = solve_levels(M, J, adiabatic=ad, check=False)
                if sJ.levels:
                    row[f"{tag}_E0{J}"] = sJ.levels[0].E
            # finite-difference mass derivatives in ln M (Richardson over h, h/2)
            try:
                d1, d2 = {}, {}
                for h in (0.02, 0.01):
                    ep, em = E0(M * math.exp(h), ad), E0(M * math.exp(-h), ad)
                    ec = E0(M, ad)
                    d1[h] = (ep - em) / (2 * h)
                    d2[h] = (ep - 2 * ec + em) / h ** 2
                row[f"{tag}_dE_dlnM"] = (4 * d1[0.01] - d1[0.02]) / 3
                row[f"{tag}_dE_dlnM_err"] = abs(d1[0.01] - d1[0.02]) / 3
                row[f"{tag}_d2E_dlnM2"] = (4 * d2[0.01] - d2[0.02]) / 3
                row[f"{tag}_d2E_dlnM2_err"] = abs(d2[0.01] - d2[0.02]) / 3 + 4 * 1e-10 / 0.01 ** 2
            except Exception as exc:  # level may disappear at M e^-h near the BO critical mass
                log(f"BO derivative failed at M={M}: {exc}")
            # radial density of v=0 on the fixed and scaled grids (DVR grid -> interpolation)
            P = s0.P[0]
            Rg = s0.R
            for grid, xs in (("fixed", R_FIXED), ("scaled", X_SCALED)):
                Rq = xs if grid == "fixed" else xs * l0.R_mean
                cs = CubicSpline(np.concatenate([[0.0], Rg]), np.concatenate([[0.0], P]))
                Pq = np.where(Rq <= Rg[-1], cs(np.minimum(Rq, Rg[-1])), 0.0)
                if grid == "scaled":
                    Pq = Pq * l0.R_mean
                for x, p in zip(xs, Pq):
                    gout.row(dict(M=M, model=tag, grid=grid, x=x, P=p))
        out.row(row)
        log(f"BO M={M:g}: E0(BO)={row.get('BO_E00', float('nan')):.12f} "
            f"E0(BO+ad)={row.get('BOad_E00', float('nan')):.12f} t={time.time() - t0:.1f}s")


# ==================================================================== exact L = 0
def settings_L0(M):
    if M < 150:
        return dict(weights=(1, 1, 1), N_opt=12, Ns=[16, 20, 24, 28])
    if M < 3000:
        return dict(weights=(1, 1, 2), N_opt=20, Ns=[28, 32, 36, 40])
    if M < 5000:
        return dict(weights=(1, 1, 2), N_opt=20, Ns=[32, 36, 40, 44])
    return dict(weights=(1, 1, 2), N_opt=20, Ns=[36, 40, 44, 48])


def descriptors(res, M):
    """Observables of root 0 at fixed basis, incl. Hellmann-Feynman and fixed-basis 2nd derivative."""
    ob = res.observables(0)
    T0 = res.dE_dinverse_mass(0, 0)
    T1 = res.dE_dinverse_mass(1, 0)
    T2 = res.dE_dinverse_mass(2, 0)
    d = dict(R=ob["r01"], R2=ob["r01^2"], lindemann=ob["lindemann01"], invR=ob["1/r01"], invR2=ob["1/r01^2"],
             rXe=ob["r02"], rXe2=ob["r02^2"], T=ob["T"], V=ob["V"], virial=ob["virial_ratio"],
             T_X=(T0 + T1) / M, T_e=T2, dE_dlnM=-(T0 + T1) / M)
    # second derivative: central difference of the exact HF first derivative, same matrices
    ms, b = res.matrices, res.basis
    D1 = {}
    E = {}
    for h in (0.02, 0.01):
        for sgn in (+1, -1):
            Mp = M * math.exp(sgn * h)
            r = solve(symmetric_family(Mp), b, nstates=1, matrices=ms)
            D1[(h, sgn)] = -(r.dE_dinverse_mass(0, 0) + r.dE_dinverse_mass(1, 0)) / Mp
            E[(h, sgn)] = r.energies[0]
    d2 = {h: (D1[(h, 1)] - D1[(h, -1)]) / (2 * h) for h in (0.02, 0.01)}
    d["d2E_dlnM2"] = (4 * d2[0.01] - d2[0.02]) / 3
    d["d2E_dlnM2_fd_err"] = abs(d2[0.01] - d2[0.02]) / 3
    # energy-only check (second difference of the eigenvalue, h = 0.01)
    d["d2E_dlnM2_Echeck"] = (E[(0.01, 1)] - 2 * res.energies[0] + E[(0.01, -1)]) / 0.01 ** 2
    # Hellmann-Feynman check of the first derivative
    d["dE_dlnM_fdcheck"] = (E[(0.01, 1)] - E[(0.01, -1)]) / 0.02
    # coalescence: g(R) at R -> 0, dimensionless <R>^3 g(0)
    Rs = np.array([1e-3, 2e-3]) * d["R"]
    g = res.g_of_r(0, 1, Rs)
    d["g0_R3"] = (2 * g[0] - g[1]) * d["R"] ** 3
    return d


def cm_density_profile(res, M):
    """rho_X about the CM on s/l in [0, 3], l = sqrt(<s^2>); small-s fit for the curvature."""
    R2 = res.expect_power(0, 1, 2)
    r022 = res.expect_power(0, 2, 2)
    s2 = R2 / 4 + (r022 - R2 / 4) / (2 * M + 1) ** 2
    ell = math.sqrt(s2)
    xs = np.concatenate([np.linspace(0, 0.1, 11), np.linspace(0.12, 0.2, 5), np.linspace(0.25, 3.0, 56)])
    rho = density_about_cm(res, M, xs * ell)
    small = xs <= 0.08 + 1e-12
    A = np.vstack([np.ones(small.sum()), xs[small] ** 2, xs[small] ** 4]).T
    coef, *_ = np.linalg.lstsq(A, rho[small], rcond=None)
    a2 = coef[1] / coef[0]
    imax = int(np.argmax(rho))
    return dict(ell=ell, s2=s2, xs=xs, rho=rho, curv_a2=a2, s_max_over_ell=xs[imax],
                contrast=rho[imax] / rho[0])


def stage_exact(masses):
    path = os.path.join(HERE, "exact_L0.csv")
    desc_keys = ["R", "R2", "lindemann", "invR", "invR2", "rXe", "rXe2", "T", "V", "virial", "T_X", "T_e",
                 "dE_dlnM", "d2E_dlnM2", "d2E_dlnM2_fd_err", "d2E_dlnM2_Echeck", "dE_dlnM_fdcheck", "g0_R3"]
    fields = (["M", "kappa", "E_th", "weights", "alpha_X", "alpha_e", "N_last", "size_last", "dropped",
               "E0", "E0_ext", "E0_err", "E0_flag", "E0_seq", "E1", "E1_ext", "E1_err", "E1_flag", "E1_seq",
               "E1_opt", "E1_opt_ext", "E1_opt_err", "E1_opt_flag", "E1_opt_seq", "alpha1_X", "alpha1_e"]
              + desc_keys + [k + "_prevN" for k in desc_keys]
              + ["cm_ell", "cm_curv_a2", "cm_curv_a2_prevN", "cm_s_max_over_ell", "cm_contrast",
                 "cm_norm_check", "E1_R", "E1_lindemann", "time_s"])
    out = Appender(path, fields)
    have = done_masses(path)
    gout = Appender(os.path.join(HERE, "gR_exact.csv"), ["M", "grid", "x", "P", "g"])
    dout = Appender(os.path.join(HERE, "density_cm.csv"), ["M", "s_over_ell", "s", "rho", "rho_over_rho0"])
    x0 = None
    for M in masses:
        if round(M, 9) in have:
            continue
        if have and x0 is None:   # warm start from the nearest mass already computed
            Mn = min(have, key=lambda m: abs(math.log(m / M)))
            x0 = (float(have[Mn]["alpha_X"]),) * 2 + (float(have[Mn]["alpha_e"]),)
        t0 = time.time()
        st = settings_L0(M)
        s = symmetric_family(M)
        if x0 is not None and M > 150 and x0[0] < 4:  # jump to the (1,1,2) regime: restart from H2+-like
            x0 = None
        guess = x0 if x0 is not None else ((8.0 * (M / 1836.15) ** 0.25,) * 2 + (1.4,) if M > 150 else None)
        al, _, nev = optimize_alphas(s, st["N_opt"], "S", x0=guess, weights=st["weights"])
        x0 = al
        seq0, seq1, sizes = [], [], []
        last = dp = cmp_ = None
        for N in st["Ns"]:
            b = PerimetricBasis(N, al, "S", st["weights"])
            r = solve(s, b, nstates=2, matrices=MatrixSet(b))
            seq0.append(float(r.energies[0]))
            seq1.append(float(r.energies[1]))
            sizes.append(b.size)
            if N == st["Ns"][-2]:   # descriptors at the second-largest basis (convergence evidence)
                dp = descriptors(r, M)
                cmp_ = cm_density_profile(r, M)
            last = r
            del r
        e0 = extrapolate(seq0)
        e1 = extrapolate(seq1)
        row = dict(M=M, kappa=M ** -0.25, E_th=family_threshold(M), weights=json.dumps(st["weights"]),
                   alpha_X=al[0], alpha_e=al[2], N_last=st["Ns"][-1], size_last=sizes[-1],
                   dropped=last.info["dropped"],
                   E0=seq0[-1], E0_ext=e0[0], E0_err=e0[1], E0_flag=e0[3], E0_seq=json.dumps(seq0),
                   E1=seq1[-1], E1_ext=e1[0], E1_err=e1[1], E1_flag=e1[3], E1_seq=json.dumps(seq1))
        d = descriptors(last, M)
        row.update(d)
        row.update({k + "_prevN": v for k, v in dp.items()})
        ob1 = last.observables(1, inverse_square=False)
        row.update(E1_R=ob1["r01"], E1_lindemann=ob1["lindemann01"])
        # v = 1 with scales optimized for root 1 (only where it can be bound and converges slowly)
        if V1_CRIT - 0.01 < M < 150:
            al1, _, _ = optimize_alphas(s, 12, "S", state=1, x0=al)
            sq = []
            for N in [16, 20, 24, 28]:
                b1 = PerimetricBasis(N, al1, "S")
                sq.append(float(solve(s, b1, nstates=2, matrices=MatrixSet(b1)).energies[1]))
            e1o = extrapolate(sq)
            row.update(E1_opt=sq[-1], E1_opt_ext=e1o[0], E1_opt_err=e1o[1], E1_opt_flag=e1o[3],
                       E1_opt_seq=json.dumps(sq), alpha1_X=al1[0], alpha1_e=al1[2])
        # pair distribution (root 0)
        for grid, xs in (("fixed", R_FIXED), ("scaled", X_SCALED)):
            Rq = xs if grid == "fixed" else xs * d["R"]
            P = last.pair_distribution(0, 1, Rq, 0)
            with np.errstate(divide="ignore", invalid="ignore"):
                g = np.where(Rq > 0, P / (4 * np.pi * Rq ** 2), np.nan)
            if grid == "scaled":
                P = P * d["R"]
            for x, p, gg in zip(xs, P, g):
                gout.row(dict(M=M, grid=grid, x=x, P=p, g=gg))
        # density of X about the centre of mass
        cm = cm_density_profile(last, M)
        # normalization check: trapezoid of 4 pi s^2 rho over s (coarse, profile to 3 ell)
        ss = cm["xs"] * cm["ell"]
        row.update(cm_ell=cm["ell"], cm_curv_a2=cm["curv_a2"], cm_curv_a2_prevN=cmp_["curv_a2"],
                   cm_s_max_over_ell=cm["s_max_over_ell"], cm_contrast=cm["contrast"],
                   cm_norm_check=float(np.trapz(4 * np.pi * ss ** 2 * cm["rho"], ss)))
        for x, sv, rv in zip(cm["xs"], ss, cm["rho"]):
            dout.row(dict(M=M, s_over_ell=x, s=sv, rho=rv, rho_over_rho0=rv / cm["rho"][0]))
        last.save(os.path.join(DATA, f"L0_M{M:.6g}.npz"), states=[0, 1])
        row["time_s"] = round(time.time() - t0, 1)
        out.row(row)
        log(f"exact L0 M={M:g}: E0={seq0[-1]:.13f} ext {e0[0]:.13f}+/-{e0[1]:.1e} {e0[3]} | "
            f"E1={seq1[-1]:.11f} | Lind={d['lindemann']:.5f} a2={cm['curv_a2']:+.4f} t={row['time_s']}s "
            f"(nev={nev}, alphas={np.round(al, 4).tolist()})")


# ==================================================================== exact L = 1, 2
def settings_rot(M):
    if M < 150:
        return dict(weights=(1, 1, 1), N_opt=12, Ns=[14, 18, 21])
    if M < 3000:
        return dict(weights=(1, 1, 2), N_opt=16, Ns=[24, 28, 32])
    if M < 5000:
        return dict(weights=(1, 1, 2), N_opt=16, Ns=[26, 30, 34])
    return dict(weights=(1, 1, 2), N_opt=20, Ns=[30, 34, 38])  # L=1 only; the L=2 row at 1e4 used [26, 30, 34]


def stage_rot(masses):
    path = os.path.join(HERE, "exact_rot.csv")
    fields = ["M", "L", "sym", "E_th", "weights", "alpha_X", "alpha_e", "N_last", "size_last", "E", "E_ext",
              "E_err", "E_flag", "E_seq", "R", "R2", "virial", "time_s"]
    out = Appender(path, fields)
    have = {}
    if os.path.exists(path):
        with open(path) as fh:
            have = {(round(float(r["M"]), 9), int(r["L"])) for r in csv.DictReader(fh)}
    x0 = None
    for M in masses:
        cases = []
        if M > L1_CRIT:
            cases.append((1, "A"))
        if M > L2_CRIT:
            cases.append((2, "S"))
        st = settings_rot(M)
        s = symmetric_family(M)
        al = None
        for L, sym in cases:
            if (round(M, 9), L) in have:
                continue
            t0 = time.time()
            if al is None:
                guess = x0
                if guess is not None and M > 150 and guess[0] < 4:
                    guess = None
                if guess is None and M > 150:
                    guess = (7.7 * (M / 1836.15) ** 0.25,) * 2 + (1.3,)
                al, _, _ = optimize_scales(s, 1, st["N_opt"], "A", x0=guess, weights=st["weights"])
                x0 = al
            tab = convergence_table_L(s, L, sym, al, st["Ns"], nstates=1, weights=st["weights"],
                                      observables=True)
            last = tab["rows"][-1]
            E, err, rho, flag = tab["extrap"][0]
            ob = tab["observables"][0]
            out.row(dict(M=M, L=L, sym=sym, E_th=family_threshold(M), weights=json.dumps(st["weights"]),
                         alpha_X=al[0], alpha_e=al[2], N_last=last["N"], size_last=last["size"],
                         E=float(last["energies"][0]), E_ext=float(E), E_err=float(err), E_flag=flag,
                         E_seq=json.dumps([float(r["energies"][0]) for r in tab["rows"]]),
                         R=ob["r01"], R2=ob["r01_2"], virial=ob["virial"], time_s=round(time.time() - t0, 1)))
            log(f"rot M={M:g} L={L}{sym}: E={last['energies'][0]:.13f} ext {E:.13f}+/-{err:.1e} {flag} "
                f"t={time.time() - t0:.0f}s")


if __name__ == "__main__":
    stage = sys.argv[1]
    ms = [float(x) for x in sys.argv[2:]] or MASSES
    if stage == "bo":
        stage_bo(ms)
    elif stage == "exact":
        stage_exact(ms)
    elif stage == "rot":
        stage_rot(ms)
    else:
        raise SystemExit(__doc__)
