"""Analysis of the entanglement scan: errors, fits, features, figures.  *Last updated: 2026-11-13*

Reads results_{exact,fine,prevN,bo}.csv and data/ent_*_M*.npz; writes table.csv, fits.txt, fig1..fig4.
"""
import csv
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.special import eval_legendre  # noqa: E402
from scipy.interpolate import CubicSpline  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
MARKERS = [(0.25, "Lindemann max"), (0.4924, "King g min"), (0.5, "E/mu max"), (0.543, "0.543"),
           (0.635, "CM density"), (1.11, "max T_X"), (1.4, "~1.4"), (4.21, "L=1 binds"),
           (9.53, "v=1 binds"), (11.49, "L=2 binds")]
S_OR_INF = 0.249338007   # electronic_l_weights(R_e = 1.99719), converged to 1e-8
R_E = 1.99719


def load(mode):
    fn = os.path.join(HERE, f"results_{mode}.csv")
    if not os.path.exists(fn):
        return {}
    out = {}
    with open(fn) as fh:
        for r in csv.DictReader(fh):
            out[round(float(r["M"]), 6)] = {k: float(v) for k, v in r.items()}
    return out


def npz(mode, M):
    fn = os.path.join(DATA, f"ent_{mode}_M{M:g}.npz")
    return np.load(fn) if os.path.exists(fn) else None


def extras(mode, M):
    """Quantities from the saved spectra: second Schmidt ratio of l=0, D(pi/2), l-tail estimate."""
    d = npz(mode, M)
    if d is None:
        return {}
    ls, w, lb = d["ls"], d["w"], d["lam_blocks"]
    b0 = lb[0] / lb[0].sum()
    D90 = float(np.sum(w * eval_legendre(ls, 0.0)))
    # tail beyond lmax with w_l ~ w_lmax (l/lmax)^-6 (measured decay, README)
    lmax = ls[-1]
    lt = np.arange(lmax + 2, 20 * lmax, 2)
    wt = w[-1] * (lt / lmax) ** -6.0
    tail_S = float(np.sum(wt * (np.log(2 * lt + 1) - np.log(wt))))
    return dict(p2_l0=float(b0[1]), p1_l0=float(b0[0]), D90=D90, tail_w=float(wt.sum()), tail_S=tail_S)


def slopes(x, y):
    lx, ly = np.log(x), np.log(np.abs(y))
    return (ly[1:] - ly[:-1]) / (lx[1:] - lx[:-1])


def fit(x, y):
    p = np.polyfit(np.log(x), np.log(np.abs(y)), 1, full=False, cov=True)
    return p[0][0], math.sqrt(p[1][0, 0]), math.exp(p[0][1])


def main():
    ex, fi, pn, bo = load("exact"), load("fine"), load("prevN"), load("bo")
    Ms = np.array(sorted(ex))
    keys = ["S", "S_lin", "S_orient", "S_rad", "S_0", "w0", "w2", "E_N", "MI"]
    lines = ["# entanglement fits and features", "*Last updated: 2026-11-13*", ""]
    # ---------------------------------------------------------------- table with errors
    rows = []
    for M in Ms:
        r = dict(ex[M])
        r.update(extras("exact", M))
        for k in keys:
            eq = abs(fi[M][k] - r[k]) if M in fi else float("nan")
            eb = abs(pn[M][k] - r[k]) if M in pn else float("nan")
            r[k + "_quad_err"] = eq
            r[k + "_basis_err"] = eb
            if M in bo:
                r[k + "_BO"] = bo[M][k]
        rows.append(r)
    cols = ["M", "kappa"] + [c for k in keys for c in (k, k + "_quad_err", k + "_basis_err", k + "_BO")] + \
        ["lam_max", "K", "p2_l0", "D90", "tail_w", "tail_S", "w_lmax", "direct_norm"]
    with open(os.path.join(HERE, "table.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        wr.writeheader()
        for r in rows:
            wr.writerow({c: (f"{r[c]:.12g}" if c in r else "") for c in cols})
    T = {k: np.array([r.get(k, np.nan) for r in rows]) for k in cols}

    # ---------------------------------------------------------------- error summary
    lines.append("## Convergence (max over M of |base - fine| quadrature, |N - (N-4)| basis, l-tail estimate)")
    for k in ["S", "S_orient", "S_rad", "E_N", "MI"]:
        lines.append(f"{k:9s} quad max {np.nanmax(T[k+'_quad_err']):.2e}  basis max {np.nanmax(T[k+'_basis_err']):.2e}")
    lines.append(f"l-tail estimate of S (w_l ~ l^-6 beyond lmax): max {np.nanmax(T['tail_S']):.2e}")
    lines.append("")

    # ---------------------------------------------------------------- large-M fits
    lines.append("## Large-M behaviour (window M >= 500 for fits; local slopes over the whole M >= 20 range)")
    kap = T["kappa"]
    big = Ms >= 20
    win = Ms >= 500
    q = {"S - S_orient_inf": T["S"] - S_OR_INF, "S_orient - S_orient_inf": T["S_orient"] - S_OR_INF,
         "S_rad": T["S_rad"], "S_0": T["S_0"], "p2 (2nd Schmidt coeff of l=0 block)": T["p2_l0"],
         "E_N (radial log-negativity)": T["E_N"], "MI(|R|;|rho|)": T["MI"], "S_lin": T["S_lin"] - (1 - np.nan)}
    q.pop("S_lin")
    for name, y in q.items():
        if np.any(~np.isfinite(y[win])):
            continue
        p, dp, a = fit(kap[win], y[win])
        sl = slopes(kap[big], y[big])
        lines.append(f"{name}: exponent in kappa {p:.3f} +- {dp:.3f}, prefactor {a:.4g};  local slopes "
                     + " ".join(f"{s:.3f}" for s in sl))
    lines.append(f"S_orient_inf (electronic sigma_g partial waves at R_e = {R_E}) = {S_OR_INF}")
    lines.append("")

    # ---------------------------------------------------------------- features
    lines.append("## Extrema / inflections on the M grid (cubic spline in ln M; grid-resolved only)")
    de = load("dense")      # fixed-basis re-solves on 0.25-0.6 (run.py dense), merged with the scan grid
    allM = sorted(set(Ms) | set(de))
    lm = np.log(allM)
    lines.append(f"(spline through {len(allM)} masses: scan grid + dense {sorted(de)})")
    for k in ["S", "S_lin", "S_orient", "S_rad", "S_0", "w0", "E_N", "MI", "D90"]:
        src = []
        for M in allM:
            if M in ex:
                src.append(T[k][list(Ms).index(M)])
            elif k in de[M]:
                src.append(de[M][k])
            else:
                src.append(extras("dense", M).get(k, np.nan))
        y = np.array(src)
        cs = CubicSpline(lm, y)
        xx = np.linspace(lm[0], lm[-1], 20001)
        d1, d2 = cs(xx, 1), cs(xx, 2)
        ext = xx[1:][np.sign(d1[1:]) != np.sign(d1[:-1])]
        inf = xx[1:][np.sign(d2[1:]) != np.sign(d2[:-1])]
        ex_s = ", ".join(f"M={math.exp(e):.3f} ({'max' if cs(e, 2) < 0 else 'min'}, {k}={float(cs(e)):.6f})"
                         for e in ext)
        lines.append(f"{k}: extrema {ex_s or 'none'}; inflections at M = "
                     + ", ".join(f"{math.exp(e):.3g}" for e in inf if Ms[0] * 1.05 < math.exp(e) < Ms[-1] / 1.05))
    open(os.path.join(HERE, "fits.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))

    # ---------------------------------------------------------------- figures
    def marks(ax):
        for m, _ in MARKERS:
            ax.axvline(m, color="0.85", lw=0.8, zorder=0)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
    ax = axs[0]
    for k, lab, c in [("S", "S (von Neumann)", "k"), ("S_orient", "S_orient", "C0"), ("S_rad", "S_rad", "C3"),
                      ("S_lin", "linear entropy", "C2"), ("E_N", "E_N radial (log-neg.)", "C1"),
                      ("MI", "I(|R|;|rho|)", "C4")]:
        ax.semilogx(Ms, T[k], "o-", ms=3, color=c, label=lab)
        if k + "_BO" in T:
            ax.semilogx(Ms, T[k + "_BO"], "--", color=c, lw=1)
    ax.axhline(S_OR_INF, color="C0", ls=":", lw=1, label="S_orient(M=inf)")
    marks(ax)
    ax.set_xlabel("M (units of m_e)")
    ax.set_ylabel("entropy (nats)")
    ax.set_title("R | rho entanglement, exact (solid) and BO product (dashed)")
    ax.legend(fontsize=7)
    ax = axs[1]
    ax.semilogx(Ms, T["w0"], "o-", ms=3, label="w_0")
    ax.semilogx(Ms, T["w2"], "o-", ms=3, label="w_2")
    ax.semilogx(Ms, T["D90"], "s-", ms=3, label="D(90 deg) orientational coherence")
    ax.semilogx(Ms, T["p2_l0"] * 10, "^-", ms=3, label="10 x p2 (l=0 block)")
    marks(ax)
    ax.set_xlabel("M")
    ax.legend(fontsize=7)
    ax.set_title("partial-wave weights and orientational coherence")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig1_entropies_vs_M.png"), dpi=150)

    fig, ax = plt.subplots(figsize=(6, 4.6))
    sel = Ms >= 3
    for name, y, c in [("S - S_or(inf)", T["S"] - S_OR_INF, "k"), ("S_orient - S_or(inf)", T["S_orient"] - S_OR_INF, "C0"),
                       ("S_rad", T["S_rad"], "C3"), ("p2 (l=0)", T["p2_l0"], "C2"), ("E_N", T["E_N"], "C1"),
                       ("MI", T["MI"], "C4")]:
        ax.loglog(kap[sel], np.abs(y[sel]), "o-", ms=3, color=c, label=name)
    kk = np.array([kap[sel].min(), kap[sel].max()])
    for p, ls_ in [(1, ":"), (2, "--")]:
        ax.loglog(kk, 0.3 * (kk / kk[1]) ** p, ls_, color="0.5", lw=0.8, label=f"kappa^{p}")
    ax.set_xlabel("kappa = M^(-1/4)")
    ax.legend(fontsize=7)
    ax.set_title("large-M approach (M >= 3)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig2_largeM_scaling.png"), dpi=150)

    fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
    ax = axs[0]
    th = np.linspace(0, np.pi, 181)
    for M in [0.05, 1.0, 9.5256, 206.768283, 1.0e4]:
        d = npz("exact", M)
        if d is None:
            continue
        D = sum(wl * eval_legendre(l, np.cos(th)) for l, wl in zip(d["ls"], d["w"]))
        ax.plot(np.degrees(th), D, label=f"M={M:g}")
    ax.set_xlabel("angle between R and R' (deg)")
    ax.set_ylabel("D = sum_l w_l P_l(cos)")
    ax.set_title("orientational coherence of the X-X axis after tracing the electron")
    ax.legend(fontsize=7)
    ax = axs[1]
    for M in [0.05, 1.0, 9.5256, 206.768283, 1.0e4]:
        d = npz("exact", M)
        if d is None:
            continue
        ax.semilogy(d["ls"], d["w"], ".-", label=f"M={M:g}")
    ax.set_xlabel("l")
    ax.set_ylabel("w_l")
    ax.legend(fontsize=7)
    ax.set_title("partial-wave weights (decay ~ l^-6, cusps)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig3_orientation.png"), dpi=150)

    fig, ax = plt.subplots(figsize=(6, 4.4))
    for k, c in [("S", "k"), ("E_N", "C1"), ("MI", "C4")]:
        ax.loglog(Ms, T[k + "_quad_err"], "o", color=c, ms=3, label=f"{k} quadrature")
        ax.loglog(Ms, T[k + "_basis_err"], "x", color=c, ms=4, label=f"{k} basis N-4")
    ax.loglog(Ms, T["tail_S"], "-", color="0.5", label="S l-tail estimate")
    ax.set_xlabel("M")
    ax.set_ylabel("|difference|")
    ax.legend(fontsize=7)
    ax.set_title("convergence")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig4_convergence.png"), dpi=150)


if __name__ == "__main__":
    main()
