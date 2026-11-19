"""Analysis and figures for paper_data/fidelity (reads results.csv, dense.csv; writes fits.txt, fig*.png)."""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
MARKERS = [(0.25, "Lindemann max"), (0.4924, "KLC g min"), (0.5, "E/mu max"), (0.543, "Kais-Shi"),
           (0.635, "CM density"), (1.11, "T_X max")]
DD_HMINUS = 1.367546381   # ||D psi||^2 of the infinite-mass H^- -like state (limit_small_M.py, N=28)


def load(fn):
    rows = list(csv.DictReader(open(os.path.join(HERE, fn))))
    out = {}
    for k in rows[0]:
        try:
            out[k] = np.array([float(r[k]) if r[k] != "" else np.nan for r in rows])
        except ValueError:
            out[k] = np.array([r[k] for r in rows])
    o = np.argsort(out["M"])
    return {k: v[o] for k, v in out.items()}


def err(d, k):
    return np.abs(d[k] - d[k + "_N4"])


def local_extrema(x, y):
    i = np.where((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:]))[0] + 1
    j = np.where((y[1:-1] < y[:-2]) & (y[1:-1] < y[2:]))[0] + 1
    res = []
    for kind, idx in (("max", i), ("min", j)):
        for n in idx:   # parabolic refinement in ln M
            xs, ys = x[n - 1:n + 2], y[n - 1:n + 2]
            p = np.polyfit(xs, ys, 2)
            res.append((kind, np.exp(-p[1] / (2 * p[0])), np.polyval(p, -p[1] / (2 * p[0])), p[0]))
    return res


def main():
    g = load("results.csv")
    out = []
    M = g["M"]; L = np.log(M)
    out.append("# fidelity susceptibilities, variable ln M; chi = <dpsi|dpsi> - |<psi|dpsi>|^2, F = 1 - chi h^2/2")
    out.append("M  chi_F  chi_inv  chi_sep  |N-(N-4)| for each  FDrel(chi_F)  scal_chi_inv  removed_inv  gram_corr")
    for i in range(len(M)):
        fdrel = (g["chiF_fd_rich"][i] / g["chiF"][i] - 1)
        out.append(f"{M[i]:<12g} {g['chiF'][i]:.9f} {g['chi_inv'][i]:.9e} {g['chi_sep'][i]:.9e}  "
                   f"{err(g,'chiF')[i]:.1e} {err(g,'chi_inv')[i]:.1e} {err(g,'chi_sep')[i]:.1e}  {fdrel:.1e} "
                   f"{g['scal_chi_inv'][i]:.1e} {g['chiF'][i]-g['chi_inv'][i]:.6f} {g['gram_corr'][i]:.4f}")
    # extrema on the grid
    for k in ("chiF", "chi_inv", "chi_sep"):
        out.append(f"grid extrema of {k}: {local_extrema(L, g[k])}")
    # large-M behaviour
    big = M >= 100
    kap = M ** -0.25
    for k in ("chiF", "chi_inv", "chi_sep"):
        dev = np.abs(g[k] - 1 / 32)
        sl = np.diff(np.log(dev[big])) / np.diff(np.log(kap[big]))
        out.append(f"large M: {k} - 1/32 = {(g[k]-1/32)[big]} ; local slopes in kappa {sl}")
        p = np.polyfit(kap[M >= 500] ** 2, g[k][M >= 500], 2)
        out.append(f"   fit chi = c0 + c1 kappa^2 + c2 kappa^4 (M>=500): c0 = {p[2]:.6f}, c1 = {p[1]:.4f}, c2 = {p[0]:.3f}")
    # small M
    sm = M <= 0.2
    out.append(f"small M: chi_F (M+1)^2 / ||D psi_inf||^2 = {g['chiF'][sm]*(M[sm]+1)**2/DD_HMINUS}")
    for k in ("chi_inv", "chi_sep"):
        out.append(f"small M: log-slope of {k} vs M: {np.diff(np.log(g[k][M <= 0.3]))/np.diff(L[M <= 0.3])}")
    # dense windows
    if os.path.exists(os.path.join(HERE, "dense.csv")):
        d = load("dense.csv")
        win = np.array([k.split("_")[0] for k in d["key"]])
        okw = [w for w in set(win) if (win == w).sum() >= 6]
        keep = np.isin(win, okw)
        d = {k: v[keep] for k, v in d.items()}
        win = win[keep]
        out.append("\n# dense windows (fixed alphas per window: analytic in M for each basis)")
        for w in sorted(set(win)):
            s = win == w
            x, Ld = d["M"][s], np.log(d["M"][s])
            for k in ("chiF", "chi_inv", "chi_sep"):
                y = d[k][s]
                e = err(d, k)[s]
                # smoothness: residual of a local quartic fit vs basis error
                p = np.polyfit(Ld, y, 6)
                res = y - np.polyval(p, Ld)
                d2 = np.diff(y, 2) / np.diff(Ld)[0] ** 2
                out.append(f"{w} {k}: M {x[0]:.3f}-{x[-1]:.3f}, n={len(x)}, extrema {local_extrema(Ld, y)}, "
                           f"max|deg-6 fit resid| {np.max(np.abs(res)):.1e}, max basis err {np.max(e):.1e}, "
                           f"d2/dlnM2 range [{d2.min():.4g}, {d2.max():.4g}]")

        # consistency between windows (different fixed alphas) and with the grid (own alphas)
        from scipy.interpolate import CubicSpline
        wins = sorted(set(win))
        for k in ("chiF", "chi_inv", "chi_sep"):
            msgs = []
            for a_, b_ in zip(wins[:-1], wins[1:]):
                sa, sb = win == a_, win == b_
                lo, hi = max(d["M"][sa].min(), d["M"][sb].min()), min(d["M"][sa].max(), d["M"][sb].max())
                if hi <= lo:
                    continue
                xx = np.geomspace(lo, hi, 20)
                fa = CubicSpline(np.log(d["M"][sa]), d[k][sa])(np.log(xx))
                fb = CubicSpline(np.log(d["M"][sb]), d[k][sb])(np.log(xx))
                msgs.append(f"{a_}/{b_} [{lo:.3f},{hi:.3f}] max rel diff {np.max(np.abs(fa/fb-1)):.1e}")
            for w in wins:
                sw = win == w
                sg = (g["M"] > d["M"][sw].min()) & (g["M"] < d["M"][sw].max())
                if sg.any():
                    f = CubicSpline(np.log(d["M"][sw]), d[k][sw])(np.log(g["M"][sg]))
                    msgs.append(f"{w} vs grid: max rel diff {np.max(np.abs(f/g[k][sg]-1)):.1e}")
            out.append(f"consistency {k}: " + "; ".join(msgs))
        # extrema from splines of each window
        for w in wins:
            sw = win == w
            for k in ("chi_inv", "chi_sep", "chiF"):
                cs = CubicSpline(np.log(d["M"][sw]), d[k][sw])
                r = [x for x in cs.derivative().roots(extrapolate=False)]
                ri = [x for x in cs.derivative(2).roots(extrapolate=False)
                      if np.log(d["M"][sw]).min() + 0.05 < x < np.log(d["M"][sw]).max() - 0.05]
                if ri:
                    out.append(f"{w} spline inflections (in lnM) of {k}: " + ", ".join(f"M={np.exp(x):.4f} val={cs(x):.8f} slope={cs(x,1):.5f}" for x in ri))
                if r:
                    out.append(f"{w} spline extrema of {k}: " + ", ".join(f"M={np.exp(x):.4f} val={cs(x):.8f} curv={cs(x,2):.4e}" for x in r))
    open(os.path.join(HERE, "fits.txt"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))
    figures(g, d if os.path.exists(os.path.join(HERE, "dense.csv")) else None, win)


def figures(g, d, win=None):
    M = g["M"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    a = ax[0]
    a.loglog(M, g["chiF"], "o-", label=r"$\chi_F$ (fixed $R,\rho$)")
    a.loglog(M, g["chi_inv"], "s-", label=r"$\chi_{\rm inv}$ (global dilation removed)")
    a.loglog(M, g["chi_sep"], "^-", label=r"$\chi_{\rm sep}$ ($R$, $\rho$ dilations removed)")
    mm = np.geomspace(0.05, 1e4, 200)
    a.loglog(mm, DD_HMINUS / (mm + 1) ** 2, "k:", lw=1, label=r"$\|D\psi_{\infty}\|^2/(M+1)^2$ (small-$M$ limit)")
    a.axhline(1 / 32, color="gray", ls="--", lw=1, label="1/32 (BO vibrational width)")
    for m, _ in MARKERS:
        a.axvline(m, color="r", lw=0.5, alpha=0.5)
    a.set_xlabel("M"); a.set_ylabel(r"$\chi$ (per $(\ln M)^2$)"); a.legend(fontsize=7); a.set_ylim(1e-4, 3)
    a.set_title("(a) fidelity susceptibilities, mass_scans grid")
    a = ax[1]
    a.semilogx(M, g["chi_inv"], "s-", label=r"$\chi_{\rm inv}$")
    a.semilogx(M, g["chi_sep"], "^-", label=r"$\chi_{\rm sep}$")
    a.semilogx(M, g["chiF"], "o-", label=r"$\chi_F$")
    a.axhline(1 / 32, color="gray", ls="--", lw=1)
    for m, lab in MARKERS:
        a.axvline(m, color="r", lw=0.5, alpha=0.5)
    a.set_ylim(0, 0.08); a.set_xlabel("M"); a.legend(fontsize=7); a.set_title("(b) linear scale (red: crossover markers)")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig1_chi_vs_M.png"), dpi=150); plt.close(fig)

    if d is not None:
        fig, ax = plt.subplots(2, 3, figsize=(14, 7))
        for j, k in enumerate(("chiF", "chi_inv", "chi_sep")):
            for w in sorted(set(win)):
                s = win == w
                x, y = d["M"][s], d[k][s]
                ax[0, j].plot(x, y, ".-", label=w)
                L = np.log(x)
                d1 = np.gradient(y, L); d2 = np.gradient(d1, L)
                ax[1, j].plot(x[2:-2], d2[2:-2], ".-", label=w)  # edge points of np.gradient dropped
            s = (g["M"] >= 0.15) & (g["M"] <= 1.4)
            ax[0, j].plot(g["M"][s], g[k][s], "kx", label="grid (own basis)")
            for m, _ in MARKERS:
                ax[0, j].axvline(m, color="r", lw=0.5, alpha=0.5); ax[1, j].axvline(m, color="r", lw=0.5, alpha=0.5)
            ax[0, j].set_xscale("log"); ax[1, j].set_xscale("log")
            ax[0, j].set_title(k); ax[1, j].set_title(f"d²{k}/d(lnM)²"); ax[1, j].set_xlabel("M")
            ax[0, j].legend(fontsize=7)
        fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig2_dense_windows.png"), dpi=150); plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    big = M >= 20
    kap = M ** -0.25
    for k, mk in (("chiF", "o"), ("chi_inv", "s"), ("chi_sep", "^")):
        ax[0].loglog(kap[big], np.abs(g[k][big] - 1 / 32), mk + "-", label=f"|{k} - 1/32|")
    ax[0].loglog(kap[big], 0.3 * kap[big] ** 2, "k:", label=r"$\propto\kappa^2$")
    ax[0].set_xlabel(r"$\kappa=M^{-1/4}$"); ax[0].legend(fontsize=7); ax[0].set_title("(a) approach to the BO value 1/32")
    for k, mk in (("chiF", "o"), ("chi_inv", "s"), ("chi_sep", "^")):
        ax[1].loglog(M, np.maximum(err(g, k) / np.abs(g[k]), 1e-16), mk + "-", label=f"{k}: |N - (N-4)|/value")
    ax[1].loglog(M, np.maximum(np.abs(g["chiF_fd_rich"] / g["chiF"] - 1), 1e-16), "x--", label="chi_F: FD Richardson vs perturbative")
    ax[1].loglog(M, np.maximum(np.abs(g["scal_chi_inv"]) / g["chi_inv"], 1e-16), "d--", label="uniform-scaling test chi_inv / chi_inv")
    ax[1].set_xlabel("M"); ax[1].legend(fontsize=7); ax[1].set_title("(b) convergence indicators (relative)")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig3_largeM_convergence.png"), dpi=150); plt.close(fig)


if __name__ == "__main__":
    main()
