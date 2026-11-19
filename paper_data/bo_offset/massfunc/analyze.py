"""Analysis of the computed mass functions (run after run.py).

    python3 paper_data/bo_offset/massfunc/analyze.py [--results FILE]

1. Plots m_v(R), m_r(R) (fig_massfunc.png).
2. Checks m_v(R_e), m_r(R_e) against the data-derived 0.3465 (kappa^6) and 0.0754 (J spacings).
3. Folds the functions with the threshold (zero-energy) weights of paper_data/bo_offset (coupling.py
   machinery; full weight densities exported to threshold_weights.csv) and predicts the
   nonadiabatic part of the BO+ad critical-mass offsets with NO fitted parameter:
     first order  (eq. 4 of THEORY.md)  and  nonperturbative (q = lam/(lam+m(r)) in the kinetic term).
4. Compares with the exact offsets (offset_table.csv) and writes analysis.txt, fig_fold.png.
"""
import argparse
import os
import sys

here = os.path.dirname(os.path.abspath(__file__))
parent = os.path.dirname(here)
sys.path.insert(0, parent)
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

import coupling as C   # parent folder: zero-energy coupling problem (dU, dD, weights, RB, RM, DHI)

R_E = 1.99719332                        # unit-mass R_e (bo_h2plus README)
REF = dict(m_v=(0.3465, 0.0003), m_r=(0.0754, 0.001))   # paper_data/bo_offset (kappa6.md, README 2b)
EFF = dict(m_v=0.4895, m_r=0.205)                         # effective values fitted to the offsets


# ------------------------------------------------------------------ interpolants
class MassFunction:
    """Cubic spline in R on the computed grid; power-law extrapolation outside."""

    def __init__(self, R, m, tail_power=None):
        o = np.argsort(R)
        self.R, self.m = R[o], m[o]
        self.sp = CubicSpline(self.R, self.m)
        r0, r1 = self.R[0], self.R[1]
        self.k0 = np.log(self.m[1] / self.m[0]) / np.log(r1 / r0)      # m ~ R^k0 below the grid
        # large-R: m - 1/2 ~ c R^-n, n from the last ~20% of the R range
        sel = self.R > 0.6 * self.R[-1]
        d = self.m[sel] - 0.5
        if tail_power is None and sel.sum() >= 3 and np.all(d != 0) and np.all(np.sign(d) == np.sign(d[-1])):
            self.n = -np.polyfit(np.log(self.R[sel]), np.log(np.abs(d)), 1)[0]
        else:
            self.n = tail_power if tail_power is not None else 3.0
        self.c = (self.m[-1] - 0.5) * self.R[-1] ** self.n

    def __call__(self, r, tail=True):
        r = np.asarray(r, float)
        out = np.empty_like(r)
        lo, hi = r < self.R[0], r > self.R[-1]
        mid = ~(lo | hi)
        out[mid] = self.sp(r[mid])
        out[lo] = self.m[0] * (r[lo] / self.R[0]) ** self.k0
        out[hi] = 0.5 + (self.c * r[hi] ** (-self.n) if tail else 0.0)
        return out


# ------------------------------------------------------------------ threshold weights
def threshold_states():
    res = pd.read_csv(os.path.join(parent, "coupling_results.csv"), comment="#")
    off = pd.read_csv(os.path.join(parent, "offset_table.csv"))
    rows = []
    for _, q in off.iterrows():
        L, v = int(q.L), int(q.v)
        ad = res[(res.L == L) & (res.c == 1)].reset_index(drop=True)
        a = ad.loc[v - 1 if L == 0 else 0]
        rows.append(dict(L=L, v=v, lam=a.lam_c, M_BOad=a.M_c, M_exact=q.M_exact, off_BOad=q.off_BOad,
                         f_rad=q.f_rad, f_cen=q.f_cen))
    return rows


def export_weights(states, path):
    """Weight densities w_rad = u'^2/(lam B), w_cen = L(L+1)u^2/(r^2 lam B) on a log grid."""
    rg = np.geomspace(1e-3, C.RM, 3000)
    out = []
    for s in states:
        I, r, u, up = C.weights(s["lam"], s["L"], 1)
        lB = s["lam"] * I["B"]
        wr = np.interp(rg, r, up ** 2 / lB)
        wc = np.interp(rg, r, s["L"] * (s["L"] + 1) * u ** 2 / r ** 2 / lB)
        out.append(pd.DataFrame(dict(L=s["L"], v=s["v"], M_c_BOad=s["M_BOad"], r=rg, w_rad=wr, w_cen=wc)))
        s["_grid"] = (r, up ** 2 / lB, s["L"] * (s["L"] + 1) * u ** 2 / r ** 2 / lB)
    df = pd.concat(out)
    with open(path, "w") as fh:
        fh.write("# zero-energy BO+ad threshold weights (paper_data/bo_offset/coupling.py); r = mu_e R (unit-mass bohr)\n"
                 "# first-order offset M_BOad - M_exact = int m_v w_rad dr + int m_r w_cen dr\n")
        df.to_csv(fh, index=False, float_format="%.10e")


def fold_first_order(s, mv, mr, tail=True):
    r, wr, wc = s["_grid"]
    return np.trapz(mv(r, tail) * wr, r) + np.trapz(mr(r, tail) * wc, r)


# ------------------------------------------------------------------ nonperturbative solve
def shoot_q(lam, L, mv, mr):
    """Zero-energy BO+ad equation with kinetic factors q = lam/(lam+m(r)); returns the growing amplitude."""
    def f(r, y):
        u, w = y
        qv = lam / (lam + mv(np.array([r]))[0])
        qr = lam / (lam + mr(np.array([r]))[0])
        return [w / qv, (qr * L * (L + 1) / r ** 2 + lam * C.dU(r) + C.dD(r)) * u]
    r0 = 1e-3
    a1 = lam / (2 * L + 2)
    y0 = [r0 ** (L + 1) * (1 + a1 * r0), (L + 1) * r0 ** L + (L + 2) * a1 * r0 ** (L + 1)]
    sol = solve_ivp(f, (r0, C.RM), y0, method="DOP853", rtol=1e-11, atol=1e-30)
    u, w = sol.y[0, -1], sol.y[1, -1]
    r = C.RM
    qinf = lam / (lam + 0.5)
    up = w / qinf
    if L == 0:
        b = np.sqrt((lam * 9 / 4 - C.DHI * C.RB ** 4) / qinf)
        th = np.arctan2(u / r, (u - r * up) / b)
        amp = np.hypot(u / r, (u - r * up) / b)
        A = amp * np.sin(th - b / r)
    else:
        A = (L * u + r * up) / ((2 * L + 1) * r ** (L + 1))
    return A / np.max(np.abs(sol.y[0]))


def lam_nonpert(s, mv, mr, guess_shift):
    lam0 = s["lam"] - guess_shift
    g = lambda x: shoot_q(x, s["L"], mv, mr)
    xs = np.linspace(lam0 - 0.15, lam0 + 0.15, 7)
    F = [g(x) for x in xs]
    roots = [brentq(g, xs[i], xs[i + 1], xtol=1e-10) for i in range(len(xs) - 1) if np.sign(F[i]) != np.sign(F[i + 1])]
    if not roots:
        return np.nan
    return min(roots, key=lambda x: abs(x - lam0))


# ------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", help="results CSV (default: results_fine.csv if present, else results.csv)")
    ap.add_argument("--no-nonpert", action="store_true", help="skip the nonperturbative check (~1 min)")
    ap.add_argument("--tag", default="", help="suffix for output files (tests)")
    a = ap.parse_args(argv)
    path = a.results
    if path is None:
        fine = os.path.join(here, "results_fine.csv")
        path = fine if os.path.exists(fine) else os.path.join(here, "results.csv")
    df = pd.read_csv(path).drop_duplicates("R").sort_values("R")
    R = df.R.values
    mv, mr = MassFunction(R, df.m_v.values), MassFunction(R, df.m_r.values)
    L = []
    P = L.append
    P(f"# massfunc analysis  (results: {os.path.relpath(path, here)}, {len(R)} R points in [{R.min()}, {R.max()}])")
    P(f"max basis-convergence indicator |level_hi - level_lo|: m_v {df.dm_v.max():.1e}, m_r {df.dm_r.max():.1e}")
    P(f"max |E0_basis - E0_exact| = {np.abs(df.dE0_basis).max():.1e}")
    P(f"extrapolation: below R={R[0]}: m ~ R^k with k_v={mv.k0:.3f}, k_r={mr.k0:.3f}; "
      f"above R={R[-1]}: m - 1/2 ~ c R^-n with n_v={mv.n:.2f} (c={mv.c:.4g}), n_r={mr.n:.2f} (c={mr.c:.4g})")
    P("")
    P("## Values at R_e = 1.99719332 (unit mass) vs data-derived values (no tuning)")
    for k, f in (("m_v", mv), ("m_r", mr)):
        val, ref, err = f(np.array([R_E]))[0], REF[k][0], REF[k][1]
        P(f"{k}(R_e) = {val:.6f}   data {ref} +- {err}   diff {val - ref:+.6f} ({(val - ref) / err:+.1f} sigma)")
    iv = np.argmax(df.m_v.values)
    P(f"m_v maximum {df.m_v.values[iv]:.5f} at R = {R[iv]:.3f};  m_r monotonic: {bool(np.all(np.diff(df.m_r.values) > -1e-12))}")
    for Rx in (0.5, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30):
        if R.min() <= Rx <= R.max():
            P(f"   R={Rx:5.1f}  m_v={mv(np.array([Rx]))[0]:.6f}  m_r={mr(np.array([Rx]))[0]:.6f}")
    P("")
    states = threshold_states()
    export_weights(states, os.path.join(here, f"threshold_weights{a.tag}.csv"))
    P("## Fold with the zero-energy BO+ad threshold weights: predicted nonadiabatic offset M_BOad - M_exact")
    P("state      lam_c     exact   pred_1st  (no tail) resid_1st  pred_nonpert resid_np  <m_v>_w  <m_r>_w")
    rows = []
    for s in states:
        r, wr, wc = s["_grid"]
        p1 = fold_first_order(s, mv, mr)
        p1nt = fold_first_order(s, mv, mr, tail=False)
        mvw = np.trapz(mv(r) * wr, r) / np.trapz(wr, r)
        mrw = np.trapz(mr(r) * wc, r) / np.trapz(wc, r) if s["L"] else np.nan
        pnp = np.nan
        if not a.no_nonpert:
            lam_np = lam_nonpert(s, mv, mr, p1)
            pnp = s["lam"] - lam_np
        name = f"L={s['L']} v={s['v']}" if s["L"] == 0 else f"L={s['L']}"
        rows.append(dict(state=name, L=s["L"], v=s["v"], lam_c=s["lam"], exact=s["off_BOad"], pred_first=p1,
                         pred_first_notail=p1nt, pred_nonpert=pnp, mv_w=mvw, mr_w=mrw))
        P(f"{name:9s} {s['lam']:8.4f}  {s['off_BOad']:.4f}   {p1:.4f}   ({p1nt:.4f})  {s['off_BOad'] - p1:+.4f}"
          f"     {pnp:.4f}     {s['off_BOad'] - pnp:+.4f}   {mvw:.4f}   {mrw:.4f}")
    P(f"(data-fitted effective values for comparison: <m_v> = {EFF['m_v']}, <m_r> = {EFF['m_r']})")
    P("")
    P("## Where the radial fold comes from (fraction of int m_v w_rad, by r bin)")
    edges = [0, 1.5, 3, 5, 8, 15, 40, 1e9]
    for s in states:
        r, wr, _ = s["_grid"]
        tot = np.trapz(mv(r) * wr, r)
        fr = [np.trapz((mv(r) * wr)[(r >= x) & (r < y)], r[(r >= x) & (r < y)]) / tot for x, y in zip(edges[:-1], edges[1:])]
        P(f"L={s['L']} v={s['v']}: " + " ".join(f"[{x},{y}):{f:.3f}" for x, y, f in zip(edges[:-1], edges[1:], fr)))
    txt = "\n".join(L) + "\n"
    print(txt)
    open(os.path.join(here, f"analysis{a.tag}.txt"), "w").write(txt)
    pd.DataFrame(rows).to_csv(os.path.join(here, f"fold_results{a.tag}.csv"), index=False)
    plots(df, mv, mr, states, rows, a.tag)


def plots(df, mv, mr, states, rows, tag):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(7, 7.5), sharex=True, gridspec_kw=dict(height_ratios=[2, 1]))
    rr = np.geomspace(max(df.R.min(), 0.05), df.R.max(), 800)
    ax[0].plot(rr, mv(rr), label="m_v(R)  (Σg response)")
    ax[0].plot(rr, mr(rr), label="m_r(R)  (Πg response)")
    ax[0].axhline(0.5, color="k", lw=0.6, ls="--")
    ax[0].errorbar([R_E], [REF["m_v"][0]], yerr=[REF["m_v"][1]], fmt="o", color="C0", ms=4, label="κ⁶ data m_v(R_e)")
    ax[0].errorbar([R_E], [REF["m_r"][0]], yerr=[REF["m_r"][1]], fmt="s", color="C1", ms=4, label="J-spacing data m_r(R_e)")
    ax[0].set_ylabel("extra mass (electron masses)")
    ax[0].legend(fontsize=8)
    ax[0].set_xscale("log")
    for s in states:
        r, wr, wc = s["_grid"]
        if s["L"] == 0:
            m = (r > 0.05) & (r < df.R.max())
            ax[1].plot(r[m], (wr * r)[m] / np.trapz(wr, r), label=f"w_rad r, L=0 v={s['v']}")
    ax[1].set_xlabel("R = μ_e R_phys (bohr)")
    ax[1].set_ylabel("weight density × r")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(here, f"fig_massfunc{tag}.png"), dpi=130)
    fig, ax = plt.subplots(figsize=(6, 3.8))
    x = np.arange(len(rows))
    ax.bar(x - 0.27, [q["exact"] for q in rows], 0.27, label="exact (BO+ad − exact)")
    ax.bar(x, [q["pred_first"] for q in rows], 0.27, label="predicted, 1st order")
    ax.bar(x + 0.27, [q["pred_nonpert"] for q in rows], 0.27, label="predicted, nonpert.")
    ax.set_xticks(x, [q["state"] for q in rows])
    ax.set_ylabel("ΔM_c (electron masses)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(here, f"fig_fold{tag}.png"), dpi=130)


if __name__ == "__main__":
    main()
