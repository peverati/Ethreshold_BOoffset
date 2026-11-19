"""Final figures for Paper 1 (PCCP).  Run:  python3 make_figures.py [F1 F2 ...]

Reads existing data files only (see README.md for sources); no physics is recomputed.
Shared style: PCCP single column 8.3 cm / double column 17.1 cm, sans-serif 7 pt, vector PDF and
300 dpi PNG, Okabe-Ito colour-blind-safe palette.  Symbol convention: exact = filled circles,
BO = dashed lines, BO+adiabatic = dotted lines.
"""
import os, sys, json, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from summary_values import (MUONIC_LEVELS, p1, ROOT, CRIT, MUONIC, LEVEL_THRESHOLDS, MARKERS, F_LEVELS, C6, C4,
                            S_INF, MV_RE, MR_RE, KUTZELNIGG, V4_PRED_MC, LV_PRED_MC, L1S, FW_11,
                            LEVEL_KEY, RUN5, RUN6, MAP_UNC, UNEQ_LENIENT, CRIT_JV)

OUT = os.path.dirname(os.path.abspath(__file__))
CM = 1 / 2.54
W1, W2 = 8.3 * CM, 17.1 * CM

# Okabe-Ito
OI = dict(black="#000000", orange="#E69F00", sky="#56B4E9", green="#009E73", yellow="#F0E442",
          blue="#0072B2", verm="#D55E00", purple="#CC79A7", grey="#999999")
LEVEL_COL = {"v0": OI["black"], "v1": OI["blue"], "v2": OI["green"], "v3": OI["orange"],
             "v4": OI["verm"], "L1": OI["purple"], "L2": OI["sky"],
             # thresholds_Lv levels (Okabe-Ito exhausted; distinct marker shapes as second cue)
             "J1v1": "#882255", "L3": "#8C7A00", "L4": "#555555", "L1S": "#9A9A9A"}
LEVEL_MK = {"J1v1": "D", "L2": "v", "L3": "s", "L4": "^", "L1S": "o"}

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7, "legend.fontsize": 6,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "lines.linewidth": 1.0, "lines.markersize": 3.2,
    "legend.frameon": False, "legend.handlelength": 1.8, "mathtext.fontset": "dejavusans",
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300,
})
EXACT = dict(marker="o", ls="none", ms=3.0, mew=0)            # exact: filled circles
BO_LS, BOAD_LS = (0, (4, 2)), (0, (1, 1.3))                    # BO dashed, BO+ad dotted


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("wrote", name)


def panel(ax, s, x=-0.16, y=1.0):
    ax.text(x, y, s, transform=ax.transAxes, fontweight="bold", fontsize=8, va="bottom", ha="left")


# ------------------------------------------------------------------ data loaders
def load_ground():
    e = pd.read_csv(p1("mass_scans", "energies.csv")).sort_values("M").drop_duplicates("M")
    return e


def load_levels():
    """Binding energies E_th - E of every bound level we have, with convergence filter
    (E_err < 0.1 binding).  Returns dict key -> DataFrame(M, Eb)."""
    out = {}
    e = load_ground()
    out["v0"] = pd.DataFrame(dict(M=e.M, Eb=e.binding, err=e.E0_err))
    # v = 1: mass-scan root 1 (prefer the re-optimised E1_opt), rotational scan L=0 root 1,
    #        near-threshold critical-mass runs
    x = pd.read_csv(p1("mass_scans", "exact_L0.csv")).sort_values("M").drop_duplicates("M")
    E1 = np.where(np.isfinite(x.E1_opt_ext), x.E1_opt_ext, x.E1_ext)
    E1e = np.where(np.isfinite(x.E1_opt_ext), x.E1_opt_err, x.E1_err)
    d1 = pd.DataFrame(dict(M=x.M, Eb=x.E_th - E1, err=E1e))
    sc = pd.read_csv(p1("rotational", "scan_L.csv"))
    s0 = sc[(sc.L == 0)]
    d1 = pd.concat([d1, pd.DataFrame(dict(M=s0.M, Eb=s0.E_th - s0.E1_ext, err=s0.E1_err))])
    ob = pd.read_csv(p1("critical_masses", "observables_state1.csv"))
    d1 = pd.concat([d1, pd.DataFrame(dict(M=ob.M, Eb=ob.binding, err=0.0))])
    out["v1"] = d1
    for k, f in (("v2", ("critical_masses", "observables_state2.csv")),
                 ("v3", ("critical_masses", "observables_state3.csv")),
                 ("v4", ("critical_masses_v4", "observables_state4.csv"))):
        ob = pd.read_csv(p1(*f))
        out[k] = pd.DataFrame(dict(M=ob.M, Eb=ob.binding, err=0.0))
    r = pd.read_csv(p1("mass_scans", "exact_rot.csv"))
    for k, L, sym, js in (("L1", 1, "A", "critical_L1A.json"), ("L2", 2, "S", "critical_L2S.json")):
        rr = r[(r.L == L) & (r.sym == sym)]
        ss = sc[(sc.L == L) & (sc.sym == sym)]
        b = json.load(open(p1("rotational", js)))["binding"]
        out[k] = pd.concat([pd.DataFrame(dict(M=rr.M, Eb=rr.E_th - rr.E_ext, err=rr.E_err)),
                            pd.DataFrame(dict(M=ss.M, Eb=ss.E_th - ss.E0_ext, err=ss.E0_err)),
                            pd.DataFrame(dict(M=[q[0] for q in b], Eb=[q[1] for q in b], err=0.0))])
    # thresholds_Lv: near-threshold binding energies (largest basis), x = M/M* - 1 <= 0.2
    for k, f in (("J1v1", "J1v1"), ("L3", "L3A"), ("L4", "L4S"), ("L1S", "L1S")):
        bd = pd.read_csv(p1("thresholds_Lv", f"{f}_binding.csv"))
        out[k] = pd.DataFrame(dict(M=bd.M, Eb=bd.binding, err=0.0))
    mc = {c["key"]: c["Mc"] for c in CRIT}
    mc["L1S"] = L1S["Mc"]
    for k, d in out.items():
        d = d[(d.Eb > 0) & (d.err < 0.1 * d.Eb)]
        if k in mc:
            d = d[d.M > mc[k]]
        d = d[d.Eb > 1e-12] if k == "v4" else d   # v=4 rows at the 1e-14 rounding floor removed
        out[k] = d.sort_values("M").drop_duplicates("M").reset_index(drop=True)
    return out


def crossover_f():
    """width.py logic: f = -4 dlnL/dlnM between neighbouring scan masses (geometric midpoints)."""
    d = pd.read_csv(p1("mass_scans", "descriptors.csv")).sort_values("M").drop_duplicates("M")
    M, L = d.M.values, d.lindemann.values
    x = np.log(M)
    xm = (x[1:] + x[:-1]) / 2
    f = -4 * np.diff(np.log(L)) / np.diff(x)
    def cross(level):
        i = np.where((f[:-1] - level) * (f[1:] - level) <= 0)[0][0]
        return float(np.exp(xm[i] + (level - f[i]) * (xm[i + 1] - xm[i]) / (f[i + 1] - f[i])))
    return np.exp(xm), f, {lev: cross(lev) for lev in (0.1, 0.5, 0.9)}, M, L


# ------------------------------------------------------------------ F1
def fig1():
    lv = load_levels()
    fig = plt.figure(figsize=(W2, 13.2 * CM))
    gsp = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.82], hspace=0.3)
    ax = fig.add_subplot(gsp[0]); mp = fig.add_subplot(gsp[1])
    names = {"v0": r"$v=0$", "v1": r"$v=1$", "v2": r"$v=2$", "v3": r"$v=3$", "v4": r"$v=4$",
             "L1": r"$L=1$ (antisym.)", "L2": r"$L=2$ (sym.)", "J1v1": r"$L=1$, $v=1$ (antisym.)",
             "L3": r"$L=3$ (antisym.)", "L4": r"$L=4$ (sym.)", "L1S": r"$L=1$ (sym., expl.)"}
    for k in ("v0", "L1", "v1", "L2", "J1v1", "L3", "v2", "L4", "v3", "v4", "L1S"):
        d = lv[k]
        ax.plot(d.M, d.Eb, "-", color=LEVEL_COL[k], lw=0.7)
        st = dict(EXACT, marker=LEVEL_MK.get(k, "o"))
        if k == "L1S":
            st.update(mfc="white", mec=LEVEL_COL[k], mew=0.6)
        ax.plot(d.M, d.Eb, color=LEVEL_COL[k], label=names[k], **st)
    # run 5 (threshold_map): 14 new levels near threshold, x = M/M* - 1 <= 0.2 (data/J*v*_binding.csv;
    # (7,0) from threshold_map_L7hp/data, double-double assembly)
    for c in RUN5:
        bd = pd.read_csv(p1("threshold_map_L7hp" if c["hp"] else "threshold_map", "data",
                            f"{c['key']}_binding.csv")).sort_values("M")
        ax.plot(bd.M, bd.binding, "-", color="0.55", lw=0.5)
        ax.plot(bd.M, bd.binding, marker="o", ms=1.6, mew=0, ls="none", color="0.55",
                label="threshold map" if c is RUN5[0] else None)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.04, 1.3e4); ax.set_ylim(1e-11, 0.3)
    ax.set_xlabel(r"mass ratio $M=m_X/m_e$")
    ax.set_ylabel(r"binding energy $E_{\rm th}-E$ ($E_h$)")
    trans = ax.get_xaxis_transform()
    for c in CRIT + [L1S]:
        expl = c["key"] == "L1S"
        ax.axvline(c["Mc"], color=LEVEL_COL[c["key"]], lw=0.5, alpha=0.8, zorder=0,
                   ls=BO_LS if expl else "-")
        left = c["key"] in ("v1", "v2", "L1S")
        ylab = {"J1v1": 0.2, "v2": 0.13, "L4": 0.2, "L1S": 0.62}.get(c["key"], 0.36)   # staggered to avoid overlaps
        ax.text(c["Mc"] * (0.95 if left else 1.05), ylab,
                f"$M_c$ = {c['Mc']:.2f}" + (" (expl.)" if expl else ""), transform=trans,
                rotation=90, fontsize=5.5, va="center", ha="right" if left else "left",
                color=LEVEL_COL[c["key"]])
    for m in MUONIC:
        ax.plot([m["M"]], [1.0], marker="v", ms=4, color=OI["grey"], transform=trans, clip_on=False)
        ax.text(m["M"], 1.035, m["name"], transform=trans, fontsize=5.5, ha="center", va="bottom",
                color="0.35")
    # 2026-11-19: legend in the empty area under the (0,0) curve, with the line/marker explanations as entries
    h, l = ax.get_legend_handles_labels()
    h += [Line2D([], [], color="0.3", lw=0.6), Line2D([], [], color="0.5", lw=0.6, ls=BO_LS),
          Line2D([], [], marker="v", ms=4, color=OI["grey"], ls="none")]
    l += ["exact $M_c$", "exact $M_c$ (expl.)", "muonic molecules"]
    ax.legend(h, l, loc="upper left", bbox_to_anchor=(0.0, 0.51), ncol=2, borderaxespad=0.3, fontsize=5.5,
              columnspacing=0.5, handletextpad=0.3, handlelength=1.4, frameon=False)
    # (b) bound-state map: every (J, v) level with M_c <= ~100 as a line starting at its threshold
    levs = [dict(J=CRIT_JV[c["key"]][0], v=CRIT_JV[c["key"]][1], Mc=c["Mc"], BO=c["BO"], BOad=c["BOad"], st="exact")
            for c in CRIT] + \
           [dict(J=c["J"], v=c["v"], Mc=c["Mc"], BO=c["BO"], BOad=c["BOad"], st="bound" if c["onesided"] else "exact")
            for c in RUN5]
    XMAX = 103
    vcol = {0: OI["black"], 1: OI["blue"], 2: OI["green"], 3: OI["orange"], 4: OI["verm"]}
    dy = 0.17
    mp.plot([0, XMAX], [0 - 0.25, 0 - 0.25], "-", color=vcol[0], lw=0.9)          # (0,0): bound at every M
    mp.text(1.0, -0.25 + 0.05, "(0,0)", fontsize=4.6, va="bottom", color=vcol[0])
    for d in levs:
        y = d["J"] + dy * d["v"] - (0.25 if d["J"] == 0 else 0)
        col = vcol[d["v"]]
        pred = d["st"] == "pred"
        mp.plot([d["Mc"], XMAX], [y, y], color=col, lw=0.9, ls=(0, (2, 1.5)) if pred else "-", alpha=0.6 if pred else 1)
        mp.plot([d["Mc"]], [y], marker="o", ms=3.0, mew=0.7 if pred else 0, mfc="white" if pred else col, mec=col,
                ls="none", zorder=4)
        mp.plot([d["BOad"]], [y], marker="|", ms=4.0, mew=0.8, color=col, ls="none", zorder=3)
        mp.plot([d["BO"]], [y], marker="|", ms=4.0, mew=0.8, color=col, alpha=0.45, ls="none", zorder=3)
        mp.text(d["Mc"] - 0.6, y, f"({d['J']},{d['v']})" + (" predicted only; exact not computed" if pred else ""),
                fontsize=4.3, ha="right", va="center", color="0.35" if pred else col)
    for m in MUONIC:
        mp.axvline(m["M"], color=OI["grey"], lw=0.5, ls=(0, (1, 1.2)), zorder=0)
        mp.text(m["M"], 7.75, m["name"], fontsize=5.0, ha="center", va="bottom", color="0.35")
    mp.set_xlim(0, XMAX); mp.set_ylim(-0.6, 7.75)
    mp.set_yticks(range(8)); mp.set_yticklabels([f"$J={J}$" for J in range(8)], fontsize=6)
    mp.tick_params(axis="y", which="both", right=False, length=0)
    mp.set_xlabel(r"mass ratio $M=m_X/m_e$")
    mp.legend(handles=[Line2D([], [], marker="o", ms=3, mew=0, color="k", ls="none", label="exact $M_c$ (level bound to the right)"),
                       Line2D([], [], marker="|", ms=4, mew=0.8, color="k", ls="none", label="BO+ad $M_c$"),
                       Line2D([], [], marker="|", ms=4, mew=0.8, color="k", alpha=0.45, ls="none", label="BO $M_c$")] +
                      [Line2D([], [], color=vcol[v], lw=1.2, label=f"$v={v}$") for v in range(5)],
              loc="upper left", ncol=3, fontsize=5.0, columnspacing=0.8, handletextpad=0.3, borderaxespad=0.3,
              frameon=True, framealpha=0.9, edgecolor="none")
    panel(ax, "a", x=-0.06); panel(mp, "b", x=-0.06)
    save(fig, "F1_levels")


# ------------------------------------------------------------------ F2 (paper Fig. 4): threshold laws
# 2026-11-19: split; the halo panel is now its own figure (F2h_halo, paper Fig. 5). Labels use (J,v) and k_b.
def _halo_data():
    out = []
    for c in CRIT[:4]:
        k = c["key"]
        ob = pd.read_csv(p1(*(("critical_masses_v4", "observables_state4.csv") if k == "v4"
                              else ("critical_masses", f"observables_state{k[1]}.csv"))))
        ob = ob[(ob.M > c["Mc"]) & (ob.binding > 1e-12)]      # near-threshold runs only
        out.append((k, c, ob))
    return out


def fig2():
    fig = plt.figure(figsize=(W1, 11.0 * CM))
    gsp = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.95], hspace=0.38)
    a = fig.add_subplot(gsp[0]); c3 = fig.add_subplot(gsp[1])
    fitinfo = {}
    for k, c, ob in _halo_data():
        x = ob.M.values - c["Mc"]
        y = np.sqrt(ob.binding.values)
        col = LEVEL_COL[k]
        a.plot(x, y, color=col, label=k.replace("v", "$v=") + "$", **EXACT)
        sel = (x > 0) & (x <= 0.31)
        A = np.vstack([x[sel], x[sel] ** 2]).T
        coef = np.linalg.lstsq(A, y[sel], rcond=None)[0]
        fitinfo[k] = (coef, int(sel.sum()))
        xx = np.logspace(np.log10(x.min()) - 0.15, np.log10(x.max()) + 0.1, 200)
        a.plot(xx, coef[0] * xx + coef[1] * xx ** 2, "-", color=col, lw=0.6)
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel(r"$x=M-M_c$")
    a.set_ylabel(r"$\sqrt{E_b}$ ($E_h^{1/2}$)")
    a.legend(loc="upper left")
    # (fit note removed 2026-11-19: the caption says it)
    # (b) J > 0 levels: E_b/(M - M_c) vs x; flat = linear threshold law E_b ~ B (M - M_c)
    by = {c["key"]: c for c in CRIT}
    nm = {"J1v1": r"(1,1)", "L2": r"(2,0)", "L3": r"(3,0)", "L4": r"(4,0)"}
    for k, f in (("J1v1", "J1v1"), ("L2", None), ("L3", "L3A"), ("L4", "L4S")):
        if f is None:   # (2,0): binding points only in the text log of the 2026-11-13 run (paper1/rotational)
            rows = [ln.split() for ln in open(p1("rotational", "critical_L2S.txt"))
                    if ln.strip().startswith("M=") and "E_th-E0=" in ln]
            bd = pd.DataFrame(dict(M=[float(r[0][2:]) for r in rows], binding=[float(r[2].split("=")[1]) for r in rows]))
        else:
            bd = pd.read_csv(p1("thresholds_Lv", f"{f}_binding.csv"))
        x = bd.M.values - by[k]["Mc"]
        y = bd.binding.values / x * 1e3
        st = dict(EXACT, marker=LEVEL_MK[k])
        c3.plot(x, y, "-", color=LEVEL_COL[k], lw=0.6)
        c3.plot(x, y, color=LEVEL_COL[k], label=nm[k], **st)
        fitinfo[k] = (float(y[0]), float(y[-1]))
    c3.set_xscale("log")
    c3.set_ylim(0, 3.8)
    c3.set_xlabel(r"$x=M-M_c$")
    c3.set_ylabel(r"$E_b/(M-M_c)$ ($10^{-3}\,E_h$)")
    c3.legend(loc="upper left", bbox_to_anchor=(0.0, 0.96), ncol=2, fontsize=5.5, columnspacing=1.0)
    # (explanatory note removed 2026-11-19: the caption says it)
    panel(a, "a", x=-0.2); panel(c3, "b", x=-0.2)
    save(fig, "F2_threshold")
    return fitinfo


# ------------------------------------------------------------------ F2h (paper Fig. 5): halo size
def fig2h():
    fig, b = plt.subplots(figsize=(W1, 5.6 * CM))
    for k, c, ob in _halo_data():
        col = LEVEL_COL[k]
        b.plot(ob.M - c["Mc"], ob.r01 / ob.one_over_2kappa, color=col, label=k.replace("v", "$v=") + "$", **EXACT)
        b.plot(ob.M - c["Mc"], ob.r01 / ob.one_over_2kappa, "-", color=col, lw=0.6)
    b.axhline(1.0, color="0.5", lw=0.5)
    b.set_xscale("log")
    b.set_ylim(0.95, 1.6)
    b.set_xlim(6e-4, 2.0)
    b.set_xlabel(r"$x=M-M_c$")
    b.set_ylabel(r"$2k_b\langle r_{XX}\rangle$")
    b.legend(loc="upper left", fontsize=5.5)
    # (halo note removed 2026-11-19: the caption says it)
    save(fig, "F2h_halo")


# ------------------------------------------------------------------ F3 (central): non-adiabatic offsets
J_COL = {0: OI["black"], 1: OI["blue"], 2: OI["green"], 3: OI["orange"], 4: OI["verm"], 5: OI["purple"],
         6: OI["sky"], 7: OI["grey"]}
ROLE_MK = {"calib": "o", "blind1": "*", "blind5": "D"}
ROLE_MS = {"calib": 3.6, "blind1": 6.5, "blind5": 3.4}


def offset_levels():
    """All symmetric-family levels with exact M_c: list of dicts (J, v, Mc, err, onesided, off, pred_off, diff, role)."""
    out = []
    for c in CRIT:
        J, v = CRIT_JV[c["key"]]
        role = "blind1" if c["prespec"] else "calib"
        out.append(dict(J=J, v=v, Mc=c["Mc"], err=c["err"], onesided=False, off=c["off"], pred_off=c["pred"],
                        diff=c["pred"] - c["off"], role=role))
    for c in RUN5:
        out.append(dict(J=c["J"], v=c["v"], Mc=c["Mc"], err=c["err"], onesided=c["onesided"], off=c["off"],
                        pred_off=c["pred_off"], diff=c["diff"], role="blind5"))
    return sorted(out, key=lambda d: d["Mc"])


def fig3():
    lv = offset_levels()
    fig = plt.figure(figsize=(W2, 11.5 * CM))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.45, 1.0], height_ratios=[1.25, 1.0], hspace=0.07, wspace=0.27)
    a = fig.add_subplot(gs[0, 0]); b = fig.add_subplot(gs[1, 0], sharex=a); c = fig.add_subplot(gs[:, 1])
    plt.setp(a.get_xticklabels(), visible=False)
    # (a) offsets vs exact M_c
    lab_dy = {(0, 1): 0.012, (1, 0): -0.03}
    for d in lv:
        col = J_COL[d["J"]]
        a.errorbar([d["Mc"] + 0.9], [d["pred_off"]], yerr=MAP_UNC, fmt="none", ecolor=col, elinewidth=0.5, capsize=1.2,
                   capthick=0.5, alpha=0.75, zorder=2)
        a.plot([d["Mc"] + 0.9], [d["pred_off"]], marker="_", ms=4, mew=0.9, color=col, ls="none", zorder=3)
        a.plot([d["Mc"]], [d["off"]], marker=ROLE_MK[d["role"]], ms=ROLE_MS[d["role"]], mew=0, color=col, ls="none",
               zorder=4)
        a.text(d["Mc"] - 0.8, d["off"] + lab_dy.get((d["J"], d["v"]), 0.0), f"({d['J']},{d['v']})", fontsize=4.3,
               ha="right", va="center", color=col)
    a.set_ylim(0.15, 0.6)
    a.set_ylabel(r"offset $M_c^{\rm BO+ad}-M_c$")
    hj = [Line2D([], [], color=J_COL[J], marker="s", ms=3.2, mew=0, ls="none", label=f"$J={J}$") for J in range(8)]
    a.legend(handles=hj, loc="upper right", ncol=4, fontsize=5.0, columnspacing=0.6, handletextpad=0.1,
             borderaxespad=0.3)
    a.text(0.015, 0.03, "filled: exact; bar: parameter-free prediction $\\pm$%.4f" % MAP_UNC, transform=a.transAxes,
           fontsize=5.2, va="bottom")
    # (b) residual exact - predicted M_c with the frozen band
    b.axhspan(-MAP_UNC, MAP_UNC, color=OI["sky"], alpha=0.18, lw=0, zorder=0)
    b.axhline(0, color="0.6", lw=0.5, zorder=0)
    for d in lv:
        b.plot([d["Mc"]], [d["diff"]], marker=ROLE_MK[d["role"]], ms=ROLE_MS[d["role"]], mew=0, color=J_COL[d["J"]],
               ls="none", zorder=3)
        if d["onesided"]:
            b.annotate("", (d["Mc"], d["diff"] - 0.008), (d["Mc"], d["diff"] - 0.0015),
                       arrowprops=dict(arrowstyle="-|>", lw=0.5, color=J_COL[d["J"]], mutation_scale=4))
    b.text(2, MAP_UNC + 0.002, r"frozen band $\pm$%.4f (max. calibration residual)" % MAP_UNC, fontsize=5.0,
           va="bottom", color="0.35")
    b.set_ylim(-0.045, 0.045)
    b.set_xlim(0, 103)
    b.set_xlabel(r"exact $M_c$")
    b.set_ylabel(r"$M_c^{\rm exact}-M_c^{\rm pred}$")
    hr = [Line2D([], [], color="k", marker=ROLE_MK["calib"], ms=ROLE_MS["calib"], mew=0, ls="none",
                 label="calibration (5)"),
          Line2D([], [], color="k", marker=ROLE_MK["blind1"], ms=ROLE_MS["blind1"], mew=0, ls="none",
                 label="pre-specified, tests 1-2 (4)"),
          Line2D([], [], color="k", marker=ROLE_MK["blind5"], ms=ROLE_MS["blind5"], mew=0, ls="none",
                 label="pre-specified, map (14)")]
    b.legend(handles=hr, loc="lower right", ncol=3, fontsize=4.8, columnspacing=0.6, handletextpad=0.1,
             borderaxespad=0.2)
    # (2026-11-19: no L = 6 arrows any more; both levels now have extrapolated values)
    # (c) unequal masses: |exact - model| on a log scale
    xs = np.arange(len(RUN6))
    c.axhspan(1e-3, 0.03, color=OI["sky"], alpha=0.18, lw=0, zorder=0)
    c.axhline(0.03, color=OI["sky"], lw=0.6)
    c.axhline(UNEQ_LENIENT, color="0.45", lw=0.6, ls=BO_LS)
    c.text(len(RUN6) - 0.45, 0.027, "frozen band 0.03", fontsize=5.0, ha="right", va="top", color="0.35")
    c.text(2.5, 0.108, "lenient 0.10", fontsize=5.0, ha="center", va="bottom", color="0.35")
    nullst = {"NS": ("s", "null NS (symmetric offset)"), "NT": ("^", "null NT (offset from exact thr.)"),
              "NP": ("v", "null NP (persistence)")}
    for i, d in enumerate(RUN6):
        c.plot([i], [abs(d["diff"])], marker="o", ms=4.0, mew=0, color=OI["black"], ls="none", zorder=4)
        for k, (mk, _) in nullst.items():
            if d[k] is not None:
                c.plot([i + 0.18], [abs(d["Mc"] - d[k])], marker=mk, ms=3.4, mfc="white", mec=OI["verm"], mew=0.7,
                       ls="none", zorder=3)
        c.text(i - 0.13, abs(d["diff"]), f"{d['null_ratio']:.0f}$\\times$", fontsize=4.6, ha="right", va="center",
               color="0.3")
    c.axvline(len(RUN6) - 1.5, color="0.7", lw=0.5)
    c.set_xticks(xs)
    c.set_xticklabels([f"({d['J']},{d['v']})" + ("\n$q$=1.5" if d["r"] != 2.0 else "") for d in RUN6], fontsize=5.5)
    c.text(2.5, 1.4e-3, r"$q=m_B/m_A=2$", fontsize=5.5, ha="center", va="bottom")
    c.set_yscale("log"); c.set_ylim(1e-3, 120); c.set_xlim(-0.6, len(RUN6) - 0.4)
    c.set_ylabel(r"$|M_c^{\rm exact}-M_c^{\rm model}|$ (light-nucleus mass)")
    c.set_xlabel(r"level $(J,v)$, unequal masses")
    c.tick_params(axis="x", which="both", top=False)
    hc = [Line2D([], [], marker="o", ms=4, mew=0, color="k", ls="none", label="theory (2ch+m), frozen")] + \
         [Line2D([], [], marker=mk, ms=3.4, mfc="white", mec=OI["verm"], mew=0.7, ls="none", label=lab)
          for mk, lab in nullst.values()]
    c.legend(handles=hc, loc="upper left", fontsize=4.8, handletextpad=0.2, borderaxespad=0.3)
    c.text(0.5, 0.2, "numbers: smallest null miss / theory miss\ntheory residuals all < 0 (overestimates $M_c$)",
           transform=c.transAxes, fontsize=4.6, ha="center", va="bottom", color="0.3")
    panel(a, "a", x=-0.13); panel(b, "b", x=-0.13); panel(c, "c", x=-0.2)
    save(fig, "F3_offsets")
    return {"n_levels": len(lv), "max|diff| prespec": max(abs(d["diff"]) for d in lv if d["role"] != "calib"),
            "null_ratio range": (min(d["null_ratio"] for d in RUN6), max(d["null_ratio"] for d in RUN6))}


# ------------------------------------------------------------------ F8: non-adiabatic mass functions (was F3a)
def fig8():
    mf = pd.read_csv(p1("bo_offset", "massfunc", "results_fine.csv"))
    fig, a = plt.subplots(figsize=(W1, 6.0 * CM))
    a.plot(mf.R, mf.m_v, "-", color=OI["blue"], label=r"$m_v(R)$ (this work)")
    a.plot(mf.R, mf.m_r, "-", color=OI["verm"], label=r"$m_r(R)$ (this work)")
    a.axhline(0.5, color="0.5", lw=0.5, ls="-")
    Re = 1.99719332
    a.axvline(Re, color="0.7", lw=0.5)
    a.text(Re + 0.2, 0.68, r"$R_e$", fontsize=6)
    a.text(16.5, 0.51, "1/2", fontsize=6, va="bottom", color="0.4")
    first = {"simple LCAO": True, "Finkelstein-Horowitz": True}
    for R, mv, mr, model, approx in KUTZELNIGG:
        mk = "s" if model == "simple LCAO" else "D"
        fc = "white"
        for val, col in ((mv, OI["blue"]), (mr, OI["verm"])):
            if val is None:
                continue
            lab = None
            if first[model] and col == OI["blue"]:
                lab = f"Kutzelnigg 2007, {model}" + (" (faded: approx.)" if model == "simple LCAO" else "")
                first[model] = False
            a.plot([R], [val], marker=mk, ms=3.4 if not approx else 3.0, mfc=fc if not approx else "none",
                   mec=col, mew=0.8, ls="none", label=lab, alpha=1.0 if not approx else 0.6)
    a.set_xlim(0, 20); a.set_ylim(0, 0.72)
    a.set_xlabel(r"$R$ (bohr)"); a.set_ylabel(r"non-adiabatic mass ($m_e$)")
    h, l = a.get_legend_handles_labels()
    a.legend(h, l, loc="lower right", fontsize=5.5)
    save(fig, "F8_massfunc")


# ------------------------------------------------------------------ F4
def fig4():
    Mm, f, lev, M, L = crossover_f()
    fig, (a, s) = plt.subplots(2, 1, figsize=(W2, 11.0 * CM), sharex=True,
                               gridspec_kw=dict(height_ratios=[1.0, 1.25]))
    gs = [m for m in MARKERS if m[4] != "thr"]
    lo_all = min(m[0] for m in gs); hi_all = max((m[1] or m[0]) for m in gs)
    for ax in (a, s):
        ax.axvspan(lo_all, hi_all, color=OI["sky"], alpha=0.12, lw=0, zorder=0)
    for m in MARKERS:
        if m[4] == "thr":
            a.axvline(m[0], color=OI["green"], lw=0.5, alpha=0.5, zorder=0)
    a.plot(Mm, f, "-", color=OI["black"], lw=0.7)
    a.plot(Mm, f, color=OI["black"], **EXACT)
    for y, Mx in lev.items():
        a.axhline(y, color="0.6", lw=0.5)
        a.plot([Mx, Mx], [-0.05, y], color="0.6", lw=0.5)
        a.text(Mx * 1.12, y - 0.025, f"{int(round(100*y))}%: $M$ = {Mx:.3g}", fontsize=5.5, va="top")
    a.set_ylabel(r"$f=-4\,d\ln\mathcal{L}/d\ln M$")
    a.set_ylim(-0.05, 1.08)
    a.text(0.012, 0.97, r"$\mathcal{L}=\Delta R/\langle R\rangle$ (exact ground state); $f\to1$: BO limit $\mathcal{L}\propto M^{-1/4}$",
           transform=a.transAxes, fontsize=5.5, va="top",
           bbox=dict(fc="white", ec="none", pad=0.5))
    # marker strip: one row per ground-state marker, thresholds on the last row
    kinds = {"lit": (OI["verm"], "s"), "ours": (OI["blue"], "o"), "thr": (OI["green"], "D")}
    short = {0.2: "Lindemann-ratio maximum", 0.357: "entanglement extrema (band)",
             0.4: "X density about CM, max to hole (Matyus 2011)", 0.491: "minimum of g(M), fitted (King 2013)",
             0.4924222: "minimum of g(M), fit-free", 0.5: r"maximum of $E/\mu_{Xe}$",
             0.543: "stability-boundary minimum (Kais-Shi 2000, converted)",
             0.635: "X density about CM, max to hole", 0.98: "marginal-distribution criterion (Salas 2022, interpolated)",
             1.11: "maximum of X kinetic energy", 1.3856: r"steepest rise of $\chi_{\rm inv}$ (fidelity)"}
    short = {k: v.replace("Matyus", r"M$\acute{\rm a}$tyus") for k, v in short.items()}
    gs_sorted = sorted(gs, key=lambda m: m[0])
    n = len(gs_sorted)
    for j, (lo, hi, desc, src, kind, status) in enumerate(gs_sorted):
        y = n - j
        col, mk = kinds[kind]
        if hi is not None:
            s.plot([lo, hi], [y, y], "-", color=col, lw=2.4, solid_capstyle="butt")
        else:
            s.plot([lo], [y], marker=mk, color=col, ms=3.6, mew=0, ls="none")
        s.text((hi or lo) * 1.12, y, short[lo], fontsize=5.3, va="center", ha="left")
    for lo, hi, desc, src, kind, status in MARKERS:
        if kind == "thr":
            s.plot([lo], [-0.6], marker="D", color=OI["green"], ms=3.4, mew=0, ls="none")
    tl = {4.2082827: "$L$=1", 9.525634: "$v$=1", 11.4894217: "$L$=2", 29.63939: "$v$=2",
          59.31866: "$v$=3", 98.26550: "$v$=4", 17.094282: "(1,1)", 21.376199: "$L$=3", 33.914349: "$L$=4"}
    for Mx, t in tl.items():
        ha = "right" if Mx in (9.525634, 29.63939) else ("left" if Mx in (11.4894217, 33.914349) else "center")
        low = Mx in (17.094282, 33.914349, 29.63939)          # second label row (crowded 11-34)
        s.text(Mx * ((0.97 if low else 0.93) if ha == "right" else 1.07 if ha == "left" else 1.0), -0.85 if low else -0.35, t,
               fontsize=5.3, ha=ha, va="top" if low else "bottom")
        s.plot([Mx, Mx], [-1.75, -0.6], color=OI["green"], lw=0.5, alpha=0.5)
    s.set_yticks([]); s.set_ylim(-1.75, n + 0.8)
    s.tick_params(axis="y", which="both", right=False, left=False)
    s.set_xscale("log"); s.set_xlim(0.04, 1.3e4)
    s.set_xlabel(r"mass ratio $M$")
    s.legend(handles=[Line2D([], [], marker="s", color=OI["verm"], ls="none", ms=3.6, label="proposed in the literature"),
                      Line2D([], [], marker="o", color=OI["blue"], ls="none", ms=3.6, label="ground-state marker, this work"),
                      Line2D([], [], marker="D", color=OI["green"], ls="none", ms=3.4, label="bound-state threshold, this work")],
             loc="upper right", fontsize=5.5)
    panel(a, "a", x=-0.07); panel(s, "b", x=-0.07)
    fig.subplots_adjust(hspace=0.05)
    save(fig, "F4_crossover")
    return lev


# ------------------------------------------------------------------ F5
def fig5():
    e = load_ground()
    s = pd.read_csv(p1("mass_scans", "spacings.csv")).sort_values("M")
    fig, (a, b) = plt.subplots(1, 2, figsize=(W2, 6.2 * CM))
    sel = e.M >= 1
    a.plot(e.kappa[sel], e.exact_minus_BO[sel], marker="s", color=OI["verm"], ls="none", ms=3.0, mew=0,
           label=r"$E-E_{\rm BO}$ ($>0$)")
    a.plot(e.kappa, -e.exact_minus_BOad, color=OI["blue"], label=r"$E_{\rm BO+ad}-E$ ($>0$)", **EXACT)
    k = np.logspace(-1.05, np.log10(2.2), 100)
    a.plot(k, C4 * k ** 4, color=OI["verm"], ls=BO_LS, lw=0.8, label=r"$0.1750\,\kappa^4$")
    a.plot(k, -C6 * k ** 6, color=OI["blue"], ls=BOAD_LS, lw=1.0, label=r"$0.03931\,\kappa^6$")
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel(r"$\kappa=M^{-1/4}$"); a.set_ylabel(r"energy difference ($E_h$)")
    a.set_xlim(0.09, 2.3); a.set_ylim(1e-8, 0.3)
    a.legend(loc="upper left", fontsize=5.5)
    ax2 = a.secondary_xaxis("top", functions=(lambda kk: kk ** -4.0, lambda MM: MM ** -0.25))
    ax2.set_xticks([1e4, 1e2, 1, 0.1]); ax2.set_xticklabels(["$10^4$", "$10^2$", "1", "0.1"])
    ax2.set_xlabel(r"$M$", labelpad=1)
    ax2.tick_params(direction="in", labelsize=6)
    a.tick_params(axis="x", which="both", top=False)
    # (b) d2/d1
    s = s.dropna(subset=["d2_over_d1"])
    b.plot(s.M, s.BO_d2_over_d1, color=OI["verm"], ls=BO_LS, lw=0.9, label="BO")
    b.plot(s.M, s.BOad_d2_over_d1, color=OI["blue"], ls=BOAD_LS, lw=1.1, label="BO+ad")
    b.plot(s.M, s.d2_over_d1, color=OI["black"], label="exact", **EXACT)
    b.axhline(3.0, color="0.5", lw=0.5)
    b.text(1.5e3, 3.005, "rigid rotor 3", fontsize=5.5, ha="right", va="bottom", color="0.4")
    b.set_xscale("log"); b.set_ylim(2.15, 3.08)
    b.set_xlabel(r"$M$"); b.set_ylabel(r"$d_2/d_1$, $d_L=E(L)-E(L{=}0)$")
    b.legend(loc="lower right", fontsize=5.5)
    ins = b.inset_axes([0.42, 0.33, 0.52, 0.36])
    ins.plot(s.M, np.abs(s.d2_over_d1 - s.BO_d2_over_d1), color=OI["verm"], ls=BO_LS, lw=0.8)
    ins.plot(s.M, np.abs(s.d2_over_d1 - s.BOad_d2_over_d1), color=OI["blue"], ls=BOAD_LS, lw=1.0)
    ins.set_xscale("log"); ins.set_yscale("log")
    ins.tick_params(labelsize=5); ins.set_ylabel(r"|exact$-$model|", fontsize=5.5, labelpad=1)
    sig = s.d2_over_d1 * np.sqrt((s.d2_err / s.d2) ** 2 + (s.d1_err / s.d1) ** 2)
    ins.plot(s.M, sig, color="0.55", lw=0.7)
    ins.text(9e3, 3e-3, "uncertainty\nof exact", fontsize=4.8, color="0.45", ha="right")
    ins.set_ylim(1e-7, 0.1)
    panel(a, "a", x=-0.17); panel(b, "b", x=-0.17)
    fig.subplots_adjust(wspace=0.32)
    save(fig, "F5_BO_approach")


# ------------------------------------------------------------------ F6
def fig6():
    t = pd.read_csv(p1("entanglement", "table.csv")).sort_values("M").drop_duplicates("M")
    fig, ax = plt.subplots(figsize=(W1, 6.0 * CM))
    for col, c, lab in (("S", OI["black"], r"$S$"), ("S_orient", OI["blue"], r"$S_{\rm orient}$"),
                        ("S_rad", OI["verm"], r"$S_{\rm rad}$")):
        ax.plot(t.M, t[col], "-", color=c, lw=0.6)
        ax.plot(t.M, t[col], color=c, label=lab + " (exact)", **EXACT)
        bo = t.dropna(subset=[col + "_BO"])
        if col == "S" and len(bo):
            ax.plot(bo.M, bo[col + "_BO"], color=c, ls=BO_LS, lw=0.8, label=r"$S$ (BO product state)")
    ax.axhline(S_INF, color=OI["green"], lw=0.7)
    ax.text(0.06, S_INF - 0.02, r"$S(\infty)=0.249338$, orientational", fontsize=5.5,
            color=OI["green"], va="top")
    ax.set_xscale("log"); ax.set_xlim(0.04, 1.3e4); ax.set_ylim(0, 1.75)
    ax.set_xlabel(r"$M$"); ax.set_ylabel(r"entropy (nats)")
    ax.legend(loc="upper right", fontsize=5.5)
    save(fig, "F6_entanglement")


# ------------------------------------------------------------------ F7
def fig7():
    levels = MUONIC_LEVELS
    fig, ax = plt.subplots(figsize=(W1, 5.4 * CM))
    nmatch = 0; nmis = 0
    for j, m in enumerate(MUONIC):
        y = len(MUONIC) - 1 - j
        for i, lv in enumerate(levels):
            thr = LEVEL_THRESHOLDS[lv]
            lit = lv in m["lit"]
            pred = None if thr is None else (m["M"] > thr)
            if pred is None:
                ax.add_patch(Rectangle((i - 0.45, y - 0.4), 0.9, 0.8, fc="0.93", ec="none"))
            else:
                match = pred == lit
                if lv != (0, 0):
                    nmatch += match; nmis += (not match)
                col = OI["blue"] if match else OI["verm"]
                ax.plot(i - 0.16, y, marker="o", ms=5.0, mfc=col if pred else "white", mec=col, mew=0.9, ls="none")
            ax.plot(i + 0.18, y, marker="s", ms=4.0, mfc="0.45" if lit else "white", mec="0.45", mew=0.8, ls="none")
    ax.set_xticks(range(len(levels)))
    tick = []
    for lv in levels:
        t = LEVEL_THRESHOLDS[lv]
        tick.append(f"({lv[0]},{lv[1]})\n" + ("always" if t == 0.0 else f"{t:.2f}"))
    ax.set_xticklabels(tick, fontsize=5.3)
    ax.set_yticks(range(len(MUONIC)))
    ax.set_yticklabels([f"{m['name']}\n$M$={m['M']:.2f}" for m in MUONIC[::-1]], fontsize=6)
    ax.set_xlim(-0.6, len(levels) - 0.4); ax.set_ylim(-0.6, len(MUONIC) - 0.4)
    ax.set_xlabel(r"level $(J,v)$ and our threshold mass $M_c$", labelpad=2)
    ax.tick_params(which="both", top=False, right=False, length=0)
    h = [Line2D([], [], marker="o", ms=5, mfc=OI["blue"], mec=OI["blue"], ls="none", label="ours: bound ($M>M_c$)"),
         Line2D([], [], marker="o", ms=5, mfc="white", mec=OI["blue"], ls="none", label="ours: not bound"),
         Line2D([], [], marker="s", ms=4, mfc="0.45", mec="0.45", ls="none", label="Frolov-Wardlaw: bound"),
         Line2D([], [], marker="s", ms=4, mfc="white", mec="0.45", ls="none", label="Frolov-Wardlaw: absent")]
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.27), ncol=2, fontsize=5.5,
              columnspacing=1.0, handletextpad=0.3)
    print(f"F7: {nmatch} agreements, {nmis} mismatches (excluding the always-bound (0,0))")
    save(fig, "F7_muonic")
    return nmatch, nmis


FIGS = dict(F1=fig1, F2=fig2, F2h=fig2h, F3=fig3, F4=fig4, F5=fig5, F6=fig6, F7=fig7, F8=fig8)

if __name__ == "__main__":
    which = sys.argv[1:] or list(FIGS)
    for w in which:
        r = FIGS[w]()
        if r is not None:
            print(w, r)
