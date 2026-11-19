"""Numbers quoted from paper1/RESULTS_SUMMARY.md (single source of truth), shared by
make_figures.py and make_tables.py.  Each entry carries the summary section it comes from and a
status label: 'rigorous' (variational bound or exact statement), 'extrapolated' (threshold-law
fits, basis extrapolation), 'fitted' (fit to data), 'derived' (parameter-free prediction from our
own theory and computed functions), 'literature' (value taken from a paper), 'converted'
(literature parameter converted by us), 'located' (extremum/zero located on our scan).

check_against_data() re-reads the raw data files and asserts agreement with these numbers to the
quoted precision; any mismatch is printed and recorded in README.md (Discrepancies).
"""
import os, re, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
P1 = os.path.normpath(os.path.join(HERE, ".."))          # paper1/
ROOT = os.path.normpath(os.path.join(P1, ".."))          # project root

def p1(*a):
    return os.path.join(P1, *a)

# ---------------------------------------------------------------- critical masses (summary 1, 6, 15)
# state, exact M_c, +-, proven upper bound, BO, BO+ad, exact non-adiabatic offset (BO+ad - exact),
# parameter-free predicted offset (non-perturbative fold of m_v(R), m_r(R)), pre-specified flag
CRIT = [
    dict(key="v1", label="$v=1$, $L=0$", Mc=9.525634, err=0.000004, ub="9.52586", BO=11.0814566,
         BOad=10.0631790, off=0.5375, pred=0.5281, prespec=False),
    dict(key="v2", label="$v=2$, $L=0$", Mc=29.63939, err=0.00001, ub="29.64016", BO=31.1677743,
         BOad=30.1455552, off=0.5062, pred=0.5124, prespec=False),
    dict(key="v3", label="$v=3$, $L=0$", Mc=59.31866, err=0.00002, ub="59.32056", BO=60.8466744,
         BOad=59.8156817, off=0.4970, pred=0.5017, prespec=False),
    dict(key="v4", label="$v=4$, $L=0$", Mc=98.26550, err=0.00011, ub="98.26962", BO=99.7973934,
         BOad=98.7585118, off=0.4930, pred=0.4958, prespec=True),
    dict(key="L1", label="$L=1$ (antisym.)", Mc=4.2082827, err=0.0000005, ub="4.2082830", BO=5.5374983,
         BOad=4.5705653, off=0.3623, pred=0.3328, prespec=False),
    dict(key="L2", label="$L=2$ (sym.)", Mc=11.4894217, err=0.0000003, ub="11.4894218", BO=12.7397228,
         BOad=11.7786164, off=0.2892, pred=0.2778, prespec=False),
    # thresholds_Lv (summary of paper1/thresholds_Lv/RESULT.md, 2026-11-15 22:45). ub = rounding-safe
    # bound (margin 1e-11 Eh, tol 1e-12); BO, BO+ad from thresholds_Lv/bo_brackets.txt (coupling.py).
    dict(key="J1v1", label="$L=1$, $v=1$ (antisym.)", Mc=17.094282, err=0.000001, ub="17.094281939",
         BO=18.4944752, BOad=17.5236378, off=0.4294, pred=0.4165, prespec=True,
         src="thresholds_Lv/RESULT.md, J1v1_crossings.csv"),
    dict(key="L3", label="$L=3$ (antisym.)", Mc=21.376199, err=0.000001, ub="21.376203385",
         BO=22.5995902, BOad=21.6380558, off=0.2619, pred=0.2562, prespec="secondary",
         src="thresholds_Lv/RESULT.md, L3A_crossings.csv"),
    dict(key="L4", label="$L=4$ (sym.)", Mc=33.914349, err=0.000002, ub="33.914359616",
         BO=35.1246258, BOad=34.1615637, off=0.2472, pred=0.2441, prespec="secondary",
         src="thresholds_Lv/RESULT.md, L4S_crossings.csv"),
]
V4_PRED_MC = (98.2627, 0.0095)      # frozen 2026-11-14 17:22, PREDICTION.md (summary 1)
# frozen 2026-11-15 20:17, thresholds_Lv/PREDICTION_J1v1.md (sha256 verified below); J1v1 primary,
# L3/L4 secondary. Caveat (thresholds_Lv/PROGRESS.md): the J1v1 scale optimisation (stage 1, N=12)
# ran 20:11-20:13, before the freeze; it computed no crossing.
LV_PRED_MC = {"J1v1": (17.1072, 0.0295), "L3": (21.3818, 0.0295), "L4": (33.9175, 0.0295)}
LV_FROZEN = "2026-11-15 20:17"
# exploratory: 2p sigma_u-like J=1 level (L=1 'S', root 0); not converged to the level of the others
L1S = dict(key="L1S", label=r"$L=1$ (sym., 2p$\sigma_u$-like)", Mc=813.54, err=0.07, ub="813.6485",
           ub_nomargin="813.6015", BO=815.558715, src="thresholds_Lv/RESULT.md, L1S_crossings.csv")
# Frolov-Wardlaw (1,1) estimate (basis-limited upper estimate; RESULT.md, arXiv:1008.3010)
FW_11 = 17.09444

# ---------------------------------------------------------------- muonic molecules (summary 1)
# M = m_nucleus / m_mu (Frolov-Wardlaw masses), literature bound (J, v) levels (Frolov & Wardlaw,
# arXiv:1008.3010 = EPJD 63, 339 (2011), via docs/literature/referee_check_2026-11-15.md Part B; same level
# lists in Frolov arXiv:2503.17389 (2025), docs/literature/frolov_2025_2026_notes.md). Existence only.
MUONIC = [
    dict(name=r"pp$\mu$", M=8.8802, lit={(0, 0), (1, 0)}),
    dict(name=r"dd$\mu$", M=17.7517, lit={(0, 0), (0, 1), (1, 0), (1, 1), (2, 0)}),
    dict(name=r"tt$\mu$", M=26.5849, lit={(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (3, 0)}),
]
# (J, v) -> our threshold mass (0.0: always bound). Odd J: exchange-antisymmetric (L 'A'), even J: 'S'.
LEVEL_THRESHOLDS = {(0, 0): 0.0, (0, 1): 9.525634, (0, 2): 29.63939, (1, 0): 4.2082827,
                    (1, 1): 17.094282, (2, 0): 11.4894217, (3, 0): 21.376199, (4, 0): 33.914349,
                    (2, 1): 27.69724706}     # (2,1): run 5 (threshold_map), the only run-5 level near the muonic M
MUONIC_LEVELS = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (4, 0)]
LEVEL_KEY = {(0, 1): "v1", (0, 2): "v2", (1, 0): "L1", (1, 1): "J1v1", (2, 0): "L2", (3, 0): "L3",
             (4, 0): "L4"}

# ---------------------------------------------------------------- crossover markers (summary 2, 7-9, 12, 13)
# (M_lo, M_hi or None, descriptor, source, kind) ; kind: 'lit' literature-proposed, 'ours' found by us,
# 'thr' bound-state threshold (ours, sharp)
MARKERS = [
    (0.2, 0.3, "Lindemann-ratio maximum", "this work (mass scan)", "ours", "located"),
    (0.357, 0.475, "electron-nuclear entanglement extrema (band)", "this work", "ours", "located"),
    (0.4, 0.8, "X density about the CM: maximum → hole", "Mátyus et al. 2011 (interval)", "lit", "literature"),
    (0.491, None, "minimum of $g(M)=(E_0-E_{\\rm th})/E_{\\rm th}$ (fitted)", "King et al. 2013", "lit", "literature"),
    (0.4924222, None, "minimum of $g(M)$, fit-free", "this work (recomputed)", "ours", "located"),
    (0.5, None, "maximum of $E/\\mu_{Xe}$ (approx.)", "this work (mass scan)", "ours", "located"),
    (0.543, None, "minimum of charge-stability boundary", "Kais & Shi 2000 (converted by us)", "lit", "converted"),
    (0.635, None, "X density about the CM: maximum → hole", "this work (same criterion as Mátyus)", "ours", "located"),
    (0.98, None, "sign change of marginal-distribution criterion", "Salas et al. 2022 (interpolated by us)", "lit", "converted"),
    (1.11, None, "maximum of X kinetic energy (zero of $d^2E/d(\\ln M)^2$)", "this work (mass scan)", "ours", "located"),
    (1.3856, None, "steepest rise of $\\chi_{\\rm inv}$ (fidelity)", "this work", "ours", "located"),
    (4.2082827, None, "$L=1$ (antisym.) level binds", "this work", "thr", "extrapolated"),
    (9.525634, None, "$v=1$ level binds", "this work", "thr", "extrapolated"),
    (11.4894217, None, "$L=2$ (sym.) level binds", "this work", "thr", "extrapolated"),
    (29.63939, None, "$v=2$ level binds", "this work", "thr", "extrapolated"),
    (59.31866, None, "$v=3$ level binds", "this work", "thr", "extrapolated"),
    (98.26550, None, "$v=4$ level binds", "this work", "thr", "extrapolated"),
    (17.094282, None, "$L=1$, $v=1$ (antisym.) level binds", "this work", "thr", "extrapolated"),
    (21.376199, None, "$L=3$ (antisym.) level binds", "this work", "thr", "extrapolated"),
    (33.914349, None, "$L=4$ (sym.) level binds", "this work", "thr", "extrapolated"),
]
# crossover width (summary 7): local Lindemann exponent levels
F_LEVELS = {0.1: 0.455, 0.5: 1.52, 0.9: 95.0}

# ---------------------------------------------------------------- asymptotics (summary 7, 13)
C6 = -0.03931          # +- 0.00003, fixed-exponent fits (kappa6.md)
C4 = 0.1750            # a of a kappa^4 + b kappa^6 fit of E - E_BO (summary 7)
S_INF = 0.249338       # nats, M -> inf (rigorous limit, summary 13)
MV_RE, MR_RE = 0.34655, 0.07542   # summary 15

# Kutzelnigg 2007 LCAO estimates (paper1/bo_offset/massfunc/LITERATURE_COMPARISON.md, table [K]).
# (R, m_v or None, m_r or None, model, approximate?)
KUTZELNIGG = [
    (2.0, 0.20, 0.06, "simple LCAO", False),
    (2.5, 0.29, 0.11, "simple LCAO", False),
    (3.0, 0.50, None, "simple LCAO", True),
    (5.0, 0.60, None, "simple LCAO", True),
    (2.0, 0.31, 0.14, "Finkelstein-Horowitz", False),
]


# ---------------------------------------------------------------- run 5: threshold map (threshold_map/RESULT.md, 2026-11-17 01:59)
# 14 new (J, v) levels of the 1s sigma_g manifold with M_c <= ~100, blind-predicted (frozen 2026-11-16 09:35,
# PREDICTION_MAP.md, sha256 checked).  Per level: exact M_c = extrapolated value of RESULT.md ('extrapolated',
# geometric in compact size) with its error; best = best proven upper bound (up to rounding, margin 0);
# safe = rounding-safe bound (margin 1e-11 Eh) = the numerically certified upper bound quoted in the tables;
# BO, BO+ad from bo_map.json; pred = frozen parameter-free prediction (prediction_map.json, M_pred_np).
# L = 6 levels ('onesided'): the compact-size sequences are not monotone (ratio 0.85 / nan) and for (6,1) the
# 'extrapolated' value lies ABOVE the best proven bound, so the extrapolation is unreliable.  For them the
# value is the best proven upper bound (rigorous up to rounding) with a one-sided (downward) error
# sig- = max(|best - extrapolated| + extrapolation error, largest step of the quoted compact-size sequence).
# (7,0) ('hp', added 2026-11-17): computed in paper1/threshold_map_L7hp with double-double (about 32-digit) assembly
# of S and H, single-pass canonical orthogonalization at tol 1e-14 and spurious-root guards (validate_L 7 passes there).
# Mc = 'M_c(J7v0) = 87.521468 +/- 1e-6' of the J7v0 note in threshold_map/RESULT.md (best bound minus the geometric
# compact-size correction); best = 87.5214689 (N = 12 + 2 diffuse); safe = margin 1e-11 Eh at tol 1e-14 (not 1e-12),
# N = 10 + 2 diffuse; crossings in threshold_map_L7hp/data/J7v0_crossings.csv.
MAP_FROZEN = "2026-11-16 09:35"
MAP_UNC = 0.0295            # frozen band = max |exact - pred| of the 9 calibration levels (level (1,0)); rms 0.0123
RUN5 = [
    dict(J=2, v=1, Mc=27.69724706, err=1.0e-9, best="27.697247055", safe="27.697247072", BO=29.019661, BOad=28.066887, pred=27.706757),
    dict(J=1, v=2, Mc=39.30862402, err=1.0e-9, best="39.308624018", safe="39.308624131", BO=40.755428, BOad=39.766531, pred=39.309346),
    dict(J=3, v=1, Mc=41.07434726, err=1.0e-9, best="41.074347265", safe="41.074347285", BO=42.355695, BOad=41.408480, pred=41.081019),
    dict(J=5, v=0, Mc=49.11537545, err=4.9e-6, best="49.115375144", safe="49.115380349", BO=50.318074, BOad=49.353402, pred=49.117289),
    dict(J=2, v=2, Mc=52.58916447, err=2.2e-9, best="52.589164468", safe="52.589164528", BO=53.968672, BOad=53.002721, pred=52.592031),
    dict(J=4, v=1, Mc=57.16494981, err=4.6e-7, best="57.164949070", safe="57.164949247", BO=58.421973, BOad=57.475700, pred=57.169673),
    # 2026-11-19: (6,0) and (6,1) recomputed with the double-double assembly (threshold_map_L7hp/data_J6, PROGRESS.md);
    # extrapolated from the two-diffuse sequences N = 8, 10, 12 (geometric ratios 0.14 and 0.04), user decision to
    # accept rows flagged only by the conditioning flag (lambda_min/lambda_max -1e-16 to -4e-16; slope guard passed).
    dict(J=6, v=0, Mc=66.983717, err=1.0e-6, dd6=True, best="66.98371834", safe="66.983798633", BO=68.181554, BOad=67.215424, pred=66.984942),
    dict(J=3, v=2, Mc=69.01388864, err=2.3e-8, best="69.013888660", safe="69.013888822", BO=70.347388, BOad=69.393851, pred=69.017053),
    dict(J=1, v=3, Mc=71.19123222, err=1.0e-9, best="71.191232217", safe="71.191232553", BO=72.664208, BOad=71.661196, pred=71.189617),
    dict(J=5, v=1, Mc=75.94865299, err=7.5e-7, best="75.948653736", safe="75.948654001", BO=77.189949, BOad=76.242851, pred=75.952068),
    dict(J=2, v=3, Mc=86.78733050, err=2.8e-8, best="86.787330525", safe="86.787330745", BO=88.206149, BOad=87.225441, pred=86.787548),
    dict(J=7, v=0, Mc=87.521468, err=1.0e-6, hp=True, best="87.521468935", safe="87.521471609", BO=88.715954,
         BOad=87.748556, pred=87.522283),
    dict(J=4, v=2, Mc=88.35419413, err=3.2e-7, best="88.354193808", safe="88.354193850", BO=89.656384, BOad=88.708986, pred=88.357074),
    dict(J=6, v=1, Mc=97.416732, err=2.0e-6, dd6=True, best="97.41673744", safe="97.416776689", BO=98.647176, BOad=97.698595, pred=97.419258),
]
for _c in RUN5:
    _c["key"] = f"J{_c['J']}v{_c['v']}"
    _c["off"] = _c["BOad"] - _c["Mc"]                     # exact non-adiabatic offset
    _c["pred_off"] = _c["BOad"] - _c["pred"]              # predicted offset
    _c["diff"] = _c["Mc"] - _c["pred"]
    _c["passed"] = abs(_c["diff"]) <= MAP_UNC + _c["err"]
    _c.setdefault("onesided", False)
    _c.setdefault("hp", False)
    _c.setdefault("dd6", False)
# (7,0): formerly 'not computed' (run 5's double-precision validate_L 7 failed); now the 'hp' entry of RUN5.
J7V0 = [c for c in RUN5 if c["key"] == "J7v0"][0]
# symmetric-family levels in the offset figure: (J, v) -> CRIT key and role
#   'calib'  = calibration only (exact known before any prediction was frozen)
#   'blind1' = pre-specified in the earlier blind tests (v=4; L=1 v=1, L=3, L=4), later part of the run-5 calibration
#   'blind5' = pre-specified in run 5 (RUN5)
CRIT_JV = {"v1": (0, 1), "v2": (0, 2), "v3": (0, 3), "v4": (0, 4), "L1": (1, 0), "L2": (2, 0),
           "J1v1": (1, 1), "L3": (3, 0), "L4": (4, 0)}

# ---------------------------------------------------------------- run 6: unequal masses {A+, B+, e-}, m_B = r m_A
# (unequal_masses/RESULT.md 2026-11-16 17:28; frozen predictions prediction_unequal.json 2026-11-16 09:48, sha256
# checked).  M_c = mass of the light nucleus A.  theory = '2ch+m' (two-channel + mass functions), frozen band unc
# (0.03; 0.0368 for J1v0) and lenient band 0.10.  Null models: NS = BO+ad minus the exact symmetric-family offset;
# NT = same from the exact threshold; NP = persistence (offset of the previous J = 0 level, J0v2-v4 only).
# best = best upper bound column, safe = certified bound (margin 1e-11) of RESULT.md.
UNEQ_FROZEN = "2026-11-16 09:48"
UNEQ_LENIENT = 0.10
RUN6 = [
    dict(t="J0v1", pk="r2_J0v1", r=2.0, J=0, v=1, Mc=11.573761, err=7.7e-6, best="11.57375731", safe="11.57393265",
         pred=11.713333, unc=0.03, NS=7.146439, NT=14.485127, NP=None, BO=8.352759, BOad=7.549564),
    dict(t="J0v2", pk="r2_J0v2", r=2.0, J=0, v=2, Mc=29.005219, err=2.9e-5, best="29.00521327", safe="29.00580969",
         pred=29.088634, unc=0.03, NS=22.230196, NT=33.814139, NP=26.6340, BO=23.417497, BOad=22.609846),
    dict(t="J0v3", pk="r2_J0v3", r=2.0, J=0, v=3, Mc=53.417413, err=9.8e-5, best="53.41737672", safe="53.41904003",
         pred=53.474144, unc=0.03, NS=44.48935, NT=60.837378, NP=51.2575, BO=45.676672, BOad=44.8621),
    dict(t="J0v4", pk="r2_J0v4", r=2.0, J=0, v=4, Mc=84.345408, err=4.4e-4, best="84.34540782", safe="84.34973515",
         pred=84.385819, unc=0.03, NS=73.699344, NT=95.51769, NP=82.6244, BO=74.889712, BOad=74.069094),
    dict(t="J1v0", pk="r2_J1v0", r=2.0, J=1, v=0, Mc=5.405970, err=2.1e-8, best="5.40597030", safe="5.40597031",
         pred=5.581935, unc=0.036796, NS=3.154713, NT=6.240383, NP=None, BO=4.19479, BOad=3.426438),
    dict(t="J2v0", pk="r2_J2v0", r=2.0, J=2, v=0, Mc=11.029758, err=9.4e-7, best="11.02975771", safe="11.02975806",
         pred=11.147797, unc=0.03, NS=8.616049, NT=11.53572, NP=None, BO=9.596459, BOad=8.832949),
    dict(t="S15J0v1", pk="r1.5_J0v1", r=1.5, J=0, v=1, Mc=11.675521, err=7.1e-6, best="11.67551689", safe="11.67568620",
         pred=11.811695, unc=0.03, NS=7.938893, NT=14.009849, NP=None, BO=9.251214, BOad=8.38681),
]
for _c in RUN6:
    _c["diff"] = _c["Mc"] - _c["pred"]
    _c["pass_unc"] = abs(_c["diff"]) <= _c["unc"] + _c["err"]
    _c["pass_len"] = abs(_c["diff"]) <= UNEQ_LENIENT + _c["err"]
    _nulls = {k: abs(_c["Mc"] - _c[k]) for k in ("NS", "NT", "NP") if _c[k] is not None}
    _c["null_best"] = min(_nulls.values())                # smallest null-model miss
    _c["null_ratio"] = _c["null_best"] / abs(_c["diff"])
# r = 1 validation of the unsymmetrized code (RESULT.md, diffuse-basis enlargement test): production basis 9.587305,
# mirror basis M 9.52565573, MB 9.52564891, vs symmetric exact 9.525634 (v = 1, L = 0).
R1_CHECK = dict(production=9.587305, mirror=9.52565573, mirror_big=9.52564891, exact_sym=9.525634)


def check_against_data(verbose=True):
    """Compare the hard-coded numbers with the raw data files. Returns list of discrepancy strings."""
    bad = []
    def chk(name, a, b, tol):
        ok = abs(a - b) <= tol
        if verbose:
            print(f"  {'ok ' if ok else 'BAD'} {name}: summary {a} vs data {b} (tol {tol:g})")
        if not ok:
            bad.append(f"{name}: summary {a} vs data {b}")
    # BO / BO+ad critical masses
    bo = pd.read_csv(p1("critical_masses", "bo_critical_masses.csv"), comment="#", header=None,
                     names=["model", "n", "Mc", "err", "a", "b"])
    for c, n in zip(CRIT[:3], (2, 3, 4)):
        chk(f"{c['key']} BO", c["BO"], float(bo[(bo.model == "BO") & (bo.n == n)].Mc.iloc[0]), 6e-8)
        chk(f"{c['key']} BO+ad", c["BOad"], float(bo[(bo.model == "BO+ad(own asymptote)") & (bo.n == n)].Mc.iloc[0]), 6e-8)
    pj = json.load(open(p1("critical_masses_v4", "prediction.json")))
    chk("v4 BO", CRIT[3]["BO"], pj["M_BO"], 6e-8)
    chk("v4 BO+ad", CRIT[3]["BOad"], pj["M_BOad"], 6e-8)
    ot = pd.read_csv(p1("bo_offset", "offset_table.csv"))
    for c, L in ((CRIT[4], 1), (CRIT[5], 2)):
        r = ot[ot.L == L].iloc[0]
        chk(f"{c['key']} BO", c["BO"], r.M_BO, 6e-8)
        chk(f"{c['key']} BO+ad", c["BOad"], r.M_BOad, 6e-8)
    # upper bounds (margin 1e-11, minimum over basis variants)
    ub = pd.read_csv(p1("critical_masses", "upper_bounds.csv"))
    ub4 = pd.read_csv(p1("critical_masses_v4", "upper_bounds.csv"))
    for c, s in zip(CRIT[:3], (1, 2, 3)):
        m = ub[(ub.state == s) & (ub.margin == 1e-11)].M_bound.min()
        chk(f"{c['key']} upper bound", float(c["ub"]), m, 6e-6)
    chk("v4 upper bound", float(CRIT[3]["ub"]), ub4[ub4.margin == 1e-11].M_bound.min(), 6e-6)
    for c, f in ((CRIT[4], "critical_L1A.json"), (CRIT[5], "critical_L2S.json")):
        d = json.load(open(p1("rotational", f)))
        chk(f"{c['key']} upper bound", float(c["ub"]), min(list(d["diffuse"].values()) + list(d.get("compact_check", {}).values())), 6e-8)
    # offsets and predicted offsets
    fr = pd.read_csv(p1("bo_offset", "massfunc", "fold_results.csv"))
    for c, st in zip([CRIT[0], CRIT[1], CRIT[2], CRIT[4], CRIT[5]],
                     ["L=0 v=1", "L=0 v=2", "L=0 v=3", "L=1", "L=2"]):
        r = fr[fr.state == st].iloc[0]
        chk(f"{c['key']} exact offset", c["off"], r.exact, 6e-5)
        chk(f"{c['key']} predicted offset", c["pred"], r.pred_nonpert, 6e-5)
    chk("v4 predicted offset", CRIT[3]["pred"], pj["fold"]["4"]["pred_nonpert"], 6e-5)
    chk("v4 predicted M_c", V4_PRED_MC[0], pj["predictions"]["pf_nonpert"], 6e-5)
    for c in CRIT:
        chk(f"{c['key']} offset = BO+ad - Mc", c["off"], c["BOad"] - c["Mc"], 6e-5)
    # m_v, m_r at R_e
    mf = pd.read_csv(p1("bo_offset", "massfunc", "results_fine.csv"))
    Re = 1.99719332
    chk("m_v(R_e)", MV_RE, float(np.interp(Re, mf.R, mf.m_v)), 2e-5)
    chk("m_r(R_e)", MR_RE, float(np.interp(Re, mf.R, mf.m_r)), 2e-5)
    # ---- thresholds_Lv: RESULT.md, *_crossings.csv, *_binding.csv, prediction_J1v1.json
    import hashlib
    tl = lambda *a: p1("thresholds_Lv", *a)
    res = open(tl("RESULT.md")).read()
    secs = {k: v for k, v in zip(["J1v1", "L3A", "L4S", "L1S"], re.split(r"\n## \(J=|\n## 2p", res)[1:5])}
    num = r"([-+0-9.eE]+)"
    def grab(sec, pat):
        return float(re.search(pat + r"\s*" + num, sec).group(1).rstrip("."))
    pj2 = json.load(open(tl("prediction_J1v1.json")))
    sha = open(tl("prediction_J1v1.sha256")).read().split()[0]
    sha_now = hashlib.sha256(open(tl("PREDICTION_J1v1.md"), "rb").read()).hexdigest()
    chk("PREDICTION_J1v1.md sha256 unchanged", 0, 0 if sha == sha_now else 1, 0)
    chk("prediction freeze date", 0, 0 if pj2["date"] == LV_FROZEN else 1, 0)
    by = {c["key"]: c for c in CRIT}
    for key, tgt, sk, pk in (("J1v1", "J1v1", "L1v1", "J1v1"), ("L3", "L3A", "L3v0", "J3v0"),
                             ("L4", "L4S", "L4v0", "J4v0")):
        c, sec = by[key], secs[tgt]
        chk(f"{key} M_c vs RESULT extrapolated", c["Mc"], grab(sec, r"extrapolated M_c ="), c["err"])
        chk(f"{key} upper bound vs RESULT", float(c["ub"]), grab(sec, r"rounding-safe bound \(margin 1e-11 Eh\):"), 1e-9)
        cr = pd.read_csv(tl(f"{tgt}_crossings.csv"))
        safe = cr[(cr.margin == 1e-11) & (cr.tol == 1e-12)].M_star.min()
        chk(f"{key} upper bound vs crossings.csv", float(c["ub"]), safe, 1e-9)
        best = cr[cr.margin == 0.0].M_star.min()
        chk(f"{key} best bound vs RESULT", best, grab(sec, r"best proven upper bound \(up to rounding\):"), 1e-9)
        chk(f"{key} M_c <= rounding-safe bound", 0, 0 if c["Mc"] - c["err"] <= safe else 1, 0)
        chk(f"{key} BO vs RESULT", c["BO"], grab(sec, r"\* BO"), 6e-7)
        chk(f"{key} BO+ad vs RESULT", c["BOad"], grab(sec, r"BO\+ad"), 6e-7)
        chk(f"{key} BO vs json", c["BO"], pj2["predictions"][pk]["M_BO"], 6e-8)
        chk(f"{key} BO+ad vs json", c["BOad"], pj2["predictions"][pk]["M_BOad"], 6e-8)
        chk(f"{key} exact offset vs RESULT", c["off"], grab(sec, r"M_c\(BO\+ad\) - M_c\(exact\) ="), 6e-5)
        chk(f"{key} predicted offset vs RESULT", c["pred"], grab(sec, r"\(predicted"), 6e-5)
        chk(f"{key} predicted offset vs json", c["pred"], pj2["states"][sk]["pred_nonpert"], 6e-5)
        chk(f"{key} frozen M_c vs json", LV_PRED_MC[key][0], pj2["predictions"][pk]["pf_nonpert"], 6e-5)
        chk(f"{key} frozen unc vs json", LV_PRED_MC[key][1], pj2["predictions"][pk]["unc"], 6e-5)
        # binding.csv: M* implied by x, and E_b/(M-M*) vs the RESULT table
        bd = pd.read_csv(tl(f"{tgt}_binding.csv"))
        Mstar = (bd.M / (1 + bd.x)).values
        chk(f"{key} binding.csv M* = best bound (max dev.)", 0.0, float(np.max(np.abs(Mstar - best))), 1e-7)
        tab = re.findall(r"\| " + num + r" \| " + num + r" \| " + num + r" \| " + num + r" \|", sec)
        rat = np.array([float(t[3]) for t in tab])
        chk(f"{key} binding.csv E_b/(M-M*) vs RESULT table (max rel. dev.)", 0.0,
            float(np.max(np.abs(bd.binding / (bd.M - Mstar) / rat - 1))), 1e-4)
    sec = secs["L1S"]
    chk("L1S M_c vs RESULT extrapolated", L1S["Mc"], grab(sec, r"extrapolated M_c ="), L1S["err"])
    chk("L1S BO vs RESULT", L1S["BO"], grab(sec, r"\* BO"), 6e-7)
    cr = pd.read_csv(tl("L1S_crossings.csv"))
    chk("L1S rounding-safe bound", float(L1S["ub"]), cr[(cr.margin == 1e-11) & (cr.tol == 1e-12)].M_star.min(), 1e-4)
    chk("L1S best bound", float(L1S["ub_nomargin"]), cr[cr.margin == 0.0].M_star.min(), 1e-4)
    chk("Frolov-Wardlaw (1,1) estimate above our bound", 0, 0 if FW_11 > float(by["J1v1"]["ub"]) else 1, 0)
    for c in CRIT[6:]:
        chk(f"{c['key']} offset = BO+ad - Mc", c["off"], c["BOad"] - c["Mc"], 6e-5)
    for lv, k in LEVEL_KEY.items():
        chk(f"level {lv} threshold = T1 {k}", LEVEL_THRESHOLDS[lv], by[k]["Mc"], 0)
    chk("level (2,1) threshold = run-5 value", LEVEL_THRESHOLDS[(2, 1)], [c for c in RUN5 if c["key"] == "J2v1"][0]["Mc"], 0)
    # ---- run 5: threshold_map (RESULT.md sections and summary table, data/*_crossings.csv, prediction_map.json,
    #      bo_map.json, sha256 of PREDICTION_MAP.md)
    tm = lambda *a: p1("threshold_map", *a)
    res5 = open(tm("RESULT.md")).read()
    sha = open(tm("prediction_map.sha256")).read().split()[0]
    chk("PREDICTION_MAP.md sha256 unchanged", 0, 0 if sha == hashlib.sha256(open(tm("PREDICTION_MAP.md"), "rb").read()).hexdigest() else 1, 0)
    pm = json.load(open(tm("prediction_map.json")))
    bm = json.load(open(tm("bo_map.json")))
    chk("map freeze date", 0, 0 if pm["date"] == MAP_FROZEN else 1, 0)
    chk("map band 0.0295 vs json unc", MAP_UNC, pm["unc"], 6e-5)
    chk("map band = max |resid| of calibration", pm["unc"], max(abs(c["resid"]) for c in pm["calibration"].values()), 1e-12)
    sec5 = {}
    for part in re.split(r"\n## \(", res5)[1:]:
        m = re.match(r"(\d),(\d)\)", part)
        if m:
            sec5[(int(m.group(1)), int(m.group(2)))] = part
    summ = {}
    for row in re.findall(r"^\| \((\d),(\d)\) \|([^\n]*)$", res5.split("## Summary: new thresholds")[1].split("## Complete map")[0], re.M):
        summ[(int(row[0]), int(row[1]))] = [x.strip() for x in row[2].split("|")]
    n_pass = 0
    for c in RUN5:
        k, sec, jv = c["key"], sec5[(c["J"], c["v"])], (c["J"], c["v"])
        if not c["dd6"]:
            chk(f"{k} best bound vs RESULT", float(c["best"]), grab(sec, r"best proven upper bound \(up to rounding\):"), 1e-9)
        chk(f"{k} certified bound vs RESULT", float(c["safe"]), grab(sec, r"rounding-safe bound \(margin 1e-11 Eh\):"), 1e-9)
        if c["hp"]:   # (7,0): double-double assembly, threshold_map_L7hp, standard tol 1e-14
            cr = pd.read_csv(p1("threshold_map_L7hp", "data", f"{k}_crossings.csv"))
            ctol = 1e-14
        else:
            cr = pd.read_csv(tm("data", f"{k}_crossings.csv"))
            ctol = 1e-12
        if c["dd6"]:   # (6,0), (6,1): double-double rerun of 2026-11-19 (threshold_map_L7hp/data_J6, J6_note.md)
            cr6 = pd.read_csv(p1("threshold_map_L7hp", "data_J6", f"{k}_crossings.csv"))
            chk(f"{k} best bound vs data_J6 crossings.csv (guard-clean)", float(c["best"]),
                cr6[(cr6.margin == 0.0) & (cr6.tol == 1e-14)].M_star.min(), 1e-8)
            note = open(p1("threshold_map_L7hp", "J6_note.md")).read()
            mm = re.search(r"M_c\(%d,%d\) = ([0-9.]+)\((\d)\)" % (c["J"], c["v"]), note)
            chk(f"{k} M_c vs J6_note.md", c["Mc"], float(mm.group(1)), 1e-12)
            chk(f"{k} error vs J6_note.md", c["err"], int(mm.group(2)) * 1e-6, 1e-12)
        else:
            chk(f"{k} best bound vs crossings.csv", float(c["best"]), cr[cr.margin == 0.0].M_star.min(), 1e-9)
        chk(f"{k} certified bound vs crossings.csv", float(c["safe"]),
            cr[(cr.margin == 1e-11) & (cr.tol == ctol)].M_star.min(), 1e-9)
        if c["hp"]:
            mm = re.search(r"M_c\(J7v0\) = ([0-9.]+) \+/- ([0-9.eE+-]+)", sec)
            ex, exe = float(mm.group(1)), float(mm.group(2))
        else:
            ex = grab(sec, r"extrapolated M_c ="); exe = float(re.search(r"extrapolated M_c = [0-9.]+ \+/- ([0-9.eE+-]+)", sec).group(1))
        if c["onesided"]:
            chk(f"{k} (L=6) value = best proven bound", c["Mc"], float(c["best"]), 0)
            chk(f"{k} (L=6) extrapolation vs RESULT", c["extrap"][0], ex, 1e-8)
            chk(f"{k} (L=6) extrapolation error vs RESULT", c["extrap"][1], exe, 1e-12)
            seq = [float(x) for x in re.search(r"N = \[[0-9, ]+\]: ([0-9., ]+);", sec).group(1).split(",")]
            chk(f"{k} (L=6) compact sequence vs RESULT", 0.0, float(np.max(np.abs(np.array(seq) - c["compact"]))), 1e-9)
            sig = max(abs(float(c["best"]) - ex) + exe, float(np.max(np.abs(np.diff(seq)))))
            chk(f"{k} (L=6) one-sided error rule", c["err"], sig, 6e-6)
        elif not c["dd6"]:
            chk(f"{k} M_c vs RESULT extrapolated", c["Mc"], ex, 1e-12)
            chk(f"{k} error vs RESULT", c["err"], exe, 1e-12)
        chk(f"{k} BO vs RESULT", c["BO"], grab(sec, r"\nBO"), 1e-6)
        chk(f"{k} BO+ad vs RESULT", c["BOad"], grab(sec, r"BO\+ad"), 1e-6)
        L = c["J"]; rts = bm["roots"]
        chk(f"{k} BO vs bo_map.json", c["BO"], rts[f"L={L} BO"][c["v"] + (1 if L == 0 else 0)], 6e-7)
        chk(f"{k} BO+ad vs bo_map.json", c["BOad"], rts[f"L={L} BO+ad"][c["v"]], 6e-7)
        pr = pm["predictions"][k]
        chk(f"{k} prediction vs json", c["pred"], pr["M_pred_np"], 6e-7)
        chk(f"{k} BO+ad vs prediction json", c["BOad"], pr["M_BOad"], 6e-7)
        chk(f"{k} predicted offset vs json", c["pred_off"], pr["off_nonpert"], 2e-6)
        row = summ[jv]
        chk(f"{k} prediction vs RESULT summary", c["pred"], float(row[4].split("+/-")[0]), 5e-5)
        if c["dd6"]:   # run 5's diff used the old value; recompute from the rerun
            miss = [float(x.replace("−", "-")) for x in re.findall(r"\*\*(−[0-9.]+)\*\*", note)]
            chk(f"{k} miss vs J6_note.md", c["diff"], miss[0 if c["v"] == 0 else 1], 1e-4)
        elif not c["onesided"]:   # RESULT's diff uses the extrapolated value
            chk(f"{k} diff vs RESULT summary", c["diff"], float(row[5]), 5e-5)
        else:
            chk(f"{k} (L=6) diff vs RESULT summary (extrapolated; differs by best - extrap)", c["diff"] - (c["Mc"] - c["extrap"][0]), float(row[5]), 5e-5)
        chk(f"{k} PASS vs RESULT summary", 0, 0 if (row[6] == "PASS") == c["passed"] else 1, 0)
        chk(f"{k} complete-map offset BO+ad - exact", c["off"], float(re.search(r"\| \(%d,%d\) \| [0-9.]+ \| [^|]+\| [0-9.]+ \| [0-9.]+ \| [0-9.]+ \| ([0-9.]+) \|" % jv, res5).group(1)), 2e-4)
        n_pass += c["passed"]
    # known inconsistency (README, discrepancies): extrapolated M_c - err above the best proven bound by more than the
    # printing precision (5e-9).  Expected: only (4,1) (7.4e-7 above, quoted error 4.6e-7).
    above = [c["key"] for c in RUN5 if not c["onesided"] and c["Mc"] - c["err"] > float(c["best"]) + 5e-9]
    chk("run 5 levels with extrapolated M_c above the best bound = ['J4v1'] (documented)", 0, 0 if above == ["J4v1"] else 1, 0)
    mt = re.search(r"Blind test: (\d+) PASS out of (\d+)", res5)
    chk("run 5: 14 PASS of 14 (RESULT)", n_pass, int(mt.group(1)), 0)
    chk("run 5: 14 targets (RESULT) = len(RUN5)", len(RUN5), int(mt.group(2)), 0)
    chk("run 5: 14 levels in RUN5", len(RUN5), 14, 0)
    s7 = sec5[(7, 0)]
    chk("J7v0 computed with double-double assembly (RESULT)", 0, 0 if "double-double" in s7 and "Not computed" not in s7 else 1, 0)
    chk("J7v0 blind test PASS stated in RESULT", 0, 0 if re.search(r"exact - pred = -0\.0008: \*\*PASS", s7) else 1, 0)
    chk("J7v0 M_c - err <= best bound", 0, 0 if J7V0["Mc"] - J7V0["err"] <= float(J7V0["best"]) else 1, 0)
    chk("J7v0 best bound <= certified bound", 0, 0 if float(J7V0["best"]) <= float(J7V0["safe"]) else 1, 0)
    cr7 = pd.read_csv(p1("threshold_map_L7hp", "data", "J7v0_crossings.csv"))
    chk("J7v0 excluded row (N = 8 + 2 diffuse) has no crossing in csv", 0,
        0 if cr7[cr7.label == "compact N=8 + 2 diffuse"].M_star.isna().all() else 1, 0)
    comp = cr7[(cr7.margin == 0.0) & cr7.label.str.fullmatch(r"compact N=\d+")].sort_values("size").M_star.values
    chk("J7v0 compact sequence monotone decreasing", 0, 0 if np.all(np.diff(comp) < 0) else 1, 0)
    chk("J7v0 in complete map: order (2,3) < (7,0) < (4,2)", 0,
        0 if "(2,3) < (7,0) < (4,2)" in res5.split("Order of appearance")[1].split("\n")[2] else 1, 0)
    # calibration (earlier levels): their predictions in prediction_map.json equal the T1 predicted offsets
    by = {c["key"]: c for c in CRIT}
    for key, (J, v) in CRIT_JV.items():
        cal = pm["calibration"][f"J{J}v{v}"]
        chk(f"{key} calibration exact vs T1", cal["exact"], by[key]["Mc"], 1e-9)
        chk(f"{key} calibration predicted offset vs T1", cal["off_nonpert"], by[key]["pred"], 6e-5)
        chk(f"{key} calibration BO+ad vs T1", cal["M_BOad"], by[key]["BOad"], 6e-7)
    # ---- run 6: unequal masses (RESULT.md blind table, prediction_unequal.json, sha256 files, crossings CSVs)
    um = lambda *a: p1("unequal_masses", *a)
    res6 = open(um("RESULT.md")).read()
    for line in open(um("prediction_unequal.sha256")):
        h, f = line.split()
        chk(f"{f} sha256 unchanged", 0, 0 if h == hashlib.sha256(open(um(f), "rb").read()).hexdigest() else 1, 0)
    pu = json.load(open(um("prediction_unequal.json")))
    chk("unequal freeze date", 0, 0 if pu["date"] == UNEQ_FROZEN else 1, 0)
    blind = res6.split("## Blind comparison")[1]
    top = res6.split("## Excluded rows")[0]
    for c in RUN6:
        t = c["t"]
        row = [x.strip() for x in re.search(r"^\| %s \|([^\n]*)$" % t, blind, re.M).group(1).split("|")]
        chk(f"run6 {t} M_c vs RESULT", c["Mc"], float(row[0].split("±")[0]), 1e-6)
        chk(f"run6 {t} error vs RESULT", c["err"], float(row[0].split("±")[1]), 1e-12)
        chk(f"run6 {t} prediction vs RESULT", c["pred"], float(row[1]), 6e-5)
        chk(f"run6 {t} unc vs RESULT", c["unc"], float(row[2]), 6e-5)
        chk(f"run6 {t} diff vs RESULT", c["diff"], float(row[3]), 6e-5)
        pa, pl = row[4].split("/")
        chk(f"run6 {t} PASS(unc) vs RESULT", 0, 0 if (pa.strip() == "PASS") == c["pass_unc"] else 1, 0)
        chk(f"run6 {t} PASS(0.10) vs RESULT", 0, 0 if (pl.strip() == "PASS") == c["pass_len"] else 1, 0)
        chk(f"run6 {t} null NS vs RESULT", c["NS"], float(row[5].split("(")[0]), 6e-5)
        chk(f"run6 {t} null NT vs RESULT", c["NT"], float(row[6].split("(")[0]), 6e-5)
        if c["NP"] is not None:
            chk(f"run6 {t} null NP vs RESULT", c["NP"], float(row[8]), 6e-5)
        else:
            chk(f"run6 {t} no NP in RESULT", 0, 0 if row[8] == "" else 1, 0)
        chk(f"run6 {t} BO vs RESULT", c["BO"], float(row[10]), 6e-5)
        chk(f"run6 {t} BO+ad vs RESULT", c["BOad"], float(row[11]), 6e-5)
        tg = pu["targets"][c["pk"]]
        chk(f"run6 {t} prediction vs json", c["pred"], tg["M_pred"], 6e-7)
        chk(f"run6 {t} unc vs json", c["unc"], tg["unc"], 6e-7)
        chk(f"run6 {t} NS vs json", c["NS"], tg["null_NS"], 6e-7)
        chk(f"run6 {t} NT vs json", c["NT"], tg["null_NT"], 6e-7)
        chk(f"run6 {t} BO vs json", c["BO"], tg["M"]["BO"], 6e-7)
        chk(f"run6 {t} BO+ad vs json", c["BOad"], tg["M"]["BOad"], 6e-7)
        chk(f"run6 {t} r, J, v vs json", 0, 0 if (tg["r"], tg["L"], tg["v"]) == (c["r"], c["J"], c["v"]) else 1, 0)
        trow = [x.strip() for x in re.search(r"^\| %s \|([^\n]*)$" % t, top, re.M).group(1).split("|")]
        chk(f"run6 {t} best bound vs RESULT", float(c["best"]), float(trow[1]), 1e-8)
        chk(f"run6 {t} certified bound vs RESULT", float(c["safe"]), float(trow[2]), 1e-8)
        cr = pd.read_csv(um(f"{t}_crossings.csv"))
        chk(f"run6 {t} certified bound occurs in crossings.csv (margin 1e-11)", 0.0,
            float(np.min(np.abs(cr[cr.margin == 1e-11].M_star - float(c["safe"])))), 1e-8)
        chk(f"run6 {t} M_c <= certified bound", 0, 0 if c["Mc"] - c["err"] <= float(c["safe"]) else 1, 0)
    chk("run6: no level passes the frozen band", 0, sum(c["pass_unc"] for c in RUN6), 0)
    rr = re.search(r"Mirror basis M: ([0-9.]+) .*?\(MB, 14000 functions\): ([0-9.]+)", res6, re.S)
    chk("r=1 mirror basis vs RESULT", R1_CHECK["mirror"], float(rr.group(1).rstrip(".")), 1e-8)
    chk("r=1 MB basis vs RESULT", R1_CHECK["mirror_big"], float(rr.group(2).rstrip(".")), 1e-8)
    chk("r=1 production vs RESULT", R1_CHECK["production"],
        float(re.search(r"\| R1J0v1 \| \d+ \| [0-9.]+ \| [0-9.]+ \| [0-9.]+ \| [^|]+\| ([0-9.]+) \|", res6).group(1)), 1e-6)
    chk("r=1 symmetric exact = T1 v=1", R1_CHECK["exact_sym"], by["v1"]["Mc"], 0)
    # crossover levels (recomputed with width.py logic in make_figures.crossover_f)
    return bad


if __name__ == "__main__":
    b = check_against_data()
    print("discrepancies:", b if b else "none")
