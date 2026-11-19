"""Tables for Paper 1 (PCCP), written as Markdown (tables.md) and LaTeX/booktabs (tables.tex).
Run: python3 make_tables.py
Numbers come from summary_values.py (= paper_data/RESULTS_SUMMARY.md), cross-checked against the raw
data files by summary_values.check_against_data() before writing.
"""
import os, re
from summary_values import (MUONIC_LEVELS, RUN5, RUN6, MAP_UNC, MAP_FROZEN, UNEQ_FROZEN, UNEQ_LENIENT,
                            R1_CHECK, CRIT, MARKERS, MUONIC, LEVEL_THRESHOLDS, LEVEL_KEY, V4_PRED_MC, LV_PRED_MC,
                            LV_FROZEN, L1S, FW_11, check_against_data)

OUT = os.path.dirname(os.path.abspath(__file__))


def fmt_pm(v, e):
    """central value with uncertainty in parentheses on the last digits, e.g. 9.525634(4), 98.26550(11)."""
    es = f"{e:.10f}".rstrip("0")
    sig = es.split(".")[1].lstrip("0")          # significant digits of the quoted uncertainty
    nd = len(es.split(".")[1])
    return f"{v:.{nd}f}({int(sig)})"


def md_to_tex(cell):
    t = cell.replace("±", r"$\pm$").replace("%", r"\%").replace("&", r"\&")
    t = t.replace("μ", r"$\mu$").replace("σ", r"$\sigma$").replace("≈", r"$\approx$").replace("→", r"$\to$")
    t = t.replace("–", "--").replace("≤", r"$\le$").replace("−", "-").replace("§", r"\S").replace("á", r"\'a")
    t = t.replace("**", "")
    t = re.sub(r"\*(.+?)\*", r"\\textit{\1}", t)
    parts = t.split("$")
    t = "$".join(p.replace("_", r"\_") if i % 2 == 0 else p for i, p in enumerate(parts))
    return t


def write(name, caption, header, rows, notes, colspec=None):
    md = [f"### {name}. {caption}", "", "| " + " | ".join(header) + " |",
          "|" + "|".join("---" for _ in header) + "|"]
    md += ["| " + " | ".join(r) + " |" for r in rows]
    md += [""] + [f"{n}" for n in notes] + [""]
    tex = [r"\begin{table*}[t]", r"\footnotesize", r"\caption{" + md_to_tex(caption) + "}",
           r"\label{tab:" + name.lower() + "}",
           r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}" + (colspec or "l" * len(header)) + "}", r"\toprule",
           " & ".join(md_to_tex(h) for h in header) + r" \\", r"\midrule"]
    tex += [" & ".join(md_to_tex(c) for c in r) + r" \\" for r in rows]
    tex += [r"\bottomrule", r"\end{tabular*}"]
    tex += [r"\par\footnotesize " + md_to_tex(n) for n in notes]
    tex += [r"\end{table*}", ""]
    return "\n".join(md), "\n".join(tex)


def table1():
    header = ["state", "$M_c$ (exact)", "upper bound", r"$M_c^{\rm BO}$", r"$M_c^{\rm BO+ad}$",
              r"$\Delta_{\rm na}$ (exact)", r"$\Delta_{\rm na}$ (pred.)"]
    rows = []
    for c in CRIT:
        lab = c["label"]
        pred = f"{c['pred']:.4f}"
        if c["prespec"]:
            fn = "$^a$" if c["key"] == "v4" else ("$^b$" if c["prespec"] is True else "$^c$")
            lab += fn; pred += fn
        rows.append([lab, fmt_pm(c["Mc"], c["err"]), c["ub"], f"{c['BO']:.7f}", f"{c['BOad']:.7f}",
                     f"{c['off']:.4f}", pred])
    rows.append([L1S["label"] + "$^d$", fmt_pm(L1S["Mc"], L1S["err"]), L1S["ub"], f"{L1S['BO']:.6f}", "–", "–", "–"])
    j, l3, l4 = (next(c for c in CRIT if c["key"] == k) for k in ("J1v1", "L3", "L4"))
    notes = [
        "Status: exact $M_c$ extrapolated (threshold-law fits for v = 1-4, L = 1, 2; geometric basis extrapolation of the "
        "crossing masses for L = 1 v = 1, L = 3, L = 4; lower side not rigorous); upper bound numerically certified at that mass "
        "(variational root below threshold, 1e-11 Eh margin for L = 0 and for L = 1 v = 1, L = 3, L = 4); BO and BO+ad from the exact BO / BO+adiabatic "
        "zero-energy problems (converged to ≤ 1e-7); offset = M_c(BO+ad) − M_c(exact); predicted offset: parameter-free "
        "(non-perturbative fold of the computed non-adiabatic masses m_v(R), m_r(R); no fitted parameter). "
        "Columns: upper bound = numerically certified upper bound to M_c (variational root below threshold at that mass; floating point, not interval arithmetic); Delta_na = non-adiabatic offset M_c(BO+ad) − M_c(exact).",
        f"$^a$ Pre-specified: prediction frozen 2026-11-14 17:22 before any exact v = 4 calculation; predicted "
        f"M_c = {V4_PRED_MC[0]} ± {V4_PRED_MC[1]}, exact {CRIT[3]['Mc']:.5f} ± {CRIT[3]['err']:.5f}.",
        f"$^b$ Pre-specified (second blind test): prediction frozen {LV_FROZEN} (sha256-recorded) before any exact crossing "
        f"was computed; predicted M_c = {LV_PRED_MC['J1v1'][0]} ± {LV_PRED_MC['J1v1'][1]}, exact {j['Mc']:.6f} "
        f"(difference {j['Mc'] - LV_PRED_MC['J1v1'][0]:+.4f}). Caveat: a scale optimisation of the J = 1, v = 1 basis (N = 12) "
        "ran 4 min before the freeze; it computed no crossing (its only information, binding at M = 17.52, was already "
        "implied by M_c < M_c(BO+ad)). Frolov & Wardlaw (arXiv:1008.3010) give a basis-limited upper estimate "
        f"M_c ≈ {FW_11} for this level, consistent with our bound.",
        f"$^c$ Secondary pre-specified predictions in the same frozen file: L = 3: {LV_PRED_MC['L3'][0]} ± {LV_PRED_MC['L3'][1]} "
        f"(exact {l3['Mc']:.6f}, difference {l3['Mc'] - LV_PRED_MC['L3'][0]:+.4f}); L = 4: {LV_PRED_MC['L4'][0]} ± "
        f"{LV_PRED_MC['L4'][1]} (exact {l4['Mc']:.6f}, difference {l4['Mc'] - LV_PRED_MC['L4'][0]:+.4f}).",
        f"$^d$ Exploratory: lowest exchange-symmetric J = 1 level (2pσu-like, bound by the long-range polarization well). "
        f"M_c from a geometric basis extrapolation that is not converged to the level of the other rows (spread "
        f"of the bases ≈ 0.1); upper bound = rounding-safe crossing of the compact N = 10 + 5 diffuse basis "
        f"(a larger basis crosses at {L1S['ub_nomargin']} without margin); BO = 2pσu BO critical mass; no BO+ad "
        "(adiabatic correction not available for 2pσu), hence no offset.",
        "Earlier estimates (finite exponential bases, no error bars; first location of these thresholds): lowest L = 1 4.219 (1989), 4.2133 (1994); "
        "lowest L = 2 11.485 (1989), bracket 11.489-11.490 (1994); v = 1 10.045, v = 2 30.70, v = 3 61.135, (1,1) 17.125 (1992); (1,1) 17.0944 (2011); "
        "Frolov and co-workers, see the main text. All are consistent with the values above. v = 4, L = 3 and L = 4 had not been computed.",
        "Sources: paper_data/critical_masses/, paper_data/critical_masses_v4/, paper_data/rotational/, "
        "paper_data/bo_offset/massfunc/fold_results.csv, paper_data/critical_masses_v4/prediction.json; summary §1, §6, §15; "
        "paper_data/thresholds_Lv/RESULT.md, *_crossings.csv, prediction_J1v1.json, bo_brackets.txt.",
    ]
    return write("T1", r"Critical masses of the excited bound states of $\{X^+,X^+,e^-\}$ ($M=m_X/m_e$)",
                 header, rows, notes)


def table2():
    header = ["$M$", "descriptor", "source", "type", "status"]
    kind = {"lit": "literature-proposed", "ours": "this work", "thr": "this work (threshold)"}
    rows = []
    for lo, hi, desc, src, k, status in sorted(MARKERS, key=lambda m: m[0]):
        Ms = f"{lo}–{hi}" if hi is not None else f"{lo}"
        if desc.startswith("minimum of g(M), fit-free"):
            Ms = "0.4924222(1)"
        if lo == 0.491:
            Ms = "0.491(3)"
        if lo == 0.635:
            Ms = "0.635(10)"
        if k == "thr":
            c = [c for c in CRIT if abs(c["Mc"] - lo) < 1e-9][0]
            Ms = fmt_pm(c["Mc"], c["err"])
        rows.append([Ms, desc, src, kind[k], status])
    notes = [
        "Status: located = extremum/zero located on our exact scan (smooth quantity; location uncertainty from the grid "
        "unless given); literature = value as published; converted = literature parameter converted to M by us, not "
        "recomputed; extrapolated = threshold-law fit (Table T1).",
        r"Crossover width (local Lindemann exponent $f=-4\,d\ln\mathcal{L}/d\ln M$, $\mathcal{L}=\Delta R/\langle R\rangle$): $f$ = 0.1/0.5/0.9 at M = 0.455/1.52/95 "
        "(2.3 decades, width 5.34 in ln M); ground-state markers 0.2–1.11 span 0.32 widths (summary §7); "
        "including the fidelity marker 1.39 they span 0.36 widths; literature-proposed 0.4–0.98 span 0.17.",
        "Sources: paper_data/mass_scans/, paper_data/king2013/, paper_data/fidelity/, paper_data/entanglement/, "
        "paper_data/crossover_width/; RESULTS_SUMMARY.md §2, §7–9, §12, §13.",
    ]
    return write("T2", "Crossover markers of the atom-to-molecule change", header, rows, notes,
                 colspec="l p{5.2cm} p{4.4cm} p{2.0cm} l")


def table3():
    levels = MUONIC_LEVELS
    header = ["molecule", "$M$", "predicted bound $(J,v)$", "literature bound $(J,v)$",
              "predicted and literature absent", "match", "closest threshold"]
    rows = []
    ntot = 0
    mpres, mabs = [], []
    for m in MUONIC:
        pred = [lv for lv in levels if m["M"] > LEVEL_THRESHOLDS[lv]]
        lit = sorted(m["lit"])
        assert all(lv in levels for lv in lit)
        absent = [lv for lv in levels if lv not in pred and lv not in m["lit"]]
        ok = sum((lv in pred) == (lv in m["lit"]) for lv in levels[1:])
        ntot += ok
        for lv in levels[1:]:
            r = m["M"] / LEVEL_THRESHOLDS[lv] - 1
            (mpres if r > 0 else mabs).append((abs(r), lv, m["name"].replace(r"$\mu$", "μ")))
        marg = min(((abs(m["M"] / LEVEL_THRESHOLDS[lv] - 1), lv) for lv in levels[1:]))
        f = lambda L: ", ".join(f"({a},{b})" for a, b in L) or "–"
        rows.append([m["name"].replace(r"$\mu$", "μ"), f"{m['M']:.4f}", f(pred), f(lit), f(absent),
                     f"{ok}/{len(levels) - 1}", f"({marg[1][0]},{marg[1][1]}): {100 * marg[0]:.1f}%"])
    thr = "; ".join(f"({a},{b}) M > {LEVEL_THRESHOLDS[(a, b)]}" for a, b in levels[1:])
    notes = [
        f"Predictions use our thresholds (Tables T1, T4): {thr}; (0,0) is bound at every M. Odd J are the "
        "exchange-antisymmetric levels (L 'A'), even J the symmetric ones (L 'S'). Every level in the three ions, and "
        "the absent (0,2) and (4,0), is now compared with a computed threshold: "
        f"{ntot}/{3 * (len(levels) - 1)} statements agree. "
        "Closest threshold: smallest relative distance |M/M_c − 1|; all distances exceed the uncertainty of M_c "
        "(≤ 1e-5 relative) by orders of magnitude.",
        "Status: level existence only (no energies compared). Presence of a level rests on the numerically certified "
        "upper bound (a bound state exists at that mass) together with the observed monotone growth of the binding "
        "energy with M; absence rests on the extrapolated M_c (lower side not rigorous). Smallest margins: present "
        f"({min(mpres)[1][0]},{min(mpres)[1][1]}) in {min(mpres)[2]} ({100 * min(mpres)[0]:.1f}% above M_c), absent ({min(mabs)[1][0]},{min(mabs)[1][1]}) in "
        f"{min(mabs)[2]} ({100 * min(mabs)[0]:.1f}% below M_c).",
        "Literature: Frolov & Wardlaw, EPJD 63, 339 (2011), arXiv:1008.3010; the same level counts (2, 5, 6 bound "
        "states) are restated in Frolov, arXiv:2503.17389 (2025). M = m_nucleus/m_μ with their masses.",
    ]
    return write("T3", "Muonic molecular ions: level existence predicted from our thresholds", header, rows, notes,
                 colspec="llp{3.3cm}p{3.3cm}p{2.6cm}lp{2.0cm}")


def fmt_mc(v, e, nd):
    """M_c at the precision printed in the source RESULT.md (nd decimals) with its quoted error."""
    return f"{v:.{nd}f} ± {e:.1e}"


def table4():
    header = ["$(J,v)$", "$M_c$ (exact)", "certified upper bound", r"$M_c^{\rm BO}$", r"$M_c^{\rm BO+ad}$",
              f"prediction (frozen, ±{MAP_UNC})", "exact − pred.", "PASS"]
    rows = []
    for c in RUN5:
        if c["onesided"]:
            mc = f"≤ {float(c['best']):.9f} (−{c['err']:.1e})$^b$"
        else:
            mc = fmt_mc(c["Mc"], c["err"], 6 if c["hp"] else 8) + ("$^c$" if c["key"] == "J4v1" else "") + \
                ("$^d$" if c["hp"] else "")
        rows.append([f"({c['J']},{c['v']})", mc, c["safe"], f"{c['BO']:.6f}", f"{c['BOad']:.6f}",
                     f"{c['pred']:.4f}", f"{c['diff']:+.4f}", "PASS" if c["passed"] else "FAIL"])
    npass = sum(c["passed"] for c in RUN5)
    md = max(abs(c["diff"]) for c in RUN5)
    notes = [
        f"$^a$ Blind test of the threshold map: all predictions frozen {MAP_FROZEN} (PREDICTION_MAP.md, sha256-recorded) "
        "before any exact calculation of these levels. Prediction = parameter-free non-perturbative fold of m_v(R), m_r(R) "
        f"into the BO+ad zero-energy problem (same code as the earlier blind tests); band ±{MAP_UNC} = largest |exact − pred.| "
        "of the nine levels known before the freeze (Table T1; level (1,0)); criterion |exact − pred.| ≤ band + exact error. "
        f"Result: {npass}/{len(RUN5)} PASS, largest |exact − pred.| = {md:.4f}.",
        "Status: exact M_c = geometric extrapolation in the compact basis size of the crossing masses (lower side not rigorous); "
        "certified upper bound = numerically certified upper bound (crossing of E = E_th − 1e-11 Eh, tol 1e-12, tol 1e-14 for (7,0); variational "
        "root below threshold at that mass; floating point, not interval arithmetic); BO and BO+ad from the exact BO / "
        "BO+adiabatic zero-energy problems (bo_map.json).",
        "$^b$ L = 6: the compact-size sequence is not monotone and the extrapolation is unreliable (for (6,1) it lies above "
        "the best crossing). Value = best crossing (rigorous upper bound up to rounding, margin 0), with a one-sided error "
        "= max(|best − extrapolated| + extrapolation error, largest compact-size step); the difference and PASS use this value "
        "(with the extrapolated value: −0.0013 and −0.0025, also PASS).",
        "$^c$ The extrapolated (4,1) value 57.16494981(46) lies 7.4e-7 above the best crossing 57.164949070 (and above the "
        "certified bound): the extrapolation error is underestimated by about a factor 2; irrelevant at the 0.03 level.",
        "$^d$ (7,0): the double-precision L = 7 basis failed validation (canonical-orthogonalization truncation at tol 1e-12, "
        "precision loss at 1e-14; threshold_map/validate_L7_diagnostics.md). Computed instead with double-double (about "
        "32-digit) assembly of the overlap and Hamiltonian matrices, single-pass canonical orthogonalization at tol 1e-14 and "
        "spurious-root guards (paper_data/threshold_map_L7hp; validate_L 7 passes; L = 6 control reproduces double precision to "
        "2e-10). Value = best crossing 87.5214689 (N = 12 + 2 diffuse) minus the geometric compact-size correction; "
        "certified bound at tol 1e-14 (not 1e-12).",
        "Sources: paper_data/threshold_map/RESULT.md (2026-11-17), data/*_crossings.csv, threshold_map_L7hp/data/J7v0_crossings.csv, "
        "prediction_map.json, bo_map.json.",
    ]
    return write("T4", r"Threshold map: blind test of the non-adiabatic offset theory on 14 new levels of $\{X^+,X^+,e^-\}$$^a$",
                 header, rows, notes, colspec="@{}l p{2.7cm} p{1.76cm} lll p{1.45cm} p{0.95cm} c@{}")


def table5():
    header = ["$q$", "$(J,v)$", "$M_c$ (exact)", "certified upper bound", "theory (frozen) ± band", "exact − theory",
              "PASS band / 0.10", "exact − null: NS / NT / NP", "null miss / theory miss"]
    rows = []
    for c in RUN6:
        nulls = " / ".join(f"{c['Mc'] - c[k]:+.2f}" if c[k] is not None else "–" for k in ("NS", "NT", "NP"))
        rows.append([f"{c['r']:g}", f"({c['J']},{c['v']})", fmt_mc(c["Mc"], c["err"], 6), c["safe"],
                     f"{c['pred']:.4f} ± {c['unc']:.3g}", f"{c['diff']:+.4f}",
                     f"{'PASS' if c['pass_unc'] else 'FAIL'} / {'PASS' if c['pass_len'] else 'FAIL'}", nulls,
                     f"{c['null_ratio']:.0f}"])
    nl = sum(c["pass_len"] for c in RUN6)
    rr = [c["null_ratio"] for c in RUN6]
    notes = [
        f"Unequal masses {{A+, B+, e−}}, m_B = q m_A; M_c = mass of the light nucleus A (in m_e). Theory = two-channel "
        f"(g/u) model with the non-adiabatic mass functions (2ch+m), predictions frozen {UNEQ_FROZEN} (sha256-recorded) "
        "before the first exact run. Band (frozen): model-variant spread + the symmetric-family residual of the same level scaled by (1+q)/(2q), floor 0.03 (0.0368 for (1,0)); lenient criterion 0.10 (also frozen).",
        f"Result: 0/{len(RUN6)} within the frozen band, {nl}/{len(RUN6)} within 0.10. Every theory residual is negative "
        f"(the theory overestimates M_c by 0.04-0.18). Null models: NS = BO+ad minus the exact symmetric-family offset of "
        "the same level; NT = the same from the exact threshold; NP = persistence (offset of the previous J = 0 level; "
        f"(0,2)-(0,4) only). The smallest null-model miss is {min(rr):.0f}-{max(rr):.0f} times the theory miss.",
        "Status: exact M_c = extrapolated (geometric in compact size); certified upper bound as in Table T4 (margin 1e-11 Eh). "
        f"Validation at q = 1 (unsymmetrized code): the production basis gives {R1_CHECK['production']} vs the symmetric "
        f"exact {R1_CHECK['exact_sym']}; with the mirror basis (diffuse ladder also in the closed channel) "
        f"{R1_CHECK['mirror']} ({R1_CHECK['mirror_big']} with a larger basis), i.e. +1.5e-5 from the exact value. "
        "The same enlargements shift the q = 2 and q = 1.5 (0,1) values by ≤ 1.2e-6 (L = 0 only; (1,0), (2,0) not covered).",
        "Sources: paper_data/unequal_masses/RESULT.md (2026-11-16 17:28), prediction_unequal.json, *_crossings.csv.",
    ]
    return write("T5", "Unequal masses: blind test of the two-channel offset theory (fails the frozen band; "
                 "much closer than null models)", header, rows, notes, colspec="ll p{2.6cm} p{1.6cm} p{1.7cm} p{1.2cm} p{1.4cm} p{2.6cm} p{1.2cm}")


if __name__ == "__main__":
    bad = check_against_data(verbose=False)
    if bad:
        print("DISCREPANCIES:", bad)
    mds, texs = [], []
    for t in (table1, table2, table3, table4, table5):
        m, x = t()
        mds.append(m); texs.append(x)
    open(os.path.join(OUT, "tables.md"), "w").write("# Paper 1 tables (generated by make_tables.py)\n\n" + "\n".join(mds))
    open(os.path.join(OUT, "tables.tex"), "w").write(
        "% Paper 1 tables, generated by make_tables.py; needs \\usepackage{booktabs}\n" + "\n".join(texs))
    print("wrote tables.md, tables.tex")
