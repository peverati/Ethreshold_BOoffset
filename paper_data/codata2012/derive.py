"""Compare CODATA-2012 physical-mass results (this folder) with the CODATA-2018 rows of ../mass_scans.
Writes results.csv (one row per quantity and mass; rebuilt from the CSVs each time, cheap) and prints it.
Spacings follow mass_scans/analyze.py (root-1 scales if lower; E01, E02 from exact_rot.csv)."""
import csv, math, os
HERE = os.path.dirname(os.path.abspath(__file__)); OLD = os.path.join(HERE, "..", "mass_scans")
HART = 219474.6313705
NEW_M = {"mu": 206.7682843, "p": 1836.15267245}; OLD_M = {"mu": 206.768283, "p": 1836.15267343}

def rd(d, f):
    p = os.path.join(d, f)
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []
def fl(x):
    try: return float(x)
    except Exception: return float("nan")
def pick(rows, M, key="M", tol=1e-6):
    for r in rows:
        if abs(fl(r[key]) - M) < tol * max(1, M): return r
def pickL(rows, M, L):
    for r in rows:
        if abs(fl(r["M"]) - M) < 1e-6 * M and int(r["L"]) == L: return r

def collect(d, M):
    out = {}
    ex = pick(rd(d, "exact_L0.csv"), M); bo = pick(rd(d, "bo.csv"), M)
    r1, r2 = pickL(rd(d, "exact_rot.csv"), M, 1), pickL(rd(d, "exact_rot.csv"), M, 2)
    if ex is None or bo is None: return None
    e0, e0err = fl(ex["E0_ext"]), fl(ex["E0_err"])
    e1, e1err = fl(ex["E1_ext"]), fl(ex["E1_err"])
    if math.isfinite(fl(ex.get("E1_opt"))) and fl(ex["E1_opt"]) < fl(ex["E1"]): e1, e1err = fl(ex["E1_opt_ext"]), fl(ex["E1_opt_err"])
    Eth = fl(ex["E_th"])
    out["E0_exact"] = (e0, e0err); out["E1_exact"] = (e1, e1err)
    for t in ("BO", "BOad"):
        out[f"E0_{t}"] = (fl(bo[f"{t}_E00"]), 0)
        out[f"E0_exact_minus_{t}"] = (e0 - fl(bo[f"{t}_E00"]), e0err)
    out["M*(E-E_BO)"] = (M * (e0 - fl(bo["BO_E00"])), M * e0err)
    out["M^1.5*(E-E_BOad)"] = (M ** 1.5 * (e0 - fl(bo["BOad_E00"])), M ** 1.5 * e0err)
    out["BO_levels_J0"] = (fl(bo["BO_n_J0"]), 0); out["BOad_levels_J0"] = (fl(bo["BOad_n_J0"]), 0)
    for k in ("R", "lindemann", "invR2"): out[f"exact_{k}"] = (fl(ex[k]), 0)
    out["exact_<R>_Egap_note"] = (float("nan"), 0)
    del out["exact_<R>_Egap_note"]
    for t in ("BO", "BOad"):
        for k in ("R", "dR", "lindemann", "invR2"):
            out[f"{t}_{k}"] = (fl(bo[f"{t}_{k}"]), 0)
    out["exact_dE_dlnM"] = (fl(ex["dE_dlnM"]), 0)
    # spacings in hartree
    sp = {}
    sp["vib"] = (e1 - e0, e1err + e0err)
    if r1: sp["d1"] = (fl(r1["E_ext"]) - e0, fl(r1["E_err"]) + e0err)
    if r2: sp["d2"] = (fl(r2["E_ext"]) - e0, fl(r2["E_err"]) + e0err)
    for t in ("BO", "BOad"):
        sp[f"{t}_vib"] = (fl(bo[f"{t}_E10"]) - fl(bo[f"{t}_E00"]), 0)
        sp[f"{t}_d1"] = (fl(bo[f"{t}_E01"]) - fl(bo[f"{t}_E00"]), 0)
        sp[f"{t}_d2"] = (fl(bo[f"{t}_E02"]) - fl(bo[f"{t}_E00"]), 0)
    for k, (v, e) in sp.items():
        out[f"{k}_hartree"] = (v, e); out[f"{k}_cm-1"] = (v * HART, e * HART)
    if "d1" in sp and "d2" in sp: out["d2/d1"] = (sp["d2"][0] / sp["d1"][0], 0); out["M*d1"] = (M * sp["d1"][0], M * sp["d1"][1])
    if "d1" in sp: out["vib/d1"] = (sp["vib"][0] / sp["d1"][0], 0)
    for t in ("BO", "BOad"):
        out[f"{t}_vib_relerr"] = (sp[f"{t}_vib"][0] / sp["vib"][0] - 1, 0)
        if "d1" in sp: out[f"{t}_d1_relerr"] = (sp[f"{t}_d1"][0] / sp["d1"][0] - 1, 0)
    out["rigid_rotor_2<1/R2>/M_hartree"] = (2 * fl(ex["invR2"]) / M, 0)
    out["rigid_rotor_2<1/R2>/M_cm-1"] = (2 * fl(ex["invR2"]) / M * HART, 0)
    out["B0_rigid_<1/R2>/(2muN)_cm-1"] = (fl(ex["invR2"]) / M * HART, 0)
    out["BO_<1/R2>/(2muN)_cm-1"] = (fl(bo["BO_invR2"]) / M * HART, 0)
    return out

rows = []
for lab in ("mu", "p"):
    new, old = collect(HERE, NEW_M[lab]), collect(OLD, OLD_M[lab])
    if new is None: continue
    for k, (v, e) in new.items():
        o = old.get(k, (float("nan"), 0))[0] if old else float("nan")
        rows.append(dict(mass=lab, M_old=OLD_M[lab], M_new=NEW_M[lab], quantity=k, old=o, new=v, diff=v - o, err_new=e))
with open(os.path.join(HERE, "results.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader()
    for r in rows: w.writerow({k: (f"{v:.15g}" if isinstance(v, float) else v) for k, v in r.items()})
for r in rows: print(f"{r['mass']:3s} {r['quantity']:34s} old {r['old']:+.12g} new {r['new']:+.12g} diff {r['diff']:+.3e} err {r['err_new']:.1e}")
