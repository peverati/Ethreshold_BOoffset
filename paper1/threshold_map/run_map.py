"""Exact (non-BO) critical masses of every (J, v) target of the threshold map of {X+, X+, e-}.

Written 2026-11-16 (threshold_map, step 3).  Generalizes paper1/thresholds_Lv/run_exact.py (J1v1, L3A, L4S: fixed
basis => mass-independent matrices; Brent root of E_k(M) - E_th(M); each fixed-basis crossing M* is a rigorous
upper bound to M_c up to rounding) to any target (J, v) = root k = v of the natural-parity L = J sector with exchange
symmetry (-1)^J.  Targets, BO and BO+ad values come from bo_map.json.  L >= 5 term lists live in a local cache
(tb_local.py; solvers/ is not modified).

Run from the project root (BLAS threads: THREEBODY_BLAS_THREADS, default 8):
    python3 paper1/threshold_map/run_map.py --list               # targets in priority order
    ... run_map.py --target J1v2                  # all stages for one target (resumable; done rows skipped)
    ... run_map.py --target J1v2 --dry            # plan: sizes, memory, time estimate from a calibration solve
    ... run_map.py --target J1v2 --test N ND M    # one E_k - E_th at mass M (short test, not stored)
    ... run_map.py --target J1v2 --stage 2 3      # some stages only
Stages (rows appended to data/<T>_crossings.csv, data/<T>_binding.csv at once):
  1 scales   optimize compact scale sets at M_opt (root 0, and root k for k > 0)      -> data/<T>_scales.json
  2 compact  nested compact sequence N = Ns (all scale sets)                          -> crossings
  3 diffuse  compact N_D + nested diffuse sets (atom x X-atom relative motion)       -> crossings
  4 checks   largest basis at tol 1e-10 and margin 1e-11 Eh (rounding-safe bound); compact N_D - 2, N_D + 2
             with all diffuse sets (geometric extrapolation in N)
  5 binding  E_th - E_k at M = M*_best (1 + x), largest basis (L > 0: binding linear in M - M_c)
Root-ordering check: at every crossing the roots below k must be bound (E_j < E_th, j < k); the gaps of roots k-1
and k+1 at M* are stored in the crossings CSV (columns gap_below, gap_above).
Log: data/<T>_log.txt.  Then analyze_map.py -> RESULT.md.
"""
import argparse
import csv
import gc
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
os.environ.setdefault("THREEBODY_BLAS_THREADS", "8")
sys.path.insert(0, HERE)
import tb_local  # noqa: E402,F401  (local angular term cache; must precede threebody computations)
import numpy as np  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

from threebody import symmetric_family, family_threshold, angular  # noqa: E402
from threebody.coupled import CoupledBasis, IdxBlock, block  # noqa: E402
from threebody.rotational import optimize_scales  # noqa: E402
from threebody.threshold import diffuse_index, total_degree_index  # noqa: E402

BOMAP = json.load(open(os.path.join(HERE, "bo_map.json"))) if os.path.exists(os.path.join(HERE, "bo_map.json")) else None
DIFFUSE = ((4, 12, 0.4), (4, 16, 0.15), (5, 16, 0.05), (5, 20, 0.02), (6, 20, 0.008))
XB = (0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2)
SIZE_CAP = 8000             # largest basis (functions) of any stage: ~12 n^2 doubles = 6.1 GB
FROLOV_FIRST = ("J1v2", "J2v1", "J2v2")   # Frolov 1992 II levels: run first


# ---------------------------------------------------------------- target configuration
def config(name):
    """Basis plan of a target, from L and k (rules justified in README/PROGRESS):
    * scale sets: k = 0: root-0 optimum + 0.6 alpha_X copy (as L3A/L4S); k = 1: root-0 and root-1 optima (as J1v1);
      k >= 2: root-0, root-k optima and their geometric mean (more extended vibrational wave function).
    * diffuse sets: L = 1: 5 (as J1v1); L = 2: 3; L = 3: 3; L >= 4: 2 (for L >= 3 they moved M* by < 1e-6 in
      thresholds_Lv: the centrifugal barrier confines the level).
    * N_D: the largest even N such that the biggest check basis (N_D + 2 with all diffuse sets) stays <= SIZE_CAP;
      compact sequence N_D - 4 .. N_D + 2 (nested).
    * scale optimization: N_opt = 10 for L <= 2, 8 for L >= 3 (cost: N = 12 for L = 2 took 1441 s for 70 evaluations on
      a loaded machine; J2v1 itself was optimized at N = 12 during the short test, not a result); M_opt = BO+ad M_c + dM, with
      dM = 0 for k <= 1 (as before) and 0.03 M for k >= 2 (root k must be bound in the small basis).
    """
    L, k = int(name[1:].split("v")[0]), int(name.split("v")[1])
    lev = {(d["J"], d["v"]): d for d in BOMAP["levels"]}[(L, k)]
    nscale = 2 if k <= 1 else 3
    nd = {0: 0, 1: 5, 2: 3, 3: 3}.get(L, 2)
    cfg = dict(name=name, L=L, k=k, sym="A" if L % 2 else "S", M_BO=lev["M_BO"], M_BOad=lev["M_BOad"],
               nscale=nscale, diffuse=DIFFUSE[:nd], N_opt=10 if L <= 2 else 8,
               M_opt=lev["M_BOad"] * (1.0 if k <= 1 else 1.03))
    ND = 8
    for N in range(8, 21, 2):
        if basis_size(cfg, [(1.0, 1.0, 1.0)] * nscale, N + 2, nd)[0] <= SIZE_CAP:
            ND = N
    cfg["ND"] = min(ND, 18)
    cfg["Ns"] = tuple(n for n in (cfg["ND"] - 4, cfg["ND"] - 2, cfg["ND"], cfg["ND"] + 2) if n >= 6)
    return cfg


def targets_in_order():
    t = list(BOMAP["targets"])
    first = [x for x in FROLOV_FIRST if x in t]
    rest = sorted([x for x in t if x not in first], key=lambda x: config(x)["M_BOad"])
    return first + rest


# ---------------------------------------------------------------- IO
class Table:
    def __init__(self, path, cols):
        self.path, self.cols = path, cols
        self.rows = list(csv.DictReader(open(path))) if os.path.exists(path) else []

    def find(self, **kw):
        return [r for r in self.rows if all(str(r.get(k)) == str(v) for k, v in kw.items())]

    def add(self, **row):
        row["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
        new = not os.path.exists(self.path)
        with open(self.path, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=self.cols)
            if new:
                w.writeheader()
            w.writerow({k: (repr(v) if isinstance(v, float) else v) for k, v in row.items()})
        self.rows.append({k: (repr(v) if isinstance(v, float) else str(v)) for k, v in row.items()})


CROSS_COLS = ["target", "label", "stage", "size", "kept", "tol", "margin", "M_star", "gap_below", "gap_above",
              "seconds", "finished"]
BIND_COLS = ["target", "label", "M", "x", "E", "E_th", "binding", "E_next", "seconds", "finished"]
_logs = {}


def path(t, kind):
    os.makedirs(DATA, exist_ok=True)
    return os.path.join(DATA, f"{t}_{kind}")


def log(t, msg):
    if t not in _logs:
        _logs[t] = open(path(t, "log.txt"), "a")
    line = time.strftime("%Y-%m-%d %H:%M:%S ") + msg
    print(line, flush=True)
    _logs[t].write(line + "\n")
    _logs[t].flush()


# ---------------------------------------------------------------- basis
def diffuse_pair(alpha_atom, alpha_d, n_c, n_d):
    idx = diffuse_index(n_c, n_d)[:, [1, 0, 2]]  # diffuse in s0: X0 far from the (X1 e) atom
    B = IdxBlock((alpha_d, alpha_atom, alpha_atom), idx, name=f"D{n_c}x{n_d}@{alpha_d:g}")
    return [B, B.mirror()]


def make_basis(cfg, scales, N, n_diffuse, mu):
    blocks = [block(N, al, name=f"N{N}s{i}") for i, al in enumerate(scales)]
    for (n_c, n_d, ad) in cfg["diffuse"][:n_diffuse]:
        blocks += diffuse_pair(mu, ad, n_c, n_d)
    return CoupledBasis(cfg["L"], blocks, cfg["sym"])


def basis_size(cfg, scales, N, n_diffuse):
    """Size without building (self-mirror compact blocks halve on the self-swap channel)."""
    chans = [c for c in angular.components(cfg["L"]) if c <= angular.swap_component(c)]
    a = total_degree_index(N, (1, 1, 1))
    n = 0
    for c in chans:
        for al in scales:
            if angular.swap_component(c) == c and al[0] == al[1]:
                n += int((a[:, 0] >= a[:, 1]).sum() if cfg["sym"] == "S" else (a[:, 0] > a[:, 1]).sum())
            else:
                n += len(a)
        for (n_c, n_d, ad) in cfg["diffuse"][:n_diffuse]:
            n += 2 * len(diffuse_index(n_c, n_d))
    return n, len(chans)


SPURIOUS = -0.15   # E_0 - E_th below this is impossible for M <= ~110 (ground-state binding < 0.1 Eh): numerical noise


def gaps(cfg, b, M, tol=1e-12, nroots=None):
    """E_j - E_th, j = 0..nroots-1.  Guard: at high L the overlap matrix has noise eigenvalues up to ~1e-12 of the
    largest (measured: lam_min = -2.3e-12 for an L = 5 H2+ basis); a noise vector kept by the projection produces a
    spurious very low root (E = -2.4 seen at tol 1e-14).  Such a root would fake a crossing, so it raises instead."""
    k = cfg["k"]
    e = b.solve(symmetric_family(M), nroots or k + 1, tol, vectors=False)[0]
    g = e - family_threshold(M)
    if g[0] < SPURIOUS:
        raise ValueError(f"spurious root: E_0 - E_th = {g[0]:.3e} at M = {M:.6f} (overlap noise; tol {tol:g})")
    return g


def crossing(cfg, b, lo, hi, tol=1e-12, margin=0.0):
    k = cfg["k"]
    f = lambda M: gaps(cfg, b, M, tol)[k] + margin
    flo, fhi = f(lo), f(hi)
    tries = 0
    while fhi > 0 and tries < 4:      # widen upwards if the basis is still too small
        lo, flo = hi, fhi
        hi *= 1.1
        fhi = f(hi)
        tries += 1
    if not (flo > 0 > fhi):
        raise ValueError(f"no bracket: f({lo:.6f})={flo:.3e} f({hi:.6f})={fhi:.3e}")
    return brentq(f, lo, hi, xtol=1e-11, rtol=1e-13)


def bracket(cfg):
    # exact M_c lay 0.25-0.54 below BO+ad for all nine known levels; a basis crossing lies above the exact M_c,
    # so lo = BO+ad - 1.5 is unbound in any basis (f(lo) > 0, checked) and hi is widened if needed
    return cfg["M_BOad"] - 1.5, cfg["M_BOad"] + 0.3


# ---------------------------------------------------------------- stages
def stage_scales(t):
    p = path(t, "scales.json")
    if os.path.exists(p):
        return json.load(open(p))
    cfg = config(t)
    fam = symmetric_family(cfg["M_opt"])
    t0 = time.perf_counter()
    a0, e0, n0 = optimize_scales(fam, cfg["L"], cfg["N_opt"], cfg["sym"], state=0)
    log(t, f"scales root 0 at M={cfg['M_opt']:.6f} N={cfg['N_opt']}: {a0} E={e0:.12f} ({n0} evals, "
           f"{time.perf_counter() - t0:.0f}s)")
    out = dict(M_opt=cfg["M_opt"], scales=[list(a0)])
    if cfg["k"] > 0:
        t0 = time.perf_counter()
        ak, ek, nk = optimize_scales(fam, cfg["L"], cfg["N_opt"], cfg["sym"], state=cfg["k"], x0=a0)
        gk = ek - family_threshold(cfg["M_opt"])
        log(t, f"scales root {cfg['k']}: {ak} E={ek:.12f} E-E_th={gk:.3e} ({nk} evals, {time.perf_counter() - t0:.0f}s)"
               + ("  WARNING: root k unbound at M_opt in the optimization basis" if gk > 0 else ""))
        out["scales"].append(list(ak))
        if cfg["nscale"] == 3:
            out["scales"].append([float(np.sqrt(a0[0] * ak[0])), float(np.sqrt(a0[1] * ak[1])), float(np.sqrt(a0[2] * ak[2]))])
    else:
        out["scales"].append([0.6 * a0[0], 0.6 * a0[1], a0[2]])
    log(t, f"scale sets: {out['scales']}")
    json.dump(out, open(p, "w"), indent=1)
    return out


def run_crossing(t, T, label, stage, N, n_diffuse, sc, lo, hi, tols=(1e-12,), margins=(0.0,)):
    cfg = config(t)
    todo = [(tol, m) for tol in tols for m in margins if not T.find(target=t, label=label, tol=tol, margin=m)]
    if not todo:
        return [float(r["M_star"]) for r in T.find(target=t, label=label, tol=1e-12, margin=0.0)]
    mu = sc["M_opt"] / (sc["M_opt"] + 1)
    b = make_basis(cfg, sc["scales"], N, n_diffuse, mu)
    out = []
    for tol, m in todo:
        t1 = time.perf_counter()
        pr = b.project(tol)
        gb = ga = float("nan")
        try:
            Ms = crossing(cfg, b, lo, hi, tol, m)
            g = gaps(cfg, b, Ms, tol, cfg["k"] + 2)
            gb = float(g[cfg["k"] - 1]) if cfg["k"] > 0 else float("nan")
            ga = float(g[cfg["k"] + 1])
            if cfg["k"] > 0 and not gb < 0:
                log(t, f"  ROOT-ORDER WARNING {label}: root {cfg['k'] - 1} not bound at M* (gap {gb:.3e})")
        except ValueError as exc:
            log(t, f"  {label} tol={tol:g} margin={m:g}: {exc}")
            Ms = float("nan")
        T.add(target=t, label=label, stage=stage, size=b.size, kept=pr["kept"], tol=tol, margin=m, M_star=Ms,
              gap_below=gb, gap_above=ga,
              seconds=round(time.perf_counter() - t1 + (b.build_time if (tol, m) == todo[0] else 0), 1))
        log(t, f"stage{stage} {label}: size={b.size} kept={pr['kept']} lam_min/lam_max={pr['lam_min'] / pr['lam_max']:.1e} "
               f"tol={tol:g} margin={m:g} M*={Ms:.10f} "
               f"gap_below={gb:.2e} ({time.perf_counter() - t1:.0f}s; build {b.build_time:.0f}s)")
        out.append(Ms)
    del b
    gc.collect()
    return out


def stage_compact(t):
    cfg, sc = config(t), stage_scales(t)
    T = Table(path(t, "crossings.csv"), CROSS_COLS)
    lo, hi = bracket(cfg)
    for N in cfg["Ns"]:
        Ms = run_crossing(t, T, f"compact N={N}", 2, N, 0, sc, lo, hi)
        if Ms and np.isfinite(Ms[0]):
            hi = Ms[0] * 1.0000001   # nested sequence: next crossing lies below


def stage_diffuse(t):
    cfg, sc = config(t), stage_scales(t)
    T = Table(path(t, "crossings.csv"), CROSS_COLS)
    lo, hi = bracket(cfg)
    r = T.find(target=t, label=f"compact N={cfg['ND']}", tol=1e-12, margin=0.0)
    if r and np.isfinite(float(r[0]["M_star"])):
        hi = float(r[0]["M_star"]) * 1.0000001
    for nd in range(1, len(cfg["diffuse"]) + 1):
        Ms = run_crossing(t, T, f"compact N={cfg['ND']} + {nd} diffuse", 3, cfg["ND"], nd, sc, lo, hi)
        if Ms and np.isfinite(Ms[0]):
            hi = Ms[0] * 1.0000001


def full_label(cfg, N):
    nd = len(cfg["diffuse"])
    return f"compact N={N}" + (f" + {nd} diffuse" if nd else "")


def stage_checks(t):
    cfg, sc = config(t), stage_scales(t)
    T = Table(path(t, "crossings.csv"), CROSS_COLS)
    lo, hi = bracket(cfg)
    nd = len(cfg["diffuse"])
    run_crossing(t, T, full_label(cfg, cfg["ND"]), 4, cfg["ND"], nd, sc, lo, hi, tols=(1e-12, 1e-10), margins=(0.0, 1e-11))
    for N in (cfg["ND"] - 2, cfg["ND"] + 2):
        run_crossing(t, T, full_label(cfg, N), 4, N, nd, sc, lo, hi)


def best_mstar(t):
    T = Table(path(t, "crossings.csv"), CROSS_COLS)
    v = [float(r["M_star"]) for r in T.rows if r["tol"] == repr(1e-12) and r["margin"] == repr(0.0)]
    v = [x for x in v if np.isfinite(x)]
    if not v:
        raise RuntimeError("no crossing yet")
    return min(v)


def stage_binding(t):
    cfg, sc = config(t), stage_scales(t)
    B = Table(path(t, "binding.csv"), BIND_COLS)
    nd = len(cfg["diffuse"])
    label = full_label(cfg, cfg["ND"])
    Ms = best_mstar(t)
    todo = [x for x in XB if not B.find(target=t, label=label, x=x)]
    if not todo:
        return
    mu = sc["M_opt"] / (sc["M_opt"] + 1)
    b = make_basis(cfg, sc["scales"], cfg["ND"], nd, mu)
    k = cfg["k"]
    for x in todo:
        t1 = time.perf_counter()
        M = Ms * (1 + x)
        e = b.solve(symmetric_family(M), k + 2, 1e-12, vectors=False)[0]
        Eth = family_threshold(M)
        B.add(target=t, label=label, M=M, x=x, E=float(e[k]), E_th=Eth, binding=float(Eth - e[k]),
              E_next=float(e[k + 1]), seconds=round(time.perf_counter() - t1, 1))
        log(t, f"stage5 M={M:.8f} E_th-E_k={Eth - e[k]:.6e} ({time.perf_counter() - t1:.0f}s)")
    del b
    gc.collect()


# ---------------------------------------------------------------- plan / dry run
def plan(cfg):
    nd = len(cfg["diffuse"])
    p = [(N, 0) for N in cfg["Ns"]] + [(cfg["ND"], i) for i in range(1, nd + 1)]
    if nd:
        p += [(cfg["ND"] - 2, nd), (cfg["ND"] + 2, nd)]
    return p


def dry(t, calibrate=True):
    cfg = config(t)
    log(t, f"=== --dry {t}: L={cfg['L']} sym={cfg['sym']} root k={cfg['k']}; BO {cfg['M_BO']:.6f}, BO+ad "
           f"{cfg['M_BOad']:.6f}; bracket {bracket(cfg)}; N_opt {cfg['N_opt']} M_opt {cfg['M_opt']:.4f}; Ns {cfg['Ns']} "
           f"ND {cfg['ND']} diffuse {len(cfg['diffuse'])}; threads {os.environ.get('THREEBODY_BLAS_THREADS')}")
    sc = json.load(open(path(t, "scales.json")))["scales"] if os.path.exists(path(t, "scales.json")) else \
        [[1.0, 1.0, 1.0], [0.6, 0.6, 1.0], [0.8, 0.8, 1.0]][:cfg["nscale"]]
    tb = te = n0 = None
    if calibrate:
        t0 = time.perf_counter()
        b = make_basis(cfg, sc, cfg["ND"], 0, 0.95)
        tb = time.perf_counter() - t0
        t1 = time.perf_counter()
        b.project(1e-12)
        gaps(cfg, b, cfg["M_BOad"])
        te = time.perf_counter() - t1
        n0 = b.size
        log(t, f"calibration: compact N={cfg['ND']}, size {n0}: build {tb:.1f}s, project+1 solve {te:.1f}s")
        del b
    total, peak = 0.0, 0
    for N, k in plan(cfg):
        n, nch = basis_size(cfg, sc, N, k)
        reps = 4 if (N, k) == (cfg["ND"], len(cfg["diffuse"])) else 1
        peak = max(peak, n)
        if calibrate:
            tt = tb * (n / n0) ** 2 + 20 * te * (n / n0) ** 3
            total += tt * reps
            log(t, f"  N={N:2d} + {k} diffuse: size {n:6d} ({nch} ch) est {tt / 60:6.1f} min x{reps}, "
                   f"memory ~{12 * n * n * 8 / 1e9:.1f} GB")
        else:
            log(t, f"  N={N:2d} + {k} diffuse: size {n:6d} ({nch} ch), memory ~{12 * n * n * 8 / 1e9:.1f} GB")
    if calibrate:
        n, _ = basis_size(cfg, sc, cfg["ND"], len(cfg["diffuse"]))
        total += len(XB) * te * (n / n0) ** 3 + tb * (n / n0) ** 2
        log(t, f"TOTAL estimated {total / 3600:.2f} h (excl. stage 1 scales); peak size {peak}, "
               f"peak memory ~{12 * peak * peak * 8 / 1e9:.1f} GB (12 n^2 doubles)")
    return peak


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target")
    ap.add_argument("--list", action="store_true", help="print the targets in run order and exit")
    ap.add_argument("--stage", type=int, nargs="*", default=[1, 2, 3, 4, 5])
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--nocal", action="store_true", help="--dry without the calibration solve")
    ap.add_argument("--test", nargs=3, metavar=("N", "NDIFF", "M"), help="one E_k - E_th (short test, not stored)")
    a = ap.parse_args()
    if a.list:
        for t in targets_in_order():
            c = config(t)
            print(f"{t}  L={c['L']} '{c['sym']}' root {c['k']}  BO {c['M_BO']:.4f}  BO+ad {c['M_BOad']:.4f}  "
                  f"ND {c['ND']} Ns {c['Ns']} diffuse {len(c['diffuse'])} scale sets {c['nscale']}")
        return
    t = a.target
    if t not in BOMAP["targets"]:
        sys.exit(f"unknown target {t}; choose from {BOMAP['targets']}")
    if a.dry:
        return dry(t, not a.nocal)
    if a.test:
        cfg = config(t)
        sc = stage_scales(t)
        N, nd, M = int(a.test[0]), int(a.test[1]), float(a.test[2])
        t0 = time.perf_counter()
        b = make_basis(cfg, sc["scales"], N, nd, sc["M_opt"] / (sc["M_opt"] + 1))
        pr = b.project(1e-12)
        g = gaps(cfg, b, M, 1e-12, cfg["k"] + 2)
        log(t, f"TEST (not a result) N={N} +{nd} diffuse M={M}: size={b.size} kept={pr['kept']} "
               f"E_j-E_th (j=0..k+1) = {', '.join(f'{x:.4e}' for x in g)} (build {b.build_time:.0f}s, "
               f"total {time.perf_counter() - t0:.0f}s)")
        return
    log(t, f"=== run_map.py --target {t} stages={a.stage} threads={os.environ.get('THREEBODY_BLAS_THREADS')}")
    funcs = {1: stage_scales, 2: stage_compact, 3: stage_diffuse, 4: stage_checks, 5: stage_binding}
    t0 = time.perf_counter()
    for s in a.stage:
        t1 = time.perf_counter()
        funcs[s](t)
        log(t, f"stage {s} finished ({time.perf_counter() - t1:.0f}s)")
    log(t, f"ALL REQUESTED STAGES DONE ({time.perf_counter() - t0:.0f}s)")


if __name__ == "__main__":
    main()
