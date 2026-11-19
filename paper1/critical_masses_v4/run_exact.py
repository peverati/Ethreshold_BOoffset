"""Exact (non-BO) critical mass of the 5th symmetric S level (v = 4, root k = 4) of {X+, X+, e-}.

Written 2026-11-15.  Same methodology as paper1/critical_masses (k = 1..3, examples/critical_mass.py),
with the basis adapted to the more diffuse v = 4 state and with a resumable, memory-lean driver.
The blind prediction (PREDICTION.md) is NOT read here.

Run from the project root:
    THREEBODY_BLAS_THREADS=8 python3 paper1/critical_masses_v4/run_exact.py          # everything
    ... run_exact.py --dry            # plan + time estimate from a one-point calibration (~1-2 min)
    ... run_exact.py --stage 1 2      # only some stages
    ... run_exact.py --test VARIANT M # one energy E_4 - E_th at mass M for a basis variant (short test)

Stages (in priority order; every result row is appended + flushed at once, done rows are skipped):
  1 crossings   fixed-basis crossing M*(alpha_s) for the ladder           -> crossings_state4.csv
                (each a rigorous upper bound to M_c up to rounding)
  2 energies    E_4 - E_th and observables at M = M*_best + DELTAS, every
                ladder basis (the fit data)                               -> energies_state4.csv, observables_all.csv
  3 checks      convergence: compact N 16 -> 20 (two most diffuse compact
                blocks, then all four), diffuse nc 4 -> 6, nd 50/30 -> 70/45, tol 1e-10   -> convergence_state4.csv (+ crossings)
  4 bounds      rounding-safe upper bounds, margins 0, 1e-11, 1e-10       -> upper_bounds.csv
Log: run_log.txt.  Then: analyze.py (fits, error bar, comparison with PREDICTION.md, RESULT.md).

Basis (justified by the short tests in PROGRESS.md):
  * 4 compact total-degree blocks, N = 16, scales (alpha_X, alpha_e):
      ground-state optimum at M = 100 (N = 12):            (4.472, 1.365)
      two intermediate sets (geometric in alpha_X):       (2.30, 1.25), (1.20, 1.15)
      root-4 optimum at M = 100 (N = 12, single scale):   (0.717, 1.104)
    The k = 3 basis had 3 compact blocks; v = 4 extends to larger R, so one more intermediate set is used.
  * 2 diffuse blocks (mu_a, alpha_s, mu_a) with a0 + a2 <= 4, a1 <= 50 (k = 1..3: 40), and
    (mu_a, 0.03, mu_a) with a1 <= 30; mu_a = M_ref/(M_ref + 1), M_ref = 99.
  * ladder alpha_s in {1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 5e-5, 3e-5} (k = 1..3 ended at 3e-5, which was resolution
    limited with a1 <= 40; here a1 <= 50; the lowest root over the ladder is used at every M).
  Near threshold the halo is set by kappa ~ c1 x sqrt(2 mu); with c1 extrapolated from k = 1..3 (about 1e-3)
  x = 0.003 gives 1/kappa ~ 3e4 bohr, so the mass grid starts at x = 0.003 (E_b ~ 1e-11) instead of 0.001.
"""
import argparse
import csv
import gc
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
os.environ.setdefault("THREEBODY_BLAS_THREADS", "8")
sys.path.insert(0, os.path.join(ROOT, "solvers", "threebody"))
import numpy as np  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

from threebody import family_threshold  # noqa: E402
from threebody.threshold import MultiBlockBasis, compact_block, diffuse_block  # noqa: E402

K = 4
OPS = ("S", "TX", "Te", "V01", "VXe", "r01", "r01^2", "rXe")
CONFIG = dict(M_ref=99.0,
              compact=[(4.472, 1.365), (2.30, 1.25), (1.20, 1.15), (0.717, 1.104)],
              M_lo=97.0, M_hi=99.5)
LADDER = [1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 5e-5, 3e-5]
DELTAS = [0.003, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0]
STD = dict(Nc=16, nc=4, nd=50, nd2=30, Nc_big=None)
CHECKS = [("compact N=20", dict(Nc_big=20)),          # N = 20 for the two most diffuse compact blocks
          ("compact N=20 all", dict(Nc=20)),        # N = 20 for all four compact blocks
          ("diffuse nc=6", dict(nc=6)),
          ("diffuse nd=70, nd2=45", dict(nd=70, nd2=45))]
CHECK_ALPHAS = (1e-3, 1e-4)
TOLS = (1e-12, 1e-10)
MARGINS = (0.0, 1e-11, 1e-10)

FILES = dict(cross="crossings_state4.csv", en="energies_state4.csv", obs="observables_all.csv",
             conv="convergence_state4.csv", ub="upper_bounds.csv")
COLS = dict(
    cross=["state", "alpha_s", "variant", "size", "kept", "dropped", "M_star", "build_s", "total_s", "finished"],
    en=["state", "M", "alpha_s", "variant", "E", "E_minus_Eth", "seconds", "finished"],
    obs=["state", "M", "alpha_s", "E", "E_th", "binding", "r01", "r01_2", "rXe", "dE_dM", "finished"],
    conv=["state", "M", "variant", "alpha_s", "size", "kept", "E_minus_Eth", "diff_from_ref", "finished"],
    ub=["state", "variant", "alpha_s", "tol", "margin", "M_bound", "slope_at_Mstar", "finished"],
)


# ---------------------------------------------------------------- basis
def basis(alpha_s, Nc=16, nc=4, nd=50, nd2=30, Nc_big=None, compact=None):
    c = CONFIG
    mu = c["M_ref"] / (c["M_ref"] + 1)
    comp = compact if compact is not None else c["compact"]
    blocks = []
    for i, (ax, ae) in enumerate(comp):
        N = Nc_big if (Nc_big and i >= len(comp) - 2) else Nc
        blocks.append(compact_block(N, ax, ae))
    blocks += [diffuse_block(mu, alpha_s, nc, nd), diffuse_block(mu, 0.03, nc, nd2)]
    return MultiBlockBasis(blocks, ops=OPS)


def slim(mb, tols=(1e-12,)):
    """Project at the needed tolerances, then drop the raw n x n matrices (memory)."""
    keep = {}
    for t in tols:
        mb._proj = None
        keep[t] = mb.project(t)
    mb.mats = {}
    gc.collect()
    mb._projs = keep
    return mb


def energies(mb, M, tol=1e-12):
    if hasattr(mb, "_projs"):
        mb._proj = mb._projs[tol]
    return mb.energies(M, K + 1, tol)


def gap(mb, M, tol=1e-12):
    return energies(mb, M, tol)[K] - family_threshold(M)


def crossing(mb, tol=1e-12, margin=0.0, lo=None, hi=None):
    lo = CONFIG["M_lo"] if lo is None else lo
    hi = CONFIG["M_hi"] if hi is None else hi
    f = lambda M: gap(mb, M, tol) + margin
    flo, fhi = f(lo), f(hi)
    if not (flo > 0 > fhi):
        raise ValueError(f"no sign change: f({lo})={flo:.3e}, f({hi})={fhi:.3e}")
    return brentq(f, lo, hi, xtol=1e-10, rtol=1e-14)


# ---------------------------------------------------------------- resumable CSV
class Table:
    def __init__(self, key):
        self.path = os.path.join(HERE, FILES[key])
        self.cols = COLS[key]
        self.rows = list(csv.DictReader(open(self.path))) if os.path.exists(self.path) else []

    def find(self, **kw):
        return [r for r in self.rows if all(_eq(r.get(k), v) for k, v in kw.items())]

    def add(self, **row):
        row["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
        new = not os.path.exists(self.path)
        with open(self.path, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=self.cols)
            if new:
                w.writeheader()
            w.writerow({k: (repr(v) if isinstance(v, float) else v) for k, v in row.items()})
        self.rows.append({k: str(v) if not isinstance(v, float) else repr(v) for k, v in row.items()})


def _eq(a, b):
    if a is None:
        return False
    if isinstance(b, float):
        try:
            return abs(float(a) - b) <= 1e-12 * max(1.0, abs(b))
        except ValueError:
            return False
    return str(a) == str(b)


_logf = None


def log(msg):
    global _logf
    if _logf is None:
        _logf = open(os.path.join(HERE, "run_log.txt"), "a")
    line = time.strftime("%Y-%m-%d %H:%M:%S ") + msg
    print(line, flush=True)
    _logf.write(line + "\n")
    _logf.flush()


# ---------------------------------------------------------------- stages
def stage_crossings():
    T = Table("cross")
    for a in LADDER:
        if T.find(alpha_s=a, variant="standard"):
            continue
        t0 = time.perf_counter()
        mb = basis(a)
        pr = mb.project(1e-12)
        slim(mb)
        try:
            Ms = crossing(mb)
        except ValueError as exc:
            Ms = float("nan")
            log(f"  alpha_s={a:g}: {exc}")
        T.add(state=K, alpha_s=a, variant="standard", size=mb.size, kept=pr["kept"], dropped=pr["dropped"],
              M_star=Ms, build_s=round(mb.build_time, 1), total_s=round(time.perf_counter() - t0, 1))
        log(f"stage1 alpha_s={a:g} size={mb.size} kept={pr['kept']} M*={Ms:.9f} ({time.perf_counter() - t0:.0f}s)")
        del mb
        gc.collect()
    return T


def best_Mstar():
    T = Table("cross")
    ms = [float(r["M_star"]) for r in T.find(variant="standard") if np.isfinite(float(r["M_star"]))]
    if len(T.find(variant="standard")) < len(LADDER) or not ms:
        raise RuntimeError("stage 1 incomplete")
    return min(ms)


def grid():
    Ms = best_Mstar()
    return [Ms + d for d in DELTAS]


def stage_energies():
    E, O = Table("en"), Table("obs")
    G = grid()
    for a in LADDER:
        todo = [M for M in G if not E.find(alpha_s=a, variant="standard", M=M)]
        if not todo:
            continue
        t0 = time.perf_counter()
        mb = slim(basis(a))
        for M in todo:
            t1 = time.perf_counter()
            ob = mb.observables(M, state=K)
            Eth = family_threshold(M)
            E.add(state=K, M=M, alpha_s=a, variant="standard", E=ob["E"], E_minus_Eth=ob["E"] - Eth,
                  seconds=round(time.perf_counter() - t1, 1))
            O.add(state=K, M=M, alpha_s=a, E=ob["E"], E_th=Eth, binding=Eth - ob["E"], r01=ob["r01"],
                  r01_2=ob["r01_2"], rXe=ob["rXe"], dE_dM=ob["dE_dM"])
        log(f"stage2 alpha_s={a:g}: {len(todo)} masses ({time.perf_counter() - t0:.0f}s)")
        del mb
        gc.collect()


def check_masses():
    Ms = best_Mstar()
    return [Ms + 0.01, Ms + 0.1, Ms + 0.5]


def stage_checks():
    C, X, E = Table("conv"), Table("cross"), Table("en")
    cm = check_masses()
    # reference energies of the standard ladder at the check masses (+ tol 1e-10 variant)
    for a in CHECK_ALPHAS:
        need = [M for M in cm if not C.find(variant="reference (tol 1e-12)", alpha_s=a, M=M)]
        if not need:
            continue
        mb = slim(basis(a), tols=TOLS)
        for M in need:
            e = gap(mb, M)
            C.add(state=K, M=M, variant="reference (tol 1e-12)", alpha_s=a, size=mb.size, kept=mb._projs[1e-12]["kept"],
                  E_minus_Eth=e, diff_from_ref=0.0)
            et = gap(mb, M, 1e-10)
            C.add(state=K, M=M, variant="tol=1e-10", alpha_s=a, size=mb.size, kept=mb._projs[1e-10]["kept"],
                  E_minus_Eth=et, diff_from_ref=et - e)
        log(f"stage3 reference alpha_s={a:g} done")
        del mb
        gc.collect()
    ref = {(float(r["alpha_s"]), float(r["M"])): float(r["E_minus_Eth"]) for r in C.find(variant="reference (tol 1e-12)")}
    for label, kw in CHECKS:
        for a in CHECK_ALPHAS:
            need = [M for M in cm if not C.find(variant=label, alpha_s=a, M=M)]
            want_cross = (a == 1e-4) and not X.find(variant=label, alpha_s=a)
            if not need and not want_cross:
                continue
            t0 = time.perf_counter()
            mb = basis(a, **{**STD, **kw})
            pr = mb.project(1e-12)
            slim(mb)
            for M in need:
                e = gap(mb, M)
                r0 = min((v for (aa, mm), v in ref.items() if aa == a and abs(mm - M) < 1e-9), default=float("nan"))
                C.add(state=K, M=M, variant=label, alpha_s=a, size=mb.size, kept=pr["kept"], E_minus_Eth=e,
                      diff_from_ref=e - r0)
            if want_cross:
                try:
                    Ms = crossing(mb)
                except ValueError:
                    Ms = float("nan")
                X.add(state=K, alpha_s=a, variant=label, size=mb.size, kept=pr["kept"], dropped=pr["dropped"],
                      M_star=Ms, build_s=round(mb.build_time, 1), total_s=round(time.perf_counter() - t0, 1))
            log(f"stage3 '{label}' alpha_s={a:g} size={mb.size} ({time.perf_counter() - t0:.0f}s)")
            del mb
            gc.collect()


def stage_bounds():
    U = Table("ub")
    for label, kw in [("standard", {}), ("diffuse nd=70, nd2=45", dict(nd=70, nd2=45))]:
        a = 1e-4
        need = [(t, m) for t in TOLS for m in MARGINS if not U.find(variant=label, tol=t, margin=m)]
        if not need:
            continue
        t0 = time.perf_counter()
        mb = slim(basis(a, **{**STD, **kw}), tols=TOLS)
        for tol, margin in need:
            try:
                Mb = crossing(mb, tol=tol, margin=margin)
            except ValueError as exc:
                log(f"stage4 {label} tol={tol:g} margin={margin:g}: {exc}")
                Mb = float("nan")
            slope = float("nan")
            if margin == 0.0 and np.isfinite(Mb):
                mb._proj = mb._projs[tol]
                ob = mb.observables(Mb, state=K, tol=tol)
                slope = ob["dE_dM"] + 1 / (2 * (Mb + 1) ** 2)
            U.add(state=K, variant=label, alpha_s=a, tol=tol, margin=margin, M_bound=Mb, slope_at_Mstar=slope)
            log(f"stage4 {label} tol={tol:g} margin={margin:g}: M_bound={Mb:.9f}")
        log(f"stage4 '{label}' done ({time.perf_counter() - t0:.0f}s)")
        del mb
        gc.collect()


# ---------------------------------------------------------------- dry run / tests
def calibrate():
    """One basis (alpha_s = 1e-4, standard): build, project, one energy. Returns timings."""
    t0 = time.perf_counter()
    mb = basis(1e-4)
    tb = time.perf_counter() - t0
    pr = mb.project(1e-12)
    tp = time.perf_counter() - t0 - tb
    t1 = time.perf_counter()
    e = mb.energies(99.0, K + 1)[K] - family_threshold(99.0)
    te = time.perf_counter() - t1
    return dict(size=mb.size, kept=pr["kept"], build=tb, project=tp, eig=te, gap_99=e)


def sizes(**kw):
    """Basis size without building matrices."""
    from threebody.threshold import Block
    mu = CONFIG["M_ref"] / (CONFIG["M_ref"] + 1)
    p = {**STD, **kw}
    n = 0
    for i, (ax, ae) in enumerate(CONFIG["compact"]):
        N = p["Nc_big"] if (p["Nc_big"] and i >= len(CONFIG["compact"]) - 2) else p["Nc"]
        n += compact_block(N, ax, ae).size
    n += diffuse_block(mu, 1e-4, p["nc"], p["nd"]).size + diffuse_block(mu, 0.03, p["nc"], p["nd2"]).size
    del Block
    return n


def dry():
    log(f"=== run_exact.py --dry  threads={os.environ.get('THREEBODY_BLAS_THREADS')}")
    log(f"config {json.dumps(CONFIG)}  ladder {LADDER}  deltas {DELTAS}")
    c = calibrate()
    log(f"calibration (alpha_s=1e-4, n={c['size']}, kept={c['kept']}): build {c['build']:.0f}s, project {c['project']:.0f}s, "
        f"one eigensolve {c['eig']:.1f}s;  E_4 - E_th at M=99: {c['gap_99']:.3e}")
    n0 = c["size"]
    b0 = c["build"] + c["project"]
    ev = c["eig"]
    scale = lambda n: (n / n0) ** 3
    est = {}
    # stage 1: per basis build+project + ~22 eigensolves (brent incl. bracket)
    est[1] = len(LADDER) * (b0 + 22 * ev)
    est[2] = len(LADDER) * (b0 + len(DELTAS) * 1.3 * ev)
    s3 = 2 * (b0 * scale(n0) * 1.0) + 2 * 3 * 2 * ev
    for label, kw in CHECKS:
        n = sizes(**kw)
        s3 += 2 * (b0 * scale(n) + 3 * ev * scale(n)) + 22 * ev * scale(n)
        log(f"  check '{label}': size {n} (x{scale(n):.1f} cost), peak memory ~{16 * n * n * 8 / 1e9:.1f} GB")
    est[3] = s3
    n_nd = sizes(nd=70, nd2=45)
    est[4] = (b0 + 6 * 22 * ev) + (b0 * scale(n_nd) + 6 * 22 * ev * scale(n_nd))
    T = Table("cross")
    log(f"done so far: crossings {len(T.rows)}, energies {len(Table('en').rows)}, checks {len(Table('conv').rows)}, "
        f"bounds {len(Table('ub').rows)} rows (done rows are skipped)")
    for s in (1, 2, 3, 4):
        log(f"stage {s}: estimated {est[s] / 60:.0f} min (before skipping done rows)")
    log(f"TOTAL estimated {sum(est.values()) / 3600:.1f} h; peak memory standard basis ~{16 * n0 * n0 * 8 / 1e9:.1f} GB")


def test(variant, M):
    variants = {"compact3": dict(compact=[CONFIG["compact"][i] for i in (0, 1, 3)]),   # k=3-like: 3 compact sets
                "compact3_noroot": dict(compact=CONFIG["compact"][:3])}               # without the root-4 set
    kw = variants.get(variant) or dict(CHECKS).get(variant, {})
    alpha = 1e-4
    t0 = time.perf_counter()
    mb = basis(alpha, **{**STD, **kw})
    pr = mb.project(1e-12)
    e = mb.energies(M, K + 1)
    log(f"test {variant} alpha_s={alpha:g} M={M}: size={mb.size} kept={pr['kept']} E_4-E_th={e[K] - family_threshold(M):.12e} "
        f"E_3-E_th={e[K - 1] - family_threshold(M):.6e} ({time.perf_counter() - t0:.0f}s)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--stage", type=int, nargs="*", default=[1, 2, 3, 4])
    ap.add_argument("--test", nargs=2, metavar=("VARIANT", "M"))
    a = ap.parse_args()
    if a.dry:
        return dry()
    if a.test:
        return test(a.test[0], float(a.test[1]))
    log(f"=== run_exact.py stages={a.stage} threads={os.environ.get('THREEBODY_BLAS_THREADS')}")
    t0 = time.perf_counter()
    funcs = {1: stage_crossings, 2: stage_energies, 3: stage_checks, 4: stage_bounds}
    for s in a.stage:
        t1 = time.perf_counter()
        funcs[s]()
        log(f"stage {s} finished ({time.perf_counter() - t1:.0f}s)")
    with open(os.path.join(HERE, "basis_state4.json"), "w") as fh:
        json.dump(dict(config=CONFIG, ladder=LADDER, deltas=DELTAS, std=STD, checks=CHECKS,
                       blocks=[dict(name=b.name, alphas=list(b.alphas), size=b.size, nmax=list(b.nmax))
                               for b in basis_blocks_only(1e-4)]), fh, indent=1)
    log(f"ALL STAGES DONE ({time.perf_counter() - t0:.0f}s)")


def basis_blocks_only(alpha_s):
    mu = CONFIG["M_ref"] / (CONFIG["M_ref"] + 1)
    bl = [compact_block(16, ax, ae) for ax, ae in CONFIG["compact"]]
    return bl + [diffuse_block(mu, alpha_s, STD["nc"], STD["nd"]), diffuse_block(mu, 0.03, STD["nc"], STD["nd2"])]


if __name__ == "__main__":
    main()
