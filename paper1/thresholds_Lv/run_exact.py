"""Exact (non-BO) critical masses of natural-parity L > 0 levels of {X+, X+, e-}: thresholds_Lv.

Written 2026-11-15.  Same method as paper1/rotational (examples/rotational_critical.py: fixed basis =>
mass-independent matrices; Brent root of E_k(M) - E_th(M); every fixed-basis crossing M* is a rigorous
upper bound to M_c up to rounding), extended to excited roots (k = 1: J = 1, v = 1), to L = 3, 4, with
two compact scale sets, nested diffuse sets, compact-size checks, rounding-safe bounds and a
binding-energy table for the L > 0 threshold law (binding linear in M - M_c).

Targets (key: L, exchange symmetry, root k):
    J1v1 : L = 1, 'A', k = 1   (1s sigma_g J=1 v=1; the (1,1) level of muonic molecules)
    L3A  : L = 3, 'A', k = 0   (1s sigma_g J=3 v=0)
    L4S  : L = 4, 'S', k = 0   (1s sigma_g J=4 v=0)
    L1S  : L = 1, 'S', k = 0   (2p sigma_u J=1; BO M_c = 815.56): fixed compact scales, the level is
           carried by 5 diffuse sets (atom x X-atom relative motion); feasibility test in PROGRESS.md.

Run from the project root (BLAS threads: THREEBODY_BLAS_THREADS, default 8):
    python3 paper1/thresholds_Lv/run_exact.py --target J1v1            # all stages
    ... --target J1v1 --dry       # sizes, memory, time estimate from a one-point calibration
    ... --target J1v1 --stage 2 3 # some stages only
Stages (results appended to <target>_crossings.csv / <target>_binding.csv at once; done rows skipped):
  1 scales     optimize two compact scale sets (root 0 and root k, N = 12) -> <target>_scales.json
  2 compact    nested compact sequence N = Ns (both scale sets)            -> crossings
  3 diffuse    compact N_D + nested diffuse sets                           -> crossings
  4 checks     compact size N_D - 2, N_D + 2 with all diffuse sets; tol 1e-10; margin 1e-11 (rounding-safe)
  5 binding    E_th - E_k at M = M*_best * (1 + x) with the largest basis  -> binding
Log: <target>_log.txt.  Then analyze.py.
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

from threebody import symmetric_family, family_threshold  # noqa: E402
from threebody.coupled import CoupledBasis, IdxBlock, block  # noqa: E402
from threebody.rotational import optimize_scales  # noqa: E402
from threebody.threshold import diffuse_index  # noqa: E402

BO = json.load(open(os.path.join(HERE, "bo_brackets.json"))) if os.path.exists(os.path.join(HERE, "bo_brackets.json")) else {}

# M_BO / M_BOad: pure BO and BO+ad critical masses (bo_brackets.json); the exact M_c lies below M_BOad
# for every level computed so far (offsets 0.29-0.54), so brackets are [0.85 M_BOad - 1, M_BOad + 0.3].
DIFFUSE = ((4, 12, 0.4), (4, 16, 0.15), (5, 16, 0.05), (5, 20, 0.02), (6, 20, 0.008))
TARGETS = {
    "J1v1": dict(L=1, sym="A", k=1, bo=("g L=1", 1), Ns=(12, 14, 16, 18, 20), ND=18, diffuse=DIFFUSE),
    # L = 3 (2 channels) and L = 4 (3 channels) multiply the size per block: compact N kept smaller so that
    # the largest check basis stays below ~8000 functions (~6 GB); see README for the sizes.
    "L3A": dict(L=3, sym="A", k=0, bo=("g L=3", 0), Ns=(8, 10, 12, 14), ND=12, diffuse=DIFFUSE[:4]),
    "L4S": dict(L=4, sym="S", k=0, bo=("g L=4", 0), Ns=(8, 10, 12), ND=10, diffuse=DIFFUSE[:3]),
    # 2p sigma_u J=1: a shallow (6.1e-5 Eh unit-mass) long-range well at R ~ 12.6 bohr; BO M_c = 815.56.
    # Fixed scales (no root-0 optimum exists in a compact basis): one short-range compact set; the level
    # itself is carried by the diffuse (atom x X-atom relative motion) sets, which are larger here.
    "L1S": dict(L=1, sym="S", k=0, bo=("u L=1", 0), Ns=(8, 10), ND=10,
                fixed_scales=[[1.0, 1.0, 1.0], [0.5, 0.5, 1.0]],
                diffuse=((4, 40, 0.5), (4, 40, 0.25), (4, 50, 0.12), (4, 50, 0.06), (4, 60, 0.03))),
}
XB = (0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2)   # binding table: M = M* (1 + x)


def bo_values(t):
    key, v = TARGETS[t]["bo"]
    if f"{key} BO+ad" not in BO:          # 2p sigma_u: BO only
        return BO[f"{key} BO"][v], None
    return BO[f"{key} BO"][v], BO[f"{key} BO+ad"][v]


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


CROSS_COLS = ["target", "label", "stage", "size", "kept", "tol", "margin", "M_star", "seconds", "finished"]
BIND_COLS = ["target", "label", "M", "x", "E", "E_th", "binding", "E_next", "seconds", "finished"]

_logs = {}


def log(t, msg):
    if t not in _logs:
        _logs[t] = open(os.path.join(HERE, f"{t}_log.txt"), "a")
    line = time.strftime("%Y-%m-%d %H:%M:%S ") + msg
    print(line, flush=True)
    _logs[t].write(line + "\n")
    _logs[t].flush()


# ---------------------------------------------------------------- basis
def diffuse_pair(alpha_atom, alpha_d, n_c, n_d):
    idx = diffuse_index(n_c, n_d)[:, [1, 0, 2]]  # diffuse in s0: X0 far from the (X1 e) atom
    B = IdxBlock((alpha_d, alpha_atom, alpha_atom), idx, name=f"D{n_c}x{n_d}@{alpha_d:g}")
    return [B, B.mirror()]


def make_basis(t, scales, N, n_diffuse, mu):
    cfg = TARGETS[t]
    blocks = [block(N, al, name=f"N{N}s{i}") for i, al in enumerate(scales)]
    for (n_c, n_d, ad) in cfg["diffuse"][:n_diffuse]:
        blocks += diffuse_pair(mu, ad, n_c, n_d)
    return CoupledBasis(cfg["L"], blocks, cfg["sym"])


def basis_size(t, scales, N, n_diffuse):
    cfg = TARGETS[t]
    from threebody import angular
    nch = len([c for c in angular.components(cfg["L"]) if c <= angular.swap_component(c)])
    n = 0
    for c in [c for c in angular.components(cfg["L"]) if c <= angular.swap_component(c)]:
        for al in scales:
            b = block(N, al)
            if angular.swap_component(c) == c:
                a = b.index
                n += int((a[:, 0] >= a[:, 1]).sum() if cfg["sym"] == "S" else (a[:, 0] > a[:, 1]).sum())
            else:
                n += b.size
        for (n_c, n_d, ad) in cfg["diffuse"][:n_diffuse]:
            n += 2 * len(diffuse_index(n_c, n_d))
    return n, nch


def gap(t, b, M, tol=1e-12):
    k = TARGETS[t]["k"]
    e = b.solve(symmetric_family(M), k + 1, tol, vectors=False)[0]
    return e[k] - family_threshold(M)


def crossing(t, b, lo, hi, tol=1e-12, margin=0.0):
    f = lambda M: gap(t, b, M, tol) + margin
    flo, fhi = f(lo), f(hi)
    tries = 0
    while fhi > 0 and tries < 4:      # widen upwards if the basis is still too small
        lo, flo = hi, fhi
        hi *= 1.15
        fhi = f(hi)
        tries += 1
    if not (flo > 0 > fhi):
        raise ValueError(f"no bracket: f({lo:.6f})={flo:.3e} f({hi:.6f})={fhi:.3e}")
    return brentq(f, lo, hi, xtol=1e-11, rtol=1e-13)


def bracket(t):
    Mbo, Mad = bo_values(t)
    if Mad is None:
        return 0.9 * Mbo, 1.05 * Mbo
    return 0.85 * Mad - 1.0, Mad + 0.3


# ---------------------------------------------------------------- stages
def scales_path(t):
    return os.path.join(HERE, f"{t}_scales.json")


def stage_scales(t):
    p = scales_path(t)
    if os.path.exists(p):
        return json.load(open(p))
    cfg = TARGETS[t]
    Mbo, Mad = bo_values(t)
    if "fixed_scales" in cfg:
        out = dict(M_opt=Mbo, scales=cfg["fixed_scales"])
        json.dump(out, open(p, "w"), indent=1)
        log(t, f"fixed scales {cfg['fixed_scales']} (M_ref {Mbo:.4f} for the atom scale mu)")
        return out
    M_opt = Mad          # above the exact M_c: the level is bound there
    fam = symmetric_family(M_opt)
    t0 = time.perf_counter()
    a0, e0, n0 = optimize_scales(fam, cfg["L"], 12, cfg["sym"], state=0)
    log(t, f"scales root 0 at M={M_opt:.6f}: {a0} E={e0:.12f} ({n0} evals, {time.perf_counter() - t0:.0f}s)")
    out = dict(M_opt=M_opt, scales=[list(a0)])
    if cfg["k"] > 0:
        t0 = time.perf_counter()
        ak, ek, nk = optimize_scales(fam, cfg["L"], 12, cfg["sym"], state=cfg["k"], x0=a0)
        log(t, f"scales root {cfg['k']} at M={M_opt:.6f}: {ak} E={ek:.12f} E-E_th={ek - family_threshold(M_opt):.3e} "
               f"({nk} evals, {time.perf_counter() - t0:.0f}s)")
        out["scales"].append(list(ak))
    else:
        # second, more diffuse set for the X-X stretch (geometric step in alpha_X), as for v = 4
        out["scales"].append([0.6 * a0[0], 0.6 * a0[1], a0[2]])
        log(t, f"second compact set (0.6 alpha_X): {out['scales'][1]}")
    json.dump(out, open(p, "w"), indent=1)
    return out


def run_crossing(t, T, label, stage, N, n_diffuse, sc, lo, hi, tols=(1e-12,), margins=(0.0,)):
    todo = [(tol, m) for tol in tols for m in margins if not T.find(target=t, label=label, tol=tol, margin=m)]
    if not todo:
        return [float(r["M_star"]) for r in T.find(target=t, label=label, tol=1e-12, margin=0.0)]
    t0 = time.perf_counter()
    mu = sc["M_opt"] / (sc["M_opt"] + 1)
    b = make_basis(t, sc["scales"], N, n_diffuse, mu)
    out = []
    for tol, m in todo:
        t1 = time.perf_counter()
        pr = b.project(tol)
        try:
            Ms = crossing(t, b, lo, hi, tol, m)
        except ValueError as exc:
            log(t, f"  {label} tol={tol:g} margin={m:g}: {exc}")
            Ms = float("nan")
        T.add(target=t, label=label, stage=stage, size=b.size, kept=pr["kept"], tol=tol, margin=m, M_star=Ms,
              seconds=round(time.perf_counter() - t1 + (b.build_time if (tol, m) == todo[0] else 0), 1))
        log(t, f"stage{stage} {label}: size={b.size} kept={pr['kept']} tol={tol:g} margin={m:g} M*={Ms:.10f} "
               f"({time.perf_counter() - t1:.0f}s; build {b.build_time:.0f}s)")
        out.append(Ms)
    del b
    gc.collect()
    return out


def stage_compact(t):
    cfg, sc = TARGETS[t], stage_scales(t)
    T = Table(os.path.join(HERE, f"{t}_crossings.csv"), CROSS_COLS)
    lo, hi = bracket(t)
    for N in cfg["Ns"]:
        Ms = run_crossing(t, T, f"compact N={N}", 2, N, 0, sc, lo, hi)
        if Ms and np.isfinite(Ms[0]):
            hi = Ms[0] * 1.0000001   # nested sequence: next crossing lies below


def stage_diffuse(t):
    cfg, sc = TARGETS[t], stage_scales(t)
    T = Table(os.path.join(HERE, f"{t}_crossings.csv"), CROSS_COLS)
    lo, hi = bracket(t)
    r = T.find(target=t, label=f"compact N={cfg['ND']}", tol=1e-12, margin=0.0)
    if r and np.isfinite(float(r[0]["M_star"])):
        hi = float(r[0]["M_star"]) * 1.0000001
    for nd in range(1, len(cfg["diffuse"]) + 1):
        Ms = run_crossing(t, T, f"compact N={cfg['ND']} + {nd} diffuse", 3, cfg["ND"], nd, sc, lo, hi)
        if Ms and np.isfinite(Ms[0]):
            hi = Ms[0] * 1.0000001


def stage_checks(t):
    cfg, sc = TARGETS[t], stage_scales(t)
    T = Table(os.path.join(HERE, f"{t}_crossings.csv"), CROSS_COLS)
    lo, hi = bracket(t)
    nd = len(cfg["diffuse"])
    # rounding-safe bound and tol check on the standard largest basis
    run_crossing(t, T, f"compact N={cfg['ND']} + {nd} diffuse", 4, cfg["ND"], nd, sc, lo, hi,
                 tols=(1e-12, 1e-10), margins=(0.0, 1e-11))
    for N in (cfg["ND"] - 2, cfg["ND"] + 2):
        run_crossing(t, T, f"compact N={N} + {nd} diffuse", 4, N, nd, sc, lo, hi)


def best_mstar(t):
    T = Table(os.path.join(HERE, f"{t}_crossings.csv"), CROSS_COLS)
    v = [float(r["M_star"]) for r in T.rows if r["tol"] == repr(1e-12) and r["margin"] == repr(0.0)]
    v = [x for x in v if np.isfinite(x)]
    if not v:
        raise RuntimeError("no crossing yet")
    return min(v)


def stage_binding(t):
    cfg, sc = TARGETS[t], stage_scales(t)
    B = Table(os.path.join(HERE, f"{t}_binding.csv"), BIND_COLS)
    nd = len(cfg["diffuse"])
    label = f"compact N={cfg['ND']} + {nd} diffuse"
    Ms = best_mstar(t)
    todo = [x for x in XB if not B.find(target=t, label=label, x=x)]
    if not todo:
        return
    mu = sc["M_opt"] / (sc["M_opt"] + 1)
    b = make_basis(t, sc["scales"], cfg["ND"], nd, mu)
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


# ---------------------------------------------------------------- dry run
def dry(t):
    cfg = TARGETS[t]
    Mbo, Mad = bo_values(t)
    log(t, f"=== --dry {t}: L={cfg['L']} sym={cfg['sym']} root k={cfg['k']}; BO {Mbo:.6f}, BO+ad {Mad}; "
           f"bracket {bracket(t)}; threads {os.environ.get('THREEBODY_BLAS_THREADS')}")
    sc = json.load(open(scales_path(t)))["scales"] if os.path.exists(scales_path(t)) else [[1.0, 1.0, 1.0]] * 2
    plan = [(N, 0) for N in cfg["Ns"]] + [(cfg["ND"], i) for i in range(1, len(cfg["diffuse"]) + 1)]
    nd = len(cfg["diffuse"])
    plan += [(cfg["ND"] - 2, nd), (cfg["ND"] + 2, nd)]
    # calibration: build + one eigensolve at a mid-size basis
    Nc = cfg["ND"]          # calibrate on the compact basis used under the diffuse sets
    t0 = time.perf_counter()
    b = make_basis(t, sc, Nc, 0, 0.95)
    tb = time.perf_counter() - t0
    t1 = time.perf_counter()
    b.project(1e-12)
    gap(t, b, Mad or Mbo)
    te = time.perf_counter() - t1
    n0 = b.size
    log(t, f"calibration: compact N={Nc}, size {n0}: build {tb:.1f}s, project+1 solve {te:.1f}s")
    total, peak = 0.0, 0
    for N, k in plan:
        n, nch = basis_size(t, sc, N, k)
        # build ~ n^2 (term evaluation), linear algebra ~ n^3; ~18 solves per crossing (measured: 17 for J1v1)
        tt = tb * (n / n0) ** 2 + 18 * te * (n / n0) ** 3
        reps = 4 if (N, k) == (cfg["ND"], nd) else 1
        total += tt * reps
        peak = max(peak, n)
        log(t, f"  N={N:2d} + {k} diffuse: size {n:6d} ({nch} channel(s)) est {tt / 60:6.1f} min x{reps}, "
               f"memory ~{12 * n * n * 8 / 1e9:.1f} GB")
    n, _ = basis_size(t, sc, cfg["ND"], nd)
    total += len(XB) * te * (n / n0) ** 3 + tb * (n / n0) ** 2
    log(t, f"TOTAL estimated {total / 3600:.2f} h (excl. stage 1 scales, ~1-5 min); peak size {peak}, "
           f"peak memory ~{12 * peak * peak * 8 / 1e9:.1f} GB (12 n^2 doubles)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True, choices=sorted(TARGETS))
    ap.add_argument("--stage", type=int, nargs="*", default=[1, 2, 3, 4, 5])
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--test", nargs=3, metavar=("N", "NDIFF", "M"), help="one E_k - E_th (short test)")
    a = ap.parse_args()
    t = a.target
    if a.dry:
        return dry(t)
    if a.test:
        sc = stage_scales(t)
        N, nd, M = int(a.test[0]), int(a.test[1]), float(a.test[2])
        t0 = time.perf_counter()
        b = make_basis(t, sc["scales"], N, nd, sc["M_opt"] / (sc["M_opt"] + 1))
        pr = b.project(1e-12)
        k = TARGETS[t]["k"]
        e = b.solve(symmetric_family(M), k + 2, 1e-12, vectors=False)[0]
        log(t, f"test N={N} +{nd} diffuse M={M}: size={b.size} kept={pr['kept']} E_k-E_th={e[k] - family_threshold(M):.6e} "
               f"E_k+1-E_th={e[k + 1] - family_threshold(M):.3e} (build {b.build_time:.0f}s, total {time.perf_counter() - t0:.0f}s)")
        return
    log(t, f"=== run_exact.py --target {t} stages={a.stage} threads={os.environ.get('THREEBODY_BLAS_THREADS')}")
    funcs = {1: stage_scales, 2: stage_compact, 3: stage_diffuse, 4: stage_checks, 5: stage_binding}
    t0 = time.perf_counter()
    for s in a.stage:
        t1 = time.perf_counter()
        funcs[s](t)
        log(t, f"stage {s} finished ({time.perf_counter() - t1:.0f}s)")
    log(t, f"ALL REQUESTED STAGES DONE ({time.perf_counter() - t0:.0f}s)")


if __name__ == "__main__":
    main()
