"""Exact (non-BO) critical masses of {A+ (M), B+ (rM), e-} levels: unequal_masses.

Written 2026-11-16.  Same methodology as paper1/critical_masses_v4 (L = 0: ladder of ultra-diffuse open-channel
scales) and paper1/thresholds_Lv (L > 0: nested compact and diffuse sequences, checks, rounding-safe bounds,
binding table), adapted to unequal masses:
  * no exchange symmetry (CoupledBasis(L, blocks, None)); particles 0 = A (mass M), 1 = B (rM), 2 = e;
  * exact threshold E_th = -mu_Be/2 = -rM/(2(rM+1)) (A+ + (B e)(1s), the lowest continuum for r > 1);
  * open-channel diffuse blocks (A far from (B e): high degree in s0, scales (alpha_d, mu_B, mu_B)) and
    closed-channel blocks (B far from (A e): high degree in s1, scales (mu_A, alpha_c, mu_A));
  * memory-lean assembly: only S, T_N = T0 + T1/r and H0 = T2 + V01 - V02 - V12 are kept
    (H(M) = T_N/M + H0), projected once (canonical orthogonalization), raw matrices dropped.
Every fixed-basis crossing M* (root k of E(M) = E_th(M)) is a rigorous upper bound to M_c (up to rounding).
The frozen prediction (PREDICTION_UNEQUAL.md) is NOT read here; brackets come from bo_brackets.json
(BO+ad own-asymptote and exact-threshold model values) with automatic widening.

Run from the project root (BLAS threads: THREEBODY_BLAS_THREADS, default 8):
    python3 paper1/unequal_masses/run_exact.py --target J0v1          # all stages
    ... --target J0v1 --dry          # sizes, memory, time estimate (one-point calibration, ~1-3 min)
    ... --target J0v1 --stage 2 3    # some stages only
    ... --target J0v1 --test N M     # one E_k - E_th at mass M with compact N and the standard diffuse sets
    ... --list                       # targets
Stages (results appended to <target>_crossings.csv / <target>_binding.csv at once; done rows are skipped):
  1 scales   optimize compact scale sets (root 0 and root k, N = 12, symmetric-pair form (aX, aX, ae))
  2 compact  nested compact sequence (no diffuse)                          -> crossings
  3 diffuse  L = 0: ladder of open-channel alpha_s (all other diffuse sets present); L > 0: nested diffuse
  4 checks   largest basis: tol 1e-10, margin 1e-11 (rounding-safe); compact N_D -/+ 2; enlarged diffuse
  5 binding  E_th - E_k at M = M*_best (1 + x) with the standard largest basis
Log: <target>_log.txt.  Then analyze.py -> RESULT.md.
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
os.environ.setdefault("OPENBLAS_NUM_THREADS", os.environ["THREEBODY_BLAS_THREADS"])
sys.path.insert(0, os.path.join(ROOT, "solvers", "threebody"))
import numpy as np  # noqa: E402
from scipy.linalg import eigh  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

from threebody.coupled import CoupledBasis, IdxBlock, block  # noqa: E402
from threebody.rotational import optimize_scales  # noqa: E402
from threebody.system import ThreeBodySystem  # noqa: E402
from threebody.threshold import diffuse_index  # noqa: E402

BRK = os.path.join(HERE, "bo_brackets.json")
LADDER = (1e-2, 3e-3, 1e-3, 3e-4, 1e-4)                      # L = 0 open-channel alpha_s
OPEN_FIXED = ((4, 30, 0.03), (4, 16, 0.15))                  # L = 0: (n_c, n_d, alpha) always present
CLOSED = ((4, 16, 0.6), (4, 20, 0.25))                       # closed channel (kappa ~ sqrt(2 mu_rel dE) ~ 0.58)
DIFF_L = ((4, 12, 0.4), (4, 16, 0.15), (5, 16, 0.05), (5, 20, 0.02), (6, 20, 0.008))   # L > 0, nested open
TARGETS = {
    # primary, r = 2
    "J0v1": dict(r=2.0, L=0, k=1, Ns=(12, 14, 16), ND=16, nd=50),
    "J0v2": dict(r=2.0, L=0, k=2, Ns=(12, 14, 16), ND=16, nd=50),
    "J0v3": dict(r=2.0, L=0, k=3, Ns=(12, 14, 16), ND=16, nd=50),
    "J0v4": dict(r=2.0, L=0, k=4, Ns=(12, 14, 16), ND=16, nd=50),
    # L > 0: 2 (L=1) or 3 (L=2) coupled components and no exchange symmetry -> sizes kept below ~8000
    # (peak ~3 GB); no zero-energy halo for L > 0, so fewer/less diffuse sets suffice (thresholds_Lv).
    "J1v0": dict(r=2.0, L=1, k=0, Ns=(8, 10, 12), ND=12, diffuse=DIFF_L[:4]),
    "J2v0": dict(r=2.0, L=2, k=0, Ns=(6, 8, 10), ND=10, diffuse=DIFF_L[:3]),
    # checks / secondary
    "R1J0v1": dict(r=1.0, L=0, k=1, Ns=(12, 14, 16), ND=16, nd=50),     # r = 1: symmetric exact 9.525634
    "S15J0v1": dict(r=1.5, L=0, k=1, Ns=(12, 14, 16), ND=16, nd=50),    # secondary ratio
}
XB = (0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2)
CROSS_COLS = ["target", "label", "stage", "size", "kept", "tol", "margin", "M_star", "seconds", "finished",
              "spurious", "reason"]
BIND_COLS = ["target", "label", "M", "x", "E", "E_th", "binding", "E_next", "seconds", "finished", "spurious", "reason"]
# Spurious-root guard (added 2026-11-16 after a fake low root at J0v1, N=16 + alpha_s=3e-4, tol 1e-12):
#  (a) hard: E_0 - E_th < GUARD_HARD at the solve mass (as in threshold_map/run_map.py);
#  (b) k >= 1 only: E_0 of the extended basis below E_0 of its compact-only part (same N, same M) by more than
#      GUARD_E0 (the ground state is compact; diffuse functions cannot lower it materially).  Not applied for
#      k = 0, where root 0 is the threshold state itself and diffuse functions lower it legitimately.
#  (c) crossings only: a genuine crossing has gap(M*(1 - d)) >= GUARD_BELOW (root k still at/above threshold just
#      below M*) and gap(M*(1 + d)) - gap(M*(1 - d)) < -GUARD_SLOPE (a real slope), d = GUARD_D.  The J0v1 failure
#      was of this kind, not (a) or (b): at tol 1e-12 the alpha_s=3e-4 basis has a near-threshold level whose
#      E - E_th is flat in M and not reproducible between builds (+6e-12 in one build, -1e-8 in another at the same
#      M; BLAS rounding in an ill-conditioned projection), so brentq found a "crossing" in the noise.  Measured at
#      d = 1e-3: genuine crossings (alpha_s = 1e-3, 1e-4) gap(+d) - gap(-d) ~ -4e-8 with gap(-d) >= 0; the
#      artefact had gap(-d) = -1e-8 .. -1e-7 and a difference of -4e-10 .. -1e-9 (check_guard.py).
# Spurious rows are written with M_star = nan, spurious = 1 and the reason; they are never used downstream.
GUARD_HARD = -0.15
GUARD_E0 = 1e-6
GUARD_D = 1e-3
GUARD_BELOW = -1e-10
GUARD_SLOPE = 1e-9


# ---------------------------------------------------------------- physics
def threshold(M, r):
    """Exact lowest threshold: -mu/2 of the heavier atom (r >= 1: B e)."""
    m = max(M, r * M)
    return -0.5 * m / (m + 1.0)


def system(M, r):
    return ThreeBodySystem((M, r * M, 1.0), (1.0, 1.0, -1.0), name=f"AB_e M={M:g} r={r:g}")


class LeanBasis:
    """Unsymmetrized natural-parity basis with only S, T_N = T0 + T1/r, H0 = T2 + V01 - V02 - V12 kept.

    Energies at mass M: eigenvalues of X^T (T_N/M + H0) X (canonical orthogonalization, tol)."""

    def __init__(self, L, blocks, r, tols=(1e-12,)):
        t0 = time.perf_counter()
        cb = CoupledBasis(L, blocks, None, ops=("S",))
        self.size, self.L, self.r = cb.size, L, r
        TN = cb._assemble("T0")
        TN += cb._assemble("T1") / r
        H0 = cb._assemble("T2")
        H0 += cb._assemble("V01")
        H0 -= cb._assemble("V02")
        H0 -= cb._assemble("V12")
        self.build_time = time.perf_counter() - t0
        S = cb.mats.pop("S")
        del cb
        gc.collect()
        self.proj = {}
        for tol in tols:
            dd = 1.0 / np.sqrt(np.diag(S))
            lam, U = np.linalg.eigh(S * np.outer(dd, dd))
            keep = lam > tol * lam[-1]
            X = (dd[:, None] * U[:, keep]) / np.sqrt(lam[keep])
            del U
            self.proj[tol] = dict(kept=int(keep.sum()), dropped=int((~keep).sum()), lam_min=float(lam[0]),
                                  TN=X.T @ TN @ X, H0=X.T @ H0 @ X)
            del X
            gc.collect()
        del S, TN, H0
        gc.collect()
        self.total_time = time.perf_counter() - t0

    def energies(self, M, n, tol=1e-12):
        p = self.proj[tol]
        return eigh(p["TN"] / M + p["H0"], eigvals_only=True, subset_by_index=[0, n - 1])


# ---------------------------------------------------------------- basis definition
def open_block(alpha_d, mu_B, n_c, n_d):
    """A far from the (B e) atom: high degree in s0 (s0 ~ 2R), atom scales mu_B in s1, s2."""
    return IdxBlock((alpha_d, mu_B, mu_B), diffuse_index(n_c, n_d)[:, [1, 0, 2]], name=f"O{n_c}x{n_d}@{alpha_d:g}")


def closed_block(alpha_c, mu_A, n_c, n_d):
    """B far from the (A e) atom: high degree in s1, atom scales mu_A in s0, s2."""
    return IdxBlock((mu_A, alpha_c, mu_A), diffuse_index(n_c, n_d), name=f"C{n_c}x{n_d}@{alpha_c:g}")


def blocks_for(t, sc, N, n_diffuse=None, alpha_s=None, big=False):
    """Block list.  L = 0: compact sets + [open alpha_s (ladder) + OPEN_FIXED + CLOSED] if alpha_s is given.
    L > 0: compact sets + first n_diffuse of DIFF_L (open) + CLOSED if n_diffuse > 0.
    big=True: enlarged diffuse (n_c + 2, n_d + 20) for the convergence check."""
    cfg = TARGETS[t]
    r, Mref = cfg["r"], sc["M_opt"]
    muA, muB = Mref / (Mref + 1), r * Mref / (r * Mref + 1)
    dc, dd = ((2, 20) if cfg["L"] == 0 else (1, 6)) if big else (0, 0)
    bl = [block(N, tuple(a), name=f"N{N}s{i}") for i, a in enumerate(sc["scales"])]
    diff = []
    if cfg["L"] == 0 and alpha_s is not None:
        diff.append(open_block(alpha_s, muB, 4 + dc, cfg["nd"] + dd))
        diff += [open_block(a, muB, nc + dc, nd + dd) for nc, nd, a in OPEN_FIXED]
    elif cfg["L"] > 0 and n_diffuse:
        diff += [open_block(a, muB, nc + dc, nd + dd) for nc, nd, a in cfg["diffuse"][:n_diffuse]]
    if diff:
        diff += [closed_block(a, muA, nc + dc, nd + dd) for nc, nd, a in (CLOSED if cfg["L"] == 0 else CLOSED[1:])]
    return bl + diff


def nfun(t, sc, N, n_diffuse=None, alpha_s=None, big=False):
    from threebody import angular
    nch = len(angular.components(TARGETS[t]["L"]))
    return nch * sum(b.size for b in blocks_for(t, sc, N, n_diffuse, alpha_s, big)), nch


# ---------------------------------------------------------------- IO
class Table:
    def __init__(self, path, cols):
        self.path, self.cols = path, cols
        self.rows = list(csv.DictReader(open(path))) if os.path.exists(path) else []
        if self.rows and any(c not in self.rows[0] for c in cols):     # old file: rewrite with the new header
            for r in self.rows:
                for c in cols:
                    r.setdefault(c, "")                                 # "" = written before the guard existed
            with open(path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=cols)
                w.writeheader()
                w.writerows({c: r[c] for c in cols} for r in self.rows)

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


_logs = {}


def log(t, msg):
    if t not in _logs:
        _logs[t] = open(os.path.join(HERE, f"{t}_log.txt"), "a")
    line = time.strftime("%Y-%m-%d %H:%M:%S ") + msg
    print(line, flush=True)
    _logs[t].write(line + "\n")
    _logs[t].flush()


def brackets(t):
    """(M_lo, M_hi, M_opt) from the model baselines in bo_brackets.json (NOT the theory value)."""
    cfg = TARGETS[t]
    if cfg["r"] == 1.0:
        return 9.0, 10.2, 10.2                       # symmetric BO+ad v=1: 10.063; exact 9.5256 (paper1)
    q = json.load(open(BRK))["targets"][f"r{cfg['r']:g}_J{cfg['L']}v{cfg['k']}"]["M"]
    lo = q["BOad"] - 0.5
    hi = q["BOad_thr"] + 0.5
    return lo, hi, hi


def gap(t, b, M, tol=1e-12):
    cfg = TARGETS[t]
    k = cfg["k"]
    return b.energies(M, k + 1, tol)[k] - threshold(M, cfg["r"])


def crossing(t, b, lo, hi, tol=1e-12, margin=0.0):
    f = lambda M: gap(t, b, M, tol) + margin
    flo, fhi = f(lo), f(hi)
    tries = 0
    while fhi > 0 and tries < 6:     # widen upwards if the basis is still too small
        lo, flo = hi, fhi
        hi *= 1.12
        fhi = f(hi)
        tries += 1
    tries = 0
    while flo < 0 and tries < 6:     # widen downwards
        hi, fhi = lo, flo
        lo *= 0.9
        flo = f(lo)
        tries += 1
    if not (flo > 0 > fhi):
        raise ValueError(f"no bracket: f({lo:.6f})={flo:.3e} f({hi:.6f})={fhi:.3e}")
    return brentq(f, lo, hi, xtol=1e-11, rtol=1e-13)


# ---------------------------------------------------------------- stages
def scales_path(t):
    return os.path.join(HERE, f"{t}_scales.json")


def stage_scales(t):
    p = scales_path(t)
    if os.path.exists(p):
        return json.load(open(p))
    cfg = TARGETS[t]
    lo, hi, M_opt = brackets(t)
    s = system(M_opt, cfg["r"])
    t0 = time.perf_counter()
    a0, e0, n0 = optimize_scales(s, cfg["L"], 12, None, state=0)
    log(t, f"scales root 0 at M={M_opt:.6f}: {a0} E={e0:.12f} ({n0} evals, {time.perf_counter() - t0:.0f}s)")
    scales = [list(a0)]
    if cfg["k"] > 0:
        t0 = time.perf_counter()
        ak, ek, nk = optimize_scales(s, cfg["L"], 12, None, state=cfg["k"], x0=a0)
        log(t, f"scales root {cfg['k']} at M={M_opt:.6f}: {ak} E-E_th={ek - threshold(M_opt, cfg['r']):.3e} "
               f"({nk} evals, {time.perf_counter() - t0:.0f}s)")
        # intermediate sets, geometric in alpha_X (v = 4 used two)
        nmid = 1 if cfg["k"] < 4 else 2
        for j in range(1, nmid + 1):
            x = j / (nmid + 1)
            scales.append([a0[0] ** (1 - x) * ak[0] ** x] * 2 + [a0[2] ** (1 - x) * ak[2] ** x])
        scales.append(list(ak))
    else:
        scales.append([0.6 * a0[0], 0.6 * a0[1], a0[2]])
    out = dict(M_opt=M_opt, scales=scales)
    json.dump(out, open(p, "w"), indent=1)
    log(t, f"compact scale sets: {scales}")
    return out


_compact = {}


def compact_basis(t, N):
    """Compact-only part (stage-1 scales, N) of the basis, cached per (target, N)."""
    if (t, N) not in _compact:
        _compact.clear()
        gc.collect()
        _compact[(t, N)] = LeanBasis(TARGETS[t]["L"], blocks_for(t, stage_scales(t), N), TARGETS[t]["r"])
    return _compact[(t, N)]


def guard(t, b, M, tol, N=None, margin=None):
    """Spurious-root check of basis b at mass M; returns "" (ok) or the reason.  N: compact size of b if b
    contains diffuse functions (enables check (b) for k >= 1); margin: M is a crossing with this margin
    (enables check (c))."""
    if not np.isfinite(M):
        return ""
    cfg = TARGETS[t]
    e0 = float(b.energies(M, 1, tol)[0])
    d = e0 - threshold(M, cfg["r"])
    if d < GUARD_HARD:
        return f"hard: E_0-E_th={d:.4f} < {GUARD_HARD} at M={M:.8f}"
    if N is not None and cfg["k"] >= 1:
        ec = float(compact_basis(t, N).energies(M, 1)[0])
        if e0 < ec - GUARD_E0:
            return f"E0: extended E_0 below compact-only E_0 by {ec - e0:.3e} > {GUARD_E0:g} at M={M:.8f}"
    if margin is not None:
        g_lo = gap(t, b, M * (1 - GUARD_D), tol) + margin
        g_hi = gap(t, b, M * (1 + GUARD_D), tol) + margin
        if g_lo < GUARD_BELOW or g_hi - g_lo > -GUARD_SLOPE:
            return (f"crossing: gap(M*(1-{GUARD_D:g}))={g_lo:+.3e}, gap(M*(1+{GUARD_D:g}))={g_hi:+.3e} "
                    f"(need >= {GUARD_BELOW:g} and difference < -{GUARD_SLOPE:g})")
    return ""


def run_crossing(t, T, label, stage, blocks, lo, hi, tols=(1e-12,), margins=(0.0,), N=None):
    todo = [(tol, m) for tol in tols for m in margins if not T.find(target=t, label=label, tol=tol, margin=m)]
    if not todo:
        return [float(r["M_star"]) for r in T.find(target=t, label=label, tol=1e-12, margin=0.0)]
    t0 = time.perf_counter()
    b = LeanBasis(TARGETS[t]["L"], blocks, TARGETS[t]["r"], tols=tuple(sorted({x for x, _ in todo})))
    out = []
    for tol, m in todo:
        t1 = time.perf_counter()
        try:
            Ms = crossing(t, b, lo, hi, tol, m)
        except ValueError as exc:
            log(t, f"  {label} tol={tol:g} margin={m:g}: {exc}")
            Ms = float("nan")
        why = guard(t, b, Ms, tol, N, margin=m)
        if why:
            log(t, f"WARNING spurious root: stage{stage} {label} tol={tol:g} margin={m:g} M*={Ms:.10f} rejected ({why})")
            why, Ms = f"M*={Ms:.10f}; {why}", float("nan")
        T.add(target=t, label=label, stage=stage, size=b.size, kept=b.proj[tol]["kept"], tol=tol, margin=m,
              M_star=Ms, seconds=round(time.perf_counter() - t1 + (b.total_time if (tol, m) == todo[0] else 0), 1),
              spurious=int(bool(why)), reason=why)
        log(t, f"stage{stage} {label}: size={b.size} kept={b.proj[tol]['kept']} tol={tol:g} margin={m:g} "
               f"M*={Ms:.10f} ({time.perf_counter() - t1:.0f}s; build+project {b.total_time:.0f}s)")
        out.append(Ms)
    del b
    gc.collect()
    log(t, f"  ({time.perf_counter() - t0:.0f}s)")
    return out


def tab(t):
    return Table(os.path.join(HERE, f"{t}_crossings.csv"), CROSS_COLS)


def stage_compact(t):
    cfg, sc = TARGETS[t], stage_scales(t)
    T = tab(t)
    lo, hi, _ = brackets(t)
    for N in cfg["Ns"]:
        Ms = run_crossing(t, T, f"compact N={N}", 2, blocks_for(t, sc, N), lo, hi)
        if Ms and np.isfinite(Ms[0]):
            hi = Ms[0] * 1.0000001


def std_label(t, alpha_s=None):
    cfg = TARGETS[t]
    if cfg["L"] == 0:
        return f"compact N={cfg['ND']} + diffuse alpha_s={alpha_s:g}"
    return f"compact N={cfg['ND']} + {len(cfg['diffuse'])} diffuse"


def stage_diffuse(t):
    cfg, sc = TARGETS[t], stage_scales(t)
    T = tab(t)
    lo, hi, _ = brackets(t)
    r = T.find(target=t, label=f"compact N={cfg['ND']}", tol=1e-12, margin=0.0)
    if r and np.isfinite(float(r[0]["M_star"])):
        hi = float(r[0]["M_star"]) * 1.0000001
    if cfg["L"] == 0:
        for a in LADDER:
            run_crossing(t, T, std_label(t, a), 3, blocks_for(t, sc, cfg["ND"], alpha_s=a), lo, hi, N=cfg["ND"])
    else:
        for nd in range(1, len(cfg["diffuse"]) + 1):
            label = f"compact N={cfg['ND']} + {nd} diffuse"
            Ms = run_crossing(t, T, label, 3, blocks_for(t, sc, cfg["ND"], n_diffuse=nd), lo, hi, N=cfg["ND"])
            if Ms and np.isfinite(Ms[0]):
                hi = Ms[0] * 1.0000001


def best_alpha(t):
    """L = 0: ladder member with the lowest crossing (used for the checks)."""
    T = tab(t)
    best, ba = np.inf, 1e-3
    for a in LADDER:
        r = T.find(target=t, label=std_label(t, a), tol=1e-12, margin=0.0)
        if r and np.isfinite(float(r[0]["M_star"])) and float(r[0]["M_star"]) < best:
            best, ba = float(r[0]["M_star"]), a
    return ba


def stage_checks(t):
    cfg, sc = TARGETS[t], stage_scales(t)
    T = tab(t)
    lo, hi, _ = brackets(t)
    if cfg["L"] == 0:
        a = best_alpha(t)
        kw = dict(alpha_s=a)
    else:
        a = None
        kw = dict(n_diffuse=len(cfg["diffuse"]))
    run_crossing(t, T, std_label(t, a), 4, blocks_for(t, sc, cfg["ND"], **kw), lo, hi,
                 tols=(1e-12, 1e-10), margins=(0.0, 1e-11), N=cfg["ND"])
    for N in (cfg["ND"] - 2, cfg["ND"] + 2):
        run_crossing(t, T, std_label(t, a).replace(f"N={cfg['ND']}", f"N={N}"), 4, blocks_for(t, sc, N, **kw), lo, hi,
                     N=N)
    run_crossing(t, T, std_label(t, a) + " (diffuse enlarged)", 4, blocks_for(t, sc, cfg["ND"], big=True, **kw), lo, hi,
                 N=cfg["ND"])


def best_mstar(t):
    T = tab(t)
    v = [float(r["M_star"]) for r in T.rows if r["tol"] == repr(1e-12) and r["margin"] == repr(0.0)
         and r.get("spurious", "") != "1"]
    v = [x for x in v if np.isfinite(x)]
    if not v:
        raise RuntimeError("no crossing yet")
    return min(v)


def stage_binding(t):
    cfg, sc = TARGETS[t], stage_scales(t)
    B = Table(os.path.join(HERE, f"{t}_binding.csv"), BIND_COLS)
    a = best_alpha(t) if cfg["L"] == 0 else None
    label = std_label(t, a)
    Ms = best_mstar(t)
    todo = [x for x in XB if not B.find(target=t, label=label, x=x)]
    if not todo:
        return
    kw = dict(alpha_s=a) if cfg["L"] == 0 else dict(n_diffuse=len(cfg["diffuse"]))
    b = LeanBasis(cfg["L"], blocks_for(t, sc, cfg["ND"], **kw), cfg["r"])
    k = cfg["k"]
    for x in todo:
        t1 = time.perf_counter()
        M = Ms * (1 + x)
        e = b.energies(M, k + 2)
        Eth = threshold(M, cfg["r"])
        why = guard(t, b, M, 1e-12, cfg["ND"])
        if why:
            log(t, f"WARNING spurious root: stage5 M={M:.8f} rejected ({why})")
        nan = float("nan")
        B.add(target=t, label=label, M=M, x=x, E=nan if why else float(e[k]), E_th=Eth,
              binding=nan if why else float(Eth - e[k]), E_next=nan if why else float(e[k + 1]),
              seconds=round(time.perf_counter() - t1, 1), spurious=int(bool(why)), reason=why)
        log(t, f"stage5 M={M:.8f} E_th-E_k={Eth - e[k]:.6e} ({time.perf_counter() - t1:.0f}s)")
    del b
    gc.collect()


# ---------------------------------------------------------------- dry run / test
def plan(t, sc):
    cfg = TARGETS[t]
    p = [(f"compact N={N}", dict(N=N)) for N in cfg["Ns"]]
    if cfg["L"] == 0:
        p += [(f"ladder {a:g}", dict(N=cfg["ND"], alpha_s=a)) for a in LADDER]
        a = 1e-3
        p += [("check tol/margin x4", dict(N=cfg["ND"], alpha_s=a, reps=4)),
              ("check N-2", dict(N=cfg["ND"] - 2, alpha_s=a)), ("check N+2", dict(N=cfg["ND"] + 2, alpha_s=a)),
              ("check diffuse enlarged", dict(N=cfg["ND"], alpha_s=a, big=True))]
    else:
        nd = len(cfg["diffuse"])
        p += [(f"+{i} diffuse", dict(N=cfg["ND"], n_diffuse=i)) for i in range(1, nd + 1)]
        p += [("check tol/margin x4", dict(N=cfg["ND"], n_diffuse=nd, reps=4)),
              ("check N-2", dict(N=cfg["ND"] - 2, n_diffuse=nd)), ("check N+2", dict(N=cfg["ND"] + 2, n_diffuse=nd)),
              ("check diffuse enlarged", dict(N=cfg["ND"], n_diffuse=nd, big=True))]
    return p


def dry(t):
    cfg = TARGETS[t]
    lo, hi, Mopt = brackets(t)
    sc = json.load(open(scales_path(t))) if os.path.exists(scales_path(t)) else \
        dict(M_opt=Mopt, scales=[[4.0, 4.0, 1.3]] * (2 if cfg["k"] < 2 else (3 if cfg["k"] < 4 else 4)))
    log(t, f"=== --dry {t}: r={cfg['r']} L={cfg['L']} root k={cfg['k']}; bracket [{lo:.4f}, {hi:.4f}]; "
           f"threads {os.environ.get('THREEBODY_BLAS_THREADS')}; scales {'(stage 1 done)' if os.path.exists(scales_path(t)) else '(placeholder)'}")
    kw0 = dict(alpha_s=1e-3) if cfg["L"] == 0 else dict(n_diffuse=2)
    Nc = cfg["ND"] - 2
    t0 = time.perf_counter()
    b = LeanBasis(cfg["L"], blocks_for(t, sc, Nc, **kw0), cfg["r"])
    tb = b.total_time
    t1 = time.perf_counter()
    gap(t, b, 0.5 * (lo + hi))
    te = time.perf_counter() - t1
    n0 = b.size
    del b
    log(t, f"calibration: N={Nc} + diffuse, size {n0}: build+project {tb:.1f}s, one eigensolve {te:.2f}s")
    total, peak = 0.0, 0
    for label, q in plan(t, sc):
        reps = q.pop("reps", 1)
        N = q.pop("N")
        n, nch = nfun(t, sc, N, **q)
        tt = tb * (n / n0) ** 2.5 + 20 * te * (n / n0) ** 3
        total += tt * (1 if reps == 1 else 2) + (reps - 1) * 20 * te * (n / n0) ** 3
        peak = max(peak, n)
        log(t, f"  {label:28s}: size {n:6d} ({nch} ch) est {tt / 60:6.1f} min x{reps}, memory ~{6 * n * n * 8 / 1e9:.1f} GB")
    log(t, f"TOTAL estimated {total / 3600:.2f} h (+ stage 1 scales ~2-10 min, stage 5 ~{len(XB) * te * 2 / 60:.0f} min); "
           f"peak size {peak}, peak memory ~{6 * peak * peak * 8 / 1e9:.1f} GB (~6 n^2 doubles)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", choices=sorted(TARGETS))
    ap.add_argument("--stage", type=int, nargs="*", default=[1, 2, 3, 4, 5])
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--test", nargs=2, metavar=("N", "M"), help="one E_k - E_th (short test; uses stage-1 scales)")
    a = ap.parse_args()
    if a.list or not a.target:
        for k, v in TARGETS.items():
            print(k, v)
        return
    t = a.target
    if a.dry:
        return dry(t)
    if a.test:
        sc = stage_scales(t)
        cfg = TARGETS[t]
        N, M = int(a.test[0]), float(a.test[1])
        kw = dict(alpha_s=1e-3) if cfg["L"] == 0 else dict(n_diffuse=len(cfg["diffuse"]))
        b = LeanBasis(cfg["L"], blocks_for(t, sc, N, **kw), cfg["r"])
        e = b.energies(M, cfg["k"] + 2)
        log(t, f"test N={N} M={M}: size={b.size} kept={b.proj[1e-12]['kept']} E_k-E_th={e[cfg['k']] - threshold(M, cfg['r']):.6e} "
               f"E_k+1-E_th={e[cfg['k'] + 1] - threshold(M, cfg['r']):.3e} (build+project {b.total_time:.0f}s)")
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
