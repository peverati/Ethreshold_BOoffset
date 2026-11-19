"""Nuclear (vibrational-rotational) problem on the BO potential.

    -(1/(2 mu_N)) u'' + [V(R) + J(J+1)/(2 mu_N R^2)] u = E u,   u(0) = 0,
    V(R) = mu_e U(mu_e R) [+ E_ad(R)],  mu_N = M/2,  mu_e = 2M/(2M+1).

Method: sinc-DVR (Colbert-Miller, radial variant with phi(0)=0) in a mapped
variable x, R = R_s sinh(x).  The map is uniform near R = 0 and logarithmic at
large R, so one grid resolves the inner wall, the well, and the long
-C4/R^4 tail.  With u = sqrt(R'(x)) phi(x):

    -phi'' + [Q(x) + J(J+1) coth^2 x + 2 mu_N R'^2 V] phi = E 2 mu_N R'^2 phi,
    Q = (3/4)(R''/R')^2 - (1/2) R'''/R' = (3/4) tanh^2 x - 1/2,

a symmetric-definite problem with diagonal metric, solved with LAPACK.  The
grid step is set from the largest local wavenumber in x (adapts to M), the box
from the length scale beta = sqrt(2 mu_N C4) of the -C4/R^4 tail.  Every call
solves on two grids (base, and finer + larger box) and returns the difference
as an error estimate per level.

Level counting (independent of the DVR): number of nodes in (0, inf) of the
zero-energy (threshold) solution, integrated with Numerov on the same map out
to R_b >> beta and continued analytically beyond R_b (J=0: exact zero-energy
solution of -C4/R^4, u = R sin(beta/R + phi); J>0: a R^{J+1} + b R^{-J}).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import eigh

try:  # cap BLAS threads: oversubscribed OpenBLAS was 40x slower on a busy machine
    from threadpoolctl import threadpool_limits
except ImportError:  # pragma: no cover
    threadpool_limits = None
import os
_NTHREADS = int(os.environ.get("BO_H2PLUS_THREADS", "4"))

from . import potential as P

__all__ = ["Level", "LevelSet", "solve_levels", "count_levels", "thresholds", "C4"]


def thresholds(M, adiabatic=False):
    """Dissociation thresholds (hartree): (BO limit, exact two-body X + (Xe) limit).

    BO: -mu_e/2 (+ mu_e^2/(4M) if the adiabatic correction is included).
    Exact: -M/(2(M+1)) (hydrogen-like atom with reduced mass M/(M+1)).
    """
    me = P.mu_e(M)
    bo = -me / 2 + (me * me / (4 * M) if adiabatic else 0.0)
    return bo, -M / (2 * (M + 1))


_CD = None
_RS = None


def C4(M, adiabatic=False):
    """Coefficient of the -C4/R^4 tail of V(R) relative to the threshold."""
    global _CD
    me = P.mu_e(M)
    c = 9.0 / (4.0 * me ** 3)
    if adiabatic:
        if _CD is None:
            _CD = (0.25 - P.D(np.array([P.R_HI]))[0]) * P.R_HI ** 4
        c += _CD / (M * me ** 2)
    return c


def _V(R, M, adiabatic):
    v = P.V_BO(R, M)
    if adiabatic:
        v = v + P.E_ad(R, M)
    return v


@dataclass
class Level:
    v: int
    J: int
    E: float
    binding: float          # threshold - E (BO threshold)
    R_mean: float
    R2_mean: float
    dR: float               # sqrt(<R^2> - <R>^2)
    inv_R: float
    inv_R2: float
    T_kin: float            # E - <V>  (nuclear kinetic incl. rotation)
    T_rad: float            # radial kinetic only
    lindemann: float        # dR/<R>
    err_E: float = float("nan")    # |E(base grid) - E(fine grid)|
    err_R: float = float("nan")    # |<R>(base) - <R>(fine)|


@dataclass
class LevelSet:
    M: float
    J: int
    adiabatic: bool
    threshold: float
    levels: list
    R: np.ndarray           # radial grid (fine grid)
    P: np.ndarray           # P[v, i] = |u_v(R_i)|^2, normalized: int P dR = 1
    n_dvr: int              # number of DVR eigenvalues below threshold (fine grid)
    grid_info: dict = field(default_factory=dict)

    def energies(self):
        return np.array([l.E for l in self.levels])


def _grid_params(M, J, adiabatic, ppw, box):
    me = P.mu_e(M)
    muN = M / 2.0
    thr = thresholds(M, adiabatic)[0]
    # Map scale: the Coulomb behaviour u = R(1 + mu_N R + ...) at R = 0 is not odd in x and
    # limits the sinc-DVR to O((R_s dx)^2) when the wavefunction reaches the origin
    # (small M).  A small R_s cures that at the price of a larger matrix norm (round-off
    # ~1e-10..1e-9), so it is used only for M < 300.
    Rs = _RS if _RS else (0.05 if M < 300 else 1.0)
    beta = math.sqrt(2 * muN * C4(M, adiabatic))
    Rmax = box * max(60.0 / me, 100.0 * beta)
    # largest local wavenumber in x over the classically allowed region at threshold
    r = np.geomspace(0.05, 60.0 / me, 4000)
    kin = thr - _V(r, M, adiabatic) - J * (J + 1) / (2 * muN * r * r)
    k = np.sqrt(2 * muN * np.clip(kin, 0, None)) * np.sqrt(r * r + Rs * Rs)
    kmax = max(k.max(), 1.0)
    dx = min(math.pi / (ppw * kmax), 0.1 / ppw)
    return dict(Rs=Rs, Rmax=Rmax, dx=dx, beta=beta, thr=thr, kmax=kmax)


def _dvr(M, J, adiabatic, Rs, Rmax, dx, thr, nmax=None):
    muN = M / 2.0
    xmax = math.asinh(Rmax / Rs)
    n = int(math.ceil(xmax / dx))
    i = np.arange(1, n + 1, dtype=float)
    x = i * dx
    R = Rs * np.sinh(x)
    fp = Rs * np.cosh(x)
    # radial sinc-DVR kinetic matrix for -d^2/dx^2 (Colbert & Miller 1992, eq. A7)
    I, Jm = np.meshgrid(i, i, indexing="ij")
    D = I - Jm
    sgn = 1.0 - 2.0 * (np.abs(D) % 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        T = sgn * (2.0 / D ** 2 - 2.0 / (I + Jm) ** 2)
    T[np.diag_indices(n)] = np.pi ** 2 / 3 - 1.0 / (2 * i * i)
    T /= dx * dx
    Q = 0.75 * np.tanh(x) ** 2 - 0.5
    B = 2 * muN * fp * fp
    diag = Q + J * (J + 1) / np.tanh(x) ** 2 + B * _V(R, M, adiabatic)
    A = T
    A[np.diag_indices(n)] += diag
    s = 1 / np.sqrt(B)
    C = A * np.outer(s, s)
    if threadpool_limits is not None:
        with threadpool_limits(limits=_NTHREADS):
            w, vec = eigh(C, subset_by_value=(-np.inf, thr), driver="evr")
    else:
        w, vec = eigh(C, subset_by_value=(-np.inf, thr), driver="evr")
    return R, fp, x, w, vec


def _props(M, J, adiabatic, R, fp, dx, w, vec, thr):
    muN = M / 2.0
    V = _V(R, M, adiabatic)
    cent = J * (J + 1) / (2 * muN * R * R)
    out = []
    Pd = []
    for v in range(len(w)):
        c2 = vec[:, v] ** 2
        c2 = c2 / c2.sum()
        Rm = np.dot(c2, R)
        R2 = np.dot(c2, R * R)
        Vm = np.dot(c2, V)
        cm = np.dot(c2, cent)
        dR = math.sqrt(max(R2 - Rm * Rm, 0.0))
        out.append(Level(v=v, J=J, E=w[v], binding=thr - w[v], R_mean=Rm, R2_mean=R2, dR=dR,
                         inv_R=np.dot(c2, 1 / R), inv_R2=np.dot(c2, 1 / (R * R)),
                         T_kin=w[v] - Vm, T_rad=w[v] - Vm - cm, lindemann=dR / Rm))
        Pd.append(c2 / (dx * fp))
    return out, np.array(Pd)


def solve_levels(M, J=0, adiabatic=False, ppw=2.0, box=1.0, check=True, auto_box=True,
                 max_box=4096.0):
    """All bound levels of the BO (or BO + adiabatic) nuclear problem for mass M, rotation J.

    ppw: grid density (dx = pi/(ppw * k_max)); box: scale of R_max.
    auto_box: enlarge the box (x4 steps) until the DVR finds as many levels as the
              independent node count (count_levels); needed for levels bound by
              much less than 1/(2 mu_N beta^2).
    check=True solves again on a grid with ppw/1.4 and a box twice as large and
    stores |difference| in Level.err_E / err_R (an honest per-level error estimate).
    """
    n_count = count_levels(M, J, adiabatic) if auto_box else None
    while True:
        g = _grid_params(M, J, adiabatic, ppw, box)
        R, fp, x, w, vec = _dvr(M, J, adiabatic, g["Rs"], g["Rmax"], g["dx"], g["thr"])
        if auto_box and len(w) < n_count and box < max_box:
            box *= 4.0
            continue
        if auto_box and len(w) > 0 and box < max_box:
            # make sure the least-bound level decays inside the box: R_max >= 25/kappa
            kappa = math.sqrt(2 * (M / 2.0) * max(g["thr"] - w[-1], 1e-300))
            need = 25.0 / kappa
            if g["Rmax"] < need:
                box = min(box * 2.0 ** math.ceil(math.log2(need / g["Rmax"])), max_box)
                continue
        break
    lv, Pd = _props(M, J, adiabatic, R, fp, g["dx"], w, vec, g["thr"])
    if check:
        g2 = _grid_params(M, J, adiabatic, ppw / 1.4, box * 2)
        R2, fp2, x2, w2, vec2 = _dvr(M, J, adiabatic, g2["Rs"], g2["Rmax"], g2["dx"], g2["thr"])
        lv2, _ = _props(M, J, adiabatic, R2, fp2, g2["dx"], w2, vec2, g2["thr"])
        for a in lv:
            if a.v < len(lv2):
                a.err_E = abs(a.E - lv2[a.v].E)
                a.err_R = abs(a.R_mean - lv2[a.v].R_mean)
    info = dict(g, npts=len(R), box=box, n_count=n_count)
    return LevelSet(M=M, J=J, adiabatic=adiabatic, threshold=g["thr"], levels=lv, R=R, P=Pd,
                    n_dvr=len(w), grid_info=info)


def count_levels(M, J=0, adiabatic=False, ppw=8.0, box=3.0, return_details=False):
    """Number of bound levels below the BO threshold by node counting at E = threshold."""
    g = _grid_params(M, J, adiabatic, ppw, 1.0)
    muN = M / 2.0
    Rs, dx, thr, beta = g["Rs"], g["dx"], g["thr"], g["beta"]
    Rb = box * max(60.0 / P.mu_e(M), 30.0 * beta)
    n = int(math.ceil(math.asinh(Rb / Rs) / dx))
    x = np.arange(n + 1) * dx
    R = Rs * np.sinh(x)
    fp = Rs * np.cosh(x)
    with np.errstate(divide="ignore", invalid="ignore"):
        gfun = (0.75 * np.tanh(x) ** 2 - 0.5 + J * (J + 1) / np.tanh(x) ** 2
                + 2 * muN * fp * fp * (_V(R, M, adiabatic) - thr))
    gfun[0] = 0.0
    h12 = dx * dx / 12.0
    # Numerov is unstable where h^2 g/12 is large (deep inside the Coulomb wall for
    # large M): start (phi = 0) at the first point where h^2 g/12 < 0.5.  The true
    # solution is exponentially negligible there, so this Dirichlet start is harmless.
    big = np.nonzero(h12 * gfun[1:] >= 0.5)[0]
    k0 = 0 if len(big) == 0 else int(big[-1]) + 1
    # (big[-1]+1 is the last offending index; all later points are fine)
    a = (1 - h12 * gfun).tolist()
    b = (2 + 10 * h12 * gfun).tolist()
    if k0 == 0:
        # start at x1, x2 from the Frobenius series of the regular solution at R = 0:
        # u = R^{J+1} (1 + c1 R + c2 R^2), V = 1/R + V0 + O(R^2) (Coulomb repulsion),
        # V0 = mu_e * E_el(0) = -2 mu_e (+ mu_e^2 D(0)/M ~ 0).  Starting with phi(0) = 0
        # instead would drop the finite limit g*phi at x = 0 and spoil the O(dx^4) order.
        me = P.mu_e(M)
        V0 = -2.0 * me
        V2 = me * (8.0 / 3.0) * me ** 2   # E_el(R) = -2 + (8/3) R^2 + ... (checked in tests)
        c = [1.0]
        for kk in range(1, 12):
            rhs = c[kk - 1] + (V0 - thr) * (c[kk - 2] if kk >= 2 else 0.0) \
                + V2 * (c[kk - 4] if kk >= 4 else 0.0)
            c.append(2 * muN * rhs / (kk * (kk + 2 * J + 1)))
        rr = R[1:3]
        us = rr ** (J + 1) * np.polynomial.polynomial.polyval(rr, c)
        phi_prev, phi = us / np.sqrt(fp[1:3])
        kstart = 2
    else:
        phi_prev, phi = 0.0, 1e-300
        kstart = k0 + 1
    nodes = 0
    phis = [phi_prev, phi]
    for k in range(kstart, n):
        new = (b[k] * phi - a[k - 1] * phi_prev) / a[k + 1]
        if new == 0.0 or (new > 0) != (phi > 0):
            nodes += 1
        phi_prev, phi = phi, new
        if abs(phi) > 1e200:
            phi_prev *= 1e-200
            phi *= 1e-200
        if k >= n - 1:
            phis = [phi_prev, phi]
    # analytic continuation beyond R_b
    Ra, Rbb = R[n - 1], R[n]
    ua, ub = math.sqrt(fp[n - 1]) * phis[0], math.sqrt(fp[n]) * phis[1]
    extra = 0
    if J == 0:
        bb = math.sqrt(2 * muN * C4(M, adiabatic))
        sa, sb = bb / Ra, bb / Rbb
        num = ub * Ra * math.sin(sa) - ua * Rbb * math.sin(sb)
        den = ua * Rbb * math.cos(sb) - ub * Ra * math.cos(sa)
        ph = math.atan2(num, den)
        # zeros of sin(theta) for theta in (ph, ph + sb)
        extra = math.floor((ph + sb) / math.pi) - math.floor(ph / math.pi)
        asym = math.sin(ph)
    else:
        Mx = np.array([[Ra ** (J + 1), Ra ** (-J)], [Rbb ** (J + 1), Rbb ** (-J)]])
        aa, bc = np.linalg.solve(Mx, [ua, ub])
        if aa * bc < 0 and (-bc / aa) ** (1.0 / (2 * J + 1)) > Rbb:
            extra = 1
        asym = aa
    total = nodes + extra
    if return_details:
        return dict(n=total, inner_nodes=nodes, extra=extra, Rb=Rbb, npts=n, dx=dx,
                    R_start=R[k0], asym=asym)
    return total
