"""Unequal nuclear masses {A+, B+, e-}: masses, thresholds, and the 1s sigma_g / 2p sigma_u coupling functions.

Added 2026-11-16 (new module; nothing else in bo_h2plus was changed).  Theory and conventions:
paper_data/unequal_masses/THEORY.md.  Summary:

Nuclear masses m_A = M, m_B = r M (r >= 1: B is the heavier nucleus), electron mass 1.  Jacobi coordinates:
R = R_B - R_A, rho = electron relative to the nuclear centre of mass (no mass-polarization term).

    mu_N = r M/(1+r),   M_eff = 2 mu_N,   mu_e = (1+r)M/((1+r)M+1),   lam = M_eff/mu_e,
    beta = (r-1)/(2(1+r))   (nuclear CM sits at +beta R from the midpoint, towards B).

The clamped-nuclei electronic problem is that of H2+ (midpoint coordinates rho_m = rho + beta R).  At fixed
rho the nuclear derivative picks up an extra electronic derivative:

    d_R|rho     = d_R|rho_m + beta d_z           (parallel)
    d_{R_x}|rho = (1/R)(x d_z - z d_x) + beta d_x  (perpendicular; R along z)

so for the electronic states i, j in {g = 1s sigma_g, u = 2p sigma_u} (unit electron mass, functions of
the scaled distance r = mu_e R):

    K_ij  = <grad_R phi_i | grad_R phi_j>_rho:  K_gg = D_g + beta^2 t_g,  K_uu = D_u + beta^2 t_u,
                                                K_gu = beta k_gu
    a_gu  = <phi_g | d_R phi_u>_rho = beta <g|d_z|u> = -a_ug
    D_i   = <d_R phi_i|d_R phi_i>_rho_m + <(L_x^2+L_y^2)>/R^2   (bo_h2plus D for i = g)
    t_i   = <p^2>_i = -2 (U_i + R U_i')                       (virial)
    k_gu  = <F_g|Z_u> + <Z_g|F_u> + 2 x (<G_g|S_u> + <S_g|G_u>)/R  (fields defined in derivative_fields)

Asymptotically (R -> inf): D -> 1/4, t -> 1, k_gu -> -1, <g|d_z|u> -> 0.

Scaled radial equations (r = mu_e R, multiplied by M_eff/mu_e^2; lam enters linearly):

    -u'' - 2 a u' - a' u + [L(L+1)/r^2 + K(r)] u + lam U(r) u = (M_eff/mu_e^2) E u      (2-channel)

The g-only BO+ad problem uses K_gg = D_g + beta^2 t_g, counted from its own asymptote 1/4 + beta^2.
"""
from __future__ import annotations

import math
import os

import numpy as np

from . import electronic as el
from .adiabatic import _nodes as _ad_nodes
from .electronic import Wavefunction

__all__ = ["Masses", "beta_of", "lam_of_M", "M_of_lam", "derivative_fields", "gu_functions", "gu_table",
           "kinetic_t_virial"]


# ---------------------------------------------------------------- masses and thresholds
def beta_of(r):
    """Offset of the nuclear centre of mass from the midpoint, in units of R (towards the heavier B)."""
    return (r - 1.0) / (2.0 * (1.0 + r))


def lam_of_M(M, r):
    """lam = M_eff/mu_e = 2 r M/(1+r) + 2 r/(1+r)^2 (exact; r = 1 gives M + 1/2)."""
    return 2.0 * r * M / (1.0 + r) + 2.0 * r / (1.0 + r) ** 2


def M_of_lam(lam, r):
    """Inverse of lam_of_M: M = (1+r) lam/(2r) - 1/(1+r)."""
    return (1.0 + r) * lam / (2.0 * r) - 1.0 / (1.0 + r)


class Masses:
    """All mass-dependent constants of {A+ (M), B+ (rM), e- (1)}."""

    def __init__(self, M, r):
        self.M, self.r = float(M), float(r)
        self.mA, self.mB = self.M, self.r * self.M
        self.mu_N = self.mA * self.mB / (self.mA + self.mB)
        self.M_eff = 2.0 * self.mu_N
        self.mu_e = (self.mA + self.mB) / (self.mA + self.mB + 1.0)
        self.lam = self.M_eff / self.mu_e
        self.beta = beta_of(self.r)
        self.mu_Ae = self.mA / (self.mA + 1.0)
        self.mu_Be = self.mB / (self.mB + 1.0)
        # exact two-body thresholds: (B e) + A+ is the lowest for r > 1
        self.E_th = -0.5 * max(self.mu_Ae, self.mu_Be)
        self.E_th_A = -0.5 * self.mu_Ae
        self.E_th_B = -0.5 * self.mu_Be
        # ion-atom reduced mass of the open channel A+ + (B e)
        self.mu_rel = self.mA * (self.mB + 1.0) / (self.mA + self.mB + 1.0)

    def scaled(self, E):
        """Energy E (hartree) in the scaled units of the radial equation (multiply by M_eff/mu_e^2)."""
        return E * self.M_eff / self.mu_e ** 2

    def thresholds_scaled(self):
        """Exact thresholds relative to the g-channel BO+ad asymptote (-lam/2 + 1/4 + beta^2), scaled units.

        To O(1) the open channel lies at -beta and the closed one at +beta (exact values returned)."""
        ref = -0.5 * self.lam + 0.25 + self.beta ** 2
        return self.scaled(self.E_th_B) - ref, self.scaled(self.E_th_A) - ref


# ---------------------------------------------------------------- electronic fields
def kinetic_t_virial(R, state="1sg"):
    """<p^2> (unit mass) from the molecular virial theorem: t = -2 (U + R dU/dR), U incl. 1/R."""
    U = el.total_energy(R, state)
    return -2.0 * (U + R * el.dEdR(R, state, total=True))


def _psi(wf, xi, eta, deriv=False):
    if not deriv:
        return wf.N * np.outer(wf.X(xi), wf.Y(eta))
    X, dX = wf.X(xi, deriv=True)
    Y, dY = wf.Y(eta, deriv=True)
    return wf.N * np.outer(X, Y), wf.N * np.outer(dX, Y), wf.N * np.outer(X, dY)


def derivative_fields(R, xi, eta, state="1sg", h=None):
    """Fields of one electronic state on the (xi, eta) grid (unit mass, eta = +1 at z = +R/2):

    psi, F = d_R psi|rho_m (midpoint coordinates fixed), Z = d_z psi, S = d_s psi (s = cylindrical radius),
    G = s d_z psi - z d_s psi (so (x d_z - z d_x) psi = cos(phi) G).  Returns (wf, dict)."""
    R = float(R)
    if h is None:
        h = 1e-3 * R
    wf = Wavefunction(R, state)
    psi, pxi, peta = _psi(wf, xi, eta, deriv=True)
    XI, ET = np.meshgrid(xi, eta, indexing="ij")
    den = XI ** 2 - ET ** 2
    vals = {s: _psi(Wavefunction(R + s * h, state), xi, eta) for s in (-2, -1, 1, 2)}
    dR = (8 * (vals[1] - vals[-1]) - (vals[2] - vals[-2])) / (12 * h)
    F = dR - (XI * (XI ** 2 - 1) * pxi + ET * (1 - ET ** 2) * peta) / (R * den)
    root = np.sqrt(np.clip((XI ** 2 - 1) * (1 - ET ** 2), 0.0, None))
    Z = (2.0 / R) * (ET * (XI ** 2 - 1) * pxi + XI * (1 - ET ** 2) * peta) / den
    S = (2.0 / R) * root * (XI * pxi - ET * peta) / den
    s, z = 0.5 * R * root, 0.5 * R * XI * ET
    G = s * Z - z * S
    return wf, dict(psi=psi, F=F, Z=Z, S=S, G=G, z=z)


def gu_functions(R, na=160, nb=160):
    """All g/u coupling functions at internuclear distance R (unit mass).  See module docstring.

    Returns dict: U_g, U_u (total energies incl. 1/R), D_g, D_u, t_g, t_u (quadrature), t_g_vir, t_u_vir,
    dz_gu = <g|d_z|u>, z_gu = <g|z|u>, dz_ug = <u|d_z|g>, k_gu, plus norms (checks)."""
    R = float(R)
    wg0 = Wavefunction(R, "1sg")
    wu0 = Wavefunction(R, "2pu")
    p = min(wg0.p, wu0.p)
    xi, wx, eta, we = _ad_nodes(p, na, nb)
    XI, ET = np.meshgrid(xi, eta, indexing="ij")
    W = 2 * math.pi * (R / 2) ** 3 * np.outer(wx, we) * (XI ** 2 - ET ** 2)   # m = 0 weight (phi integrated)
    wg, g = derivative_fields(R, xi, eta, "1sg")
    wu, u = derivative_fields(R, xi, eta, "2pu")
    I = lambda A, B: float(np.sum(W * A * B))
    out = dict(R=R, U_g=wg.E + 1.0 / R, U_u=wu.E + 1.0 / R,
               norm_g=I(g["psi"], g["psi"]), norm_u=I(u["psi"], u["psi"]), overlap_gu=I(g["psi"], u["psi"]))
    for n, f in (("g", g), ("u", u)):
        out[f"D_{n}"] = I(f["F"], f["F"]) + I(f["G"], f["G"]) / R ** 2
        out[f"Drad_{n}"] = I(f["F"], f["F"])
        out[f"t_{n}"] = I(f["Z"], f["Z"]) + I(f["S"], f["S"])
    out["dz_gu"] = I(g["psi"], u["Z"])
    out["dz_ug"] = I(u["psi"], g["Z"])
    out["z_gu"] = I(g["psi"], u["z"] * u["psi"])
    out["k_gu"] = I(g["F"], u["Z"]) + I(g["Z"], u["F"]) + (I(g["G"], u["S"]) + I(g["S"], u["G"])) / R
    out["k_gu_rad"] = I(g["F"], u["Z"]) + I(g["Z"], u["F"])
    out["F_overlap_gu"] = I(g["F"], u["F"])          # 0 by parity (check)
    return out


# ---------------------------------------------------------------- cached table
_CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_cache")


def default_grid():
    """R grid for the coupling tables: 0.02 .. 60 bohr (denser below 20)."""
    return np.unique(np.concatenate([np.geomspace(0.02, 1.0, 40), np.arange(1.05, 20.0, 0.1),
                                     np.arange(20.25, 40.0, 0.25), np.arange(40.5, 60.01, 0.5)]))


def gu_table(grid=None, path=None, rebuild=False, verbose=False):
    """Table of gu_functions on a grid (cached in _cache/gu_v1.npz).  About 1-3 min to build."""
    save = path if path is not None else (os.path.join(_CACHE, "gu_v1.npz") if grid is None else None)
    if grid is None and save and os.path.exists(save) and not rebuild:
        path = save
        d = np.load(path)
        return {k: d[k] for k in d.files}
    grid = default_grid() if grid is None else np.asarray(grid, float)
    rows = []
    for i, R in enumerate(grid):
        rows.append(gu_functions(R))
        if verbose and i % 25 == 0:
            print(f"  gu_table {i}/{len(grid)} R={R:.3f}", flush=True)
    keys = [k for k in rows[0] if k != "R"]
    tab = dict(R=grid, **{k: np.array([q[k] for q in rows]) for k in keys})
    tab["t_g_vir"] = np.array([kinetic_t_virial(R, "1sg") for R in grid])
    tab["t_u_vir"] = np.array([kinetic_t_virial(R, "2pu") for R in grid])
    if save:
        os.makedirs(os.path.dirname(save), exist_ok=True)
        np.savez(save, **tab)
    return tab
