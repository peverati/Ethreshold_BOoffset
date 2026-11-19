"""Diagonal Born-Oppenheimer (adiabatic) correction for {X+, X+, e-}.

In Jacobi coordinates (R = X-X vector, rho = electron relative to the X-X
midpoint) the kinetic operator is -(1/2mu_N) Lap_R - (1/2mu_e) Lap_rho with no
cross (mass-polarization) term, because the X-X midpoint is the nuclear centre
of mass for equal masses.  The diagonal correction is therefore only

    E_ad(R) = (1/(2 mu_N)) < grad_R psi | grad_R psi >_rho      (rho held fixed)

For electron mass mu the electronic function is psi_mu(rho;R) = mu^{3/2}
psi_1(mu rho; mu R), so

    E_ad(R; M) = mu_e^2 D(mu_e R) / (2 mu_N) = mu_e^2 D(mu_e R) / M,
    D(R) = <d_R psi|d_R psi>_rho + <L_x^2 + L_y^2>/R^2        (unit mass),

(the second term comes from rotating the R direction at fixed rho).  In
prolate spheroidal coordinates

    d_R psi|_rho = d_R psi|_{xi,eta}
                   - (1/R) [xi(xi^2-1) psi_xi + eta(1-eta^2) psi_eta]/(xi^2-eta^2)
    d_theta psi  = sqrt((xi^2-1)(1-eta^2)) (eta psi_xi - xi psi_eta)/(xi^2-eta^2)

The (xi, eta) integrals use xi = cosh(a), eta = cos(b) with Gauss-Legendre
nodes in a and b; this tames the mild singularity of the integrand at the nuclei.
Limits: D(R -> inf) = 1/4  (each atom: (1/4)<p^2> = 1/4), giving
E_ad(inf) = mu_e^2/(4M).
"""
from __future__ import annotations

import math

import numpy as np

from .electronic import Wavefunction

__all__ = ["D_of_R", "adiabatic_correction"]


def _nodes(p, na, nb):
    amax = math.acosh(1 + 45.0 / p)
    xa, wa = np.polynomial.legendre.leggauss(na)
    a = 0.5 * amax * (xa + 1)
    wa = 0.5 * amax * wa
    xb, wb = np.polynomial.legendre.leggauss(nb)
    b = 0.5 * math.pi * (xb + 1)
    wb = 0.5 * math.pi * wb
    xi, eta = np.cosh(a), np.cos(b)
    return xi, np.sinh(a) * wa, eta, np.sin(b) * wb


def _psi_on(wf, xi, eta, deriv=False):
    if not deriv:
        return wf.N * np.outer(wf.X(xi), wf.Y(eta))
    X, dX = wf.X(xi, deriv=True)
    Y, dY = wf.Y(eta, deriv=True)
    return wf.N * np.outer(X, Y), wf.N * np.outer(dX, Y), wf.N * np.outer(X, dY)


def D_of_R(R, na=160, nb=160, h=None, method="spheroidal", details=False):
    """D(R) = <grad_R psi|grad_R psi> (unit electron mass, 1s sigma_g).

    method='spheroidal': analytic chain rule (formula in module docstring).
    method='cartesian' : independent check, d_R psi at fixed Cartesian rho by
                         finite differences of psi evaluated at R +- h.
    """
    R = float(R)
    if h is None:
        h = 1e-3 * R
    wf = Wavefunction(R)
    xi, wx, eta, we = _nodes(wf.p, na, nb)
    XI, ET = np.meshgrid(xi, eta, indexing="ij")
    W = 2 * math.pi * (R / 2) ** 3 * np.outer(wx, we) * (XI ** 2 - ET ** 2)
    psi, pxi, peta = _psi_on(wf, xi, eta, deriv=True)
    den = XI ** 2 - ET ** 2

    wfs = {s: Wavefunction(R + s * h) for s in (-2, -1, 1, 2)}
    if method == "spheroidal":
        vals = {s: _psi_on(w, xi, eta) for s, w in wfs.items()}
        dR_fixed = (8 * (vals[1] - vals[-1]) - (vals[2] - vals[-2])) / (12 * h)
        F = dR_fixed - (XI * (XI ** 2 - 1) * pxi + ET * (1 - ET ** 2) * peta) / (R * den)
    elif method == "cartesian":
        z = 0.5 * R * XI * ET
        s2 = 0.25 * R * R * (XI ** 2 - 1) * (1 - ET ** 2)
        vals = {}
        for s, w in wfs.items():
            Rs = R + s * h
            r1 = np.sqrt((z + Rs / 2) ** 2 + s2)
            r2 = np.sqrt((z - Rs / 2) ** 2 + s2)
            xs = (r1 + r2) / Rs
            es = np.clip((r1 - r2) / Rs, -1, 1)
            vals[s] = w.N * w.X(xs) * w.Y(es)
        F = (8 * (vals[1] - vals[-1]) - (vals[2] - vals[-2])) / (12 * h)
    else:
        raise ValueError(method)
    G = np.sqrt((XI ** 2 - 1) * (1 - ET ** 2)) * (ET * pxi - XI * peta) / den
    radial = np.sum(W * F * F)
    angular = np.sum(W * G * G) / R ** 2
    if details:
        return dict(D=radial + angular, radial=radial, angular=angular,
                    norm=np.sum(W * psi * psi), overlap=np.sum(W * psi * F))
    return radial + angular


def adiabatic_correction(R, M, **kw):
    """E_ad(R) in hartree for two particles of mass M (electron mass 1)."""
    mu_e = 2 * M / (2 * M + 1)
    return mu_e ** 2 * D_of_R(mu_e * R, **kw) / M
