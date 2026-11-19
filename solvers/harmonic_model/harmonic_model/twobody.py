"""Two-body references: masses M and 1, reduced mass mu = M/(M+1), relative coordinate r.

Harmonic, V = (k/2) r^2:  psi = phi_a(r), a = mu omega = sqrt(k mu), E = (3/2) sqrt(k/mu).
Coulomb,  V = -Z/r:       psi = (beta^3/pi)^(1/2) exp(-beta r), beta = Z mu, E = -Z^2 mu/2.

dln mu/dln M = 1/(M+1).

Both ground states depend on M only through ONE length scale, so a global dilation
maps psi(M) exactly onto psi(M'): F_inv = 1 and chi_inv = 0 identically.
"""
import numpy as np
from . import common as C


def mu(M):
    M = np.asarray(M, float)
    return M / (M + 1)


# ---------------------------------------------------------------- harmonic
def harmonic_energy(M, k=1.0):
    return 1.5 * np.sqrt(k / mu(M))


def harmonic_fidelity(M1, M2, k=1.0):
    """(2 sqrt(a1 a2)/(a1+a2))^(3/2) = sech(ln(mu2/mu1)/4)^(3/2)."""
    return C.gauss3d_overlap(np.sqrt(k * mu(M1)), np.sqrt(k * mu(M2)))


def harmonic_chi_fixed(M):
    """(3/8)(dln a/dlnM)^2 = 3/(32 (M+1)^2)."""
    M = np.asarray(M, float)
    return 3 / (32 * (M + 1) ** 2)


# ----------------------------------------------------------------- Coulomb
def coulomb_energy(M, Z=1.0):
    return -0.5 * Z ** 2 * mu(M)


def coulomb_overlap(b1, b2):
    """<e^{-b1 r}|e^{-b2 r}> normalized = 8 (b1 b2)^(3/2)/(b1+b2)^3 = sech(ln(b2/b1)/2)^3."""
    return 8 * (b1 * b2) ** 1.5 / (b1 + b2) ** 3


def coulomb_fidelity(M1, M2, Z=1.0):
    return coulomb_overlap(Z * mu(M1), Z * mu(M2))


def coulomb_chi_fixed(M):
    """(3/4)(dln mu/dlnM)^2 = (3/4)/(M+1)^2."""
    M = np.asarray(M, float)
    return 0.75 / (M + 1) ** 2


def fidelity_inv(M1, M2):
    """Global dilation by s = mu1/mu2 maps psi(M2) onto psi(M1) exactly: 1."""
    return np.ones_like(np.asarray(M1, float) * np.asarray(M2, float))


def chi_inv(M):
    return 0 * np.asarray(M, float)
