"""Shared closed-form building blocks.

Conventions
-----------
* 3D Gaussian  phi_b(x) = (b/pi)^(3/4) exp(-b x^2 / 2)            (normalized)
* "chi" function f_{lam,a}(R) = N R^lam exp(-a R^2/2), normalized in 3D,
  4 pi N^2 Gamma(lam+3/2) / (2 a^(lam+3/2)) = 1.  lam = 0 is phi_a.
* Fidelity susceptibility in a parameter t (here t = ln M):
  F(t, t+d) = 1 - chi_F d^2 / 2 + O(d^3),  chi_F = <d psi|d psi> - |<psi|d psi>|^2.
* Global dilation of n internal coordinates: (U_s psi)(x) = s^(n/2) psi(s x).
  It maps phi_b -> phi_{b s^2} and f_{lam,a} -> f_{lam, a s^2}.
"""
import numpy as np
from scipy.special import gammaln, polygamma, ive


# ---------------------------------------------------------------- overlaps
def gauss3d_overlap(b1, b2):
    """<phi_b1|phi_b2> = (2 sqrt(b1 b2)/(b1+b2))^(3/2) = sech(ln(b2/b1)/2)^(3/2)."""
    b1, b2 = np.asarray(b1, float), np.asarray(b2, float)
    return (2.0 * np.sqrt(b1 * b2) / (b1 + b2)) ** 1.5


def chi_overlap(lam1, a1, lam2, a2):
    """<f_{lam1,a1}|f_{lam2,a2}> (both real, normalized, 3D).

    = Gamma(km) ((a1+a2)/2)^(-km) / sqrt(Gamma(k1) a1^-k1 Gamma(k2) a2^-k2),
    k_i = lam_i + 3/2, km = (k1+k2)/2.
    """
    k1, k2 = np.asarray(lam1, float) + 1.5, np.asarray(lam2, float) + 1.5
    km = 0.5 * (k1 + k2)
    ln = (0.5 * (k1 * np.log(a1) + k2 * np.log(a2))
          - 0.5 * (gammaln(k1) + gammaln(k2))
          + gammaln(km) - km * np.log(0.5 * (np.asarray(a1) + np.asarray(a2))))
    return np.exp(ln)


def chi_moment(lam, a, p):
    """<R^p> for |f_{lam,a}|^2 d^3R:  Gamma(k + p/2) / (Gamma(k) a^(p/2)), k = lam+3/2."""
    k = np.asarray(lam, float) + 1.5
    return np.exp(gammaln(k + 0.5 * p) - gammaln(k)) / np.asarray(a, float) ** (0.5 * p)


def chi_value(lam, a, R):
    """f_{lam,a}(R) (normalized)."""
    k = lam + 1.5
    lnN2 = np.log(2.0) + k * np.log(a) - gammaln(k) - np.log(4 * np.pi)
    R = np.asarray(R, float)
    with np.errstate(divide="ignore"):
        lnR = np.where(R > 0, np.log(np.where(R > 0, R, 1.0)), -np.inf)
    out = np.exp(0.5 * lnN2 + lam * lnR - 0.5 * a * R ** 2)
    if lam == 0:
        out = np.exp(0.5 * lnN2 - 0.5 * a * R ** 2)
    return out


# ------------------------------------------------- fidelity metric of f_{lam,a}
def chi_fisher(lam):
    """Fisher matrix of |f|^2 in (lam, ln a): [[psi1(k), -1], [-1, k]], k = lam+3/2.

    For a real positive state, <d psi|d psi> = (1/4) Fisher, so
    chi_F = (1/4) (lam', t') I (lam', t')^T with t = ln a.
    """
    k = np.asarray(lam, float) + 1.5
    return polygamma(1, k), -np.ones_like(k), k


# ------------------------------------------------------------ entanglement
def entropy_from_nu(nu):
    """von Neumann entropy (nats) of a one-mode Gaussian with symplectic eigenvalue nu >= 1/2."""
    nu = np.asarray(nu, float)
    p, m = nu + 0.5, nu - 0.5
    with np.errstate(divide="ignore", invalid="ignore"):
        t2 = np.where(m > 1e-300, m * np.log(np.where(m > 1e-300, m, 1.0)), 0.0)
    return p * np.log(p) - t2


def gaussian_bipartition(A11, A12, A22, ncomp=3):
    """Entanglement of psi ~ exp(-(A11 u^2 + 2 A12 u v + A22 v^2)/2) per Cartesian
    component, for ncomp identical components (3D isotropic: ncomp = 3).

    Single-component purity P1 = sqrt(1 - A12^2/(A11 A22)); nu = 1/(2 P1).
    Returns dict(purity, S_L = 1 - P1^ncomp, S_vN = ncomp * S(nu), nu).
    """
    A11, A12, A22 = (np.asarray(x, float) for x in (A11, A12, A22))
    P1 = np.sqrt(1.0 - A12 ** 2 / (A11 * A22))
    nu = 0.5 / P1
    return dict(purity=P1 ** ncomp, S_L=1.0 - P1 ** ncomp,
                S_vN=ncomp * entropy_from_nu(nu), nu=nu)


def partial_wave_gaussian_kernel(f, c, rmax, n_r=160, lmax=None, tol=1e-14):
    """Spectrum of the rotationally invariant reduced density matrix

        gamma(u, u') = f(|u|) f(|u'|) exp(-c |u - u'|^2)      on L^2(R^3),

    with f radial and normalized (int f^2 d^3u = 1). Uses
    exp(2c u.u') = sum_l (2l+1) i_l(2c u u') P_l(cos), so each (l, m) block is
    the radial kernel 4 pi f f' exp(-c(u^2+u'^2)) i_l(2 c u u') on u^2 du.
    Eigenvalues of block l have multiplicity 2l+1.

    Returns dict(trace, purity, S_vN, S_L, lmax_used, eigs=list of (l, array)).
    """
    x, w = np.polynomial.legendre.leggauss(n_r)
    r = 0.5 * rmax * (x + 1.0)
    w = 0.5 * rmax * w
    fr = f(r)
    sw = np.sqrt(w) * r * fr
    U, V = np.meshgrid(r, r, indexing="ij")
    z = 2.0 * c * U * V
    pref = 4.0 * np.pi * np.exp(-c * (U - V) ** 2)   # e^{-c(u^2+u'^2)} e^{z} rescaled
    trace = 0.0
    pur = 0.0
    S = 0.0
    eigs = []
    l = 0
    while True:
        # i_l(z) e^{-z} = sqrt(pi/(2z)) I_{l+1/2}(z) e^{-z}
        il = np.sqrt(np.pi / (2.0 * z)) * ive(l + 0.5, z)
        Kmat = sw[:, None] * pref * il * sw[None, :]
        ev = np.linalg.eigvalsh(Kmat)
        ev = ev[ev > 1e-300]
        eigs.append((l, ev[::-1]))
        tr_l = (2 * l + 1) * ev.sum()
        trace += tr_l
        pur += (2 * l + 1) * np.sum(ev ** 2)
        S -= (2 * l + 1) * np.sum(ev * np.log(ev))
        l += 1
        if lmax is not None and l > lmax:
            break
        if lmax is None and (tr_l < tol or l > 400):
            break
    return dict(trace=trace, purity=pur, S_L=1.0 - pur, S_vN=S,
                lmax_used=l - 1, eigs=eigs)
