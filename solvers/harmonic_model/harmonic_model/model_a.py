"""Model A: all-harmonic three-body system (Moshinsky type).

    H = p1^2/2M + p2^2/2M + p3^2/2 + (k/2)(r13^2 + r23^2) + (k'/2) r12^2

Jacobi coordinates R = r2 - r1, rho = r3 - (r1+r2)/2 (mass independent):
    r13^2 + r23^2 = 2 rho^2 + R^2/2,
    H_int = -(1/M) Lap_R - (2M+1)/(4M) Lap_rho + (K/4) R^2 + k rho^2,   K = k + 2k'.
It separates EXACTLY into two isotropic oscillators:
    mu_R = M/2,          omega_R = sqrt(K/M),          a = mu_R omega_R = sqrt(M K)/2
    mu_rho = 2M/(2M+1),  omega_rho = sqrt(k(2M+1)/M),  b = mu_rho omega_rho = 2 sqrt(k M/(2M+1))
    psi_0 = phi_a(R) phi_b(rho),  E_0 = (3/2)(omega_R + omega_rho).
Bound (normalizable ground state) iff k > 0 and K = k + 2k' > 0, i.e. k' > -k/2.

All functions accept numpy arrays for M (broadcasting).
"""
import numpy as np
from . import common as C

SQ = np.sqrt


def check_params(k, kp):
    if not (k > 0 and k + 2 * kp > 0):
        raise ValueError("unbound: need k > 0 and k + 2k' > 0 (k' > -k/2)")


def params(M, k=1.0, kp=0.0):
    check_params(k, kp)
    M = np.asarray(M, float)
    K = k + 2 * kp
    return dict(M=M, K=K, k=k, muR=M / 2, murho=2 * M / (2 * M + 1),
                omR=SQ(K / M), omrho=SQ(k * (2 * M + 1) / M),
                a=SQ(M * K) / 2, b=2 * SQ(k * M / (2 * M + 1)))


# ------------------------------------------------------------------ energy
def energy(M, k=1.0, kp=0.0):
    p = params(M, k, kp)
    return 1.5 * (p["omR"] + p["omrho"])


def dE_dlnM(M, k=1.0, kp=0.0):
    """dE/dlnM = -(3/4) omega_R - (3/4) omega_rho/(2M+1)."""
    p = params(M, k, kp)
    return -0.75 * p["omR"] - 0.75 * p["omrho"] / (2 * M + 1)


def d2E_dlnM2(M, k=1.0, kp=0.0):
    """d2E/dlnM^2 = (3/8) omega_R + (3/8) omega_rho (4M+1)/(2M+1)^2."""
    p = params(M, k, kp)
    M = p["M"]
    return 0.375 * p["omR"] + 0.375 * p["omrho"] * (4 * M + 1) / (2 * M + 1) ** 2


def kinetic(M, k=1.0, kp=0.0):
    """<T_R> = 3 omega_R/4, <T_rho> = 3 omega_rho/4 (virial)."""
    p = params(M, k, kp)
    return 0.75 * p["omR"], 0.75 * p["omrho"]


def hellmann_feynman(M, k=1.0, kp=0.0):
    """-<T_R> - <T_rho>/(2M+1): must equal dE_dlnM."""
    TR, Tr = kinetic(M, k, kp)
    return -TR - Tr / (2 * np.asarray(M, float) + 1)


# ---------------------------------------------------- distributions, moments
def pair_moments(M, k=1.0, kp=0.0):
    """Like-pair distance: Maxwell distribution P(R) = 4 pi R^2 phi_a(R)^2.

    <R> = 2/sqrt(pi a), <R^2> = 3/(2a), Delta R/<R> = sqrt(3 pi/8 - 1) (M independent).
    """
    a = params(M, k, kp)["a"]
    R1 = 2 / SQ(np.pi * a)
    R2 = 1.5 / a
    dR = SQ(R2 - R1 ** 2)
    return dict(R=R1, R2=R2, dR=dR, eta=dR / R1)


LINDEMANN = np.sqrt(3 * np.pi / 8 - 1)


def pair_distribution(R, M, k=1.0, kp=0.0):
    """Radial pair distribution 4 pi R^2 |phi_a(R)|^2 (normalized in dR); peak at 1/sqrt(a)."""
    a = params(M, k, kp)["a"]
    return 4 * np.pi * R ** 2 * C.chi_value(0.0, a, R) ** 2


def unlike_pair_sigma2(M, k=1.0, kp=0.0):
    """r13 = rho + R/2 is Gaussian; variance per Cartesian component 1/(2b) + 1/(8a)."""
    p = params(M, k, kp)
    return 1 / (2 * p["b"]) + 1 / (8 * p["a"])


def heavy_density_sigma2(M, k=1.0, kp=0.0, ref="cm"):
    """One-particle density of particle 1 (or 2) is an isotropic Gaussian.

    About the centre of mass: r1 - R_cm = -R/2 - rho/(2M+1):
        sigma^2 = 1/(8a) + 1/(2 b (2M+1)^2)     per Cartesian component.
    About the pair midpoint: r1 - (r1+r2)/2 = -R/2:  sigma^2 = 1/(8a).
    Always maximal at the origin: no hole for any M, k, k'.
    """
    p = params(M, k, kp)
    if ref == "cm":
        return 1 / (8 * p["a"]) + 1 / (2 * p["b"] * (2 * p["M"] + 1) ** 2)
    if ref == "mid":
        return 1 / (8 * p["a"])
    raise ValueError(ref)


def heavy_density(x, M, k=1.0, kp=0.0, ref="cm"):
    s2 = heavy_density_sigma2(M, k, kp, ref)
    return (2 * np.pi * s2) ** -1.5 * np.exp(-x ** 2 / (2 * s2))


# ---------------------------------------------------------------- fidelity
def log_rates(M):
    """x_a = dln a/dlnM = 1/2, x_b = dln b/dlnM = 1/(2(2M+1)) (independent of k, k')."""
    M = np.asarray(M, float)
    return 0.5 + 0 * M, 0.5 / (2 * M + 1)


def fidelity(M1, M2, k=1.0, kp=0.0):
    """|<psi(M1)|psi(M2)>| in fixed Jacobi coordinates
    = sech(ea/2)^(3/2) sech(eb/2)^(3/2), ea = ln(a2/a1), eb = ln(b2/b1)."""
    p1, p2 = params(M1, k, kp), params(M2, k, kp)
    return C.gauss3d_overlap(p1["a"], p2["a"]) * C.gauss3d_overlap(p1["b"], p2["b"])


def fidelity_dilated(M1, M2, s, k=1.0, kp=0.0):
    """|<psi(M1)|U_s psi(M2)>|, U_s the global dilation of (R, rho)."""
    p1, p2 = params(M1, k, kp), params(M2, k, kp)
    return (C.gauss3d_overlap(p1["a"], p2["a"] * s ** 2)
            * C.gauss3d_overlap(p1["b"], p2["b"] * s ** 2))


def fidelity_inv(M1, M2, k=1.0, kp=0.0):
    """max_s |<psi(M1)|U_s psi(M2)>| = sech(Delta/4)^3,  Delta = ln(a2/a1) - ln(b2/b1).
    (optimum at ln s^2 = -(ea+eb)/2; log sech is concave and even)."""
    p1, p2 = params(M1, k, kp), params(M2, k, kp)
    D = np.log(p2["a"] / p1["a"]) - np.log(p2["b"] / p1["b"])
    return 1 / np.cosh(D / 4) ** 3


def fidelity_sep(M1, M2, k=1.0, kp=0.0):
    """Maximum over SEPARATE dilations of R and rho: identically 1."""
    return np.ones_like(np.asarray(M1, float) * np.asarray(M2, float))


def chi_fixed(M):
    """chi_F(ln M) = (3/8)(x_a^2 + x_b^2) = 3/32 + 3/(32 (2M+1)^2)  (k, k' drop out)."""
    M = np.asarray(M, float)
    return 3 / 32 + 3 / (32 * (2 * M + 1) ** 2)


def chi_inv(M):
    """Dilation-invariant: (3/16)(x_a - x_b)^2 = (3/16) M^2/(2M+1)^2."""
    M = np.asarray(M, float)
    return 3 / 16 * M ** 2 / (2 * M + 1) ** 2


def chi_removed(M):
    """Part removed by the global dilation: (3/16)(x_a + x_b)^2 = (3/16)((M+1)/(2M+1))^2."""
    M = np.asarray(M, float)
    return 3 / 16 * ((M + 1) / (2 * M + 1)) ** 2


# ------------------------------------------------------------ entanglement
def entanglement(M, k=1.0, kp=0.0, which="R|rho"):
    """Entanglement of the internal ground state for several bipartitions.

    'R|rho'  : Jacobi bipartition (= particle 3 vs pair internal motion): EXACT PRODUCT, 0.
    'R|s'    : s = r3 - r1 (light particle referred to heavy particle 1, the
               ATOM-MOL-nonBO coordinates): A = [[a + b/4, -b/2], [-b/2, b]],
               P1 = sqrt(4a/(4a+b)).
    'x1|x2'  : x_i = r_i - R_cm (heavy-heavy, CM frame; x3 = -M(x1+x2)):
               R = x2-x1, rho = -(2M+1)(x1+x2)/2, B = b(2M+1)^2/4,
               A = [[a+B, B-a], [B-a, a+B]], P1 = 2 sqrt(aB)/(a+B).
    'x1|x3'  : heavy 1 vs light particle, CM frame: R = -2 x1 - x3/M,
               rho = (2M+1) x3/(2M): A = [[4a, 2a/M], [2a/M, a/M^2 + B/M^2]],
               P1 = sqrt(B/(a+B)).
    """
    p = params(M, k, kp)
    a, b, M = p["a"], p["b"], p["M"]
    B = b * (2 * M + 1) ** 2 / 4
    z = 0 * M
    if which == "R|rho":
        return C.gaussian_bipartition(a, z, b)
    if which == "R|s":
        return C.gaussian_bipartition(a + b / 4, -b / 2, b)
    if which == "x1|x2":
        return C.gaussian_bipartition(a + B, B - a, a + B)
    if which == "x1|x3":
        return C.gaussian_bipartition(4 * a, 2 * a / M, (a + B) / M ** 2)
    raise ValueError(which)


def summary(M, k=1.0, kp=0.0):
    """Dict of all scalar quantities at mass(es) M."""
    M = np.asarray(M, float)
    pm = pair_moments(M, k, kp)
    out = dict(M=M, E=energy(M, k, kp), dE_dlnM=dE_dlnM(M, k, kp),
               d2E_dlnM2=d2E_dlnM2(M, k, kp), R_mean=pm["R"], R2_mean=pm["R2"],
               dR=pm["dR"], eta=pm["eta"],
               sig2_cm=heavy_density_sigma2(M, k, kp, "cm"),
               sig2_mid=heavy_density_sigma2(M, k, kp, "mid"),
               chi_fixed=chi_fixed(M), chi_inv=chi_inv(M), chi_sep=0 * M)
    for w in ("R|rho", "R|s", "x1|x2", "x1|x3"):
        e = entanglement(M, k, kp, w)
        tag = w.replace("|", "_")
        out["SL_" + tag] = e["S_L"]
        out["SvN_" + tag] = e["S_vN"]
    return out
