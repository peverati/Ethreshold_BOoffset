"""Model B: harmonic unlike pairs + inverse-square repulsion between the like pair.

    H = p1^2/2M + p2^2/2M + p3^2/2 + (k/2)(r13^2 + r23^2) + (k'/2) r12^2 + g / r12^2

(k' optional, default 0; K = k + 2k'). In Jacobi coordinates
    H_int = [-(1/M) Lap_R + (K/4) R^2 + g/R^2] + [-(2M+1)/(4M) Lap_rho + k rho^2],
again EXACTLY separable in R | rho for any like-pair potential V(R).
The R part with angular momentum J: u = R f,
    -(1/M) u'' + [(J(J+1) + g M)/(M R^2) + K R^2/4] u = E u,
i.e. a 3D oscillator with non-integer "angular momentum" lam_J,
    lam_J (lam_J + 1) = J(J+1) + g M,   lam_J = (-1 + sqrt((2J+1)^2 + 4 g M))/2,
    E_{nJ}^R = omega_R (2n + lam_J + 3/2),   f_0(R) = N R^lam exp(-a R^2/2),
with omega_R = sqrt(K/M), a = sqrt(M K)/2 (as in Model A). lam = lam_0.
Ground state: psi_0 = f_{lam,a}(R) phi_b(rho), total L = 0, closed form.
Requires g > -1/(4M) (lam real; we use g >= 0) and k > 0, K > 0.

g = 0 reproduces Model A.
"""
import numpy as np
import mpmath
from scipy.special import gammaln, polygamma
from scipy.optimize import brentq, minimize_scalar
from . import common as C

SQ = np.sqrt


def params(M, g=1.0, k=1.0, kp=0.0):
    if not (k > 0 and k + 2 * kp > 0):
        raise ValueError("unbound: need k > 0 and k + 2k' > 0")
    if g < 0:
        raise ValueError("this module assumes g >= 0")
    M = np.asarray(M, float)
    K = k + 2 * kp
    u = g * M
    lam = 0.5 * (SQ(1 + 4 * u) - 1)
    return dict(M=M, g=g, k=k, K=K, lam=lam, kk=lam + 1.5,
                omR=SQ(K / M), omrho=SQ(k * (2 * M + 1) / M),
                a=SQ(M * K) / 2, b=2 * SQ(k * M / (2 * M + 1)),
                lam_dot=u / SQ(1 + 4 * u),                  # d lam / d ln M
                lam_ddot=u * (1 + 2 * u) / (1 + 4 * u) ** 1.5)


def bond_length(g=1.0, k=1.0, kp=0.0):
    """Classical minimum of g/R^2 + K R^2/4: R_e = (4g/K)^(1/4) (M independent)."""
    return (4 * g / (k + 2 * kp)) ** 0.25


# ------------------------------------------------------------------ energy
def energy(M, g=1.0, k=1.0, kp=0.0):
    """E_0 = omega_R (lam + 3/2) + (3/2) omega_rho."""
    p = params(M, g, k, kp)
    return p["omR"] * p["kk"] + 1.5 * p["omrho"]


def dE_dlnM(M, g=1.0, k=1.0, kp=0.0):
    p = params(M, g, k, kp)
    M = p["M"]
    return (-0.5 * p["omR"] * p["kk"] + p["omR"] * p["lam_dot"]
            - 0.75 * p["omrho"] / (2 * M + 1))


def d2E_dlnM2(M, g=1.0, k=1.0, kp=0.0):
    p = params(M, g, k, kp)
    M = p["M"]
    return (0.25 * p["omR"] * p["kk"] - p["omR"] * p["lam_dot"]
            + p["omR"] * p["lam_ddot"]
            + 0.375 * p["omrho"] * (4 * M + 1) / (2 * M + 1) ** 2)


def expectation_R(M, g=1.0, k=1.0, kp=0.0):
    """<R^2> = kk/a, <1/R^2> = a/(lam+1/2), <T_R> = E_R - K<R^2>/4 - g<1/R^2>."""
    p = params(M, g, k, kp)
    R2 = p["kk"] / p["a"]
    Rm2 = p["a"] / (p["lam"] + 0.5)
    TR = p["omR"] * p["kk"] - p["K"] * R2 / 4 - g * Rm2
    return dict(R2=R2, Rm2=Rm2, T_R=TR, T_rho=0.75 * p["omrho"])


def hellmann_feynman(M, g=1.0, k=1.0, kp=0.0):
    e = expectation_R(M, g, k, kp)
    return -e["T_R"] - e["T_rho"] / (2 * np.asarray(M, float) + 1)


def levels(M, n=0, J=0, g=1.0, k=1.0, kp=0.0):
    """Exact E_{n,J} (rho in its ground state, total L = J)."""
    p = params(M, g, k, kp)
    lamJ = 0.5 * (SQ((2 * J + 1) ** 2 + 4 * g * p["M"]) - 1)
    return p["omR"] * (2 * n + lamJ + 1.5) + 1.5 * p["omrho"]


# ---------------------------------------------------- distributions, moments
def pair_moments(M, g=1.0, k=1.0, kp=0.0):
    """<R^p> = Gamma(kk + p/2)/(Gamma(kk) a^(p/2)); eta = Delta R/<R> depends on gM only:
    eta^2 = kk Gamma(kk)^2/Gamma(kk+1/2)^2 - 1."""
    p = params(M, g, k, kp)
    R1 = C.chi_moment(p["lam"], p["a"], 1)
    R2 = C.chi_moment(p["lam"], p["a"], 2)
    dR = SQ(R2 - R1 ** 2)
    kk = p["kk"]
    eta = SQ(kk * np.exp(2 * (gammaln(kk) - gammaln(kk + 0.5))) - 1)
    return dict(R=R1, R2=R2, dR=dR, eta=eta,
                R_peak=SQ((p["lam"] + 1) / p["a"]))


def pair_distribution(R, M, g=1.0, k=1.0, kp=0.0):
    """4 pi R^2 f(R)^2 ~ R^(2 lam + 2) exp(-a R^2); peak at sqrt((lam+1)/a)."""
    p = params(M, g, k, kp)
    return 4 * np.pi * R ** 2 * C.chi_value(float(p["lam"]), float(p["a"]), R) ** 2


def heavy_density_mid(y, M, g=1.0, k=1.0, kp=0.0):
    """Density of particle 1 about the pair midpoint (y = -R/2): 8 f(2y)^2.
    Vanishes at y = 0 for every M > 0 when g > 0 (hole at all masses)."""
    p = params(M, g, k, kp)
    return 8 * C.chi_value(float(p["lam"]), float(p["a"]), 2 * np.asarray(y, float)) ** 2


def _cm_pars(p):
    M = p["M"]
    s2 = 1 / (2 * p["b"] * (2 * M + 1) ** 2)   # variance/component of rho/(2M+1)
    P = 4 * p["a"] + 1 / (2 * s2)
    return s2, P


def heavy_density_cm(x, M, g=1.0, k=1.0, kp=0.0):
    """Density of particle 1 about the centre of mass, x = -R/2 - rho/(2M+1):

    n(x) = A exp(-x^2/(2 s2)) 1F1(kk; 3/2; x^2/(4 P s2^2)),
    s2 = 1/(2 b (2M+1)^2), P = 4a + 1/(2 s2),
    A = (2 pi s2)^(-3/2) * 8 * 4pi N^2 2^(2lam) Gamma(kk)/(2 P^kk),  N^2 = 2 a^kk/(4 pi Gamma(kk)).
    """
    p = params(M, g, k, kp)
    s2, P = _cm_pars(p)
    kk, lam, a = p["kk"], p["lam"], p["a"]
    lnA = (-1.5 * np.log(2 * np.pi * s2) + np.log(8.0) + np.log(2.0) + kk * np.log(a)
           + 2 * lam * np.log(2.0) - kk * np.log(P) - np.log(2.0))
    x = np.asarray(x, float)
    z = x ** 2 / (4 * P * s2 ** 2)
    kkf, lnAf, s2f = float(kk), float(lnA), float(s2)
    # log-space evaluation with mpmath (1F1 overflows in double precision at large M)
    out = [float(mpmath.exp(lnAf - xx ** 2 / (2 * s2f) + mpmath.log(mpmath.hyp1f1(kkf, 1.5, zz))))
           for xx, zz in zip(np.atleast_1d(x).ravel(), np.atleast_1d(z).ravel())]
    return np.array(out).reshape(x.shape) if x.ndim else out[0]


def cm_curvature_coefficient(M, g=1.0, k=1.0, kp=0.0):
    """alpha/n(0) where n_cm(x) = n(0)[1 + (alpha/n0) x^2 + ...]:
    alpha/n0 = kk/(3/2) /(4 P s2^2) - 1/(2 s2).  Hole (alpha > 0) iff lam > 12 a s2,
    i.e. lam(gM) > (3/2) sqrt(K/k) (2M+1)^(-3/2)."""
    p = params(M, g, k, kp)
    s2, P = _cm_pars(p)
    return p["kk"] / 1.5 / (4 * P * s2 ** 2) - 1 / (2 * s2)


def cm_curvature_dimensionless(M, g=1.0, k=1.0, kp=0.0):
    """2 s2 alpha/n0 = (lam - 12 a s2)/(12 a s2 + 3/2): curvature at x = 0 in units of the
    width of the rho/(2M+1) smearing. Positive = hole."""
    p = params(M, g, k, kp)
    s2, P = _cm_pars(p)
    c = 12 * p["a"] * s2
    return (p["lam"] - c) / (c + 1.5)


def critical_mass_cm(g=1.0, k=1.0, kp=0.0):
    """Mass M_c at which the CM-referenced heavy density changes curvature at x = 0.
    Root of lam(gM) = (3/2) sqrt(K/k) (2M+1)^(-3/2); unique (LHS increasing, RHS
    decreasing). Returns nan for g = 0 (never a hole)."""
    if g <= 0:
        return np.nan
    K = k + 2 * kp
    f = lambda M: 0.5 * (np.sqrt(1 + 4 * g * M) - 1) - 1.5 * np.sqrt(K / k) * (2 * M + 1) ** -1.5
    hi = 1.0
    while f(hi) < 0:
        hi *= 10
    return brentq(f, 1e-14, hi, xtol=1e-15, rtol=1e-14)


def cm_shell_radius(M, g=1.0, k=1.0, kp=0.0):
    """Radius of the maximum of n_cm(x) (0 if the maximum is at the origin)."""
    if cm_curvature_coefficient(M, g, k, kp) <= 0:
        return 0.0
    p = params(M, g, k, kp)
    xmax = 2 * float(SQ(p["kk"] / p["a"]))       # = 2 <R^2>^(1/2) > R_peak/2
    xs = np.linspace(0, xmax, 400)
    ns = heavy_density_cm(xs, M, g, k, kp)
    i = int(np.argmax(ns))
    lo, hi = xs[max(i - 1, 0)], xs[min(i + 1, len(xs) - 1)]
    r = minimize_scalar(lambda x: -heavy_density_cm(x, M, g, k, kp),
                        bounds=(lo, hi), method="bounded", options=dict(xatol=1e-10))
    return r.x


# ---------------------------------------------------------------- fidelity
def overlap(M1, M2, g=1.0, k=1.0, kp=0.0, s=1.0, sR=None, srho=None):
    """<psi(M1)| U psi(M2)> in fixed Jacobi coordinates; U = global dilation by s,
    or separate dilations (sR on R, srho on rho) if given."""
    p1, p2 = params(M1, g, k, kp), params(M2, g, k, kp)
    sR = s if sR is None else sR
    srho = s if srho is None else srho
    return (C.chi_overlap(p1["lam"], p1["a"], p2["lam"], p2["a"] * sR ** 2)
            * C.gauss3d_overlap(p1["b"], p2["b"] * srho ** 2))


def fidelity(M1, M2, g=1.0, k=1.0, kp=0.0):
    return overlap(M1, M2, g, k, kp)


def fidelity_inv(M1, M2, g=1.0, k=1.0, kp=0.0):
    """max over a global dilation s (1D numerical maximization of the closed-form overlap)."""
    r = minimize_scalar(lambda t: -overlap(M1, M2, g, k, kp, s=np.exp(t)),
                        bounds=(-20, 20), method="bounded", options=dict(xatol=1e-12))
    return -r.fun


def fidelity_sep(M1, M2, g=1.0, k=1.0, kp=0.0):
    """max over separate dilations of R and rho (rho factor then equals 1)."""
    p1, p2 = params(M1, g, k, kp), params(M2, g, k, kp)
    r = minimize_scalar(lambda t: -C.chi_overlap(p1["lam"], p1["a"], p2["lam"], p2["a"] * np.exp(t)),
                        bounds=(-40, 40), method="bounded", options=dict(xatol=1e-12))
    return -r.fun


def chi_parts(M, g=1.0, k=1.0, kp=0.0):
    """Fidelity susceptibilities in ln M (closed form).

    With kk = lam + 3/2, ld = dlam/dlnM, td = dln a/dlnM = 1/2, xb = dln b/dlnM = 1/(2(2M+1)):
      chi_R     = (1/4)[psi1(kk) ld^2 - 2 ld td + kk td^2]
      chi_rho   = (3/8) xb^2
      chi_fixed = chi_R + chi_rho
      <d psi|D psi> = (1/2)(kk td - ld) + (3/4) xb,   ||D psi||^2 = kk + 3/2
      chi_inv   = chi_fixed - <d psi|D psi>^2 / ||D psi||^2          (global dilation)
      chi_sep   = (1/4) ld^2 [psi1(kk) - 1/kk]                       (separate dilations)
    """
    p = params(M, g, k, kp)
    M = p["M"]
    kk, ld, td = p["kk"], p["lam_dot"], 0.5
    xb = 0.5 / (2 * M + 1)
    t1 = polygamma(1, kk)
    chiR = 0.25 * (t1 * ld ** 2 - 2 * ld * td + kk * td ** 2)
    chirho = 0.375 * xb ** 2
    proj = 0.5 * (kk * td - ld) + 0.75 * xb
    normD2 = kk + 1.5
    chi_fixed = chiR + chirho
    chi_inv = chi_fixed - proj ** 2 / normD2
    chi_sep = 0.25 * ld ** 2 * (t1 - 1 / kk)
    return dict(chi_R=chiR, chi_rho=chirho, chi_fixed=chi_fixed, chi_inv=chi_inv,
                chi_removed=proj ** 2 / normD2, chi_sep=chi_sep)


# ------------------------------------------------------------ entanglement
def entanglement_R_s(M, g=1.0, k=1.0, kp=0.0, n_r=160):
    """Bipartition R | s with s = r3 - r1 (light particle referred to heavy 1).

    psi(R, s) = f(R) phi_b(s - R/2)  =>  gamma_R(R, R') = f(R) f(R') exp(-c|R-R'|^2),
    c = b/16. Spectrum by partial waves (numerical radial quadrature; see
    common.partial_wave_gaussian_kernel). The R | rho bipartition is an exact product.
    """
    p = params(M, g, k, kp)
    lam, a, b = float(p["lam"]), float(p["a"]), float(p["b"])
    c = b / 16
    rmax = float(SQ(p["kk"] / a)) + 12 / np.sqrt(a)
    return C.partial_wave_gaussian_kernel(lambda r: C.chi_value(lam, a, r), c, rmax, n_r=n_r)


def entanglement_R_s_large_M(g=1.0, k=1.0, kp=0.0, lmax=200):
    """M -> infinity limit of R|s: thin shell at R_e; eigenvalues p_l = e^-x i_l(x)
    (multiplicity 2l+1), x = 2 c R_e^2 = b_inf R_e^2/8, b_inf = sqrt(2k)."""
    from scipy.special import ive
    x = np.sqrt(2 * k) * bond_length(g, k, kp) ** 2 / 8
    l = np.arange(lmax + 1)
    pl = np.sqrt(np.pi / (2 * x)) * ive(l + 0.5, x)
    m = pl > 0
    S = -np.sum((2 * l[m] + 1) * pl[m] * np.log(pl[m]))
    P = np.sum((2 * l + 1) * pl ** 2)
    return dict(x=x, S_vN=S, S_L=1 - P, trace=np.sum((2 * l + 1) * pl))


def summary(M, g=1.0, k=1.0, kp=0.0):
    M = float(M)
    pm = pair_moments(M, g, k, kp)
    ch = chi_parts(M, g, k, kp)
    ent = entanglement_R_s(M, g, k, kp)
    p = params(M, g, k, kp)
    return dict(M=M, lam=float(p["lam"]), E=float(energy(M, g, k, kp)),
                dE_dlnM=float(dE_dlnM(M, g, k, kp)), d2E_dlnM2=float(d2E_dlnM2(M, g, k, kp)),
                R_mean=float(pm["R"]), dR=float(pm["dR"]), eta=float(pm["eta"]),
                R_peak=float(pm["R_peak"]),
                cm_curv=float(cm_curvature_dimensionless(M, g, k, kp)),
                cm_shell=float(cm_shell_radius(M, g, k, kp)),
                chi_fixed=float(ch["chi_fixed"]), chi_inv=float(ch["chi_inv"]),
                chi_sep=float(ch["chi_sep"]),
                SL_R_rho=0.0, SL_R_s=float(ent["S_L"]), SvN_R_s=float(ent["S_vN"]),
                E10_minus_E00=float(levels(M, 0, 1, g, k, kp) - levels(M, 0, 0, g, k, kp)),
                B_from_invR2=float(expectation_R(M, g, k, kp)["Rm2"] / M))
