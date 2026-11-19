"""Independent numerics used ONLY to validate the closed forms.

Nothing here uses the closed-form wave functions: radial equations are solved
by finite differences directly from the potentials, derivatives with respect
to ln M and to the dilation parameter are taken by central differences of
eigenvectors on a fixed grid, and integrals are done by quadrature.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal


class RadialGrid:
    """Uniform grid r_i = i h, i = 1..n, h = rmax/(n+1); Dirichlet u(0) = u(rmax) = 0."""

    def __init__(self, rmax, n=4000):
        self.h = rmax / (n + 1)
        self.r = self.h * np.arange(1, n + 1)

    def ground(self, ckin, V, nstate=0):
        """Lowest eigenpair of -ckin u'' + V(r) u (second-order FD).
        Returns E and u normalized to sum u^2 h = 1 with u > 0 near its maximum."""
        h, r = self.h, self.r
        d = 2 * ckin / h ** 2 + V(r)
        e = -ckin / h ** 2 * np.ones(len(r) - 1)
        w, v = eigh_tridiagonal(d, e, select="i", select_range=(nstate, nstate))
        u = v[:, 0] / np.sqrt(h)
        if u[np.argmax(np.abs(u))] < 0:
            u = -u
        return w[0], u

    def dot(self, u1, u2):
        return self.h * np.dot(u1, u2)


def energy_richardson(ckin, V, rmax, n=4000):
    """Richardson-extrapolated FD energy (h and h/2)."""
    E1, _ = RadialGrid(rmax, n).ground(ckin, V)
    E2, _ = RadialGrid(rmax, 2 * n + 1).ground(ckin, V)
    return (4 * E2 - E1) / 3


# ------------------------------------------------- the two radial problems
def R_problem(M, g, k=1.0, kp=0.0):
    """-(1/M) Lap_R + (K/4) R^2 + g/R^2, L = 0  ->  ckin = 1/M, V = K R^2/4 + g/R^2."""
    K = k + 2 * kp
    return 1.0 / M, (lambda r: K * r ** 2 / 4 + g / r ** 2)


def rho_problem(M, k=1.0):
    """-(2M+1)/(4M) Lap_rho + k rho^2."""
    return (2 * M + 1) / (4 * M), (lambda r: k * r ** 2)


def dilated(ckin, V, s):
    """If H psi = E psi, U_s psi (U_s psi(x) = s^{3/2} psi(s x)) is the ground state of
    -(ckin/s^2) Lap + V(s x)."""
    return ckin / s ** 2, (lambda r: V(s * r))


def states_and_derivs(grid, problem, M, dl=1e-3, dt=1e-3):
    """u(M), d u/d lnM and D u (dilation generator) by central differences on a fixed grid."""
    def st(MM, s=1.0):
        c, V = problem(MM)
        c, V = dilated(c, V, s)
        return grid.ground(c, V)[1]
    u0 = st(M)
    du = (st(M * np.exp(dl)) - st(M * np.exp(-dl))) / (2 * dl)
    Du = (st(M, np.exp(dt)) - st(M, np.exp(-dt))) / (2 * dt)
    return u0, du, Du


def chi_numeric(M, g, k=1.0, kp=0.0, n=6000):
    """Fixed, global-dilation and separate-dilation susceptibilities from FD eigenvectors.

    Uses the exact product structure psi = f(R) phi(rho) (Jacobi separation, verified
    symbolically elsewhere): with <f|df> = <f|Df> = 0 for real normalized states,
      chi = ||df||^2 + ||dphi||^2,  <dpsi|Dpsi> = <df|Df> + <dphi|Dphi>,
      ||Dpsi||^2 = ||Df||^2 + ||Dphi||^2.
    """
    a = np.sqrt(M * (k + 2 * kp)) / 2
    lam = 0.5 * (np.sqrt(1 + 4 * g * M) - 1)
    RmaxR = np.sqrt((lam + 1.5) / a) + 14 / np.sqrt(a)
    b = 2 * np.sqrt(k * M / (2 * M + 1))
    gR = RadialGrid(RmaxR, n)
    gr = RadialGrid(14 / np.sqrt(b), n)
    f, df, Df = states_and_derivs(gR, lambda MM: R_problem(MM, g, k, kp), M)
    p, dp, Dp = states_and_derivs(gr, lambda MM: rho_problem(MM, k), M)
    # project out residual components along the state (numerical noise only)
    df -= gR.dot(f, df) * f
    dp -= gr.dot(p, dp) * p
    chiR = gR.dot(df, df)
    chirho = gr.dot(dp, dp)
    cross = gR.dot(df, Df) + gr.dot(dp, Dp)
    nD = gR.dot(Df, Df) + gr.dot(Dp, Dp)
    chi = chiR + chirho
    chi_inv = chi - cross ** 2 / nD
    chi_sep = (chiR - gR.dot(df, Df) ** 2 / gR.dot(Df, Df)
               + chirho - gr.dot(dp, Dp) ** 2 / gr.dot(Dp, Dp))
    return dict(chi_R=chiR, chi_rho=chirho, chi_fixed=chi, chi_inv=chi_inv, chi_sep=chi_sep,
                normD2=nD)


def overlap_numeric(M1, M2, g, k=1.0, kp=0.0, s=1.0, n=6000):
    """<psi(M1)|U_s psi(M2)> from FD states on a common grid."""
    Mmin = min(M1, M2)
    a = np.sqrt(Mmin * (k + 2 * kp)) / 2 * min(1.0, s ** 2)
    lam = 0.5 * (np.sqrt(1 + 4 * g * max(M1, M2)) - 1)
    b = 2 * np.sqrt(k * Mmin / (2 * Mmin + 1)) * min(1.0, s ** 2)
    gR = RadialGrid(np.sqrt((lam + 1.5) / a) + 14 / np.sqrt(a), n)
    gr = RadialGrid(14 / np.sqrt(b), n)
    f1 = gR.ground(*R_problem(M1, g, k, kp))[1]
    f2 = gR.ground(*dilated(*R_problem(M2, g, k, kp), s))[1]
    p1 = gr.ground(*rho_problem(M1, k))[1]
    p2 = gr.ground(*dilated(*rho_problem(M2, k), s))[1]
    return gR.dot(f1, f2) * gr.dot(p1, p2)


def expectations_numeric(M, g, k=1.0, kp=0.0, n=8000):
    """E, <R>, <R^2>, <T_R>, <T_rho> from FD states (radial quadrature)."""
    a = np.sqrt(M * (k + 2 * kp)) / 2
    lam = 0.5 * (np.sqrt(1 + 4 * g * M) - 1)
    gR = RadialGrid(np.sqrt((lam + 1.5) / a) + 14 / np.sqrt(a), n)
    c, V = R_problem(M, g, k, kp)
    ER, u = gR.ground(c, V)
    r, h = gR.r, gR.h
    P = u ** 2 * h
    VR = np.sum(P * V(r))
    b = 2 * np.sqrt(k * M / (2 * M + 1))
    gr = RadialGrid(14 / np.sqrt(b), n)
    cr, Vr = rho_problem(M, k)
    Er, w = gr.ground(cr, Vr)
    Vrho = np.sum(w ** 2 * gr.h * Vr(gr.r))
    return dict(E=ER + Er, R=np.sum(P * r), R2=np.sum(P * r ** 2), T_R=ER - VR,
                T_rho=Er - Vrho, r=r, u=u)


def normal_mode_frequencies(M, k=1.0, kp=0.0):
    """Model A in LAB coordinates, one Cartesian component: sqrt of eigenvalues of
    the mass-weighted Hessian. Expect {0, omega_R, omega_rho}."""
    m = np.array([M, M, 1.0])
    # V = k/2 (x3-x1)^2 + k/2 (x3-x2)^2 + kp/2 (x1-x2)^2
    H = np.array([[k + kp, -kp, -k], [-kp, k + kp, -k], [-k, -k, 2 * k]], float)
    W = H / np.sqrt(np.outer(m, m))
    ev = np.linalg.eigvalsh(W)
    return np.sqrt(np.clip(ev, 0, None))


# ------------------------------------------------------- entanglement checks
def purity_R_s_quadrature(fR, r, h, c):
    """Tr gamma^2 = int int f(R)^2 f(R')^2 exp(-2c|R-R'|^2) d^3R d^3R'
    with the angular integral done by Gauss-Legendre (no closed-form kernel)."""
    stride = max(1, len(r) // 800)
    r, fR, h = r[::stride], fR[::stride], h * stride
    x, w = np.polynomial.legendre.leggauss(64)
    U, Vv = np.meshgrid(r, r, indexing="ij")
    ang = np.zeros_like(U)
    for xi, wi in zip(x, w):
        ang += wi * np.exp(-2 * c * (U ** 2 + Vv ** 2 - 2 * U * Vv * xi))
    ang *= 2 * np.pi * 4 * np.pi          # int dOmega dOmega' = 4pi * 2pi int dcos
    q = fR ** 2 * r ** 2 * h
    return q @ ang @ q


def gaussian_1d_purity_grid(A11, A12, A22, L=None, n=400):
    """Single Cartesian component: purity of the u-marginal of
    exp(-(A11 u^2 + 2 A12 u v + A22 v^2)/2) on a grid, via SVD."""
    L = L or 10 / np.sqrt(min(A11 - A12 ** 2 / A22, A22 - A12 ** 2 / A11))
    x = np.linspace(-L, L, n)
    U, V = np.meshgrid(x, x, indexing="ij")
    psi = np.exp(-(A11 * U ** 2 + 2 * A12 * U * V + A22 * V ** 2) / 2)
    s = np.linalg.svd(psi, compute_uv=False) ** 2
    p = s / s.sum()
    return np.sum(p ** 2), -np.sum(p[p > 0] * np.log(p[p > 0]))


def cm_density_numeric(x, y, n_mid, s2):
    """n_cm(x) = int n_mid(y) G_{s2}(x - y) d^3y, angular integral by Gauss-Legendre,
    radial by trapezoid on the given y grid."""
    t, w = np.polynomial.legendre.leggauss(80)
    out = []
    for xx in np.atleast_1d(x):
        d2 = xx ** 2 + y[:, None] ** 2 - 2 * xx * y[:, None] * t[None, :]
        ang = 2 * np.pi * (np.exp(-d2 / (2 * s2)) @ w)
        out.append(np.trapz(y ** 2 * n_mid * ang, y) * (2 * np.pi * s2) ** -1.5)
    return np.array(out)
