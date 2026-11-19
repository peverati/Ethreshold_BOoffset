"""Compare Muller-Herold critical relation with solvers/harmonic_model Model B (critical_mass_cm),
plus independent quadrature of the pair-density curvature. Run: python3 compare.py"""
import sys, os
import numpy as np
from scipy.optimize import brentq
from scipy import integrate
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "solvers", "harmonic_model"))
from harmonic_model import model_b as B

def t_of(me): return 0.5*(-1+np.sqrt(1+4*me))
def mu_MH(me):                     # corrected Eq. 3.13: m_u(m_e)
    return 2*me/((1.5/t_of(me))**(2/3)-1)
def me_MH(mu):                     # inverse by root finding
    return brentq(lambda m: mu_MH(m)-mu, 1e-12, 3.75*(1-1e-13), xtol=1e-15, rtol=1e-14)

# ---- mapping wrapper: MH(m_e,m_u) == Model B with M=m_e/m_u, g=m_u, k=K=4 m_u, lengths unchanged
def mapB(mu): return dict(g=mu, k=4*mu, kp=0.0)
def me_modelB(mu):
    return mu*B.critical_mass_cm(**mapB(mu))

# ---- independent quadrature: curvature of n(z) on the z axis, directly from MH wavefunctions
def density_z(z, me, mu, nq=200):
    """unnormalised pair density section n(0,0,z), particle 1 about CM (MH Eqs 3.6-3.9, re-derived):
       q1 = -r1/2 - (m_u/M) r2 ; |psi1|^2 = r^{2t} exp(-sqrt(me) r^2); |psi2|^2 = exp(-sqrt(8 mu) r2^2)."""
    t = t_of(me); M = 2*me+mu; c1 = np.sqrt(me); c2 = np.sqrt(8*(2*me*mu/M))   # sqrt(8 mu_red)
    f = lambda cth, r: r**2*r**(2*t)*np.exp(-c1*r*r)*np.exp(-c2*(M/mu)**2*(r*r/4+z*z+r*z*cth))
    # use scaled exponent to avoid underflow: integrate with Gauss-Legendre in cos, adaptive in r
    val = integrate.dblquad(f, 0, 12/np.sqrt(c1)+10, -1, 1, epsabs=1e-14, epsrel=1e-12)[0]
    return val
def curvature_quad(me, mu, h=None):
    # second derivative at 0 via Richardson on symmetric differences, normalised by n(0)
    eps2 = np.sqrt(16*me*mu/(2*me+mu))*((2*me+mu)/mu)**2
    h = h or 0.15/np.sqrt(eps2)
    n0 = density_z(0, me, mu)
    d = lambda hh: (density_z(hh, me, mu)-2*n0+density_z(-hh, me, mu))/hh**2/n0
    return (4*d(h/2)-d(h))/3
def me_quad(mu):
    hi = me_MH(mu)
    return brentq(lambda m: curvature_quad(m, mu), 0.6*hi, 1.4*hi if 1.4*hi<3.75 else 3.74999, xtol=1e-9)

if __name__ == "__main__":
    print("%8s %14s %14s %14s %12s" % ("m_u","m_e MH(3.13)","m_e ModelB","m_e quadrature","rel diff B"))
    rows = []
    for mu in [0.5, 1.0, 5.0, 20.0, 200.0, 1836.15, 1e4]:
        a = me_MH(mu); b = me_modelB(mu)
        q = me_quad(mu) if mu <= 200 else float('nan')
        rows.append((mu, a, b, q))
        print("%8g %14.10f %14.10f %14.8f %12.2e" % (mu, a, b, q, abs(a-b)/a))
    # pointwise density check: package heavy_density_cm (mapped) vs direct quadrature, normalised on z axis
    print("\npointwise density check (me=2.0, mu=10): package vs quadrature, ratio n(z)/n(0)")
    me, mu = 2.0, 10.0; Mx = me/mu; p = mapB(mu)
    for z in [0.1, 0.3, 0.6]:
        qa = density_z(z, me, mu)/density_z(0, me, mu)
        pb = B.heavy_density_cm(z, Mx, **p)/B.heavy_density_cm(0.0, Mx, **p)
        print("z=%.2f quad %.10f package %.10f" % (z, qa, pb))
    print("\ncurvature coefficient sign at mapped M (me=3.0,mu=200 -> unimodal expected, <0):",
          B.cm_curvature_coefficient(3.0/200, **mapB(200.0)) )
    print("curvature coefficient (me=4.0 ... beyond asymptote, mu=200 -> hole expected, >0):",
          B.cm_curvature_coefficient(4.0/200, **mapB(200.0)))
