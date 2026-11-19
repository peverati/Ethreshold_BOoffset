"""Re-derivation of Muller-Herold (JCP 124, 014105) critical mass relation, sympy + numerics."""
import sympy as sp, numpy as np
from scipy.optimize import brentq
me, mu, t, x, r, c = sp.symbols('m_e m_u t x r c', positive=True)
M = 2*me + mu
muR = 2*me*mu/M                      # reduced mass of third particle
eps2 = sp.simplify(sp.sqrt(8*muR)*(M/mu)**2)  # |psi2|^2 = exp(-sqrt(8 mu) r2^2), r2=-(M/m_u)(x+r1/2)
print("eps^2 =", eps2)
# curvature at x=0 of n(x) = int d3r r^{2t} e^{-c r^2} e^{-eps2 (x^2 + r^2/4 + x.r)}
# angular average of exp(-eps2 x r cos) expanded: 1 + eps2^2 x^2 r^2/6 ; x^2 coefficient:
e2, cc, tt = sp.symbols('e2 c tt', positive=True)
w = r**(2*tt+2)*sp.exp(-(cc+e2/4)*r**2)
num = sp.integrate(w*(-e2 + e2**2*r**2/6), (r, 0, sp.oo), conds='none')
num = sp.simplify(num)
print("curvature numerator ~", num)
# proportional to  (tt/6)*e2 - cc  (times positive factors): verify
chk = sp.simplify(sp.gammasimp((num/(sp.integrate(w,(r,0,sp.oo),conds='none')) - (-e2 + e2**2*(tt+sp.Rational(3,2))/(6*(cc+e2/4)))).rewrite(sp.gamma)))
print("check =", chk)
# critical condition: eps2 * t/6 = c with c = sqrt(t(t+1)) = sqrt(me), t = (-1+sqrt(1+4me))/2
tme = (-1+sp.sqrt(1+4*me))/2
cond = sp.Eq(eps2, 6*sp.sqrt(me)/tme)
mu_c = 2*me/((sp.Rational(3,2)/tme)**sp.Rational(2,3) - 1)
print("closed form m_u = 2 m_e/((3/(2t))^(2/3)-1):", mu_c)
res = sp.simplify((eps2 - 6*sp.sqrt(me)/tme).subs(mu, mu_c))
print("residual of eps2 - 6 sqrt(me)/t at closed form (numeric @ me=2):", sp.N(res.subs(me,2)), sp.N(res.subs(me,0.7)))
f = sp.lambdify(me, mu_c, 'numpy')
# asymptote m_u -> inf: (M/m_u)->1 => t = 3/2 => m_e = t(t+1) = 15/4
print("asymptotic m_e = 15/4 =", 15/4, " m_u at m_e=3.74:", f(3.74))
for m in [0.1,0.5,1,2,3,3.46,3.7]:
    print(m, f(m))
h = lambda m: f(m)-200
print("m_e(m_u=200) =", brentq(h, 1e-3, 3.74999))
print("m_u at m_e=15/16:", f(15/16), "(finite, so 15/16 is not an asymptote of this relation)")
