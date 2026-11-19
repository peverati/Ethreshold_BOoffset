# harmonic_model: exactly solvable three-body testbeds for the (M, M, 1) family

*Last updated: 2026-11-12 18:00*

Two models of two identical particles (mass $M$, labels 1, 2) and one particle of
mass 1 (label 3), in which **every** quantity planned for the Coulomb study
(energy derivatives, pair distribution, Lindemann ratio, one-particle densities
about two reference points, fidelity susceptibilities with and without a dilation,
entanglement) is available in closed form as a function of $M$. Purpose: (i) test the
*definitions* (which measures are meaningful, which are conventions), (ii) validate
numerical code.

Units $\hbar=1$, light mass $=1$. Python `python3`
(numpy 1.26, scipy 1.13, sympy 1.13, mpmath; nothing installed). BLAS threads
are limited to 2 (`HARMONIC_BLAS_THREADS`, enforced with threadpoolctl).

```
harmonic_model/
  common.py     overlaps of Gaussians and R^lam e^{-aR^2/2} functions, Fisher metric,
                Gaussian-mode entanglement, partial-wave spectrum of gamma = f f' e^{-c|u-u'|^2}
  model_a.py    Model A closed forms (all harmonic)
  model_b.py    Model B closed forms (harmonic + g/R^2); g = 0 reproduces Model A
  twobody.py    two-body harmonic and Coulomb fidelities
  numerics.py   INDEPENDENT numerics used only by the tests (FD radial solver, quadrature, SVD)
  __main__.py   CLI: CSV over a list of masses
tests/          49 tests, ~5-10 s:  python3 -m pytest -q tests
examples/make_tables.py   regenerates results/*.csv and results/tables.md (~30 s, mpmath densities)
```

CLI: `python3 -m harmonic_model --model B --g 1 --masses 0.1 1 10 100 1836 -o out.csv`
(`--model A|B|twobody`, `--k`, `--kp`, `--g`).

---------------------------------------------------------------------------

## 1. Coordinates (both models)

$\mathbf R=\mathbf r_2-\mathbf r_1$, $\boldsymbol\rho=\mathbf r_3-\tfrac12(\mathbf r_1+\mathbf r_2)$, $M_t=2M+1$:
$$\hat T_{\rm int}=-\frac1M\nabla_R^2-\frac{2M+1}{4M}\nabla_\rho^2,\qquad
\mu_R=\frac M2,\ \ \mu_\rho=\frac{2M}{2M+1},$$
$$r_{13}^2+r_{23}^2=|\boldsymbol\rho+\tfrac12\mathbf R|^2+|\boldsymbol\rho-\tfrac12\mathbf R|^2=2\rho^2+\tfrac12R^2 .$$
**Consequence:** with equal harmonic springs between the unlike pairs, *any* like-pair
interaction $V(R)$ gives a Hamiltonian that separates **exactly** into an $R$ part and a
$\rho$ part. Both models below are therefore exact products $\psi=f(R)\,\phi(\rho)$, and
**the $R|\rho$ entanglement is identically zero**. These models cannot test $R|\rho$
(electron-nuclear) entanglement or any adiabatic coupling; see §5.

Positions relative to the centre of mass ($x_i=r_i-R_{\rm cm}$):
$x_1=-\tfrac12\mathbf R-\boldsymbol\rho/M_t$, $x_2=\tfrac12\mathbf R-\boldsymbol\rho/M_t$, $x_3=2M\boldsymbol\rho/M_t$.
(Symbolic checks: `tests/test_symbolic.py`.)

Notation: $\phi_b(x)=(b/\pi)^{3/4}e^{-bx^2/2}$; $f_{\lambda,a}(R)=N R^\lambda e^{-aR^2/2}$ (normalized
in 3D, $4\pi N^2=2a^{\kappa}/\Gamma(\kappa)$, $\kappa\equiv\lambda+\tfrac32$). Dots are $d/d\ln M$.

## 2. Model A (all harmonic)

$V=\tfrac k2(r_{13}^2+r_{23}^2)+\tfrac{k'}2 r_{12}^2=k\rho^2+\tfrac K4R^2$, $K=k+2k'$.

* **Stability:** bound iff $k>0$ and $K>0$, i.e. $k'>-k/2$ (a repulsive like-pair spring is
  allowed up to that point). Normal-mode frequencies in lab coordinates $\{0,\omega_R,\omega_\rho\}$ checked.
* $\omega_R=\sqrt{K/M}$, $a=\mu_R\omega_R=\tfrac12\sqrt{MK}$; $\omega_\rho=\sqrt{k(2M+1)/M}$, $b=\mu_\rho\omega_\rho=2\sqrt{kM/(2M+1)}$.
* $\psi_0=\phi_a(R)\phi_b(\rho)$, $E_0=\tfrac32(\omega_R+\omega_\rho)$.
* $\dot E=-\tfrac34\omega_R-\tfrac34\omega_\rho/(2M+1)$ $=-\langle T_R\rangle-\langle T_\rho\rangle/(2M+1)$ (Hellmann–Feynman, with virial $\langle T\rangle=\tfrac34\omega$);
  $\ddot E=\tfrac38\omega_R+\tfrac38\omega_\rho\,(4M+1)/(2M+1)^2$.
* Pair distance: Maxwell law $P(R)=4\pi R^2\phi_a^2$, $\langle R\rangle=2/\sqrt{\pi a}$, $\langle R^2\rangle=3/(2a)$,
  **$\Delta R/\langle R\rangle=\sqrt{3\pi/8-1}=0.42202$ for all $M,k,k'$** (no localization ever).
* Heavy-particle density: isotropic Gaussian, variance per component
  $\sigma^2_{\rm cm}=1/(8a)+1/[2b(2M+1)^2]$ about the CM, $\sigma^2_{\rm mid}=1/(8a)$ about the pair midpoint.
  Maximum at the origin for every $M$: **no hole in any reference frame**.
* Large $M$: $E=\tfrac32\sqrt{2k}+\tfrac32\sqrt{K/M}+\tfrac38\sqrt{2k}/M+\dots$; all lengths of $R$ shrink as $M^{-1/4}$, $\rho$ fixed.

## 3. Model B (harmonic + inverse-square like-pair repulsion)

$V=k\rho^2+\tfrac K4R^2+g/R^2$, $g\ge0$. With $u=Rf$ and angular momentum $J$ of $\mathbf R$:
$$-\frac1M u''+\Big[\frac{J(J+1)+gM}{MR^2}+\frac K4R^2\Big]u=Eu,$$
a 3D oscillator with non-integer "angular momentum" $\lambda_J(\lambda_J+1)=J(J+1)+gM$. **Fully closed form**
(verified by substitution in sympy and by FD diagonalization):
$$\lambda=\tfrac12(\sqrt{1+4gM}-1),\quad \psi_0=f_{\lambda,a}(R)\,\phi_b(\rho),\quad E_{nJ}=\omega_R(2n+\lambda_J+\tfrac32)+\tfrac32\omega_\rho ,$$
same $a,b,\omega$ as Model A. Total $L=J$ (the $\rho$ factor is $s$-wave). $g=0$ gives Model A.

* $\dot\lambda=gM/\sqrt{1+4gM}$, $\ddot\lambda=gM(1+2gM)/(1+4gM)^{3/2}$;
  $\dot E=-\tfrac12\omega_R\kappa+\omega_R\dot\lambda-\tfrac34\omega_\rho/(2M+1)$;
  $\ddot E=\tfrac14\omega_R\kappa-\omega_R\dot\lambda+\omega_R\ddot\lambda+\tfrac38\omega_\rho(4M+1)/(2M+1)^2$.
  HF check: $\langle T_R\rangle=\omega_R\kappa-\tfrac K4\langle R^2\rangle-g\langle R^{-2}\rangle$, $\langle R^2\rangle=\kappa/a$, $\langle R^{-2}\rangle=a/(\lambda+\tfrac12)$.
* **M-independent bond length**: classical minimum of $g/R^2+KR^2/4$, $R_e=(4g/K)^{1/4}$; $\langle R^2\rangle\to R_e^2$,
  vibrational quantum $2\omega_R=2\sqrt{K/M}$ (equals $\sqrt{V''(R_e)/\mu_R}$).
  $E=\tfrac32\sqrt{2k}+\sqrt{gK}+\sqrt{K/M}+O(M^{-1})$.
* Pair distribution $4\pi R^2f^2\propto R^{2\lambda+2}e^{-aR^2}$ (chi law), peak $\sqrt{(\lambda+1)/a}$;
  $\langle R^p\rangle=\Gamma(\kappa+p/2)/[\Gamma(\kappa)a^{p/2}]$;
  $\eta^2=(\Delta R/\langle R\rangle)^2=\kappa\,\Gamma(\kappa)^2/\Gamma(\kappa+\tfrac12)^2-1$, a function of $gM$ only;
  $\eta\simeq\tfrac12(gM)^{-1/4}$ at large $M$ (BO-like). $g_{XX}(R)=f^2\propto R^{2\lambda}$ vanishes at $R=0$ for **every** $M>0$.
* **Rotational tower** (bonus, exact): $E_{0J}-E_{00}=\omega_R(\lambda_J-\lambda_0)\to J(J+1)/(MR_e^2)$. The rigid-rotor estimate
  $2B=2\langle R^{-2}\rangle/M$ is 23% too large at $M=1$ (0.894 vs 0.685) and 0.03% at $M=1836$.

### Heavy-particle densities and the reference point
* About the pair midpoint: $n_{\rm mid}(y)=8f(2y)^2\propto y^{2\lambda}e^{-4ay^2}$: **zero at $y=0$ for all $M>0$** (hole at every mass).
* About the CM, $x_1=-\tfrac12\mathbf R-\boldsymbol\rho/M_t$ is the convolution of $n_{\rm mid}$ with a Gaussian of
  variance $s^2=1/[2bM_t^2]$:
$$n_{\rm cm}(x)=\frac{8\cdot 2^{2\lambda}a^\kappa}{(2\pi s^2)^{3/2}P^\kappa}\,e^{-x^2/2s^2}\,{}_1F_1\!\Big(\kappa;\tfrac32;\frac{x^2}{4Ps^4}\Big),\qquad P=4a+\frac1{2s^2}.$$
  Curvature at the origin: $2s^2\,n''/n|_0=(\lambda-12as^2)/(12as^2+\tfrac32)$. **Hole iff**
$$\lambda(gM)>\tfrac32\sqrt{K/k}\,(2M+1)^{-3/2},$$
  a unique crossing $M_c(g)$ (left side increasing, right side decreasing). $g=1,k=1$: $M_c=0.63306$; the shell radius grows as
  $(M-M_c)^{1/2}$ (generic pitchfork) and tends to $R_e/2$. $M_c$ versus $g$ (k=1):

| g | 0.01 | 0.1 | 1 | 10 | 100 |
|---|---|---|---|---|---|
| $M_c$ | 4.689 | 1.779 | 0.6331 | 0.1825 | 0.03230 |
| $\eta(M_c)$ | 0.416 | 0.401 | 0.369 | 0.327 | 0.301 |

  The criterion is set by a competition between $\lambda$ and the *kinematic* smearing of the CM by the light particle.
  The midpoint frame has a hole at all $M$; the CM frame only above $M_c$; Model A never.

## 4. Fidelity and dilation

Convention: $F(M,M')=|\langle\psi(M)|\psi(M')\rangle|$ in fixed Jacobi coordinates,
$F(M,Me^{\delta})=1-\tfrac12\chi\delta^2+\dots$, $\chi=\langle\dot\psi|\dot\psi\rangle-|\langle\psi|\dot\psi\rangle|^2$.
Dilation $(U_s\psi)(x)=s^{3}\psi(sx)$ on all 6 internal coordinates; $U_s\phi_b=\phi_{bs^2}$, $U_sf_{\lambda,a}=f_{\lambda,as^2}$.
$F_{\rm inv}=\max_s|\langle\psi(M)|U_s\psi(M')\rangle|$. To second order,
$$\chi_{\rm inv}=\chi-\frac{\langle\dot\psi|D\psi\rangle^2}{\|D\psi\|^2},\qquad D=\mathbf x\cdot\nabla+\tfrac n2$$
(real states): the dilation removes **only the component of $\dot\psi$ along the dilation generator**, i.e. the part of
the change that is a change of the overall length unit; everything that changes a *ratio* of lengths or the *shape* of a
distribution survives. A further variant, $\chi_{\rm sep}$, maximizes over separate dilations of $\mathbf R$ and $\boldsymbol\rho$.

Building blocks: $\langle\phi_{b_1}|\phi_{b_2}\rangle={\rm sech}^{3/2}(\epsilon/2)$, $\epsilon=\ln(b_2/b_1)$;
$\langle f_{\lambda_1a_1}|f_{\lambda_2a_2}\rangle=\Gamma(\bar\kappa)(\tfrac{a_1+a_2}2)^{-\bar\kappa}/\sqrt{\Gamma(\kappa_1)a_1^{-\kappa_1}\Gamma(\kappa_2)a_2^{-\kappa_2}}$.
Fisher metric of $|f|^2$ in $(\lambda,\ln a)$: $\begin{pmatrix}\psi_1(\kappa)&-1\\-1&\kappa\end{pmatrix}$ ($\psi_1$ = trigamma), and $\chi=\tfrac14\dot\theta^TI\dot\theta$.
Log rates: $\dot{\ln a}=\tfrac12$, $\dot{\ln b}=1/[2(2M+1)]$ (independent of $k,k',g$).

**Model A** (independent of $k,k'$):
$$F={\rm sech}^{3/2}\tfrac{\epsilon_a}2\,{\rm sech}^{3/2}\tfrac{\epsilon_b}2,\qquad
F_{\rm inv}={\rm sech}^3\frac{\epsilon_a-\epsilon_b}4,\qquad F_{\rm sep}\equiv1,$$
$$\chi=\tfrac38(\dot{\ln a}^2+\dot{\ln b}^2)=\frac3{32}+\frac3{32(2M+1)^2},\quad
\chi_{\rm inv}=\tfrac3{16}(\dot{\ln a}-\dot{\ln b})^2=\frac{3}{16}\frac{M^2}{(2M+1)^2},\quad
\text{removed}=\frac3{16}\Big(\frac{M+1}{2M+1}\Big)^2 .$$

| | $M\to0$ | $M=1$ | $M\to\infty$ |
|---|---|---|---|
| $\chi$ fixed | 3/16 | 5/48 | 3/32 |
| removed by dilation | 3/16 (all) | 1/12 | 3/64 |
| $\chi_{\rm inv}$ | 0 | 1/48 | **3/64** |
| $\chi_{\rm sep}$ | 0 | 0 | 0 |

At $M\to0$ the whole state is a pure dilation (both lengths $\propto M^{-1/4}$). At large $M$ the $R$ lengths shrink
as $M^{-1/4}$ while $\rho$ is fixed: a change of the *ratio* of two lengths, which a global dilation cannot absorb.

**Model B**, with $\kappa=\lambda+\tfrac32$, $x_b=1/[2(2M+1)]$:
$$\chi=\tfrac14\big[\psi_1(\kappa)\dot\lambda^2-\dot\lambda+\tfrac14\kappa\big]+\tfrac38x_b^2,\quad
\langle\dot\psi|D\psi\rangle=\tfrac12(\tfrac12\kappa-\dot\lambda)+\tfrac34x_b,\quad \|D\psi\|^2=\kappa+\tfrac32,$$
$$\chi_{\rm sep}=\tfrac14\dot\lambda^2\big[\psi_1(\kappa)-1/\kappa\big].$$
$F_{\rm inv}$, $F_{\rm sep}$ at finite separation: 1D maximization of the closed-form overlap. Large $M$ (fitted
coefficients, agree to 4-5 digits for g = 0.1, 1, 10, *not derived analytically*):
$\chi\simeq\tfrac1{32}+\tfrac{1}{96}(gM)^{-1/2}$, $\chi_{\rm inv}\simeq\chi_{\rm sep}\simeq\tfrac1{32}-\tfrac5{96}(gM)^{-1/2}$.
Small $M$: $\chi\to3/16$, $\chi_{\rm inv},\chi_{\rm sep}\to0$. For $g=1$ all three are **monotone** in $M$ (no maximum).
The plateau 1/32 is $\tfrac12(d\ln\eta/d\ln M)^2$ with $\eta\propto M^{-1/4}$: the fidelity of a 1D Gaussian whose
width relative to a fixed $R_e$ shrinks as $M^{-1/4}$. **It survives both the global and the separate dilation.**

**Two body** (masses $M,1$, $\mu=M/(M+1)$, $\dot{\ln\mu}=1/(M+1)$):
harmonic $F={\rm sech}^{3/2}(\tfrac14\ln\tfrac{\mu_2}{\mu_1})$, $\chi=3/[32(M+1)^2]$; Coulomb
$F=8(\mu_1\mu_2)^{3/2}/(\mu_1+\mu_2)^3={\rm sech}^3(\tfrac12\ln\tfrac{\mu_2}{\mu_1})$, $\chi=\tfrac34/(M+1)^2$ (**verified**).
Both: one length scale, so $U_{s=\mu_1/\mu_2}$ maps $\psi(M_2)$ exactly onto $\psi(M_1)$: $F_{\rm inv}\equiv1$, $\chi_{\rm inv}\equiv0$ (**verified** by quadrature).

## 5. Entanglement

Gaussian modes (Model A): for $\psi\propto e^{-(A_{11}u^2+2A_{12}uv+A_{22}v^2)/2}$ per Cartesian component,
single-component purity $P_1=\sqrt{1-A_{12}^2/(A_{11}A_{22})}$, symplectic eigenvalue $\nu=1/(2P_1)$,
$S_{\rm vN}=3[(\nu+\tfrac12)\ln(\nu+\tfrac12)-(\nu-\tfrac12)\ln(\nu-\tfrac12)]$, $S_L=1-P_1^3$.

| bipartition | Model A $P_1$ | large-$M$ | Model B |
|---|---|---|---|
| $R\,|\,\rho$ (= particle 3 vs pair internal) | 1 | **exact product** | **exact product** |
| $R\,|\,s$, $s=r_3-r_1$ (ATOM-MOL-nonBO coords) | $\sqrt{4a/(4a+b)}$ | $S_L\simeq\tfrac38 b/a\propto M^{-1/2}$ | $\to$ const: $S_{\rm vN}=1.0315$, $S_L=0.4648$ ($g=k=1$) |
| $x_1\,|\,x_2$ (heavy-heavy, CM frame) | $2\sqrt{aB}/(a+B)$, $B=bM_t^2/4$ | $S_{\rm vN}\simeq\tfrac94\ln M$ diverges | not computed |
| $x_1\,|\,x_3$ (heavy-light, CM frame) | $\sqrt{B/(a+B)}$ | $S_L\simeq\tfrac32a/B\propto M^{-3/2}$ | not computed |

Model B, $R|s$: $\psi(R,s)=f(R)\phi_b(s-R/2)$ gives $\gamma_R(\mathbf R,\mathbf R')=f(R)f(R')e^{-b|\mathbf R-\mathbf R'|^2/16}$, diagonalized
by partial waves (radial kernel $4\pi ff'e^{-c(R^2+R'^2)}i_l(2cRR')$, multiplicity $2l+1$; Gauss–Legendre radial
quadrature, i.e. *numerical*, converged to ~1e-12). Large-$M$ limit (thin shell at $R_e$): eigenvalues
$p_l=e^{-x}i_l(x)$, $x=bR_e^2/8$: a closed-form orientational entanglement that does **not** vanish. Note that a
subsystem "particle 1" is not defined after CM removal: its conjugate momentum depends on which second coordinate
completes it ($x_1|x_2$ and $x_1|x_3$ give completely different answers for the "same" $x_1$).

## 6. Tables (k = 1, k' = 0, Model B g = 1; full precision in `results/*.csv`)

Model A:

| M | E | dE/dlnM | d2E/dlnM2 | ⟨R⟩ | η | σ²_cm | σ²_mid | χ | χ_inv | S_vN(R\|s) | S_vN(x1\|x2) | S_vN(x1\|x3) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 9.9396 | -4.5368 | 2.4488 | 2.8377 | 0.42202 | 1.392 | 0.79057 | 0.15885 | 0.0013021 | 1.5761 | 0.089419 | 1.4166 |
| 1 | 4.0981 | -1.183 | 0.73584 | 1.5958 | 0.42202 | 0.29811 | 0.25 | 0.10417 | 0.020833 | 1.1967 | 1.5094 | 0.56601 |
| 10 | 2.648 | -0.28893 | 0.16911 | 0.89737 | 0.42202 | 0.079878 | 0.079057 | 0.093963 | 0.042517 | 0.61998 | 5.7017 | 0.054079 |
| 100 | 2.2766 | -0.08029 | 0.042777 | 0.50463 | 0.42202 | 0.025009 | 0.025 | 0.093752 | 0.04641 | 0.26331 | 10.774 | 0.0027215 |
| 1836 | 2.1566 | -0.017792 | 0.0090405 | 0.24378 | 0.42202 | 0.0058345 | 0.0058345 | 0.09375 | 0.046849 | 0.08007 | 17.311 | 4.9526e-05 |

Model B ($R_e=\sqrt2$, $M_c=0.63306$; curv = $2s^2n''/n$ at the CM, >0 means hole; shell = radius of the CM-density maximum):

| M | λ | E | dE/dlnM | d2E/dlnM2 | ⟨R⟩ | R_peak | η | curv | shell | χ | χ_inv | χ_sep | S_L(R\|s) | S_vN(R\|s) | E₁−E₀ (J) | ⟨R⁻²⟩/M |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 0.091608 | 10.229 | -4.4144 | 2.4831 | 2.9365 | 2.6275 | 0.40908 | -0.39737 | 0 | 0.14499 | 0.00059803 | 0.00042157 | 0.64023 | 1.6302 | 2.9769 | 2.6726 |
| 1 | 0.61803 | 4.7161 | -1.0448 | 0.71147 | 1.9412 | 1.7989 | 0.35231 | 0.18414 | 0.46707 | 0.061013 | 0.0096932 | 0.0064156 | 0.60178 | 1.4734 | 0.68474 | 0.44721 |
| 10 | 2.7016 | 3.5024 | -0.22222 | 0.14178 | 1.5824 | 1.5301 | 0.24731 | 1.7722 | 0.65293 | 0.036128 | 0.021643 | 0.018626 | 0.53084 | 1.2308 | 0.094374 | 0.049386 |
| 100 | 9.5125 | 3.2279 | -0.056539 | 0.031524 | 1.4673 | 1.45 | 0.1515 | 6.3391 | 0.68965 | 0.032445 | 0.027121 | 0.02648 | 0.48869 | 1.1023 | 0.0099382 | 0.0049938 |
| 1836 | 42.351 | 3.145 | -0.012026 | 0.0061913 | 1.4266 | 1.4225 | 0.075612 | 28.234 | 0.70299 | 0.031502 | 0.030106 | 0.03006 | 0.47064 | 1.0491 | 0.00054448 | 0.00027231 |

Two body: $\chi_{\rm harm}=0.077479,\,0.023438,\,7.7479\times10^{-4},\,9.19\times10^{-6},\,2.78\times10^{-8}$;
$\chi_{\rm Coul}=0.61983,\,0.1875,\,6.198\times10^{-3},\,7.35\times10^{-5},\,2.22\times10^{-7}$ at $M=0.1,1,10,100,1836$; $\chi_{\rm inv}=0$.

## 7. Lessons for the Coulomb problem

1. **Meaningful, reference free:** the pair distribution $g_{XX}(R)$ and the Lindemann ratio $\eta$. In Model A
   $\eta=0.422$ at every $M$ (correctly: no structure); in Model B $\eta$ decreases smoothly as $\tfrac12(gM)^{-1/4}$
   (localization at fixed $R_e$). $\eta$ separates the two models; no sharp mass is singled out.
2. **Convention dependent:** one-particle densities. The same Model B state has a hole about the pair midpoint at
   *every* $M>0$, a hole about the CM only for $M>M_c$ (closed-form criterion, controlled by the kinematic CM
   smearing $\propto(2M+1)^{-3/2}$; $M_c$ ranges 0.03 to 4.7 as $g$ varies, with $\eta(M_c)$ = 0.30 to 0.42), and Model A
   never has one. The shell radius $\propto(M-M_c)^{1/2}$ is the generic pitchfork (referee M3, M11 confirmed).
3. **Fixed-coordinate $\chi$** is dominated by dilation: at small $M$ it is entirely dilation (3/16 in both models,
   removed completely); the two-body values (including $\tfrac34/(M+1)^2$) are entirely dilation.
4. **Global-dilation $\chi_{\rm inv}$** removes only the overall length change. What survives is any *differential*
   scaling of lengths. Crucially, Model A, which has **no** structure, has $\chi_{\rm inv}\to3/64$, *larger* than Model B's
   $1/32$: in Model A the $R$ lengths shrink as $M^{-1/4}$ relative to the fixed $\rho$ cloud. So a plateau of
   $\chi_{\rm inv}$ is **not** by itself a signature of molecular structure; it signals that two length scales have different
   power laws in $M$. In Coulomb $\{X^+,X^+,e^-\}$ this happens too ($R_e$ and the electron cloud fixed, $\Delta R\propto M^{-1/4}$).
5. **Separate dilation of $R$ and $\rho$** ($\chi_{\rm sep}$) discriminates: 0 identically for Model A (and for two bodies),
   $\to1/32$ for Model B. The surviving 1/32 is $\tfrac12(d\ln\eta/d\ln M)^2$ with $\eta\propto M^{-1/4}$: it is the BO
   exponent of the relative width, generic for any harmonic-like localization at a fixed bond length, not a critical
   signature. In the BO variable $\kappa=M^{-1/4}$ it becomes $1/(2\kappa^2)$ (change of variable only). Recommendation:
   report $\chi$, $\chi_{\rm inv}$ and $\chi_{\rm sep}$, and interpret $\chi_{\rm sep}$ as "$d\ln\eta/d\ln M$ squared".
   For the Coulomb case it is natural because $R,\rho$ are mass-independent Jacobi coordinates for $(M,M,1)$.
6. **Entanglement:** only $R|\rho$ has the BO meaning (it should vanish as $M\to\infty$). Here it is exactly zero, which makes
   both models a strict null test for code (any numerical $R|\rho$ entanglement routine must return 0 on them) but gives
   no benchmark for its value. Other bipartitions are dominated by kinematics: $R|s$ (the program's particle-1-referenced
   coordinates) stays O(1) at $M\to\infty$ in Model B (orientational correlation, closed-form limit 1.0315), $x_1|x_2$
   diverges as $\tfrac94\ln M$ from the CM constraint. Do not report single-particle or particle-referenced bipartitions as
   "electron–nuclear entanglement". A benchmark for nonzero $R|\rho$ entanglement needs a model with $R$-dependent
   electron coupling (e.g. unequal or anharmonic unlike interactions), which is no longer separable; not built here.
7. **Energy derivatives** are smooth and featureless in both models (HF identity checked); they carry no structural signal.
8. **Rotational tower:** Model B gives exact $J>0$ levels, a ready benchmark for the planned $J>0$ extension of
   `solvers/threebody` (if harmonic/inverse-square potentials are added) and a measure of how wrong the $\langle R^{-2}\rangle$
   rotor estimate is (23% at $M=1$).
9. **What these models cannot say:** no dissociation threshold (confining potential, infinitely many bound states at every
   $M$), hence nothing about the threshold critical masses of referee M6; no adiabatic/non-adiabatic coupling; no Coulomb
   cusps; and at small $M$ lengths scale as $M^{-1/4}$ rather than the Coulomb $M^{-1}$.

## 8. Verification (tests/)

* `test_symbolic.py`: Jacobi kinetic matrix and spring identity, CM-frame coordinate maps, Model B radial solution
  satisfies the ODE exactly, $\dot E,\ddot E$ for general $(g,k,k')$ against sympy differentiation.
* `test_numeric.py`: FD diagonalization of the radial equations (directly from the potentials) for energies (≤2e-5),
  HF with numerical kinetic energies, moments, $J=1,2$ levels; lab-frame normal modes; overlaps by quadrature;
  fidelities and global-dilation fidelities from FD eigenvectors maximized numerically over $s$; $\chi,\chi_{\rm inv},\chi_{\rm sep}$
  from central differences of FD eigenvectors in $\ln M$ and in the dilation parameter; limits; two-body susceptibilities and
  exact dilation mapping; Gaussian entanglement against grid SVD; partial-wave $R|s$ spectrum against Model A closed form and
  against an independent double-integral purity; CM density ($_1F_1$ form) against numerical 3D convolution;
  $M_c$ against the sign of a numerical curvature; shell radius.
* FD accuracy is limited near $R=0$ when $\lambda$ is small (non-smooth $R^{\lambda+1}$): relative ~3e-5 at $M=0.1$.

Not verified / not done: the large-$M$ *subleading* coefficients for Model B $\chi$ (1/96, 5/96) are fitted numerically,
not derived; Model B $x_1|x_2$ and $x_1|x_3$ entanglement not computed; $R|s$ von Neumann entropy for Model B is numerical
(partial waves), only its purity is cross-checked independently.
