"""bo_h2plus: exact Born-Oppenheimer reference for the {X+, X+, e-} family.

Modules
-------
electronic : exact two-center one-electron energies (prolate spheroidal, Jaffe/Legendre)
adiabatic  : diagonal BO correction D(R) (validated, optional)
potential  : cached piecewise-Chebyshev interpolants U(R), D(R) and mass scaling
nuclear    : all bound vibrational levels for mass M and rotation J (mapped sinc-DVR),
             level counting by threshold node counting
scan       : mass scans (CSV) and BO critical masses
"""
from .electronic import electronic_energy, total_energy, dEdR, equilibrium, Wavefunction
from .nuclear import solve_levels, count_levels, thresholds
from .scan import mass_scan, critical_masses, critical_mass, m_eff, HARTREE_CM, M_H2PLUS

__version__ = "1.0.0"
