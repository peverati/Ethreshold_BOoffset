"""Physical definition of a three-particle Coulomb system."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

from .operators import PAIRS

M_PROTON = 1836.15267343
M_MUON = 206.7682830


@dataclass
class ThreeBodySystem:
    """Three particles with masses m_i (inf allowed) and charges q_i.

    couplings: optional dict {(i,j): Z_ij} overriding the Coulomb strength
    q_i q_j of pair (i, j) (e.g. Z_01 = 0 to switch off one interaction).
    identical: (0, 1) if particles 0 and 1 are identical (same mass and
    charge); the solver then works in a basis of definite exchange symmetry.
    """

    masses: Sequence[float]
    charges: Sequence[float]
    couplings: Optional[dict] = None
    identical: Optional[tuple] = None
    name: str = ""
    extra: dict = field(default_factory=dict)

    def __post_init__(self):
        self.masses = tuple(float(m) for m in self.masses)
        self.charges = tuple(float(q) for q in self.charges)
        if self.identical is not None:
            if tuple(self.identical) != (0, 1):
                raise ValueError("identical particles must be labelled 0 and 1")
            if self.masses[0] != self.masses[1] or self.charges[0] != self.charges[1]:
                raise ValueError("particles 0 and 1 declared identical but differ")
            c = self.coupling_dict()
            if c[(0, 2)] != c[(1, 2)]:
                raise ValueError("couplings break the 0<->1 exchange symmetry")

    @property
    def inverse_masses(self):
        return tuple(0.0 if np.isinf(m) else 1.0 / m for m in self.masses)

    def coupling_dict(self):
        c = {p: self.charges[p[0]] * self.charges[p[1]] for p in PAIRS}
        if self.couplings:
            for (i, j), v in self.couplings.items():
                c[(min(i, j), max(i, j))] = float(v)
        return c

    def reduced_mass(self, i, j):
        im = self.inverse_masses
        s = im[i] + im[j]
        return np.inf if s == 0 else 1.0 / s

    def two_body_thresholds(self):
        """Ground-state energies of each attractive pair with the third particle at rest at infinity."""
        out = {}
        for (i, j), z in self.coupling_dict().items():
            mu = self.reduced_mass(i, j)
            out[(i, j)] = -0.5 * z * z * mu if z < 0 else 0.0
        return out

    def threshold(self):
        """Lowest two-body dissociation threshold (continuum edge for S states)."""
        return min(min(self.two_body_thresholds().values()), 0.0)


def symmetric_family(M: float) -> ThreeBodySystem:
    """{X+, X+, e-}: particles 0,1 of mass M, charge +1; particle 2 of mass 1, charge -1."""
    return ThreeBodySystem((M, M, 1.0), (1.0, 1.0, -1.0), identical=(0, 1), name=f"XXe_M={M:g}")


def helium_like(Z: float = 2.0, nuclear_mass: float = np.inf) -> ThreeBodySystem:
    """Two electrons (0, 1) and a nucleus (2) of charge Z."""
    return ThreeBodySystem((1.0, 1.0, nuclear_mass), (-1.0, -1.0, Z), identical=(0, 1), name=f"Z={Z:g}")


def family_threshold(M: float) -> float:
    return -M / (2.0 * (M + 1.0))
