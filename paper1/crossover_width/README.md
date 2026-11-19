# Crossover width of the ground-state atom-to-molecule change

*Last updated: 2026-11-15 17:17*

**Why:** referee report 2026-11-15 (M5): the spread of the "emergence" markers must be compared with the width of the crossover, not with numerical error.

**What:** `width.py` (analysis of existing data only, `../mass_scans/descriptors.csv`, 37 masses) defines two widths from the exact ground-state Lindemann ratio L(M) and expresses the marker spread in units of them. Output: `results.txt`.

**Result (2026-11-15):**
* (a) local exponent f = -4 dlnL/dlnM (0 at the Lindemann maximum, 1 in the BO limit L ~ M^-1/4): f = 0.1 at M = 0.455, 0.5 at M = 1.52, 0.9 at M = 95. **Width 5.3 in ln M (2.3 decades).**
* (b) L from 90% to 10% of its maximum: M = 1.66 to about 3e4 (10% point from the M^-1/4 tail, an extrapolation). Width 9.9 in ln M (4.3 decades).
* Ground-state markers 0.2 to 1.11: spread 1.71 in ln M = **0.32 of width (a)**, 0.17 of (b). Literature-proposed markers only (0.4 to 0.98): 0.17 of (a).
* Sharp thresholds (4.21 to 98.27): spread 1.16 of width (a): they span the whole crossover.

**Reading:** the ground-state markers do not coincide, but they **cluster at the onset** of a broad crossover (all below its midpoint M = 1.5 of definition (a)), within about a third of its width. "No unique emergence point" is therefore true only in a weak sense; the stronger statements are (i) the crossover is broad (2 to 4 decades in M), (ii) the descriptor markers are onset indicators that disagree by a factor of about 5, (iii) the only sharp points are the bound-state thresholds, which lie inside and across the crossover.

**Caveats:** both widths depend on the choice of descriptor (Lindemann ratio) and on the level convention (10/90%); definition (b) needs a tail extrapolation. Finite differences on the 37-point grid; f levels move by a few percent with the grid.
