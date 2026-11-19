# CODATA-2022 copy of examples/mass_scan.py (made 2026-11-17 from the CODATA-2012 copy); writes into this folder.
"""Mass scan used for the README tables (results/scan_ground.csv, results/scan_excited1.csv).

Usage: python examples/mass_scan.py    (~10 minutes)
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
M_PROTON, M_MUON = 1836.152673426, 206.7682827  # CODATA 2022 (Paper 1)
from threebody.scan import scan_masses, write_csv  # noqa: E402

here = os.path.dirname(os.path.abspath(__file__))
Ms = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, M_MUON, 500.0, 1000.0, M_PROTON]
rows = scan_masses(Ms, [16, 20, 24, 28], nstates=4, N_opt=12)
write_csv(os.path.join(here, "scan_ground.csv"), rows)
# scales optimized for the first excited root: probes whether a 2nd S state is bound
rows = scan_masses([2.0, 5.0, 10.0, 20.0, 50.0], [16, 20, 24], nstates=3, N_opt=12, opt_state=1, observables=False)
write_csv(os.path.join(here, "scan_excited1.csv"), rows)
