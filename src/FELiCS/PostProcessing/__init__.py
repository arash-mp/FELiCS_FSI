"""
FELiCS.PostProcessing -- post-processing utilities.

These modules deliberately do NOT import FELiCS itself (no dolfinx / PETSc /
SLEPc).  They work on the CSV files written by a run, so they can be used on a
laptop against results copied off the HPC, where the FE stack is not installed.

Entry points
------------
    python -m FELiCS.PostProcessing.plot_fsi_modes   fsi_structural_modes.csv
    python -m FELiCS.PostProcessing.plot_fsi_growth  fsi_structural_modes.csv

or, after running install.sh, simply::

    fsi-modes   fsi_structural_modes.csv
    fsi-growth  fsi_structural_modes.csv --t-end 150
"""

__all__ = ["fsi_common"]
