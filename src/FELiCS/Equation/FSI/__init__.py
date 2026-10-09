"""
FELiCS FSI subpackage.

Rigid-body-motion fluid-structure-interaction coupling for linear stability
analysis, following Negi, Hanifi & Henningson (J. Fluid Mech. 903, 2020, A35).

Modal analysis      -> RigidBodyMotionFSI      (augmented operators, size N + 2n)
Resolvent analysis  -> FSIResolventBorder      (borders the weights/norms/restrictors)
                       resolvent_shift         (standard vs discounted shift)
"""
from FELiCS.Equation.FSI.RigidBodyFSI import RigidBodyMotionFSI
from FELiCS.Equation.FSI.ResolventFSI import FSIResolventBorder, resolvent_shift

__all__ = [
    "RigidBodyMotionFSI",
    "FSIResolventBorder",
    "resolvent_shift",
]