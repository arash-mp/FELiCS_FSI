#  ___________________________________   _______________________________________________________
# /-----------------------------------\ /-------------------------------------------------------\
# |   (         (                (     |  This source code is part of FELiCS                     |
# |   )\ )     ) )        (      )\ )  |  (F)inite (E)lement (Li)nearized (C)ombustion (S)olver  |
# |  (()/(  (  (()/( (    )\   (()/(   |                                                         |
# |  /(_)) )\  /(_)))\  (((_)  /(_))   |  Licensed under the GNU GPLv3                           |
# |  (_)_)((_) (_)) ((_) )\___ (_))    |                                                         |
# |  | __|| __|| |   (_)((/ __|/ __|   |  (C) 2018-2025: The FELiCS Developers (www.felics.eu)   |
# |  | _| | _| | |__ | | | (__ \__ \   |  Visit          www.felics.eu                           |
# |  |_|  |___||____||_|  \___||___/   |  Contact        info@felics.eu                          |
# \___________________________________/ \_______________________________________________________/
#
"""
Rigid-body-motion fluid-structure-interaction (FSI) coupling for FELiCS.

Implements the linearised FSI formulation of Negi, Hanifi & Henningson,
J. Fluid Mech. 903 (2020) A35, for the DIRECT modal analysis of a rigid body
that is free to undergo one OR several rigid-body degrees of freedom (DOFs)
simultaneously -- any mix of translations along fixed directions and rotations
about fixed axes (e.g. coupled transverse-translation + rotation of a square
cylinder or an inclined ellipse).

Generalised model (n DOFs)
--------------------------
Generalised coordinates  eta = (eta_1, ..., eta_n)  (rigid-body amplitudes),
velocities  phi_k = d(eta_k)/dt.  The rigid-body displacement of a material
point x on the body is the superposition

    Delta x(x) = sum_k eta_k d_k(x),

with the per-DOF displacement field

    translation along e_k :  d_k(x) = e_k                       (constant)
    rotation about c_k     :  d_k(x) = (-(y - c_ky), (x - c_kx))

The fluid sees this through the transpiration boundary condition on Gamma
(paper Eq. 2.28), summed over all DOFs:

    u'_i + sum_k [ eta_k (d_k . grad)U0_i - phi_k (d_k)_i ] = 0   on Gamma

and each DOF feels its own generalized fluid force / moment (paper Eq. 3.2i)

    F'_k = int_Gamma (d_k)_i [ p' delta_ij - tau'_ij ] n_j dS .

The structure obeys the coupled second-order system

    sum_j M_kj d(phi_j)/dt + sum_j D_kj phi_j + sum_j K_kj eta_j - F'_k = 0
    d(eta_k)/dt - phi_k = 0 ,

with the GENERALISED MASS MATRIX obtained from the kinetic energy
T = 1/2 phi^T M phi,

    M_jk = int_body rho  d_j(x) . d_k(x)  dA .

Evaluated analytically from the body mass m, centre of mass (xcm, ycm) and the
moment of inertia, this gives

    translation-translation :  M_jk = m (e_j . e_k)
    rotation-rotation       :  M_kk = I_rot   (inertia about the rotation centre)
    translation-rotation    :  M_jk = m [ e_y (xcm - cx) - e_x (ycm - cy) ]
                                          (the "static unbalance" coupling)

The translation-rotation term is non-zero exactly when the centre of mass is
offset from the rotation centre -- the mechanism behind pitch/plunge (bending/
torsion) flutter.  D and K are diagonal (each DOF carries its own damper /
spring); full-matrix overrides are accepted for the rare coupled case.

Discrete implementation (matrix bordering)
------------------------------------------
The 2n structural DOFs cannot live in the dolfinx mixed space, so the fluid
operator A and weight matrix B are bordered with 2n extra rows/columns
(eta_k, phi_k for each k) and the interface-velocity rows are replaced by the
(summed) transpiration constraint.  The augmented GEVP

    A_aug q = lambda B_aug q ,   q = (u', p', eta_1, phi_1, ..., eta_n, phi_n)

is solved by the unchanged LinearSolver / SLEPc machinery.

Sign / scaling convention (unchanged from the single-DOF version)
-----------------------------------------------------------------
FELiCS builds, per transported equation,  B = (d/dt coefficient),
A = -1j*(algebraic terms).  Hence, with structural index layout
eta_k -> N+2k, phi_k -> N+2k+1:

    eta_k row :  B[eta_k, eta_k] = 1 ,   A[eta_k, phi_k] = +1j
    phi_k row :  B[phi_k, phi_j] = M_kj ,
                 A[phi_k, phi_j] = -1j*D_kj ,  A[phi_k, eta_j] = -1j*K_kj ,
                 A[phi_k, fluid] = +1j * m^(k)   (m^(k) = force/moment of DOF k)
    constraint:  u'_i + sum_k [ b^(k)_i eta_k - a^(k)_i phi_k ] = 0
                 (algebraic; B-row = 0; a^(k)=d_k, b^(k)=(d_k.grad)U0)

For n = 1 this is identical to the original single-DOF implementation.

Mean-flow base states (option FSI.MovingBaseState)
--------------------------------------------------
The formulation above assumes that the base state solves the steady equations
the operator is linearised from.  A time-averaged turbulent flow does not: it
leaves a residual  R0 = (U0.grad)U0 + grad P0 - div(2 mu S0)  (the part of the
Reynolds-stress force the frozen eddy viscosity mu = mu_lam + mu_turb does not
represent).  R0 is a force field fixed in the laboratory frame, and when the
body moves through it the linear system sees a spurious force proportional to
the displacement (a static heave produces lift).  With MovingBaseState the
residual force and the eddy viscosity are carried rigidly with the body, which
restores the exact result that a rigidly displaced body with its displaced mean
flow is a steady solution.  For every DOF j (displacement field d_j) this adds

  momentum :  eta_j * [ -(d_j.grad)R0  + r_j z x R0  - div(2 (d_j.grad mu) S0) ]
  force_k  :  eta_j * [ int_G 2 (d_j.grad mu) d_k.(S0 n) dS
                        + int_G (d_j.n)(d_k.R0) dS ]

(r_j = 1 for a rotation, 0 for a translation; n = outward normal of the fluid).
The last term is the base-stress ("added stiffness") contribution discussed by
Negi et al. (2020); it vanishes only for R0 = 0 at the wall.  For an exact
steady base flow (R0 = 0, uniform mu) all terms vanish and the original
formulation is recovered.
Option FSI.MovingBaseStateWake changes how the fields are carried downstream
of the trailing edge: rigid motion upstream, blended (smoothstep over
MovingBaseStateWakeLength) into a wake transport downstream.  The wake
transport is chosen per DOF type:
  MovingBaseStateWakeRotation    "Translate" (default): the trailing-edge
                                 displacement (wake shifted, not swung);
                                 "Frozen": none (the wake does not follow the
                                 pitch); "Rigid": rigid rotation everywhere.
  MovingBaseStateWakeTranslation "Rigid" (default): the wake follows the
                                 translation; "Frozen": it does not.
E.g. a wake that follows the heave but not the pitch: Wake = true,
WakeRotation = "Frozen", WakeTranslation = "Rigid".  The transport field xi_j
then replaces d_j in the momentum terms (with div(xi_j) and the local rotation
of xi_j); on the body xi_j = d_j, so the force rows are unchanged.
Option FSI.MovingBaseStateForceTerms ("All", "NoGeo", "NoG", "None") keeps or
drops the two force-row terms (G: moving eddy viscosity, Geo: base stress at the
displaced wall) for sensitivity tests; the momentum forcing is always kept.
Option FSI.MovingBaseStateResidual selects R0: "Discrete" (default) is the
residual of the P2 base flow under FELiCS's own momentum operator, computed
here; "File" reads C, V, T from the mean-flow file.
Derivation: documentation/GoverningEquations/fsi_moving_base_state.md

Limitations
-----------
* DIRECT modal analysis only (no adjoint yet).
* Serial execution is assumed for the bordering (a warning is issued on MPI>1).
"""

# Third party libraries
import numpy as np
from petsc4py import PETSc
from dolfinx.fem import Function, Expression, form, locate_dofs_topological
from dolfinx.fem.petsc import assemble_matrix, assemble_vector
from dolfinx.mesh import compute_midpoints
from ufl import (
    FacetNormal,
    Measure,
    TestFunction,
    TrialFunction,
    as_vector,
    dx,
    grad,
    inner,
    max_value,
    min_value,
)

# Local libraries and methods
from FELiCS.Equation.UflDecorator import UflDecorator
from FELiCS.Misc.tensorUtils import Tensor, i_dot, i_grad, i_inner, i_t, i_conj
from FELiCS.Misc.logging import Logger, log_and_raise
from FELiCS.Parameters.Config import MOVING_BASE_STATE_FIELDS

# Get the logger
logger = Logger.get_logger("felics")

# FSI.MovingBaseStateForceTerms: which moving-base-state terms enter the force rows
MBS_FORCE_TERMS = {
    "All":   "G + Geo",
    "NoGeo": "G only (no base stress at the displaced wall)",
    "NoG":   "Geo only (no moving eddy-viscosity traction)",
    "None":  "neither G nor Geo (momentum forcing only)",
}


# FSI.MovingBaseStateWakeRotation / ...WakeTranslation: transport downstream of the TE
MBS_WAKE_ROTATION = {
    "Translate": "translated by the trailing-edge displacement",
    "Frozen":    "not moved (the wake does not follow the rotation)",
    "Rigid":     "rotated rigidly, as near the body",
}
MBS_WAKE_TRANSLATION = {
    "Rigid":  "translated with the body",
    "Frozen": "not moved (the wake does not follow the translation)",
}


def _interpolation_points(V):
    """Interpolation points of the element of V (robust to dolfinx 0.8/0.9
    exposing ``interpolation_points`` as method or property)."""
    ip = V.element.interpolation_points
    return ip() if callable(ip) else ip


class _RigidBodyDOF:
    """
    A single rigid-body degree of freedom (one translation direction or one
    rotation axis), holding its geometry and its diagonal structural
    properties (generalised inertia, damping, stiffness).
    """

    def __init__(self, motion_type, inertia, damping, stiffness,
                 direction=None, center=None, name=None,
                 inertia_about="center"):
        self.motion_type   = motion_type      # "translation" or "rotation"
        self.inertia       = float(inertia)   # diagonal generalised inertia (mass or I)
        self.damping       = float(damping)   # diagonal damping
        self.stiffness     = float(stiffness) # diagonal stiffness
        self.direction     = direction        # (ex, ey) unit vector  (translation)
        self.center        = center           # (cx, cy)              (rotation)
        self.name          = name or motion_type
        # for a rotation DOF: is `inertia` given about the rotation "center"
        # or about the "cm" (centre of mass)?  In the latter case the parallel
        # axis theorem is applied later (needs body mass + centre of mass).
        self.inertia_about = inertia_about

    @property
    def label(self):
        """Short human-readable tag identifying this DOF in output files,
        e.g. 'transY', 'transX', 'trans(0.6;0.8)', 'rot(0;0)'.

        NOTE: must never contain a comma -- these labels go into CSV column
        headers, and a comma would split the header into bogus extra columns.
        Coordinates are therefore separated by ';'.
        """
        if self.motion_type == "rotation":
            cx, cy = self.center
            return f"rot({cx:g};{cy:g})"
        ex, ey = self.direction
        if abs(ex) < 1e-12 and abs(ey - 1.0) < 1e-12:
            return "transY"
        if abs(ey) < 1e-12 and abs(ex - 1.0) < 1e-12:
            return "transX"
        return f"trans({ex:.3g};{ey:.3g})"

    # -- displacement field d_k(x) --------------------------------------
    def displacement_ufl(self, x):
        """UFL vector d_k(x)."""
        if self.motion_type == "rotation":
            cx, cy = self.center
            return as_vector((-(x[1] - cy), (x[0] - cx)))
        return as_vector((self.direction[0], self.direction[1]))

    def displacement_component_callable(self, c):
        """Python callable f(x)->values for the c-component of d_k (robust for
        the constant translation direction as well as the rotation arm)."""
        if self.motion_type == "rotation":
            cx, cy = self.center
            if c == 0:
                return lambda x: -(x[1] - cy)
            return lambda x: (x[0] - cx)
        val = self.direction[c]
        return lambda x: np.full(x.shape[1], val, dtype=np.complex128)


class RigidBodyMotionFSI:
    """
    Build the FSI-augmented operators (A_aug, B_aug) for a rigid body with one
    or more rigid-body degrees of freedom, for direct modal analysis.

    Attributes
    ----------
    n_dof : int
        Number of rigid-body DOFs.
    n_struct : int
        Number of appended structural rows/columns ( = 2 * n_dof ).
    """

    def __init__(self, param, equation, meanFlow):
        self.param    = param
        self.equation = equation
        self.meanFlow = meanFlow

        # FELiCS objects
        self.FEMSpaces = equation._FEMSpaces
        self.mesh      = equation._mesh
        self.cs        = equation._coordinateSystem
        self.x         = equation.x
        self.n         = equation.n
        self.X         = equation.X

        fsi = param.FSI
        self._check_geometry(param)
        self.interface_id = int(fsi.InterfaceBoundaryID)

        # centre of mass (optional; enables inertial coupling between DOFs).
        # An empty list counts as "not supplied" (Config may not allow None).
        self.center_of_mass = None
        com_raw = getattr(fsi, "CenterOfMass", None)
        if com_raw is not None and len(list(com_raw)) >= 2:
            com = list(com_raw)
            self.center_of_mass = (float(com[0]), float(com[1]))

        # explicit body mass (used for the translation diagonal inertia and for
        # the static-unbalance coupling).  Resolve it BEFORE parsing the DOFs so
        # that a translation DOF can inherit it when it omits its own 'Mass'.
        self.body_mass = self._resolve_body_mass(fsi)

        # ---- parse the list of DOFs ---------------------------------------
        self.dofs = self._parse_dofs(fsi)
        self.n_dof    = len(self.dofs)
        self.n_struct = 2 * self.n_dof

        # final fallback for the body mass (e.g. pure-rotation body without an
        # explicit Mass): take it from a translation DOF if one exists.
        if self.body_mass is None:
            for dof in self.dofs:
                if dof.motion_type == "translation":
                    self.body_mass = dof.inertia
                    break

        # indices in the mixed state vector
        self.u_index = param.Case.SolutionList.index("u")
        self.p_index = param.Case.SolutionList.index("p")
        self.n_comp  = param.BoundaryCondition.nVelocityComponents

        # structural matrices M, D, K (n_dof x n_dof)
        self.M_mat, self.D_mat, self.K_mat = self._build_structural_matrices(fsi)

        # ---- sanity checks -------------------------------------------------
        tags = np.unique(self.mesh.facet_tags.values)
        if self.interface_id not in tags:
            log_and_raise(
                logger,
                f"FSI interface boundary ID {self.interface_id} not found in "
                f"the mesh facet tags {tags}.",
                ValueError,
            )
        if self.mesh.dolfinxMesh.comm.size > 1:
            logger.warning(
                "RigidBodyMotionFSI bordering currently assumes serial "
                "execution; MPI > 1 rank is not supported yet for the "
                "FSI-augmented system."
            )

        logger.info(
            f"FSI enabled with {self.n_dof} rigid-body DOF(s), "
            f"interface boundary = {self.interface_id}, "
            f"centre of mass = {self.center_of_mass}."
        )
        for k, dof in enumerate(self.dofs):
            if dof.motion_type == "rotation":
                logger.info(
                    f"  DOF {k}: rotation about {dof.center}, "
                    f"I = {dof.inertia}, D = {dof.damping}, K = {dof.stiffness}."
                )
            else:
                logger.info(
                    f"  DOF {k}: translation along "
                    f"({dof.direction[0]:.4f}, {dof.direction[1]:.4f}), "
                    f"m = {dof.inertia}, D = {dof.damping}, K = {dof.stiffness}."
                )
        with np.printoptions(precision=4, suppress=True):
            logger.info(f"  generalised mass matrix M =\n{self.M_mat}")
            if np.any(self.D_mat):
                logger.info(f"  damping matrix D =\n{self.D_mat}")
            if np.any(self.K_mat):
                logger.info(f"  stiffness matrix K =\n{self.K_mat}")

        # ---- mean-flow base states -----------------------------------------
        self.moving_base_state = bool(getattr(fsi, "MovingBaseState", False))
        self.static_check      = bool(getattr(fsi, "StaticCheck", False))
        self.mbs_residual      = str(getattr(fsi, "MovingBaseStateResidual", "Discrete"))
        self.mbs_wake          = bool(getattr(fsi, "MovingBaseStateWake", False))
        self.mbs_wake_length   = float(getattr(fsi, "MovingBaseStateWakeLength", 0.5))
        self.mbs_wake_rotation    = str(getattr(fsi, "MovingBaseStateWakeRotation", "Translate"))
        self.mbs_wake_translation = str(getattr(fsi, "MovingBaseStateWakeTranslation", "Rigid"))
        self.mbs_force_terms   = str(getattr(fsi, "MovingBaseStateForceTerms", "All"))
        self._trailing_edge    = None
        self._mbs_columns      = None
        self._mbs_force_matrix = np.zeros((self.n_dof, self.n_dof), dtype=complex)
        self._force_arrays     = None
        if self.mbs_wake and not self.moving_base_state:
            logger.warning(
                "  MovingBaseStateWake = true has no effect because "
                "MovingBaseState = false."
            )
        if self.moving_base_state and self.mbs_force_terms not in MBS_FORCE_TERMS:
            log_and_raise(
                logger,
                "FSI.MovingBaseStateForceTerms must be one of "
                f"{list(MBS_FORCE_TERMS)}, got '{self.mbs_force_terms}'.",
                ValueError,
            )
        if self.moving_base_state and self.mbs_force_terms != "All":
            logger.info(
                f"  MovingBaseStateForceTerms = {self.mbs_force_terms}: force rows "
                f"use {MBS_FORCE_TERMS[self.mbs_force_terms]} (sensitivity test)."
            )
        if self.moving_base_state and self.mbs_residual not in ("Discrete", "File"):
            log_and_raise(
                logger,
                "FSI.MovingBaseStateResidual must be 'Discrete' or 'File', got "
                f"'{self.mbs_residual}'.",
                ValueError,
            )
        if self.moving_base_state and self.mbs_wake:
            if self.mbs_wake_length <= 0.0:
                log_and_raise(
                    logger,
                    "FSI.MovingBaseStateWakeLength must be positive, got "
                    f"{self.mbs_wake_length}.",
                    ValueError,
                )
            for key, value, allowed in (
                    ("MovingBaseStateWakeRotation", self.mbs_wake_rotation, MBS_WAKE_ROTATION),
                    ("MovingBaseStateWakeTranslation", self.mbs_wake_translation, MBS_WAKE_TRANSLATION)):
                if value not in allowed:
                    log_and_raise(
                        logger,
                        f"FSI.{key} must be one of {list(allowed)}, got '{value}'.",
                        ValueError,
                    )
            self._trailing_edge = self._find_trailing_edge()
        if self.moving_base_state:
            logger.info(
                "  MovingBaseState = true: the base-flow residual and the eddy "
                "viscosity move rigidly with the body (mean-flow base state)."
            )
            if self.mbs_wake:
                x_te, y_te = self._trailing_edge
                logger.info(
                    "  MovingBaseStateWake = true: the fields move rigidly with "
                    f"the body for x <= {x_te:.5f} (trailing edge at "
                    f"({x_te:.5f}, {y_te:.5f})); for x >= "
                    f"{x_te + self.mbs_wake_length:.5f} (smoothstep blend) they are, "
                    f"for rotation DOFs, {MBS_WAKE_ROTATION[self.mbs_wake_rotation]} "
                    f"(WakeRotation = {self.mbs_wake_rotation}), and for translation "
                    f"DOFs, {MBS_WAKE_TRANSLATION[self.mbs_wake_translation]} "
                    f"(WakeTranslation = {self.mbs_wake_translation})."
                )
            else:
                logger.info(
                    "  MovingBaseStateWake = false: rotation DOFs rotate the "
                    "fields rigidly everywhere, including the wake."
                )
        else:
            logger.info(
                "  MovingBaseState = false: original formulation (the base "
                "state is assumed to solve the steady equations)."
            )

    # ----------------------------------------------------------------------
    # Configuration parsing
    # ----------------------------------------------------------------------
    @staticmethod
    def _check_geometry(param):
        """The displacement fields d_k(x) (in particular the rotation arm
        (-(y-cy), x-cx)) are written for planar 2-D Cartesian geometry.  In an
        axisymmetric / cylindrical or 3-D case they would be silently wrong."""
        cs   = str(getattr(param.Case, "CoordinateSystem", "Cartesian"))
        ndim = int(getattr(param.Case, "nDim", 2))
        if cs.lower() != "cartesian" or ndim != 2:
            log_and_raise(
                logger,
                f"Rigid-body FSI is implemented for 2-D Cartesian cases only "
                f"(got CoordinateSystem = '{cs}', nDim = {ndim}). The rigid-body "
                f"displacement fields (rotation arm, translation directions) are "
                f"not defined for this geometry.",
                NotImplementedError,
            )

    @staticmethod
    def _dof_get(d):
        """Return a dict-like getter for a DOF entry (plain dict or Dotdict)."""
        return (d.get if hasattr(d, "get")
                else (lambda k, default=None: getattr(d, k, default)))

    def _resolve_body_mass(self, fsi):
        """Resolve the body mass m used for the translation diagonal and the
        static-unbalance coupling, with priority:
          1. top-level FSI.Mass,
          2. the 'Mass' of the first translation DOF that states one,
          3. None (filled later, ultimately 1.0).
        """
        # Config registers FSI.Mass with the sentinel default -1.0, meaning
        # "not given in the case file" (a real mass is >= 0).  Only a value
        # the user actually wrote counts as the top-level body mass.
        top = getattr(fsi, "Mass", None)
        if top is not None and float(top) >= 0.0:
            return float(top)
        raw = getattr(fsi, "DOFs", None)
        translation_aliases = ("translation", "translationx", "translationy",
                               "horizontal", "vertical")
        if raw:
            for d in raw:
                get = self._dof_get(d)
                if not bool(get("Enabled", True)):
                    continue          # disabled DOFs contribute nothing
                if str(get("MotionType", "")).lower() in translation_aliases:
                    mval = get("Mass", None)
                    if mval is not None:
                        return float(mval)
        else:  # legacy single-DOF translation at the FSI top level
            if str(getattr(fsi, "MotionType", "")).lower() in translation_aliases:
                return 1.0          # legacy default (top-level Mass not given)
        return None

    def _parse_dofs(self, fsi):
        """Build the list of _RigidBodyDOF objects from the FSI config, accepting
        both the new multi-DOF 'DOFs' list and the legacy single-DOF style.

        Each entry of 'DOFs' may carry an "Enabled" flag (default true).  DOFs
        with "Enabled": false are skipped entirely -- they contribute no rows,
        no columns and no coupling -- so a single case file can hold all the
        candidate motions and you switch between 1-DOF / 2-DOF / 3-DOF runs by
        flipping those flags.
        """
        raw = getattr(fsi, "DOFs", None)
        if raw:
            active  = []
            skipped = []
            for d in raw:
                get = self._dof_get(d)
                if bool(get("Enabled", True)):
                    active.append(self._parse_single_dof(d))
                else:
                    skipped.append(str(get("MotionType", "?")))
            if skipped:
                logger.info(
                    f"FSI: {len(skipped)} DOF(s) disabled in the case file and "
                    f"skipped: {', '.join(skipped)}."
                )
            if not active:
                log_and_raise(
                    logger,
                    "FSI.Enabled is true but every entry of FSI.DOFs has "
                    '"Enabled": false. Enable at least one DOF, or set '
                    "FSI.Enabled = false to run the fluid-only (rigid body) case.",
                    ValueError,
                )
            return active
        # legacy single-DOF style: motion described at the FSI top level.
        # This path is also what you land on if Config.py has NOT registered the
        # 'DOFs' key -- in which case a multi-DOF case file would be silently
        # reduced to one default DOF.  Warn loudly so that can never go unnoticed.
        logger.warning(
            "FSI: no 'DOFs' list was found in the parsed configuration -- "
            "falling back to the legacy single-DOF style (MotionType / "
            "MomentOfInertia / Mass at the FSI top level). If your case file "
            "DOES define an 'FSI.DOFs' list, then Config.py has not registered "
            "the key and is silently dropping it: add 'DOFs', 'CenterOfMass' and "
            "'MomentOfInertiaCM' to the FSI block of get_all_settings_dict."
        )
        legacy = {
            "MotionType":      fsi.MotionType,
            "MotionDirection": getattr(fsi, "MotionDirection", [0.0, 1.0]),
            "RotationCenter":  getattr(fsi, "RotationCenter", [0.0, 0.0]),
            "MomentOfInertia": getattr(fsi, "MomentOfInertia", 1.0),
            "Mass":            (float(fsi.Mass)
                                if float(getattr(fsi, "Mass", -1.0)) >= 0.0
                                else 1.0),
            "Damping":         getattr(fsi, "Damping", 0.0),
            "Stiffness":       getattr(fsi, "Stiffness", 0.0),
        }
        return [self._parse_single_dof(legacy)]

    def _parse_single_dof(self, d):
        """Parse one DOF dict into a _RigidBodyDOF."""
        # dict-like access that works for plain dicts and Dotdict
        get = (d.get if hasattr(d, "get")
               else (lambda k, default=None: getattr(d, k, default)))
        motion    = str(get("MotionType", "Rotation")).lower()
        damping   = float(get("Damping", 0.0))
        stiffness = float(get("Stiffness", 0.0))

        if motion == "rotation":
            center  = get("RotationCenter", [0.0, 0.0])
            # inertia may be given about the rotation centre (MomentOfInertia)
            # or about the centre of mass (MomentOfInertiaCM); the latter is
            # converted with the parallel-axis theorem once mass + CM are known.
            i_cm = get("MomentOfInertiaCM", None)
            if i_cm is not None and float(i_cm) > 0.0:
                inertia, about = float(i_cm), "cm"
            else:
                inertia, about = float(get("MomentOfInertia", 1.0)), "center"
            return _RigidBodyDOF(
                "rotation", inertia, damping, stiffness,
                center=(float(center[0]), float(center[1])),
                name="rotation", inertia_about=about,
            )
        elif motion in ("translation", "translationx", "translationy",
                        "horizontal", "vertical"):
            if motion in ("vertical", "translationy"):
                direction = [0.0, 1.0]
            elif motion in ("horizontal", "translationx"):
                direction = [1.0, 0.0]
            else:
                direction = list(get("MotionDirection", [0.0, 1.0]))
            nrm = float(np.hypot(direction[0], direction[1]))
            if nrm == 0.0:
                log_and_raise(logger,
                              "FSI MotionDirection must be a non-zero vector.",
                              ValueError)
            direction = [direction[0] / nrm, direction[1] / nrm]
            # a translation DOF's generalised inertia IS the body mass; default
            # to the resolved body mass when this DOF omits its own 'Mass'.
            mass_raw = get("Mass", None)
            if mass_raw is not None:
                inertia = float(mass_raw)
            elif self.body_mass is not None:
                inertia = float(self.body_mass)
            else:
                inertia = 1.0
            return _RigidBodyDOF(
                "translation", inertia, damping, stiffness,
                direction=direction, name="translation",
            )
        else:
            log_and_raise(
                logger,
                f"Unknown FSI MotionType '{motion}'. Use 'Rotation' or "
                f"'Translation' (or 'Vertical'/'Horizontal').",
                NotImplementedError,
            )

    def _build_structural_matrices(self, fsi):
        """Assemble the n_dof x n_dof mass (M), damping (D) and stiffness (K)
        matrices.  Diagonals come from each DOF; the mass off-diagonals are the
        analytic generalised-mass couplings (static unbalance) when a centre of
        mass is given.  Full-matrix overrides take precedence if supplied."""
        n = self.n_dof
        M = np.zeros((n, n))
        D = np.zeros((n, n))
        K = np.zeros((n, n))

        # Normalise every rotation DOF's diagonal inertia to be about ITS OWN
        # rotation centre.  If the user supplied it about the centre of mass
        # (MomentOfInertiaCM), apply the parallel-axis theorem
        #     I_rot = I_cm + m * |x_cm - x_center|^2 .
        # This makes the centre of mass a meaningful, separate input even for a
        # single rotation DOF.
        for dof in self.dofs:
            if dof.motion_type == "rotation" and dof.inertia_about == "cm":
                if (self.center_of_mass is not None
                        and self.center_of_mass == tuple(dof.center)):
                    pass   # CM on the rotation axis: I_rot = I_cm, no mass needed
                elif self.center_of_mass is None or self.body_mass is None:
                    logger.warning(
                        "MomentOfInertiaCM given but CenterOfMass / Mass is "
                        "missing; cannot apply the parallel-axis theorem. "
                        "Using the value as if it were about the rotation centre."
                    )
                else:
                    xcm, ycm = self.center_of_mass
                    cx, cy   = dof.center
                    d2       = (xcm - cx) ** 2 + (ycm - cy) ** 2
                    dof.inertia = dof.inertia + self.body_mass * d2
                dof.inertia_about = "center"   # now expressed about the centre

        for k, dof in enumerate(self.dofs):
            M[k, k] = dof.inertia
            D[k, k] = dof.damping
            K[k, k] = dof.stiffness

        # a rigid body has ONE mass: every translation DOF must carry it
        if self.body_mass is not None:
            for dof in self.dofs:
                if (dof.motion_type == "translation"
                        and abs(dof.inertia - self.body_mass) > 1e-12 * max(1.0, abs(self.body_mass))):
                    logger.warning(
                        f"FSI: translation DOF '{dof.label}' has Mass = {dof.inertia} "
                        f"but the body mass is {self.body_mass}. A rigid body has one "
                        f"mass; the inertial couplings use {self.body_mass}. Give the "
                        f"mass once, at FSI.Mass, unless this is intentional.")

        # off-diagonal generalised mass  M_jk = int rho d_j.d_k dA
        #  * translation-translation: m (e_j . e_k)  -- needs NO centre of mass
        #  * anything involving a rotation: needs the centre of mass
        if n > 1:
            needs_com = []
            for j in range(n):
                for k in range(j + 1, n):
                    dj, dk = self.dofs[j], self.dofs[k]
                    both_trans = (dj.motion_type == "translation"
                                  and dk.motion_type == "translation")
                    if both_trans:
                        Mjk = self.body_mass * (dj.direction[0] * dk.direction[0]
                                                + dj.direction[1] * dk.direction[1])
                    elif self.center_of_mass is not None:
                        Mjk = self._mass_coupling(dj, dk)
                    else:
                        needs_com.append(f"{dj.label}-{dk.label}")
                        continue
                    M[j, k] = Mjk
                    M[k, j] = Mjk
            if needs_com:
                logger.warning(
                    "FSI: no CenterOfMass given, so the inertial coupling(s) "
                    f"{', '.join(needs_com)} are set to ZERO. That is exact only "
                    "if the centre of mass lies on the rotation axis. Set "
                    "FSI.CenterOfMass to include the static-unbalance coupling.")

        # explicit full-matrix overrides (advanced / exotic coupling)
        M = self._matrix_override(getattr(fsi, "MassMatrix",      None), M, "MassMatrix")
        D = self._matrix_override(getattr(fsi, "DampingMatrix",   None), D, "DampingMatrix")
        K = self._matrix_override(getattr(fsi, "StiffnessMatrix", None), K, "StiffnessMatrix")

        # M must be symmetric positive-definite
        if not np.allclose(M, M.T, atol=1e-10):
            logger.warning("Generalised mass matrix M is not symmetric; "
                           "symmetrising.")
            M = 0.5 * (M + M.T)
        eigM = np.linalg.eigvalsh(M)
        if np.any(eigM <= 0.0):
            logger.warning(
                f"Generalised mass matrix M is not positive-definite "
                f"(eigenvalues {eigM}). Check the inertia / centre-of-mass "
                f"inputs (e.g. I_rot must exceed m*offset^2)."
            )
        return M, D, K

    def _mass_coupling(self, dj, dk):
        """Off-diagonal generalised-mass entry M_jk = int rho d_j.d_k dA,
        evaluated analytically from the body mass and centre of mass."""
        xcm, ycm = self.center_of_mass
        m = self.body_mass
        tj, tk = dj.motion_type, dk.motion_type

        if tj == "translation" and tk == "translation":
            return m * (dj.direction[0] * dk.direction[0]
                        + dj.direction[1] * dk.direction[1])

        if tj == "translation" and tk == "rotation":
            e = dj.direction
            cx, cy = dk.center
            return m * (e[1] * (xcm - cx) - e[0] * (ycm - cy))

        if tj == "rotation" and tk == "translation":
            e = dk.direction
            cx, cy = dj.center
            return m * (e[1] * (xcm - cx) - e[0] * (ycm - cy))

        # rotation-rotation about (possibly) different centres:
        #   M_jk = I_cm + m[(xcm-cxj)(xcm-cxk) + (ycm-cyj)(ycm-cyk)]
        # with I_cm backed out from DOF j's inertia about its own centre.
        if m is None:
            logger.warning("Centre-of-mass coupling between two rotation DOFs "
                           "needs the body mass; set FSI.Mass. Using 0 coupling.")
            return 0.0
        cxj, cyj = dj.center
        cxk, cyk = dk.center
        I_cm = dj.inertia - m * ((xcm - cxj) ** 2 + (ycm - cyj) ** 2)
        return I_cm + m * ((xcm - cxj) * (xcm - cxk) + (ycm - cyj) * (ycm - cyk))

    def _matrix_override(self, raw, default, name):
        # "not supplied" may arrive as None or as an empty list/array: Config.py
        # registers these keys with a [] default (some validators reject None).
        if raw is None:
            return default
        try:
            if len(raw) == 0:
                return default
        except TypeError:
            return default          # not sized -> treat as absent
        arr = np.array(raw, dtype=float)
        if arr.shape != (self.n_dof, self.n_dof):
            log_and_raise(
                logger,
                f"FSI {name} must be {self.n_dof}x{self.n_dof} for "
                f"{self.n_dof} DOF(s); got shape {arr.shape}.",
                ValueError,
            )
        logger.info(f"Using user-supplied {name}.")
        return arr

    # ----------------------------------------------------------------------
    # Fluid operator without the interface BC (keeps interface columns)
    # ----------------------------------------------------------------------
    def _build_fluid_linear_operator(self):
        """Assemble fluid A WITHOUT the velocity BC on the FSI interface, so the
        interface-velocity columns are retained (replicates
        EquationCollection.get_linear_operator with a filtered BC list)."""
        bcs = self.equation.boundaryHandler.get_list_of_dirichlet_b_cs_for_dolfinx(
            self.FEMSpaces.VMixed,
            exclude_ids=[self.interface_id],
        )
        A_ufl = UflDecorator()
        for eq in self.equation.equationList:
            eq.add_linear_expression(A_ufl, self.meanFlow)
        return A_ufl.get_assembled_matrix(self.mesh, bcs=bcs)

    # ----------------------------------------------------------------------
    # Generalized force / moment of one DOF
    # ----------------------------------------------------------------------
    def _build_force_vector(self, dof):
        """Assemble F'_k = int_Gamma (d_k)_i (p' d_ij - tau'_ij) n_j dS as a
        functional of the fluid state, returned as a PETSc vector (size N) with
        m_l = F'_k(phi_l).  A moment for a rotation DOF, a force component for a
        translation DOF."""
        cs      = self.cs
        X_u     = i_conj(self.X[self.u_index])   # conjugated velocity test
        X_p     = i_conj(self.X[self.p_index])   # conjugated pressure test
        n       = self.n
        mean_mu = self.meanFlow.mu_tot      # dynamic viscosity (upstream: nu -> mu)
        J_hat   = cs.J_hat
        ds_if   = self.equation.ds(self.interface_id)

        d_tens  = Tensor(dof.displacement_ufl(self.x), cs)

        gradXu        = i_grad(X_u)
        strain2       = gradXu + i_t(gradXu)
        visc_traction = i_dot(strain2, n)
        moment_press  = X_p * i_dot(d_tens, n)
        moment_visc   = i_dot(d_tens, visc_traction)
        moment_integrand = moment_press - mean_mu * moment_visc

        force_ufl = UflDecorator()
        force_ufl.add(moment_integrand.ufl_tens * J_hat * ds_if)
        return force_ufl.get_assembled_vector(self.mesh)

    # ----------------------------------------------------------------------
    # Transpiration couplings (summed over DOFs) at the interface velocity DOFs
    # ----------------------------------------------------------------------
    def _build_transpiration_couplings(self):
        """For every interface velocity DOF return its constraint couplings to
        all structural DOFs:  (mixed_dof, b_list, a_list) where, for DOF k,
        b_list[k] = (d_k.grad)U0 component  and  a_list[k] = (d_k) component at
        that DOF.  The constraint row is

            u'_i + sum_k [ b_list[k] eta_k - a_list[k] phi_k ] = 0 .
        """
        V_mixed = self.FEMSpaces.VMixed
        u_index = self.u_index
        fdim    = self.mesh.gdim - 1

        facet_tags       = self.mesh.facet_tags
        interface_facets = facet_tags.indices[facet_tags.values == self.interface_id]

        # per-DOF UFL for b_k = (d_k.grad)U0
        b_ufls = []
        for dof in self.dofs:
            d_tens = Tensor(dof.displacement_ufl(self.x), self.cs)
            b_tens = i_dot(i_grad(self.meanFlow.u), d_tens)
            b_ufls.append(b_tens.ufl_tens)

        # accumulate per mixed-dof
        a_by_dof = {}
        b_by_dof = {}
        order    = []          # preserve discovery order of mixed dofs
        for c in range(self.n_comp):
            Vc, _ = V_mixed.sub(u_index).sub(c).collapse()
            ip    = _interpolation_points(Vc)

            mixed_dofs, coll_dofs = locate_dofs_topological(
                (V_mixed.sub(u_index).sub(c), Vc), fdim, interface_facets)

            # interpolate the c-component of a_k and b_k for every DOF
            a_vals = np.empty((self.n_dof, len(mixed_dofs)), dtype=complex)
            b_vals = np.empty((self.n_dof, len(mixed_dofs)), dtype=complex)
            for k, dof in enumerate(self.dofs):
                fa = Function(Vc)
                fa.interpolate(dof.displacement_component_callable(c))
                fb = Function(Vc)
                fb.interpolate(Expression(b_ufls[k][c], ip))
                a_arr = fa.x.array
                b_arr = fb.x.array
                for idx, cd in enumerate(coll_dofs):
                    a_vals[k, idx] = a_arr[cd]
                    b_vals[k, idx] = b_arr[cd]

            for idx, md in enumerate(mixed_dofs):
                md = int(md)
                a_by_dof[md] = a_vals[:, idx].copy()
                b_by_dof[md] = b_vals[:, idx].copy()
                order.append(md)

        couplings = [(md, b_by_dof[md], a_by_dof[md]) for md in order]
        return couplings

    # ----------------------------------------------------------------------
    # Public entry point
    # ----------------------------------------------------------------------
    def assemble_augmented_operators(self):
        """Return (A_aug, B_aug) of size (N + 2*n_dof), ready for
        LinearSolver.solve_general_eigenproblem."""
        A_ff = self._build_fluid_linear_operator()
        B_ff = self.equation.get_weight_matrix(self.meanFlow)
        A_ff.assemble()
        B_ff.assemble()

        # keep the fluid mass matrix and size for the energy-ratio diagnostic
        self.B_fluid = B_ff
        self.N_fluid = A_ff.getSize()[0]

        # per-DOF force/moment vectors
        force_vecs = [self._build_force_vector(dof) for dof in self.dofs]
        self._force_arrays = [v.getArray().copy() for v in force_vecs]

        # summed transpiration couplings
        couplings = self._build_transpiration_couplings()

        # mean-flow base state carried with the body
        if self.moving_base_state:
            self._build_moving_base_state(couplings)
            self._set_direct_solver_ordering()

        A_aug = self._border_A(A_ff, force_vecs, couplings)
        B_aug = self._border_B(B_ff, couplings)

        if self.static_check:
            self._static_check(A_aug, couplings)

        logger.info(
            f"FSI-augmented operators assembled: fluid size N = {self.N_fluid}, "
            f"augmented size = {A_aug.getSize()[0]} "
            f"({self.n_dof} DOF -> {self.n_struct} structural rows), "
            f"interface velocity DOFs constrained = {len(couplings)}."
        )
        return A_aug, B_aug

    def _eta_col(self, k):
        return self.N_fluid + 2 * k

    @property
    def dof_labels(self):
        """Labels of the ACTIVE DOFs, in the same order as the 'DOF k' entries in
        the log and as the (eta_k, phi_k) columns of the augmented system."""
        return [dof.label for dof in self.dofs]

    def _phi_col(self, k):
        return self.N_fluid + 2 * k + 1

    # ----------------------------------------------------------------------
    # Bordering helper: build the augmented matrix in ONE shot
    # ----------------------------------------------------------------------
    def _bordered_matrix(self, M_ff, zero_rows, diag, extra_rows, extra_cols,
                         extra_vals):
        """Return the (N + n_struct) square bordered copy of M_ff, with

          * the rows `zero_rows` emptied and, if diag != 0, a `diag` put on
            their diagonal (what Mat.zeroRows(rows, diag) does), then
          * the entries (extra_rows, extra_cols) = extra_vals INSERTED
            (a later entry overwrites an earlier one at the same position,
            like Mat.setValue with INSERT_VALUES),

        assembled ONCE from a CSR built with numpy.  Setting the entries one
        by one into an already assembled PETSc matrix (the previous code)
        makes PETSc reallocate - i.e. copy - the whole matrix for nearly
        every new nonzero, which takes hours for N ~ 1e6.

        Zeroed rows and structural rows are rebuilt from the inserted entries
        alone; kept fluid rows may receive entries in structural columns only
        (c >= N), appended after their own entries; all other rows are copied
        block-wise from M_ff (no global sort, ~1 extra copy of the matrix in
        memory)."""
        N     = M_ff.getSize()[0]
        Nnew  = N + self.n_struct
        comm  = M_ff.getComm()
        if comm.getSize() > 1:
            log_and_raise(logger, "FSI bordering is serial only (MPI size 1).",
                          NotImplementedError)

        ai, aj, av = M_ff.getValuesCSR()
        ai = np.asarray(ai, dtype=np.int64)

        zero_rows = np.unique(np.asarray(zero_rows, dtype=np.int64))
        is_zero   = np.zeros(Nnew, dtype=bool)
        is_zero[zero_rows] = True

        # the small set of entries that rows get rebuilt from: the zeroRows
        # diagonal first, then the inserted entries (later ones win)
        er = [zero_rows] if diag != 0.0 else []
        ec = [zero_rows] if diag != 0.0 else []
        ev = ([np.full(len(zero_rows), diag, dtype=PETSc.ScalarType)]
              if diag != 0.0 else [])
        er.append(np.asarray(extra_rows, dtype=np.int64))
        ec.append(np.asarray(extra_cols, dtype=np.int64))
        ev.append(np.asarray(extra_vals, dtype=PETSc.ScalarType))
        r = np.concatenate(er); c = np.concatenate(ec); v = np.concatenate(ev)
        order = np.lexsort((np.arange(len(r)), c, r))
        r, c, v = r[order], c[order], v[order]
        last = np.ones(len(r), dtype=bool)
        last[:-1] = (r[1:] != r[:-1]) | (c[1:] != c[:-1])
        r, c, v = r[last], c[last], v[last]

        # An inserted entry in a fluid row that is kept must lie in a
        # structural column (c >= N): it is then appended after the row's
        # own entries, which all have c < N, and the row stays sorted.
        kept = ~(is_zero[r] | (r >= N))
        if np.any(c[kept] < N):
            log_and_raise(
                logger,
                "FSI bordering: an inserted entry falls in a fluid column of a "
                "fluid row that is not zeroed - not supported by the fast "
                "bordering.",
                RuntimeError,
            )

        # row lengths of the result
        counts = np.zeros(Nnew, dtype=np.int64)
        counts[:N] = np.diff(ai)
        counts[is_zero] = 0
        counts += np.bincount(r, minlength=Nnew)
        indptr = np.zeros(Nnew + 1, dtype=np.int64)
        np.cumsum(counts, out=indptr[1:])
        nnz = int(indptr[-1])
        indices = np.empty(nnz, dtype=PETSc.IntType)
        data    = np.empty(nnz, dtype=PETSc.ScalarType)

        # rebuilt rows (zeroed or structural): sorted, from (r, c, v)
        special = np.union1d(zero_rows, np.unique(r))
        rstart  = np.searchsorted(r, special, side="left")
        rend    = np.searchsorted(r, special, side="right")

        # copy the untouched fluid rows block-wise, write the special rows
        # (kept fluid rows with appended structural columns: own entries first)
        prev = 0
        for s, a0, a1 in zip(special.tolist(), rstart.tolist(), rend.tolist()):
            if s > prev:                            # untouched rows prev..s-1
                hi = min(s, N)
                if hi > prev:
                    indices[indptr[prev]:indptr[hi]] = aj[ai[prev]:ai[hi]]
                    data[indptr[prev]:indptr[hi]] = av[ai[prev]:ai[hi]]
            start = indptr[s]
            if s < N and not is_zero[s]:
                nown = ai[s + 1] - ai[s]
                indices[start:start + nown] = aj[ai[s]:ai[s + 1]]
                data[start:start + nown] = av[ai[s]:ai[s + 1]]
                start += nown
            indices[start:indptr[s + 1]] = c[a0:a1]
            data[start:indptr[s + 1]] = v[a0:a1]
            prev = s + 1
        if prev < N:
            indices[indptr[prev]:indptr[N]] = aj[ai[prev]:ai[N]]
            data[indptr[prev]:indptr[N]] = av[ai[prev]:ai[N]]
        del ai, aj, av

        M = PETSc.Mat().createAIJ(
            (Nnew, Nnew),
            csr=(indptr.astype(PETSc.IntType), indices, data),
            comm=comm,
        )
        M.assemble()
        return M

    # ----------------------------------------------------------------------
    # Bordering of the linear operator A
    # ----------------------------------------------------------------------
    def _border_A(self, A_ff, force_vecs, couplings):
        er, ec, ev = [], [], []

        # --- transpiration constraint: overwrite interface velocity rows ----
        #     (rows emptied, 1 on the diagonal, + structural couplings)
        interface_dofs = np.array([cp[0] for cp in couplings], dtype=np.int64)
        for (row, b_list, a_list) in couplings:
            for k in range(self.n_dof):
                er.append(row); ec.append(self._eta_col(k)); ev.append(b_list[k])    # + b^(k) eta_k
                er.append(row); ec.append(self._phi_col(k)); ev.append(-a_list[k])   # - a^(k) phi_k

        # --- structural ODE rows (A = -1j * algebraic terms) ----------------
        fr, fc, fv = [], [], []
        for k in range(self.n_dof):
            # eta_k row :  d(eta_k)/dt - phi_k = 0
            er.append(self._eta_col(k)); ec.append(self._phi_col(k)); ev.append(1j)
            # phi_k row :  sum_j M_kj d(phi_j)/dt + D_kj phi_j + K_kj eta_j - F'_k = 0
            # F'_k contains  sum_j C_kj eta_j  from the moving base state
            # (C = 0 otherwise), hence one combined entry per (k, j)
            for j in range(self.n_dof):
                a_eta = -1j * self.K_mat[k, j] + 1j * self._mbs_force_matrix[k, j]
                if a_eta != 0.0:
                    er.append(self._phi_col(k)); ec.append(self._eta_col(j))
                    ev.append(a_eta)
                if self.D_mat[k, j] != 0.0:
                    er.append(self._phi_col(k)); ec.append(self._phi_col(j))
                    ev.append(-1j * self.D_mat[k, j])
            # force/moment coupling to the fluid:  -1j*(-F'_k) = +1j*m^(k)
            m_arr = force_vecs[k].getArray()
            nz = np.nonzero(m_arr)[0]
            fr.append(np.full(len(nz), self._phi_col(k), dtype=np.int64))
            fc.append(nz.astype(np.int64))
            fv.append(1j * m_arr[nz])

        # --- moving base state: forcing  eta_j * f_j  in the momentum rows
        #     (A = -1j * algebraic terms, the forcing is on the right-hand
        #     side, hence +1j).  Interface rows are excluded: they hold the
        #     transpiration constraint.
        if self._mbs_columns is not None:
            for j, g in enumerate(self._mbs_columns):
                nz = np.flatnonzero(self._mbs_rows_mask(g, interface_dofs))
                fr.append(nz.astype(np.int64))
                fc.append(np.full(len(nz), self._eta_col(j), dtype=np.int64))
                fv.append(1j * g[nz])

        rows = np.concatenate([np.asarray(er, dtype=np.int64)] + fr)
        cols = np.concatenate([np.asarray(ec, dtype=np.int64)] + fc)
        vals = np.concatenate([np.asarray(ev, dtype=PETSc.ScalarType)]
                              + [np.asarray(x, dtype=PETSc.ScalarType) for x in fv])
        return self._bordered_matrix(A_ff, interface_dofs, 1.0, rows, cols, vals)

    # ----------------------------------------------------------------------
    # Mean-flow base state carried with the body (FSI.MovingBaseState)
    # ----------------------------------------------------------------------
    def _base_residual(self):
        """Residual R0 of the base flow under the frozen-viscosity momentum
        equations, as a Tensor, from the source set by
        FSI.MovingBaseStateResidual ("Discrete" or "File")."""
        if self.mbs_residual == "File":
            return self._base_residual_from_file()
        return self._base_residual_discrete()

    def _base_residual_discrete(self):
        """Residual of the P2 base flow under FELiCS's own momentum operator,

            (R, v) = int [(grad U0)(rho U0)] . v + tau0 : grad v - P0 div v dx
                     + int_dOmega (P0 n - tau0 n) . v dS ,   tau0 = mu (grad U0 + grad U0^T),

        for all P2 vector test functions v (no boundary conditions), i.e. the
        L2 projection of the strong residual including the inter-element flux
        jumps.  This is the force the discrete linear operator "sees", so moving
        it keeps the discrete displaced base state in equilibrium; the residual
        of the continuous mean flow (e.g. from DNS statistics) can differ from
        it by O(1) in an under-resolved boundary layer, where it is the small
        difference of large convective and viscous terms."""
        velocity = self.meanFlow.u.ufl_tens
        space    = velocity.function_space
        domain   = space.mesh
        pressure = self.meanFlow.p.ufl_tens
        density  = self.meanFlow.rho.ufl_tens
        mu_total = self.meanFlow.mu_tot.ufl_tens
        normal   = FacetNormal(domain)
        dx_base  = Measure("dx", domain=domain)
        ds_base  = Measure("ds", domain=domain)
        trial, test = TrialFunction(space), TestFunction(space)
        grad_u = grad(velocity)
        tau = mu_total * (grad_u + grad_u.T)
        weak_residual = (
            inner(grad_u * (density * velocity), test)
            + inner(tau, grad(test))
            - inner(pressure, test[0].dx(0) + test[1].dx(1))
        ) * dx_base + inner(pressure * normal - tau * normal, test) * ds_base
        mass = assemble_matrix(form(inner(trial, test) * dx_base))
        mass.assemble()
        rhs = assemble_vector(form(weak_residual))
        rhs.assemble()
        ksp = PETSc.KSP().create(mass.getComm())
        ksp.setOperators(mass)
        ksp.setType(PETSc.KSP.Type.PREONLY)
        ksp.getPC().setType(PETSc.PC.Type.LU)
        try:
            ksp.getPC().setFactorSolverType("mumps")
            ksp.setUp()
        except PETSc.Error:
            ksp.getPC().setFactorSolverType("petsc")
            ksp.setUp()
        residual = Function(space)
        ksp.solve(rhs, residual.x.petsc_vec)
        residual.x.scatter_forward()
        ksp.destroy()
        mass.destroy()
        rhs.destroy()
        values = residual.x.array
        logger.info(
            "  base-flow residual R0: discrete residual of the P2 base flow "
            f"(MovingBaseStateResidual = Discrete), nodal rms "
            f"{np.sqrt(np.mean(np.abs(values)**2)):.3e}, max {np.abs(values).max():.3e}"
        )
        return Tensor(residual, self.cs)

    def _base_residual_from_file(self):
        """Residual of the base flow under the frozen-viscosity equations,
        R0 = C - mu_lam V - T (T only with TurbulenceModel "File"), as a
        Tensor.  The pieces come from the mean-flow file (computed there from
        the DNS/LES statistics with the same eddy viscosity)."""
        missing = [name for name in MOVING_BASE_STATE_FIELDS
                   if not self.meanFlow.was_read_from_file(name)]
        if missing:
            log_and_raise(
                logger,
                "FSI.MovingBaseStateResidual = 'File' needs the base-flow "
                f"residual fields {MOVING_BASE_STATE_FIELDS} in the mean-flow "
                f"file, but {missing} are missing.  Write them with the "
                "mean-flow pre-processing (make_felics_baseflow.py, "
                "WRITE_FSI_RESIDUAL) or use MovingBaseStateResidual = "
                "'Discrete'.",
                ValueError,
            )
        if self.meanFlow.was_read_from_file("rho"):
            logger.warning(
                "  FSI.MovingBaseState: the mean-flow file has a 'rho' field, "
                "but the residual fields assume rho = 1 (incompressible, "
                "non-dimensional).  Make sure they were computed with the "
                "same density."
            )
        f = {name: self.meanFlow.get_field_tensor(name).ufl_tens
             for name in MOVING_BASE_STATE_FIELDS}
        mu_lam = self.meanFlow.get_field_tensor("mulam").ufl_tens
        use_turb = self.param.Case.TurbulenceModel == "File"
        r_x = f["fsiCx"] - mu_lam * f["fsiVx"]
        r_y = f["fsiCy"] - mu_lam * f["fsiVy"]
        if use_turb:
            r_x = r_x - f["fsiTx"]
            r_y = r_y - f["fsiTy"]
        logger.info(
            "  base-flow residual R0 from the mean-flow file = C - mu_lam V"
            + (" - T  (TurbulenceModel File)" if use_turb else
               "  (TurbulenceModel not File: no eddy-viscosity part)")
        )
        return Tensor(as_vector((r_x, r_y)), self.cs)

    def _base_strain(self):
        """Mean strain-rate tensor S0 = (grad U0 + grad U0^T)/2."""
        grad_u = i_grad(self.meanFlow.u)
        return (grad_u + i_t(grad_u)) * 0.5

    def _build_moving_base_state(self, couplings):
        """Momentum columns f_j and force-row matrix C_kj of the moving base
        state; see the module docstring."""
        residual = self._base_residual()
        self._mbs_columns = [
            self._moving_base_state_vector(dof, residual) for dof in self.dofs
        ]
        G, Geo = self._moving_base_state_force_matrices(residual)
        self._mbs_force_matrix = self._selected_mbs_force_matrix(G, Geo)
        interface_dofs = np.array([cp[0] for cp in couplings], dtype=np.int64)
        for j, g in enumerate(self._mbs_columns):
            kept = g[self._mbs_rows_mask(g, interface_dofs)]
            logger.info(
                f"  moving base state, DOF {j}: momentum forcing column "
                f"||f_j|| = {np.linalg.norm(kept):.4e} "
                f"({np.count_nonzero(kept)} entries)"
            )
        with np.printoptions(precision=5, suppress=False):
            logger.info(
                "  moving base state, force rows (per unit displacement):\n"
                f"    eddy-viscosity transport  G   =\n{G.real}\n"
                f"    wall residual (added stiffness) Geo =\n{Geo.real}\n"
                f"    used in the force rows: {MBS_FORCE_TERMS[self.mbs_force_terms]}\n"
                f"    (structural stiffness K =\n{self.K_mat})"
            )

    def _find_trailing_edge(self):
        """Trailing edge = midpoint of the FSI-interface facet with the
        largest x (to within half a wall-cell size)."""
        mesh  = getattr(self.mesh, "dolfinxMesh", self.mesh)
        fdim  = mesh.topology.dim - 1
        facet_tags = self.mesh.facet_tags
        facets = facet_tags.indices[facet_tags.values == self.interface_id]
        midpoints = compute_midpoints(mesh, fdim, facets)
        index = int(np.argmax(midpoints[:, 0]))
        return float(midpoints[index, 0]), float(midpoints[index, 1])

    def _transport_field(self, dof):
        """Field xi_j(x) that carries the base-flow residual and the eddy
        viscosity for DOF j, its local rotation a_j = (curl xi_j)_z / 2 and its
        divergence (None where identically zero).

        Rigid (no MovingBaseStateWake, or the wake mode "Rigid"):
        translation xi = d_j (a = None, div = None), rotation xi = d_j, a = 1.
        With MovingBaseStateWake:

            xi = w d_j + (1 - w) d_wake,
            w  = 1 - smoothstep((x - x_TE) / L),   L = MovingBaseStateWakeLength,

        i.e. rigid motion upstream of the trailing edge and d_wake downstream:
        d_j(x_TE) for a rotation with WakeRotation = "Translate", 0 for
        WakeRotation / WakeTranslation = "Frozen".  xi = d_j on the body."""
        d = dof.displacement_ufl(self.x)
        is_rotation = dof.motion_type == "rotation"
        rigid = (d, 1.0, None) if is_rotation else (d, None, None)
        if not (self.moving_base_state and self.mbs_wake):
            return rigid
        mode = self.mbs_wake_rotation if is_rotation else self.mbs_wake_translation
        if mode == "Rigid":
            return rigid
        if mode == "Translate":
            x_te, y_te = self._trailing_edge
            cx, cy = dof.center
            d_wake = as_vector((-(y_te - cy), x_te - cx))
        else:                                    # "Frozen"
            d_wake = as_vector((0.0, 0.0))
        x_te = self._trailing_edge[0]
        ramp = (self.x[0] - x_te) / self.mbs_wake_length
        ramp = max_value(0.0, min_value(1.0, ramp))
        weight = 1.0 - (3.0 * ramp**2 - 2.0 * ramp**3)
        xi = weight * d + (1.0 - weight) * d_wake
        rotation = 0.5 * (xi[1].dx(0) - xi[0].dx(1))
        divergence = xi[0].dx(0) + xi[1].dx(1)
        return xi, rotation, divergence

    def _set_direct_solver_ordering(self):
        """The moving-base-state columns are dense (the base-flow residual is
        non-zero almost everywhere).  MUMPS' automatic ordering choice may
        fail on such a matrix (PORD crashed with MUMPS 5.6), so use QAMD,
        which treats quasi-dense rows/columns separately, for the shift-and-
        invert factorisation - unless the user chose an ordering with
        -st_mat_mumps_icntl_7 (e.g. via PETSC_OPTIONS)."""
        options = PETSc.Options()
        key = "st_mat_mumps_icntl_7"
        if options.hasName(key):
            logger.info(
                f"  MUMPS ordering for the eigensolver set by the user: "
                f"-{key} {options.getString(key)}"
            )
            return
        options.setValue(key, 6)
        logger.info(
            "  MUMPS ordering for the eigensolver set to QAMD "
            f"(-{key} 6) because of the dense moving-base-state columns; "
            "override with PETSC_OPTIONS if needed."
        )

    def _mbs_rows_mask(self, g, interface_dofs):
        """Rows of a momentum forcing column that enter the matrix: nonzero
        and not an interface (transpiration-constraint) row."""
        mask = g != 0.0
        mask[interface_dofs] = False
        return mask

    def _selected_mbs_force_matrix(self, G, Geo):
        """Force-row matrix of the moving base state for the terms selected by
        FSI.MovingBaseStateForceTerms (default "All": G + Geo)."""
        terms = getattr(self, "mbs_force_terms", "All")
        use_g = terms in ("All", "NoGeo")
        use_geo = terms in ("All", "NoG")
        return (G if use_g else np.zeros_like(G)) + (Geo if use_geo else np.zeros_like(Geo))

    def _moving_base_state_force_matrices(self, residual):
        """G_kj = int_G 2 (d_j.grad mu) d_k.(S0 n) dS  (perturbation of the
        viscous traction by the moving eddy viscosity) and
        Geo_kj = int_G (d_j.n)(d_k.R0) dS  (base stress at the displaced wall,
        zero for an exact base flow).  n: outward normal of the fluid."""
        n_dof   = self.n_dof
        G       = np.zeros((n_dof, n_dof), dtype=complex)
        Geo     = np.zeros((n_dof, n_dof), dtype=complex)
        J_hat   = self.cs.J_hat
        ds_if   = self.equation.ds(self.interface_id)
        mean_mu = self.meanFlow.mu_tot
        strain  = self._base_strain()
        traction = i_dot(strain, self.n)
        for k, dof_k in enumerate(self.dofs):
            d_k = Tensor(dof_k.displacement_ufl(self.x), self.cs)
            for j, dof_j in enumerate(self.dofs):
                d_j = Tensor(dof_j.displacement_ufl(self.x), self.cs)
                mu_shift = i_dot(i_grad(mean_mu), d_j)
                g_form = UflDecorator()
                g_form.add((mu_shift * i_dot(d_k, traction) * 2.0).ufl_tens
                           * J_hat * ds_if)
                G[k, j] = g_form.get_assembled_scalar(self.mesh)
                geo_form = UflDecorator()
                geo_form.add((i_dot(d_j, self.n) * i_dot(d_k, residual)).ufl_tens
                             * J_hat * ds_if)
                Geo[k, j] = geo_form.get_assembled_scalar(self.mesh)
        return G, Geo

    def _moving_base_state_vector(self, dof, residual):
        """Weak form of the momentum forcing per unit eta_j,
        int f_j . conj(v) dx  with transport field xi_j (= d_j unless
        MovingBaseStateWake changes it, see _transport_field):

          f_j = -(xi_j.grad)R0 + a_j z x R0 - div(2 (xi_j.grad mu) S0),

        integrated by parts as

          int R0 . (xi_j.grad)conj(v) dx + int div(xi_j) R0 . conj(v) dx
          - int_dOmega (xi_j.n)(R0 . conj(v)) dS
          + int a_j (z x R0) . conj(v) dx + int 2 (xi_j.grad mu) S0 : grad conj(v) dx

        (the viscous boundary term is dropped like in the momentum equation
        itself).  Dirichlet rows are zero."""
        J_hat    = self.cs.J_hat
        test_u   = i_conj(self.X[self.u_index])
        grad_v   = i_grad(test_u)

        # ufl inner() conjugates its second argument itself, hence the
        # unconjugated test function there (as in the momentum equation)
        grad_v_plain = i_grad(self.X[self.u_index])
        xi_ufl, rotation, divergence = self._transport_field(dof)
        xi       = Tensor(xi_ufl, self.cs)
        mu_shift = i_dot(i_grad(self.meanFlow.mu_tot), xi)
        strain   = self._base_strain()

        volume = i_dot(residual, i_dot(grad_v, xi))
        volume = volume + i_inner(strain * mu_shift * 2.0, grad_v_plain)
        if divergence is not None:
            volume = volume + Tensor(divergence, self.cs) * i_dot(residual, test_u)
        if rotation is not None:
            z_cross_r = Tensor(
                as_vector((-residual.ufl_tens[1], residual.ufl_tens[0])),
                self.cs,
            )
            rotation_term = i_dot(z_cross_r, test_u)
            if isinstance(rotation, float):
                volume = volume + rotation_term * rotation
            else:
                volume = volume + Tensor(rotation, self.cs) * rotation_term
        boundary = i_dot(xi, self.n) * i_dot(residual, test_u) * (-1.0)

        g_form = UflDecorator()
        g_form.add(volume.ufl_tens * J_hat * dx)
        g_form.add(boundary.ufl_tens * J_hat * self.equation.all_ds)
        bcs = self.equation.boundaryHandler.get_list_of_dirichlet_b_cs_for_dolfinx(
            self.FEMSpaces.VMixed,
            exclude_ids=[self.interface_id],
        )
        return g_form.get_assembled_vector(self.mesh, bcs=bcs).getArray().copy()

    # ----------------------------------------------------------------------
    # Diagnostics
    # ----------------------------------------------------------------------
    def _static_check(self, A_aug, couplings):
        """Log the static aerodynamic force matrices of
        _static_force_matrices: for a pure translation in a uniform stream
        the physical answer is (close to) zero."""
        for label, force in self._static_force_matrices(A_aug, couplings).items():
            with np.printoptions(precision=5, suppress=False):
                logger.info(
                    f"  STATIC CHECK, {label}: aerodynamic force on DOF k per "
                    f"unit static displacement of DOF j (rows k, columns j; "
                    f"an aerodynamic stiffness is minus this):\n{force.real}"
                    f"\n    imaginary part (should be ~0): "
                    f"{np.abs(force.imag).max():.2e}"
                )

    def _static_force_matrices(self, A_aug, couplings):
        """Solve the steady (omega = 0) fluid problem for a unit static
        displacement of each DOF and return the aerodynamic force on every
        DOF, for the frozen base state (original formulation) and, if
        enabled, the moving base state.  One extra LU factorisation.

        Returns
        -------
        dict
            label -> (n_dof x n_dof) complex array, entry [k, j] = force on
            DOF k for eta_j = 1.
        """
        N    = self.N_fluid
        comm = A_aug.getComm()
        index_set = PETSc.IS().createGeneral(
            np.arange(N, dtype=PETSc.IntType), comm=comm)
        A_fluid = A_aug.createSubMatrix(index_set, index_set)
        ksp = PETSc.KSP().create(comm)
        ksp.setOperators(A_fluid)
        ksp.setType(PETSc.KSP.Type.PREONLY)
        ksp.getPC().setType(PETSc.PC.Type.LU)
        try:
            ksp.getPC().setFactorSolverType("mumps")
            ksp.setUp()
        except PETSc.Error:
            logger.warning(
                "  STATIC CHECK: MUMPS not available, using the PETSc LU "
                "(slow for large meshes)."
            )
            ksp.getPC().setFactorSolverType("petsc")
            ksp.setUp()
        interface_dofs = np.array([cp[0] for cp in couplings], dtype=np.int64)

        cases = [("frozen base state (original formulation)", False)]
        if self.moving_base_state:
            cases.append(("moving base state", True))
        result = {}
        for label, moving in cases:
            force = np.zeros((self.n_dof, self.n_dof), dtype=complex)
            for j in range(self.n_dof):
                column = np.zeros(N, dtype=complex)
                for (row, b_list, a_list) in couplings:
                    column[row] = b_list[j]
                if moving:
                    g = self._mbs_columns[j]
                    mask = self._mbs_rows_mask(g, interface_dofs)
                    column[mask] += 1j * g[mask]
                rhs = A_fluid.createVecLeft()
                rhs.setArray(-column)
                sol = A_fluid.createVecRight()
                ksp.solve(rhs, sol)
                q_fluid = sol.getArray()
                for k in range(self.n_dof):
                    force[k, j] = self._force_arrays[k] @ q_fluid
                    if moving:
                        force[k, j] += self._mbs_force_matrix[k, j]
            result[label] = force
        ksp.destroy()
        A_fluid.destroy()
        index_set.destroy()
        return result

    def aerodynamic_force_diagnostics(self, eigen_value, q_aug):
        """Fluid force on every DOF in one mode, and its split into an
        apparent aerodynamic stiffness and damping.

        With perturbations ~ exp(-1j lambda t), s = -1j lambda and phi = s eta.
        F_k = m^(k) . q_fluid (+ sum_j C_kj eta_j for the moving base state).
        The structural equation  sum_j (M_kj s phi_j + D_kj phi_j + K_kj eta_j)
        = F_k  must hold to solver accuracy; its relative residual is returned.
        Writing  F_k / eta_k = -k_a - c_a s  with real k_a, c_a  gives the
        apparent (per-DOF, cross-coupling included) aerodynamic stiffness and
        damping; for a single DOF k_a = m |lambda|^2 - K, c_a = -2 m Im(lambda) - D.

        Returns
        -------
        dict
            ``force``, ``ratio`` (F_k/eta_k), ``k_aero``, ``c_aero`` (lists,
            one entry per DOF), ``residual`` and ``phi_minus_s_eta``.
        """
        s       = -1j * complex(eigen_value)
        q_fluid = np.asarray(q_aug[:self.N_fluid])
        eta = np.array([q_aug[self._eta_col(k)] for k in range(self.n_dof)],
                       dtype=complex)
        phi = np.array([q_aug[self._phi_col(k)] for k in range(self.n_dof)],
                       dtype=complex)
        force = np.array([
            self._force_arrays[k] @ q_fluid + self._mbs_force_matrix[k] @ eta
            for k in range(self.n_dof)
        ])
        lhs = self.M_mat @ (s * phi) + self.D_mat @ phi + self.K_mat @ eta
        scale = max(np.linalg.norm(force), np.linalg.norm(lhs), 1e-300)
        residual = float(np.linalg.norm(lhs - force) / scale)
        phi_err = float(np.linalg.norm(phi - s * eta)
                        / max(np.linalg.norm(phi), 1e-300))
        ratio, k_aero, c_aero = [], [], []
        for k in range(self.n_dof):
            if abs(eta[k]) == 0.0 or s.imag == 0.0:
                ratio.append(complex("nan")); k_aero.append(float("nan"))
                c_aero.append(float("nan"))
                continue
            z = force[k] / eta[k]
            c = -z.imag / s.imag
            ratio.append(z)
            c_aero.append(float(c))
            k_aero.append(float(-z.real - c * s.real))
        return dict(force=list(force), ratio=ratio, k_aero=k_aero,
                    c_aero=c_aero, residual=residual, phi_minus_s_eta=phi_err)

    # ----------------------------------------------------------------------
    # Bordering of the weight matrix B
    # ----------------------------------------------------------------------
    def _border_B(self, B_ff, couplings):
        # interface-velocity rows are algebraic constraints -> zero in B
        interface_dofs = np.array([cp[0] for cp in couplings], dtype=np.int64)

        # structural mass rows
        er, ec, ev = [], [], []
        for k in range(self.n_dof):
            er.append(self._eta_col(k)); ec.append(self._eta_col(k)); ev.append(1.0)   # d(eta_k)/dt
            for j in range(self.n_dof):                                               # M_kj d(phi_j)/dt
                if self.M_mat[k, j] != 0.0:
                    er.append(self._phi_col(k)); ec.append(self._phi_col(j))
                    ev.append(self.M_mat[k, j])
        return self._bordered_matrix(B_ff, interface_dofs, 0.0, er, ec, ev)

    # ----------------------------------------------------------------------
    # Energy-weighted structural weight of a mode
    # ----------------------------------------------------------------------
    def dof_participation(self, q_aug):
        """Per-DOF energy participation fractions of a mode.

        Returns (p, phase_deg) with

            p_k = [ M_kk |phi_k|^2 + K_kk |eta_k|^2 ]
                  / sum_j [ M_jj |phi_j|^2 + K_jj |eta_j|^2 ]

        i.e. the share of the mode's structural energy carried by DOF k.  Unlike a
        raw amplitude ratio |eta_1|/|eta_0| -- which mixes a length with an angle
        and is therefore dimensionally inconsistent -- p_k is dimensionless, lies
        in [0,1], sums to 1, and is independent of the eigenvector normalisation.
        It is the meaningful way to say "this mode is plunge-dominated" vs
        "rotation-dominated" vs "hybrid":

            p_k ~ 1     -> the mode is essentially DOF k alone
            p_k ~ 1/n   -> the mode is genuinely hybrid  (flutter lives here)

        Only the DIAGONAL entries of M and K are used, so each p_k >= 0 and the
        fractions partition unity.  (The off-diagonal static-unbalance energy is a
        shared, cross-DOF quantity and cannot be attributed to one DOF alone.)

        phase_deg[k] is arg(eta_k) - arg(eta_0) in degrees, wrapped to [-180, 180).
        It is the phase lead/lag of DOF k relative to the first DOF, and is the
        quantity that governs the direction of energy transfer between the DOFs
        (for a free body, plunge/pitch phase near +/-90 deg is the classic
        flutter-favourable condition).  phase_deg[0] is 0 by construction.

        If the mode has no structural content at all (a pure fluid mode), the
        fractions are returned as zeros.
        """
        eta = np.array([q_aug[self._eta_col(k)] for k in range(self.n_dof)],
                       dtype=complex)
        phi = np.array([q_aug[self._phi_col(k)] for k in range(self.n_dof)],
                       dtype=complex)

        e_k = np.array(
            [float(self.M_mat[k, k] * abs(phi[k])**2
                   + self.K_mat[k, k] * abs(eta[k])**2)
             for k in range(self.n_dof)]
        )
        tot = float(e_k.sum())
        p = (e_k / tot) if tot > 0.0 else np.zeros(self.n_dof)

        # phase of each DOF relative to DOF 0
        phase = np.zeros(self.n_dof)
        if abs(eta[0]) > 0.0:
            ref = np.angle(eta[0])
            for k in range(self.n_dof):
                if abs(eta[k]) > 0.0:
                    phase[k] = np.degrees(
                        (np.angle(eta[k]) - ref + np.pi) % (2*np.pi) - np.pi
                    )
        return p.tolist(), phase.tolist()

    def structural_energy_ratio(self, q_aug):
        """Energy-based structural weight  sqrt(E_s)/sqrt(E_f + E_s) with

            E_f = u^H M_u u                          (fluid kinetic energy, x2)
            E_s = phi^H M phi + eta^H K eta          (struct. kinetic + potential, x2)

        using the full coupled M and K.  Ratio in [0,1], normalisation-independent.
        """
        N   = self.N_fluid
        uf  = np.ascontiguousarray(q_aug[:N], dtype=complex)
        eta = np.array([q_aug[self._eta_col(k)] for k in range(self.n_dof)], dtype=complex)
        phi = np.array([q_aug[self._phi_col(k)] for k in range(self.n_dof)], dtype=complex)

        # fluid kinetic energy via the velocity mass matrix
        v  = self.B_fluid.createVecRight()
        v.setArray(uf)
        Bv = self.B_fluid.createVecLeft()
        self.B_fluid.mult(v, Bv)
        e_fluid = max(float(np.real(np.vdot(uf, Bv.getArray()))), 0.0)

        # structural energy with the coupled matrices
        e_struct = float(np.real(np.vdot(phi, self.M_mat @ phi)
                                 + np.vdot(eta, self.K_mat @ eta)))
        e_struct = max(e_struct, 0.0)

        e_total = e_fluid + e_struct
        return float(np.sqrt(e_struct / e_total)) if e_total > 0.0 else 0.0
