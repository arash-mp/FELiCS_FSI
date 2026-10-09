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

Limitations
-----------
* DIRECT modal analysis only (no adjoint yet).
* Serial execution is assumed for the bordering (a warning is issued on MPI>1).
"""

# Third party libraries
import numpy as np
from petsc4py import PETSc
from dolfinx.fem import Function, Expression, locate_dofs_topological
from ufl import as_vector

# Local libraries and methods
from FELiCS.Equation.UflDecorator import UflDecorator
from FELiCS.Misc.tensorUtils import Tensor, i_dot, i_grad, i_t, i_conj
from FELiCS.Misc.logging import Logger, log_and_raise

# Get the logger
logger = Logger.get_logger("felics")


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

        # summed transpiration couplings
        couplings = self._build_transpiration_couplings()

        A_aug = self._border_A(A_ff, force_vecs, couplings)
        B_aug = self._border_B(B_ff, couplings)

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
    # Bordering of the linear operator A
    # ----------------------------------------------------------------------
    def _border_A(self, A_ff, force_vecs, couplings):
        N     = A_ff.getSize()[0]
        ns    = self.n_struct
        Nnew  = N + ns
        comm  = A_ff.getComm()

        ai, aj, av = A_ff.getValuesCSR()
        d_nnz          = np.empty(Nnew, dtype=PETSc.IntType)
        d_nnz[:N]      = np.diff(ai) + ns          # room for the structural cols
        d_nnz[N:]      = Nnew                       # generous for structural rows
        ai_pad         = np.concatenate(
            [ai, np.full(ns, ai[-1], dtype=ai.dtype)])

        A = PETSc.Mat().createAIJ((Nnew, Nnew), nnz=d_nnz, comm=comm)
        A.setUp()
        A.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
        A.setValuesCSR(ai_pad, aj, av)
        A.assemble()

        # --- transpiration constraint: overwrite interface velocity rows ----
        interface_dofs = np.array([c[0] for c in couplings], dtype=PETSc.IntType)
        A.zeroRows(interface_dofs, diag=1.0)
        for (row, b_list, a_list) in couplings:
            for k in range(self.n_dof):
                A.setValue(row, self._eta_col(k),  b_list[k])   # + b^(k) eta_k
                A.setValue(row, self._phi_col(k), -a_list[k])   # - a^(k) phi_k

        # --- structural ODE rows (A = -1j * algebraic terms) ----------------
        for k in range(self.n_dof):
            # eta_k row :  d(eta_k)/dt - phi_k = 0
            A.setValue(self._eta_col(k), self._phi_col(k), 1j)
            # phi_k row :  sum_j M_kj d(phi_j)/dt + D_kj phi_j + K_kj eta_j - F'_k = 0
            for j in range(self.n_dof):
                if self.K_mat[k, j] != 0.0:
                    A.setValue(self._phi_col(k), self._eta_col(j),
                               -1j * self.K_mat[k, j])
                if self.D_mat[k, j] != 0.0:
                    A.setValue(self._phi_col(k), self._phi_col(j),
                               -1j * self.D_mat[k, j])
            # force/moment coupling to the fluid:  -1j*(-F'_k) = +1j*m^(k)
            m_arr = force_vecs[k].getArray()
            for l in np.nonzero(m_arr)[0]:
                A.setValue(self._phi_col(k), int(l), 1j * m_arr[l])

        A.assemble()
        return A

    # ----------------------------------------------------------------------
    # Bordering of the weight matrix B
    # ----------------------------------------------------------------------
    def _border_B(self, B_ff, couplings):
        N     = B_ff.getSize()[0]
        ns    = self.n_struct
        Nnew  = N + ns
        comm  = B_ff.getComm()

        bi, bj, bv = B_ff.getValuesCSR()
        d_nnz          = np.empty(Nnew, dtype=PETSc.IntType)
        d_nnz[:N]      = np.diff(bi) + 1
        d_nnz[N:]      = ns + 1
        bi_pad         = np.concatenate(
            [bi, np.full(ns, bi[-1], dtype=bi.dtype)])

        B = PETSc.Mat().createAIJ((Nnew, Nnew), nnz=d_nnz, comm=comm)
        B.setUp()
        B.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
        B.setValuesCSR(bi_pad, bj, bv)
        B.assemble()

        # interface-velocity rows are algebraic constraints -> zero in B
        interface_dofs = np.array([c[0] for c in couplings], dtype=PETSc.IntType)
        B.zeroRows(interface_dofs, diag=0.0)

        # structural mass rows
        for k in range(self.n_dof):
            B.setValue(self._eta_col(k), self._eta_col(k), 1.0)   # d(eta_k)/dt
            for j in range(self.n_dof):                            # M_kj d(phi_j)/dt
                if self.M_mat[k, j] != 0.0:
                    B.setValue(self._phi_col(k), self._phi_col(j),
                               self.M_mat[k, j])
        B.assemble()
        return B

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
