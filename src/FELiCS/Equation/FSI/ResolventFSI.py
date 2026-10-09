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
ResolventFSI.py -- border the resolvent machinery for rigid-body FSI.

The resolvent operator itself needs no new derivation.  FELiCS solves
``A q = lambda B q`` with ``q ~ exp(-i*lambda*t)``, and run_resolvent forms

    R = A - omega*B

so R is exactly the shifted operator that shift-invert already factorises. The
FSI module already produces the augmented ``A_aug`` / ``B_aug`` of size N + 2n.
Feeding those to ``ResolventOperator`` therefore "just works" -- PROVIDED the
five auxiliary matrices are bordered consistently:

    W_FEM        (N+2n) x (N+2n)      how a forcing coefficient enters the eqs
    P_forcing    (N+2n) x n_f         which rows may be forced
    W_forcing     n_f    x n_f        metric on the forcing  (MUST be invertible)
    P_response    n_r    x (N+2n)     which rows are measured
    W_response    n_r    x n_r        metric on the response

This module builds those.  Everything else (SVD, MUMPS solves, gains.csv) is the
existing FELiCS code, untouched.

Four input/output combinations
------------------------------
Selected from the case file with ForceFluid / ForceStructure / MeasureFluid /
MeasureStructure:

    force    measure     question answered
    -----    -------     -----------------
    fluid    fluid       wake receptivity (FSI only modifies the operator)
    fluid    body        GUST RESPONSE: flow disturbance -> body motion
    body     fluid       ACTUATOR AUTHORITY: applied force/moment -> wake
    both     both        full coupled input-output

Which rows may be forced
------------------------
The ``eta`` rows are NEVER forced.  ``eta_dot = phi`` is a kinematic definition,
not a force balance; putting a source in it is meaningless and would silently
corrupt the gain.  (The pressure rows are excluded for the same reason, which
FELiCS already does: for the TKE norm ``resolventForcingIndices`` is the
collapsed velocity sub-space.)  Structural forcing therefore acts on the ``phi``
rows only -- a physical external force / moment on the body.

Because the eta rows are unforced, the eta row of ``(A - omega B) q = rhs`` reads
``i*phi_k - omega*eta_k = 0``, i.e. ``|phi_k| = |omega| * |eta_k|`` exactly.  The
kinetic and displacement response norms therefore differ only by a factor
omega^2 at each frequency -- both are legitimate measures of "how much the body
moves", they just scale the gain curve differently.

Norms (the part that is a modelling choice, not a coding detail)
---------------------------------------------------------------
StructuralForcingNorm:
    'Inertia'  (default)  W_f,struct = M^-1
        The generalised forces are a force [M L T^-2] for a translation DOF and a
        moment [M L^2 T^-2] for a rotation DOF -- NOT commensurable.  The metric
        f^H M^-1 f has units M L^2 T^-4 for BOTH, so M^-1 is the dimensionally
        consistent metric.  It is also SPD (M is SPD), so W_forcing stays
        invertible, which ResolventOperator requires (it LU-factorises it).
    'Identity'            W_f,struct = I / StructuralForcingScale^2
        Only sensible for a single DOF, or if the user knows what they are doing.

StructuralResponseNorm:
    'Kinetic'  (default)  measure phi;   W = M      (structural kinetic energy)
    'Displacement'        measure eta;   W = M      (M as the metric fixes units)
    'Energy'              measure eta+phi; W = diag(K, M)
        NOTE with K = 0 (free body) the eta block is ZERO, so displacement
        contributes nothing and 'Energy' silently degenerates to 'Kinetic'.
        A warning is issued.
"""

import numpy as np
from petsc4py import PETSc

from FELiCS.Misc.logging import Logger

logger = Logger.get_logger("felics")


# ----------------------------------------------------------------------
def _embed(mat, n_rows, n_cols, extra=None, diag_fill=None):
    """Copy `mat` into the top-left corner of a new (n_rows x n_cols) PETSc matrix.

    extra     : list of (i, j, v) entries to add outside the copied block
    diag_fill : (start, stop, value) -> put `value` on the diagonal for
                rows/cols in [start, stop)

    The copy is done through the CSR arrays, so it is cheap even for the full
    FE mass matrix.
    """
    mat.assemble()
    ai, aj, av = mat.getValuesCSR()
    m_rows = mat.getSize()[0]

    extra = extra or []
    # nonzeros per row: original + anything we add
    nnz = np.zeros(n_rows, dtype=np.int32)
    nnz[:m_rows] = np.diff(ai)
    for (i, j, v) in extra:
        nnz[i] += 1
    if diag_fill is not None:
        s, e, _ = diag_fill
        for i in range(s, e):
            nnz[i] += 1

    out = PETSc.Mat().createAIJ([n_rows, n_cols], nnz=nnz)
    out.setUp()
    out.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)

    # copy the original block row by row (columns keep their indices)
    for i in range(m_rows):
        lo, hi = ai[i], ai[i + 1]
        if hi > lo:
            out.setValues(i, aj[lo:hi].astype(np.int32), av[lo:hi])

    if diag_fill is not None:
        s, e, val = diag_fill
        for i in range(s, e):
            out.setValue(i, i, val)

    for (i, j, v) in extra:
        out.setValue(i, j, v)

    out.assemble()
    return out


# ----------------------------------------------------------------------
class FSIResolventBorder:
    """Borders W_FEM / P_forcing / W_forcing / P_response / W_response so that
    ResolventOperator can be driven with the augmented FSI operator."""

    def __init__(self, fsi, param):
        self.fsi = fsi
        self.n   = fsi.n_dof
        self.N   = fsi.N_fluid                 # fluid size
        self.Na  = fsi.N_fluid + 2 * fsi.n_dof  # augmented size

        io = param.IOResolvent
        self.force_fluid       = bool(getattr(io, 'ForceFluid',        True))
        self.force_structure   = bool(getattr(io, 'ForceStructure',    False))
        self.measure_fluid     = bool(getattr(io, 'MeasureFluid',      True))
        self.measure_structure = bool(getattr(io, 'MeasureStructure',  False))
        self.struct_f_norm     = str(getattr(io, 'StructuralForcingNorm',  'Inertia'))
        self.struct_r_norm     = str(getattr(io, 'StructuralResponseNorm', 'Kinetic'))
        self.struct_f_scale    = float(getattr(io, 'StructuralForcingScale', 1.0))

        if not (self.force_fluid or self.force_structure):
            raise ValueError(
                "FSI resolvent: ForceFluid and ForceStructure are both false -- "
                "there is nothing to force.")
        if not (self.measure_fluid or self.measure_structure):
            raise ValueError(
                "FSI resolvent: MeasureFluid and MeasureStructure are both false -- "
                "there is nothing to measure.")

        if self.struct_r_norm == 'Energy' and not np.any(fsi.K_mat):
            logger.warning(
                "FSI resolvent: StructuralResponseNorm='Energy' but the stiffness "
                "matrix K is zero (free body), so the eta block of the response "
                "norm is zero and displacement contributes NOTHING to the gain. "
                "'Energy' has degenerated to 'Kinetic'. Use 'Displacement' if you "
                "want to measure eta.")

        logger.info(
            f"FSI resolvent I/O: force = "
            f"{'fluid' if self.force_fluid else ''}"
            f"{'+' if self.force_fluid and self.force_structure else ''}"
            f"{'structure' if self.force_structure else ''}"
            f" | measure = "
            f"{'fluid' if self.measure_fluid else ''}"
            f"{'+' if self.measure_fluid and self.measure_structure else ''}"
            f"{'structure' if self.measure_structure else ''}")
        if self.force_structure:
            logger.info(f"  structural forcing norm  : {self.struct_f_norm} "
                        f"(applied to the phi rows only)")
        if self.measure_structure:
            logger.info(f"  structural response norm : {self.struct_r_norm}")

    # -- index helpers (must match RigidBodyFSI's layout) -----------------
    def _eta(self, k):
        return self.N + 2 * k

    def _phi(self, k):
        return self.N + 2 * k + 1

    # ------------------------------------------------------------------
    def border_fem_weighting(self, W_FEM):
        """W_FEM_aug = blkdiag(W_FEM, I_2n).

        W_FEM turns a forcing coefficient vector into the right-hand side of the
        equations.  For the fluid that is the FE mass matrix (int f.v dx).  For a
        structural DOF the generalised force IS the right-hand side of the phi
        row, so the weighting is the identity.

        (Strictly the phi row of A carries the A = -i*alg convention, so the
        physical force enters as -i*f.  That is a unimodular factor: it does not
        change any singular value, only the phase of the reported optimal
        structural forcing.)
        """
        return _embed(W_FEM, self.Na, self.Na,
                      diag_fill=(self.N, self.Na, 1.0))

    # ------------------------------------------------------------------
    def border_forcing(self, P_forcing, W_forcing):
        """Return (P_forcing_aug, W_forcing_aug).

        Structural forcing acts on the phi rows ONLY (see module docstring).
        """
        nf_fluid = P_forcing.getSize()[1] if self.force_fluid else 0
        nf_struct = self.n if self.force_structure else 0
        nf = nf_fluid + nf_struct
        if nf == 0:
            raise ValueError("FSI resolvent: empty forcing space.")

        # ---- P_forcing_aug : (N + 2n) x nf --------------------------------
        extra = []
        for k in range(nf_struct):
            extra.append((self._phi(k), nf_fluid + k, 1.0))

        if self.force_fluid:
            P_f_aug = _embed(P_forcing, self.Na, nf, extra=extra)
        else:
            P_f_aug = PETSc.Mat().createAIJ([self.Na, nf], nnz=1)
            P_f_aug.setUp()
            P_f_aug.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
            for (i, j, v) in extra:
                P_f_aug.setValue(i, j, v)
            P_f_aug.assemble()

        # ---- W_forcing_aug : nf x nf --------------------------------------
        Wf_struct = self._structural_forcing_norm()   # n x n (dense, small)
        extra_w = []
        for a in range(nf_struct):
            for b in range(nf_struct):
                if Wf_struct[a, b] != 0.0:
                    extra_w.append((nf_fluid + a, nf_fluid + b,
                                    complex(Wf_struct[a, b])))

        if self.force_fluid:
            W_f_aug = _embed(W_forcing, nf, nf, extra=extra_w)
        else:
            W_f_aug = PETSc.Mat().createAIJ([nf, nf], nnz=max(nf_struct, 1))
            W_f_aug.setUp()
            W_f_aug.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
            for (i, j, v) in extra_w:
                W_f_aug.setValue(i, j, v)
            W_f_aug.assemble()

        return P_f_aug, W_f_aug

    def _structural_forcing_norm(self):
        """Metric on the generalised forces.  MUST be SPD: ResolventOperator
        LU-factorises W_forcing."""
        if self.struct_f_norm == 'Inertia':
            # f^H M^-1 f : the only dimensionally consistent choice when the
            # body has both translation (force) and rotation (moment) DOFs.
            W = np.linalg.inv(self.fsi.M_mat)
        elif self.struct_f_norm == 'Identity':
            W = np.eye(self.n) / (self.struct_f_scale ** 2)
            if self.n > 1:
                logger.warning(
                    "FSI resolvent: StructuralForcingNorm='Identity' with more "
                    "than one DOF mixes forces and moments in one norm, which is "
                    "dimensionally inconsistent. 'Inertia' is recommended.")
        else:
            raise ValueError(
                f"FSI resolvent: unknown StructuralForcingNorm "
                f"'{self.struct_f_norm}' (use 'Inertia' or 'Identity').")

        ev = np.linalg.eigvalsh(0.5 * (W + W.T))
        if np.any(ev <= 0):
            raise ValueError(
                f"FSI resolvent: structural forcing norm is not positive definite "
                f"(eigenvalues {ev}); ResolventOperator cannot factorise it.")
        return W

    # ------------------------------------------------------------------
    def border_response(self, P_response, W_response):
        """Return (P_response_aug, W_response_aug)."""
        nr_fluid = P_response.getSize()[0] if self.measure_fluid else 0

        if not self.measure_structure:
            rows, W_struct = [], np.zeros((0, 0))
        elif self.struct_r_norm == 'Kinetic':
            rows = [self._phi(k) for k in range(self.n)]
            W_struct = self.fsi.M_mat.copy()
        elif self.struct_r_norm == 'Displacement':
            rows = [self._eta(k) for k in range(self.n)]
            W_struct = self.fsi.M_mat.copy()      # M as the metric fixes units
        elif self.struct_r_norm == 'Energy':
            rows = ([self._eta(k) for k in range(self.n)] +
                    [self._phi(k) for k in range(self.n)])
            W_struct = np.zeros((2 * self.n, 2 * self.n))
            W_struct[:self.n, :self.n] = self.fsi.K_mat     # potential energy
            W_struct[self.n:, self.n:] = self.fsi.M_mat     # kinetic energy
        else:
            raise ValueError(
                f"FSI resolvent: unknown StructuralResponseNorm "
                f"'{self.struct_r_norm}' (use 'Kinetic', 'Displacement' or "
                f"'Energy').")

        nr_struct = len(rows)
        nr = nr_fluid + nr_struct
        if nr == 0:
            raise ValueError("FSI resolvent: empty response space.")

        # ---- P_response_aug : nr x (N + 2n) -------------------------------
        extra = [(nr_fluid + a, col, 1.0) for a, col in enumerate(rows)]
        if self.measure_fluid:
            P_r_aug = _embed(P_response, nr, self.Na, extra=extra)
        else:
            P_r_aug = PETSc.Mat().createAIJ([nr, self.Na], nnz=1)
            P_r_aug.setUp()
            P_r_aug.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
            for (i, j, v) in extra:
                P_r_aug.setValue(i, j, v)
            P_r_aug.assemble()

        # ---- W_response_aug : nr x nr -------------------------------------
        extra_w = []
        for a in range(nr_struct):
            for b in range(nr_struct):
                if W_struct[a, b] != 0.0:
                    extra_w.append((nr_fluid + a, nr_fluid + b,
                                    complex(W_struct[a, b])))
        if self.measure_fluid:
            W_r_aug = _embed(W_response, nr, nr, extra=extra_w)
        else:
            W_r_aug = PETSc.Mat().createAIJ([nr, nr], nnz=max(nr_struct, 1))
            W_r_aug.setUp()
            W_r_aug.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
            for (i, j, v) in extra_w:
                W_r_aug.setValue(i, j, v)
            W_r_aug.assemble()

        return P_r_aug, W_r_aug

    # ------------------------------------------------------------------
    def layout(self, P_forcing, P_response):
        """Sizes of the forcing / response spaces, so the driver can slice the
        singular vectors into fluid and structural parts."""
        nf_fluid  = P_forcing.getSize()[1] if self.force_fluid else 0
        nf_struct = self.n if self.force_structure else 0
        nr_fluid  = P_response.getSize()[0] if self.measure_fluid else 0
        if not self.measure_structure:
            nr_struct = 0
        elif self.struct_r_norm == 'Energy':
            nr_struct = 2 * self.n
        else:
            nr_struct = self.n
        return dict(nf_fluid=nf_fluid, nf_struct=nf_struct,
                    nr_fluid=nr_fluid, nr_struct=nr_struct)


# ----------------------------------------------------------------------
def resolvent_shift(omega, param):
    """Return the (possibly complex) shift used to form R = A - shift*B.

    Standard resolvent
    ------------------
        shift = omega   (real)
    Valid only when the base flow is STABLE.  If an eigenvalue sits above the
    real axis the harmonic 'steady state' does not exist: the impulse response
    grows, the gain near that frequency is dominated by the unstable pole, and
    the whole curve is meaningless in exactly the band you care about.

    Discounted resolvent
    --------------------
        shift = omega + i*beta,   beta > max_j Im(lambda_j)

    Derivation (in FELiCS's convention).  With B dq/dt = L q + f and A = i*L, the
    discounted state q~ = exp(-beta t) q obeys

        B dq~/dt = (L - beta B) q~ + f~

    whose growth rates are those of the original shifted DOWN by beta.  In terms
    of A this is A - i*beta*B, so the resolvent of the discounted system at real
    omega is

        (A - i*beta*B - omega*B)^-1 = (A - (omega + i*beta) B)^-1

    Choosing beta larger than the largest growth rate renders the discounted
    operator stable and the gains finite and meaningful.  The gains are then
    those of the DISCOUNTED system: they answer "how much does the flow amplify
    forcing over a time horizon ~1/beta", not "what is the asymptotic steady
    harmonic response" (which does not exist).  Report beta alongside the gains;
    curves at different beta are not comparable.
    """
    if not bool(getattr(param.IOResolvent, 'Discounted', False)):
        return complex(omega), 0.0

    beta = float(getattr(param.IOResolvent, 'DiscountBeta', 0.0))
    if beta <= 0.0:
        raise ValueError(
            "IOResolvent.Discounted is true but DiscountBeta <= 0. Set "
            "DiscountBeta strictly larger than the largest growth rate "
            "(max omega_i) of your spectrum, which you can read from the modal "
            "run's spectrum.csv.")
    return complex(omega, beta), beta