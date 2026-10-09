# Standard libraries
import  os
import  time

# Third party libraries
import  numpy as np

# Local Libraries and methods
from    FELiCS.Equation.EquationCollection  import EquationCollectionClass
from    FELiCS.Fields.MeanFlowClass         import MeanFlowClass
from    FELiCS.Fields.ModeCollection        import ModeCollection
from    FELiCS.IO.Writer                    import Writer
from 	FELiCS.Misc.logging                 import Logger
from    FELiCS.Solvers.LinearSolver         import LinearSolver

# Get the logger
logger = Logger.get_logger("felics")


def _export_fsi_structural_modes(solution, writer, dof_labels=None):
    '''Write the structural amplitudes (eta_k, phi_k for each DOF) of every direct
    FSI mode to a CSV next to spectrum.csv.

    The rows are written in the same order as the direct modes in spectrum.csv,
    so row i of fsi_structural_modes.csv corresponds to row i of spectrum.csv.
    For n_dof degrees of freedom the columns are

        omega_r, omega_i,
        eta0_r, eta0_i, phi0_r, phi0_i, [eta1_r, ...],
        struct_ratio, struct_ratio_energy, error

    The eigenvalue convention is identical to spectrum.csv (omega_i is the growth
    rate).  eta_hat / phi_hat are stored as lists (one entry per DOF).
    '''
    rows  = []
    n_dof = 0
    for mode in solution.modeList:
        if getattr(mode, "isAdjoint", False):
            continue
        eta = getattr(mode, "eta_hat", None)
        phi = getattr(mode, "phi_hat", None)
        if eta is None or phi is None:
            continue
        # normalise to lists (single-DOF legacy may have stored scalars)
        if np.ndim(eta) == 0:
            eta = [eta]
            phi = [phi]
        n_dof  = max(n_dof, len(eta))
        lam    = mode.eigen_value
        ratio  = getattr(mode, "struct_ratio", float("nan"))
        eratio = getattr(mode, "struct_ratio_energy", float("nan"))
        err    = getattr(mode, "_error", float("nan"))   # SLEPc relative residual
        frac   = getattr(mode, "dof_frac",  [float("nan")] * len(eta))
        phase  = getattr(mode, "dof_phase", [float("nan")] * len(eta))

        # per-DOF blocks: raw complex amplitudes, then |eta|, participation, phase
        row = [lam.real, lam.imag]
        for k in range(len(eta)):
            row += [eta[k].real, eta[k].imag, phi[k].real, phi[k].imag]
            row += [abs(eta[k]),
                    frac[k]  if k < len(frac)  else float("nan"),
                    phase[k] if k < len(phase) else float("nan")]

        # which DOF carries most of this mode's structural energy?
        try:
            dominant = int(np.argmax(frac)) if len(frac) else -1
        except (TypeError, ValueError):
            dominant = -1
        row += [dominant, ratio, eratio, err]
        rows.append(row)
    if not rows:
        return

    # Build the header.  Each DOF gets its own block, TAGGED with the motion so the
    # file is self-documenting.  Per DOF k:
    #     eta{k}_{tag}_r/_i, phi{k}_{tag}_r/_i   raw complex amplitudes
    #     eta{k}_{tag}_abs                       |eta_k|
    #     dof{k}_{tag}_frac                      energy participation in [0,1]
    #     dof{k}_{tag}_phase_deg                 phase of eta_k relative to eta_0
    # then a global 'dominant_dof' (index of the largest participation).
    # The index k is the position among the ENABLED DOFs, i.e. exactly the "DOF k"
    # numbering printed in the log at start-up.
    cols = ["omega_r", "omega_i"]
    for k in range(n_dof):
        if dof_labels is not None and k < len(dof_labels):
            # commas would split the CSV header into bogus columns -> strip them
            safe = str(dof_labels[k]).replace(",", ";")
            tag  = f"{k}_{safe}"
        else:
            tag = f"{k}"
        cols += [f"eta{tag}_r", f"eta{tag}_i", f"phi{tag}_r", f"phi{tag}_i",
                 f"eta{tag}_abs", f"dof{tag}_frac", f"dof{tag}_phase_deg"]
    cols += ["dominant_dof", "struct_ratio", "struct_ratio_energy", "error"]

    # pad rows to equal length (in case some modes had fewer DOFs recorded)
    width = len(cols)
    rows  = [r + [float("nan")] * (width - len(r)) for r in rows]

    filePath = os.path.join(writer.exportFolder, "fsi_structural_modes.csv")
    np.savetxt(filePath, np.array(rows), delimiter=',',
               header=",".join(cols), comments='')
    label_str = (", ".join(f"DOF {k} = {l}" for k, l in enumerate(dof_labels))
                 if dof_labels else f"{n_dof} DOF")
    logger.info(f"FSI structural amplitudes ({len(rows)} modes; {label_str}) "
                f"exported to {filePath}.")


def _fsi_is_enabled(param):
    '''Return True if the case file requests rigid-body FSI.'''
    fsi = getattr(param, "FSI", None)
    return bool(fsi is not None and getattr(fsi, "Enabled", False))


def run_modal(param):
    '''This function runs the calculations preset in param
    Input:
        param: Parameter objects (see parameters.py), defining the case
    '''

    logger.info("Running Modal analysis")
    #-----------------------------------------------------------------------
    ## INITIALIZATION
    #-----------------------------------------------------------------------
    # Get FELiCS objects required for analysis
    # mesh
    mesh            = param.get_mesh()
    from    FELiCS.SpaceDisc.FEMSpaces          import FEMSpaces
    
    # FEMSpaces
    FEMSpaces       = FEMSpaces(
        param,
        mesh,
    )

    # writer to export the results in files
    writer          = Writer(
        mesh,
        param.Export.ExportFolder,
    )

    # read in mean flow and export to h5-file
    meanFlow        = MeanFlowClass(
        param,
        FEMSpaces,
        mesh,
    )
    meanFlow.import_data_from_file_and_export_to_h5(writer)

    # equation
    equation        = EquationCollectionClass(
        param,
        FEMSpaces,
        meanFlow,
        mesh,
    )

    #-----------------------------------------------------------------------
    ## MAIN PART
    #-----------------------------------------------------------------------
    # get matrices for eigenproblem
    fsi_enabled = _fsi_is_enabled(param)
    if fsi_enabled:
        # ---- FSI-augmented system (fluid + rigid-body rotation DOFs) -------
        from FELiCS.Equation.FSI.RigidBodyFSI import RigidBodyMotionFSI

        logger.info("FSI is enabled: building the augmented (fluid + structure) operators.")
        fsi      = RigidBodyMotionFSI(param, equation, meanFlow)
        A, B     = fsi.assemble_augmented_operators()
        n_struct = fsi.n_struct
    else:
        # ---- standard fluid-only system (unchanged behaviour) --------------
        A        = equation.get_linear_operator(meanFlow)
        B        = equation.get_weight_matrix  (meanFlow)
        n_struct = 0

    # get parameters for eigenproblem
    guesses         = param.Numerics.EigenValueGuess
    nSol            = param.Numerics.nSolut
    adjoint         = param.Case.CalculateAdjoint

    # eigensolver settings from the case file (fall back to the solver defaults
    # if the case file does not specify them)
    eig_tol    = getattr(param.Numerics, "EigenSolverTolerance",     1.e-12)
    eig_max_it = getattr(param.Numerics, "EigenSolverMaxIterations", 200)
    logger.info(
        f"Eigensolver settings: tolerance = {eig_tol}, max iterations = {eig_max_it}."
    )

    # FSI adjoint is not implemented yet -> direct modal only
    if fsi_enabled and adjoint:
        logger.warning(
            "Adjoint modal analysis is not implemented for the FSI-augmented "
            "system; computing only the DIRECT spectrum."
        )
        adjoint = False

    # track time
    start           = time.time()

    # solve eigenproblem for each guess
    solution        = ModeCollection(
        FEMSpaces.VMixed,
        mesh,
        analysisType = "modal",
    )
    for guess in guesses:

        logger.info("Solving direct GEVP for guess: omega = " + str(guess))
        tmp         = LinearSolver.solve_general_eigenproblem(
            A,
            B,
            guess,
            nSol,
            tol=eig_tol,
            max_it=eig_max_it,
        )

        if fsi_enabled:
            # Split the augmented eigenvectors into the fluid part (first N
            # entries, stored as the mode) and the structural amplitudes
            # (eta_k, phi_k for each DOF) in the last n_struct entries.
            eig_vals, eig_vecs_aug, error = tmp
            N             = A.getSize()[0] - n_struct
            n_dof         = fsi.n_dof
            eig_vecs_flu  = eig_vecs_aug[:, :N]
            struct_amp    = eig_vecs_aug[:, N:]            # shape (nev, 2*n_dof)
            tmp           = (eig_vals, eig_vecs_flu, error)

            solution.append_solution_of_eigen_problem(
                tmp,
                guess,
                m       = param.Case.m
            )

            # attach the structural amplitudes to the freshly added modes
            new_modes = solution.modeList[-len(eig_vals):]
            for i, mode in enumerate(new_modes):
                # per-DOF amplitudes: eta_k = struct[2k], phi_k = struct[2k+1]
                mode.eta_hat = [complex(struct_amp[i, 2 * k])     for k in range(n_dof)]
                mode.phi_hat = [complex(struct_amp[i, 2 * k + 1]) for k in range(n_dof)]

                # Relative weight of the structural part in the eigenvector:
                #   ||q_s|| / ||q|| ,  q_s = (0, 0, eta_1, phi_1, ...).
                # A ratio, hence independent of the eigenvector normalisation.
                full_norm   = np.linalg.norm(eig_vecs_aug[i, :])
                struct_norm = np.linalg.norm(struct_amp[i, :])
                mode.struct_ratio = (
                    float(struct_norm / full_norm) if full_norm > 0.0 else 0.0
                )
                # energy-weighted version: sqrt(E_s)/sqrt(E_f + E_s), using the
                # velocity mass matrix for the fluid kinetic energy and the full
                # coupled (M, K) for the structural kinetic + potential energy.
                mode.struct_ratio_energy = fsi.structural_energy_ratio(
                    eig_vecs_aug[i, :]
                )
                # per-DOF energy participation + phase relative to DOF 0.
                # p[k] answers "how much of THIS mode is DOF k?" directly.
                mode.dof_frac, mode.dof_phase = fsi.dof_participation(
                    eig_vecs_aug[i, :]
                )
                logger.debug(
                    f"  FSI mode lambda = {mode.eigen_value}: "
                    f"eta_hat = {mode.eta_hat}, phi_hat = {mode.phi_hat}, "
                    f"||q_s||/||q|| = {mode.struct_ratio:.3e}, "
                    f"energy ratio = {mode.struct_ratio_energy:.3e}, "
                    f"participation = {mode.dof_frac}"
                )
        else:
            solution.append_solution_of_eigen_problem(
                tmp,
                guess,
                m       = param.Case.m
            )

        n_calculated = nSol

        if adjoint:

            logger.info("Solving adjoint GEVP for guess: omega = " + str(np.conj(guess)))
            tmp = LinearSolver.solve_general_eigenproblem(
                A,
                B,
                guess,
                nSol,
                tol=eig_tol,
                max_it=eig_max_it,
                adjoint = True,
            )

            solution.append_solution_of_eigen_problem(
                tmp,
                guess,
                adjoint = True,
                m = param.Case.m
            )

            n_calculated += nSol

        # Exporting the (temporary) spectrum to a file
        solution.export_spectrum_to_csv(writer)

        # Exporting the structural (eta, phi) amplitudes of the FSI modes
        if fsi_enabled:
            _export_fsi_structural_modes(solution, writer, fsi.dof_labels)

        # Exporting only the newly calculated modes to files (for this guess)
        solution.export_modes(
            writer,
            onlyNewN = n_calculated,
        )

    # End tracking time
    end             = time.time() - start
    logger.info('Solving the general eigenproblem took %4g s' % end)
    residuum_max    = solution.get_maximum_error()
    logger.debug('Maximum residuum of all solutions:  %12g' % (residuum_max))