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
# Standard libraries
import  time

# Local Libraries and methods
from    FELiCS.Equation.EquationCollection  import EquationCollectionClass
from    FELiCS.Fields.MeanFlowClass         import MeanFlowClass
from    FELiCS.Fields.Mode                  import Mode
from    FELiCS.Fields.ModeCollection        import ModeCollection
from    FELiCS.IO.Writer                    import Writer
from 	FELiCS.Misc.logging                 import Logger
from    FELiCS.Solvers.LinearSolver         import LinearSolver, ResolventOperator
from    FELiCS.Equation.FSI.RigidBodyFSI    import RigidBodyMotionFSI
from    FELiCS.Equation.FSI.ResolventFSI    import FSIResolventBorder, resolvent_shift

import  csv
import  os
import  numpy   as  np

# Get the logger
logger                          = Logger.get_logger("felics")

def run_resolvent(param):
    '''This function runs the calculations preset in param
    Input:
        param: Parameter objects (see parameters.py), defining the case
    '''

    logger.warning("Running Resolvent analysis")
    
    #-----------------------------------------------------------------------
    ## INITIALIZATION
    #-----------------------------------------------------------------------
    # mesh
    mesh      = param.get_mesh()
    from    FELiCS.SpaceDisc.FEMSpaces          import FEMSpaces

    # FEMSpaces
    FEMSpaces = FEMSpaces(
        param,
        mesh,
    )
     
    # writer to export the results in files
    writer    = Writer(
        mesh,
        param.Export.ExportFolder,
    )

    # read in mean flow and export to h5-file
    meanFlow = MeanFlowClass(
        param,
        FEMSpaces,
        mesh,
    )
    meanFlow.import_data_from_file_and_export_to_h5(writer)

    # equation
    equation = EquationCollectionClass(
        param,
        FEMSpaces,
        meanFlow,
        mesh,
    )

    #-----------------------------------------------------------------------
    ## MAIN PART
    #-----------------------------------------------------------------------
    # Is rigid-body FSI switched on in the case file?
    fsi_enabled                 = bool(getattr(getattr(param, 'FSI', None), 'Enabled', False))

    # ---- fluid-side operators (unchanged) ------------------------------
    W_FEM                       = equation.get_resolvent_weighting_fem(meanFlow)  #=> including tensorial stuff

    P_shrink_f                  = equation.get_shrinker_mat_forcing()   
    P_forcing                   = equation.get_restrictor_mat_forcing().matMult(P_shrink_f)

    P_shrink_r                  = equation.get_shrinker_mat_response()
    P_response                  = P_shrink_r.transposeMatMult(equation.get_restrictor_mat_response())

    W_f_big                     = equation.get_resolvent_norm_forcing (meanFlow)
    W_r_big                     = equation.get_resolvent_norm_response(meanFlow)

    W_forcing                   = P_shrink_f.transposeMatMult(W_f_big.matMult(P_shrink_f))
    W_response                  = P_shrink_r.transposeMatMult(W_r_big.matMult(P_shrink_r))

    fsi                         = None
    border                      = None
    if not fsi_enabled:
        # ---- plain fluid resolvent, exactly as before -------------------
        A                       = equation.get_linear_operator(meanFlow)
        B                       = equation.get_weight_matrix  (meanFlow)
    else:
        # ---- rigid-body FSI: border everything to size N + 2*n_dof ------
        fsi                     = RigidBodyMotionFSI(
            param,
            equation,
            meanFlow,
        )
        A, B                    = fsi.assemble_augmented_operators()

        border                  = FSIResolventBorder(fsi, param)
        layout                  = border.layout(P_forcing, P_response)

        W_FEM                   = border.border_fem_weighting(W_FEM)
        P_forcing, W_forcing    = border.border_forcing (P_forcing,  W_forcing)
        P_response, W_response  = border.border_response(P_response, W_response)

        logger.info(
            f"FSI resolvent spaces: forcing = {layout['nf_fluid']} fluid + "
            f"{layout['nf_struct']} structural, response = {layout['nr_fluid']} "
            f"fluid + {layout['nr_struct']} structural.")

    # get parameters for resolvent analysis
    omegas                      = param.IOResolvent.Omegas
    nSol                        = param.Numerics.nSolut
    
    # Start tracking time
    start                       = time.time()
    
    # Prepare solution object
    solution                    = ModeCollection(
        FEMSpaces.VMixed,
        mesh,
        analysisType            = 'Resolvent',
    )
    
    gain_records                = []

    # Run analysis at each frequency
    for i, omega in enumerate(omegas):

        logger.info("Solving resolvent SVD for omega = " + str(omega))
        
        # Catch bug for omega = 1.0
        if omega == 1.0:
            omega               += 1.e-4
            logger.warning("Omega was equal to 1.0, which can lead to numerical issues. Added 1.e-4 to omega.")

        # Standard resolvent : shift = omega            (valid for a STABLE base flow)
        # Discounted         : shift = omega + i*beta   (beta > max growth rate)
        # See ResolventFSI.resolvent_shift for the derivation.
        shift, beta             = resolvent_shift(omega, param)
        if beta > 0.0 and i == 0:
            logger.warning(
                f"DISCOUNTED resolvent active (beta = {beta:g}). The gains are "
                f"those of the discounted system -- they are NOT comparable with "
                f"gains computed at a different beta, or with beta = 0.")

        # initialize the resolvent operator, which is a class that imitates 
        # a matrix to use matrix-free methods
        R                       = A.copy()         
        R.axpy(
            -shift,
            B,
        )       # R = A - (omega + i*beta)*B
        resolventOperator       = ResolventOperator(
            R,
            W_FEM,
            W_forcing,
            W_response,
            P_forcing,
            P_response,
        )
        
        # Perform eigenvalue decomposition of the linear operator defined in the class "ResolventOperator"         
        # via the matrix vector multiplation "mult"
        gains, eigenvectors_c   = LinearSolver.solve_svd_of_resolvent(
            resolventOperator,
            nev                 = nSol,
            tol                 = 1.e-16,
            max_it              = 200,
        )

        if not fsi_enabled:
            solution.append_solution_of_svd_problem(
                eigenvectors_c,
                omega,
                gains,
                resolventOperator,
                m                   = param.Case.m
            )
            # export spectrum and newly calculated modes
            solution.export_spectrum_to_csv(writer)
            solution.export_modes(
                writer,
                onlyNewN            = nSol*2,
            )
        else:
            # ModeCollection expects vectors of the FLUID size; the augmented
            # singular vectors carry the structural entries as well.  Build the
            # Mode objects from the fluid slice (so gains.csv and the ParaView
            # mode export keep working) and record the structural amplitudes
            # alongside -- in a single KSP solve.
            _append_and_record_fsi(
                solution,
                gains,
                eigenvectors_c,
                omega,
                beta,
                resolventOperator,
                fsi,
                border,
                layout,
                gain_records,
                m = param.Case.m,
            )
            solution.export_spectrum_to_csv(writer)
            solution.export_modes(
                writer,
                onlyNewN            = nSol*2,
            )
            # Rewrite the FSI csv after EVERY frequency (it is tiny, so this is
            # free) so that partial results are usable immediately and survive a
            # crash or a wall-clock kill.
            _export_fsi_resolvent_csv(
                gain_records,
                writer,
                fsi.dof_labels,
                border,
                quiet = (i < len(omegas) - 1),
            )

        resolventOperator.destroy_self()
        
    # End tracking time
    logger.info(f"Solving the SVD(s) took {time.time()-start:.4g} s")




# ----------------------------------------------------------------------
# FSI-specific bookkeeping
# ----------------------------------------------------------------------
def _append_and_record_fsi(solution, gains, forcing_vecs, omega, beta,
                           resolventOperator, fsi, border, layout, records, m=0):
    """Append the FLUID part of the forcing/response modes to the ModeCollection,
    and record the STRUCTURAL amplitudes, for one frequency.

    ModeCollection.append_solution_of_svd_problem cannot be reused directly: it
    calls set_coefficient_array on the full vector, which for FSI is of size
    N + 2n and does not fit the FE space.  So we reproduce its chain here and
    slice off the fluid part -- reusing the KSP that ResolventOperator already
    built, so there is no extra factorisation and only ONE solve per mode.

    The response is  q = -1j * R^-1 * W_FEM * P_f * f  (the -1j is FELiCS's
    convention: ModeCollection does Y2.scale(-1j)).
    """
    ksp     = resolventOperator.getKSP()
    W_FEM   = resolventOperator._W_FEM
    P_f     = resolventOperator._P_forcing
    W_f     = resolventOperator._W_forcing

    nff, nfs = layout['nf_fluid'], layout['nf_struct']
    n        = fsi.n_dof
    N        = fsi.N_fluid
    nf       = W_f.getSize()[0]           # size of the forcing space

    # SLEPc returns the singular vectors as COLUMNS (ModeCollection uses
    # forcingArray[:, i]).  Be defensive about the orientation anyway.
    fv = np.asarray(forcing_vecs)
    if fv.ndim != 2:
        raise ValueError(f"unexpected singular-vector array shape {fv.shape}")
    if fv.shape[0] == nf:
        get_vec = lambda i: fv[:, i]
    elif fv.shape[1] == nf:
        get_vec = lambda i: fv[i, :]
    else:
        raise ValueError(
            f"singular-vector array {fv.shape} matches neither axis of the "
            f"forcing space (n_f = {nf}).")

    x_f, _  = W_f.getVecs()          # forcing-space vector (size n_f)
    o1, o2  = W_FEM.getVecs()        # full augmented-space vectors (size N + 2n)

    for i in range(gains.shape[0]):
        f_hat = np.asarray(get_vec(i)).ravel()

        # --- response to this forcing:  q = -1j * R^-1 W_FEM P_f f ----------
        # The -1j is FELiCS's convention (ModeCollection does Y2.scale(-1j)):
        # with A = i*L, the forced solution is q = -i (A - omega B)^-1 f.
        # It is unimodular, so it does not affect the gains -- only the phase.
        x_f.setArray(f_hat)
        P_f.mult(x_f, o1)
        W_FEM.mult(o1, o2)
        o2.scale(-1j)
        f_full = o1.getArray().copy()      # P_f * f, full augmented size
        ksp.solve(o2, o1)
        q = o1.getArray().copy()

        eta = np.array([q[N + 2 * k]     for k in range(n)])
        phi = np.array([q[N + 2 * k + 1] for k in range(n)])

        # --- fluid slice -> ModeCollection (keeps gains.csv + mode export) ---
        for arr, mtype in ((f_full[:N], 'forcing'), (q[:N], 'response')):
            md = Mode(
                solution.femSpace,
                solution.mesh,
                isStateVector = True,
                analysisType  = 'Resolvent',
                modeType      = mtype,
                m             = m,
            )
            md.gain        = float(np.real(gains[i]))   # NOTE: gains are sigma^2
            md.gain_number = i
            md.frequency   = float(np.real(omega))
            md.set_coefficient_array(arr)
            solution.modeList.append(md)

        # structural part of the optimal forcing (zero if not forcing the body)
        f_struct = (np.array(f_hat[nff:nff + nfs]) if nfs else np.zeros(0, complex))

        # how much of the optimal forcing sits on the body vs in the fluid
        e_fluid  = float(np.sum(np.abs(f_hat[:nff]) ** 2)) if nff else 0.0
        e_struct = float(np.sum(np.abs(f_struct) ** 2))    if nfs else 0.0
        tot      = e_fluid + e_struct
        frac_struct_forcing = (e_struct / tot) if tot > 0 else 0.0

        # structural vs fluid energy of the RESPONSE (physical, M-weighted)
        Es = float(np.real(np.conj(phi) @ fsi.M_mat @ phi))
        if np.any(fsi.K_mat):
            Es += float(np.real(np.conj(eta) @ fsi.K_mat @ eta))
        Bf = fsi.B_fluid
        v1, v2 = Bf.getVecs()
        v1.setArray(q[:N])
        Bf.mult(v1, v2)
        Ef = float(np.real(np.vdot(q[:N], v2.getArray())))
        ratio = Es / (Es + Ef) if (Es + Ef) > 0 else 0.0
        v1.destroy(); v2.destroy()

        records.append(dict(
            omega=float(np.real(omega)),
            beta=float(beta),
            gain_number=i,
            gain=float(np.sqrt(np.abs(np.real(gains[i])))),   # gains are sigma^2
            gain_sq=float(np.abs(np.real(gains[i]))),
            eta=eta, phi=phi, f_struct=f_struct,
            frac_struct_forcing=frac_struct_forcing,
            struct_ratio_energy=ratio,
        ))

    x_f.destroy(); o1.destroy(); o2.destroy()


def _export_fsi_resolvent_csv(records, writer, labels, border, quiet=False):
    """Write fsi_resolvent_gains.csv -- one row per (omega, gain number)."""
    if not records:
        return
    path = os.path.join(writer.exportFolder, "fsi_resolvent_gains.csv")

    tags = [f"{k}_{str(l).replace(',', ';')}" for k, l in enumerate(labels)]
    header = ["omega", "beta", "gain_number", "gain", "gain_sq"]
    for t in tags:
        header += [f"eta{t}_r", f"eta{t}_i", f"eta{t}_abs",
                   f"phi{t}_r", f"phi{t}_i", f"phi{t}_abs"]
    if border.force_structure:
        for t in tags:
            header += [f"fstruct{t}_r", f"fstruct{t}_i", f"fstruct{t}_abs"]
        header += ["frac_struct_forcing"]
    header += ["struct_ratio_energy"]

    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in records:
            row = [f"{r['omega']:.10e}", f"{r['beta']:.10e}",
                   r['gain_number'], f"{r['gain']:.10e}", f"{r['gain_sq']:.10e}"]
            for k in range(len(labels)):
                e, p = r['eta'][k], r['phi'][k]
                row += [f"{e.real:.10e}", f"{e.imag:.10e}", f"{abs(e):.10e}",
                        f"{p.real:.10e}", f"{p.imag:.10e}", f"{abs(p):.10e}"]
            if border.force_structure:
                for k in range(len(labels)):
                    f = r['f_struct'][k]
                    row += [f"{f.real:.10e}", f"{f.imag:.10e}", f"{abs(f):.10e}"]
                row += [f"{r['frac_struct_forcing']:.10e}"]
            row += [f"{r['struct_ratio_energy']:.10e}"]
            w.writerow(row)

    if not quiet:
        logger.info(f"FSI resolvent gains exported to {path}")