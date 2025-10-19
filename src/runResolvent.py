import  time
from 	FELiCS.Misc.logging     import Logger

# Get the logger
logger                          = Logger.get_logger("felics")

def runResolvent(param):
    '''This function runs the calculations preset in param
    Input:
        param: Parameter objects (see parameters.py), defining the case
    '''

    from    FELiCS.IO.ExportSolution            import ExportFromFile 
    from    FELiCS.SpaceDisc.FEMSpaces          import FEMSpaces
    from    FELiCS.Fields.meanFlowClass         import meanFlowClass
    from    FELiCS.Equation.EquationCollection  import EquationCollectionClass
    from    FELiCS.Solvers.LinearSolver         import LinearSolver, ResolventOperator
    from    FELiCS.Fields.ModeCollection        import ModeCollection

    logger.warning("Running Resolvent analysis")
    
    # Get FELiCS objects required for analysis
    mesh                        = param.getMesh()
    FEMSpaces                   = FEMSpaces(param, mesh)
    meanFlow                    = meanFlowClass(param, FEMSpaces, mesh)
    meanFlow.importDataFromFileAndExportToH5()
    equation                    = EquationCollectionClass(param,FEMSpaces,meanFlow, mesh)

    # Get operators required for resolvent analysis
    A                           = equation.getLinearOperator(meanFlow)
    B                           = equation.getWeightMatrix  (meanFlow)
    W_FEM                       = equation.getResolventWeighting_FEM(meanFlow)  #=> including tensorial stuff

    P_shrink_f                  = equation.getShrinkerMatForcing()   
    P_forcing                   = equation.getRestrictorMatForcing().matMult(P_shrink_f)

    P_shrink_r                  = equation.getShrinkerMatResponse()
    P_response                  = P_shrink_r.transposeMatMult(equation.getRestrictorMatResponse())

    W_f_big                     = equation.getResolventNorm_forcing (meanFlow)
    W_r_big                     = equation.getResolventNorm_response(meanFlow)

    W_forcing                   = P_shrink_f.transposeMatMult(W_f_big.matMult(P_shrink_f))
    W_response                  = P_shrink_r.transposeMatMult(W_r_big.matMult(P_shrink_r))

    # get parameters for resolvent analysis
    omegas                      = param.IOResolvent.Omegas
    nSol                        = param.Numerics.nSolut
    
    # Start tracking time
    start                       = time.time()
    
    # Run analysis at each frequency
    solution                    = ModeCollection(
        FEMSpaces.VMixed, 
        mesh, 
        analysis='Resolvent'
    )
    for i, omega in enumerate(omegas):

        logger.info("Solving resolvent SVD for omega = " + str(omega))
        # R = A-omega*B         
        R                       = A.copy()         
        R.axpy(-omega, B)            
        resolventOperator       = ResolventOperator(
            R,                                         
            W_FEM,                                         
            W_forcing,                                         
            W_response,                                         
            P_forcing,                                         
            P_response
        )
        
        # Perform eigenvalue decomposition of the linear operator defined in the class "ResolventOperator"         
        # via the matrix vector multiplation "mult"
        gains, eigenvectors_c   = LinearSolver.solveSVDOfResolvent(
                resolventOperator,
                nev             = nSol,
                tol             = 1.e-16,
                max_it          = 200,
                )
        solution.appendSolutionOfSVDProblem( # NOTE: this also compute the response modes
            eigenvectors_c, 
            omega, 
            gains, 
            resolventOperator, 
        )
        
    # End tracking time
    logger.info(f"Solving the SVD(s) took {time.time()-start:.4g} s")

    # Convert solution to fluctuationSolution objects for export
    fluctSolutList              = solution.getOldSolutionObject(
        meanFlow, 
        param, 
        FEMSpaces
    )
    ExportFromFile(
        param, 
        FEMSpaces, 
        fluctSolutList, 
        meanFlow
    )
