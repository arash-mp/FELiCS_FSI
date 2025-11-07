import  time
from 	FELiCS.Misc.logging     import Logger
from    FELiCS.IO.Reader        import Reader
from    FELiCS.Fields.Mode      import Mode

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
    
    # ====================== Test mode collection for resolvent ======================
    # Get a mode to check structure
    # exampleMode     = solution.modeList[0]
    
    # Get the spectrum from the solution
    spectrum, header    = solution.getSpectrum()
    spectrumFile        = param.Export.ExportFolder + "/gains_check.csv"
    solution.exportSpectrumToCSV(spectrumFile)
    
    # Loop over all modes and print their attributes to verify successful conversion to mode collection
    modeList = solution.modeList
    for i, mode in enumerate(solution.modeList):
        logger.info(f"Mode {i}: Frequency = {mode.getFrequency()}, Gain = {mode.getGain()}, Gain Number = {mode.getGainNumber()}")
        
    # ======= Test importing a single resolvent mode =======
    logger.info("Testing importing a single resolvent mode from file")
    importMode      = Mode(
        FEMSpaces.VMixed, 
        mesh,
        name='q_hat_import',
        isStateVector=True,
        analysis='Resolvent',
    )
    
    # Instantiate a reader and use to load mode
    r               = Reader(
        sourceDir           = param.Export.ExportFolder,
        needInterpolation   = False,
        felicsMeshFilePath  = None, 
        isComplex           = True,
        cacheData           = True,
    )
    r.bind_env(param, FEMSpaces)
    importMode.importData(
        r,
        param.Export.ExportFolder,
    )
    
    # Set the gain number of the imported mode to distinguish it
    # importMode.setGainNumber(importMode.getGainNumber()+100.0)
    
    
    # ======= TEST IMPORTING A MODE COLLECTION =======
    # Create empty mode collection to import into
    # logger.info("Testing importing a mode collection from directory")
    # importSolution     = ModeCollection(FEMSpaces.VMixed, mesh, analysis='Resolvent')
    
    # # Instantiate a reader and use to load mode collection
    # r               = Reader(
    #     sourceDir           = param.Export.ExportFolder, 
    #     needInterpolation   = False,
    #     felicsMeshFilePath  = None, 
    #     isComplex           = True,
    #     cacheData           = True,
    # )
    # # TEMPORARY: Bind the params and FEM spaces to the reader
    # r.bind_env(param, FEMSpaces)
    # # Import test for collection
    # importSolution.importData(
    #     r,
    #     param.Export.ExportFolder,
    # )
    
    # # NOTE: for now the gains are not stored in fles, so we cannot check them here
    # # Loop over all modes and print their attributes to verify successful import
    # for i, mode in enumerate(importSolution.modeList):
    #     logger.info(f"Imported Mode {i}: Frequency = {mode.getFrequency()}, Gain = {mode.getGain()}, Gain Number = {mode.getGainNumber()}")
    
    # # Change the gain number in new collection to distinguish imported modes
    # for i, mode in enumerate(importSolution.modeList):
    #     importSolution.modeList[i].setGainNumber(mode.getGainNumber()+10)
    
    # # Export again the mode collection so we can check the imported modes
    # fluctSolutList_new      = importSolution.getOldSolutionObject(meanFlow, param, FEMSpaces)
    # ExportFromFile(param,FEMSpaces,fluctSolutList_new,meanFlow)
