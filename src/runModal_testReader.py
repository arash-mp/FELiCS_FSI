import  time
from 	FELiCS.Misc.logging     import Logger
from    FELiCS.IO.Reader        import Reader
from    FELiCS.Fields.Mode      import Mode

# Get the logger
logger = Logger.get_logger("felics")

def runModal(param):
    '''This function runs the calculations preset in param
    Input:
        param: Parameter objects (see parameters.py), defining the case
    '''
    # import  FELiCS.IO.Import as Import
    from    FELiCS.IO.Writer                    import Writer
    from    FELiCS.SpaceDisc.FEMSpaces          import FEMSpaces
    from    FELiCS.Fields.meanFlowClass         import meanFlowClass
    from    FELiCS.Equation.EquationCollection  import EquationCollectionClass
    from    FELiCS.Solvers.LinearSolver         import LinearSolver 
    from    FELiCS.Fields.ModeCollection        import ModeCollection
    
    logger.info("Running Modal analysis")
    #-----------------------------------------------------------------------
    ## INITIALIZATION
    #-----------------------------------------------------------------------
    # Mesh
    mesh            = param.getMesh()
    
    # FEMSpaces
    FEMSpaces       = FEMSpaces(
        param,
        mesh,
    )
    
    # Writer to export the results in files
    writer          = Writer(mesh, param.Export.ExportFolder)
     
    # Read in mean flow and export it to h5-file
    meanFlow        = meanFlowClass(
        param, 
        FEMSpaces, 
        mesh
    )
    meanFlow.importDataFromFileAndExportToH5(writer)

    # Equation
    equation        = EquationCollectionClass(
        param,
        FEMSpaces,
        meanFlow,
        mesh
    )

    #-----------------------------------------------------------------------
    ## MAIN PART
    #-----------------------------------------------------------------------
    # Get matrices for eigenproblem
    A               = equation.getLinearOperator(meanFlow)
    B               = equation.getWeightMatrix(meanFlow)

    # Get parameters for eigenproblem
    guesses         = param.Numerics.EigenValueGuess
    nSol            = param.Numerics.nSolut
    adjoint         = param.Case.CalculateAdjoint

    # track time
    start           = time.time()

    # solve eigenproblem for each guess
    solution        = ModeCollection(FEMSpaces.VMixed, mesh, analysisType = "modal")
    for guess in guesses:
        
        logger.info("Solving direct GEVP for guess: omega = " + str(guess))
        tmp         = LinearSolver.solveGeneralEigenproblem(
            A,
            B,
            guess,
            nSol,
        )
        solution.appendSolutionOfEigenProblem(tmp, guess)
        n_calculated = nSol

        if adjoint:
            logger.info("Solving adjoint GEVP for guess: omega = " + str(guess))
            tmp     = LinearSolver.solveGeneralEigenproblem(
                A,
                B,
                guess,
                nSol,
                adjoint=True
            )
            solution.appendSolutionOfEigenProblem(tmp, guess, adjoint=True)
            n_calculated += nSol

        # Exporting the (temporary) spectrum to a file
        solution.exportSpectrumToCSV(writer)
        # Exporting only the newly calculated modes to files (for this guess)
        solution.exportModes(writer, onlyNewN = n_calculated)

    # end tracking time
    end             = time.time() - start
    logger.info('Solving the general eigenproblem took %4g s' % end)
    residuum_max    = solution.getMaximumError()
    logger.debug('Maximum residuum of all solutions:  %12g' % (residuum_max))
    
    #-----------------------------------------------------------------------
    ## TEST IMPORTING A MODE COLLECTION
    #-----------------------------------------------------------------------
    # logger.info("Testing importing a mode collection from directory")
    
    # # Create empty mode collection to import into
    # importSolution  = ModeCollection(
    #     FEMSpaces.VMixed, 
    #     mesh
    # )
    
    # # Instantiate a reader and use to load mode collection
    # reader          = Reader(
    #     sourceDir           = param.Export.ExportFolder,
    #     needInterpolation   = False,
    #     felicsMeshFilePath  = None, 
    #     isComplex           = True,
    #     cacheData           = True,
    # )

    # # Import test for collection
    # importSolution.importData(
    #     reader,
    #     param.Export.ExportFolder,
    #     modeType            = 'Direct'
    # )
    
    # # Adding +100 to eigenvalues to distinguish imported modes
    # for mode in importSolution.modeList:
    #     mode.setEigenValue(mode.getEigenValue()+100.0)
        
    # # Export again the mode collection so we can check the imported modes
    # fluctSolutList_new      = importSolution.getOldSolutionObject(meanFlow, param, FEMSpaces)
    # ExportFromFile(param,FEMSpaces,fluctSolutList_new,meanFlow)
    
    #-----------------------------------------------------------------------
    ## TEST IMPORTING A SINGLE MODE
    #-----------------------------------------------------------------------
    logger.info("Testing importing a single mode from file")
    
    # Create an empty mode to import into 
    # NOTE set isStateVector to False to test setting subnames manually
    importMode      = Mode(
        FEMSpaces.VMixed, 
        mesh,
        name            = 'q_hat_import',
        isStateVector   = False,
        modeType        = 'Direct',
    )

    # Test if the setSubFieldNames method works
    importMode.setNamesOfSubFields(solution.modeList[0].getNamesOfSubFields())
    
    # Set the eigenvalue of the mode we want to import
    importMode.eigenValue = solution.modeList[0].eigenValue

    # Instantiate a reader and use to load mode
    r               = Reader(
        sourceDir           = param.Export.ExportFolder,
        needInterpolation   = False,
        felicsMeshFilePath  = None,
        isComplex           = True,
        cacheData           = True,
    )
    importMode.importData(
        r,
        param.Export.ExportFolder,
    )
    
    # Set a fake eigenvalue and append this mode to the collection
    importMode.eigenValue   = 999.999
    solution.appendMode(importMode)
    
    # Export the new mode
    solution.exportModes(writer, onlyNewN = 1)
    
    
    
