import  time
from 	FELiCS.Misc.logging import Logger

# Get the logger
logger = Logger.get_logger("felics")

def runModal(param):
    '''This function runs the calculations preset in param
    Input:
        param: Parameter objects (see parameters.py), defining the case
    '''
    # import  FELiCS.IO.Import as Import
    from    FELiCS.IO.ExportSolution            import ExportFromFile 
    from    FELiCS.SpaceDisc.FEMSpaces          import FEMSpaces
    from    FELiCS.Fields.meanFlowClass         import meanFlowClass
    from    FELiCS.Equation.EquationCollection  import EquationCollectionClass
    from    FELiCS.Solvers.LinearSolver         import LinearSolver 
    from    FELiCS.Fields.ModeCollection        import ModeCollection
    from    FELiCS.Fields.Mode                  import Mode
    from    FELiCS.IO.Reader                    import Reader

    logger.info("Running Modal analysis")
    #-----------------------------------------------------------------------
    ## INITIALIZATION
    #-----------------------------------------------------------------------
    # Mesh
    mesh            = param.getMesh()
    mesh.saveInFELiCSFormat(f'{param.Export.ExportFolder}/mesh.h5')
    
    # FEMSpaces
    FEMSpaces       = FEMSpaces(
        param,
        mesh,
    )
     
    # read in mean flow and export to h5-file
    meanFlow        = meanFlowClass(
        param, 
        FEMSpaces, 
        mesh
    )
    meanFlow.importDataFromFileAndExportToH5()

    # equation
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
    solution       = ModeCollection(FEMSpaces.VMixed, mesh)
    for guess in guesses:
        
        logger.info("Solving direct GEVP for guess: omega = " + str(guess))
        tmp         = LinearSolver.solveGeneralEigenproblem(
            A,
            B,
            guess,
            nSol,
        )
        solution.appendSolutionOfEigenProblem(tmp, guess)

        if adjoint:
            logger.info("Solving adjoint GEVP for guess: omega = " + str(guess))
            tmp     = LinearSolver.solveGeneralEigenproblem(
                A,
                B,
                guess,
                nSol,
                adjoint = True
            )
            solution.appendSolutionOfEigenProblem(tmp, guess, adjoint=True)

    # Exporting the eigenvalue spectrum to a file
    spectrumFile    = param.Export.ExportFolder + "/spectrum.csv"
    solution.exportSpectrumToCSV(spectrumFile)

    # end tracking time
    end             = time.time() - start
    logger.info('Solving the general eigenproblem took %4g s' % end)
    residuum_max    = solution.getMaximumError()
    logger.debug('Maximum residuum of all solutions:  %12g' % (residuum_max))
    
    #-----------------------------------------------------------------------
    ## EXPORT SOLUTION
    #-----------------------------------------------------------------------
    fluctSolutList  = solution.getOldSolutionObject(meanFlow, param, FEMSpaces)
    ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)
    
    #-----------------------------------------------------------------------
    ## TEST IMPORTING A MODE COLLECTION
    #-----------------------------------------------------------------------
    logger.info("Testing importing a mode collection from directory")
    
    # Create empty mode collection to import into
    importSolution  = ModeCollection(
        FEMSpaces.VMixed, 
        mesh
    )
    
    # Instantiate a reader and use to load mode collection
    reader          = Reader(
        sourceDir           = param.Export.ExportFolder,
        needInterpolation   = False,
        felicsMeshFilePath  = None, 
        isComplex           = True,
        cacheData           = True,
    )

    # Import test for collection
    importSolution.importData(
        reader,
        param.Export.ExportFolder,
        modeType            = 'Direct'
    )
    
    # Adding +100 to eigenvalues to distinguish imported modes
    for mode in importSolution.modeList:
        mode.setEigenValue(mode.getEigenValue()+100.0)
        
    # Export again the mode collection so we can check the imported modes
    fluctSolutList_new      = importSolution.getOldSolutionObject(meanFlow, param, FEMSpaces)
    ExportFromFile(param,FEMSpaces,fluctSolutList_new,meanFlow)
    
    #-----------------------------------------------------------------------
    ## TEST IMPORTING A SINGLE MODE
    #-----------------------------------------------------------------------
    # logger.info("Testing importing a single mode from file")
    
    # # Create an empty mode to import into 
    # # NOTE set isStateVector to False to test setting subnames manually
    # importMode      = Mode(
    #     FEMSpaces.VMixed, 
    #     mesh,
    #     name='q_hat_import',
    #     isStateVector=False
    # )
    # # Test if the setSubFieldNames method works
    # importMode.setNamesOfSubFields(['u', 'p'])
    
    # # Set the eigenvalue of the mode we want to import
    # importMode.setEigenValue(999.)

    # # Instantiate a reader and use to load mode
    # r               = Reader(
    #     sourceDir           = param.Export.ExportFolder,
    #     needInterpolation   = False,
    #     felicsMeshFilePath  = None,
    #     isComplex           = True,
    #     cacheData           = True,
    # )
    # # TEMPORARY: Bind the params and FEM spaces to the reader
    # r.bind_env(param, FEMSpaces)
    # importMode.importData(
    #     r,
    #     param.Export.ExportFolder,
    # )
    
    # # Set a fake eigenvalue and append this mode to the collection
    # importMode.setEigenValue(999.999)
    # solution.appendMode(importMode)
    
    # # Export again the mode collection so we can check the imported mode
    # fluctSolutList_new      = solution.getOldSolutionObject(meanFlow, param, FEMSpaces)
    # ExportFromFile(param,FEMSpaces,fluctSolutList_new,meanFlow)
    
    
