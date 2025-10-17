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
    from    FELiCS.IO.reader                    import Reader

    logger.info("Running Modal analysis")
    
    # Get FELiCS objects required for analysis
    mesh            = param.getMesh()
    FEMSpaces       = FEMSpaces(param, mesh)
    meanFlow        = meanFlowClass(param, FEMSpaces, mesh)
    meanFlow.importDataFromFileAndExportToH5()
    equation        = EquationCollectionClass(param,FEMSpaces,meanFlow,mesh)
    A               = equation.getLinearOperator(meanFlow)
    B               = equation.getWeightMatrix(meanFlow)

    # Parameters for eigenproblem
    guesses         = param.Numerics.EigenValueGuess
    nSol            = param.Numerics.nSolut
    adjoint         = param.Case.CalculateAdjoint

    # Solve eigenproblem for each guess
    solution        = ModeCollection(FEMSpaces.VMixed, mesh)
    for guess in guesses:
        logger.info("Solving direct GEVP for guess: omega = " + str(guess))
        tmp         = LinearSolver.solveGeneralEigenproblem(A,B,guess,nSol)
        solution.appendSolutionOfEigenProblem(tmp, guess)
        if adjoint:
            logger.info("Solving adjoint GEVP for guess: omega = " + str(guess))
            tmp     = LinearSolver.solveGeneralEigenproblem(A,B,guess,nSol,adjoint=True)
            solution.appendSolutionOfEigenProblem(tmp,guess,adjoint=True)

    # Export solutions to file
    fluctSolutList  = solution.getOldSolutionObject(meanFlow, param, FEMSpaces)
    ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)
    
    # Get a mode to check structure
    exampleMode     = solution.modeList[0]
    
    # ======= TEST IMPORTING A MODE COLLECTION =======
    # Create empty mode collection to import into
    logger.info("Testing importing a mode collection from directory")
    importSolution     = ModeCollection(FEMSpaces.VMixed, mesh)
    
    # Instantiate a reader and use to load mode collection
    r               = Reader(
        needInterpolation   = False,
        originalMeshFile    = None, 
        isComplex           = True,
        cacheData           = True,
    )
    # TEMPORARY: Bind the params and FEM spaces to the reader
    r.bind_env(param, FEMSpaces)
    # Import test for collection
    importSolution.importData(
        r,
        param.Export.ExportFolder,
        modeType='Direct'
    )
    # Adding +100 to eigenvalues to distinguish imported modes
    for mode in importSolution.modeList:
        mode.setEigenValue(mode.getEigenValue()+100.0)
    # Export again the mode collection so we can check the imported modes
    fluctSolutList_new      = importSolution.getOldSolutionObject(meanFlow, param, FEMSpaces)
    ExportFromFile(param,FEMSpaces,fluctSolutList_new,meanFlow)
    
    # ======= TEST IMPORTING A MODE FROM FILE =======
    # logger.info("Testing importing a mode from file")
    
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
    # importMode.setEigenValue(exampleMode.getEigenValue())

    # # Instantiate a reader and use to load mode
    # r               = Reader(
    #     needInterpolation   = False,
    #     originalMeshFile    = None,
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
    
    
