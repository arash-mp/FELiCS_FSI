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
        if not adjoint:
            logger.info("Solving direct GEVP for guess: omega = " + str(guess))
            tmp     = LinearSolver.solveGeneralEigenproblem(A,B,guess,nSol)
            solution.appendSolutionOfEigenProblem(tmp, guess)
        elif adjoint:
            logger.info("Solving adjoint GEVP for guess: omega = " + str(guess))
            tmp     = LinearSolver.solveGeneralEigenproblem(A,B,guess,nSol,adjoint=True)
            solution.appendSolutionOfEigenProblem(tmp,guess,adjoint=True)

    # Export solutions to file
    fluctSolutList  = solution.getOldSolutionObject(meanFlow, param, FEMSpaces)
    ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)
    
    # Get a mode to check structure
    exampleMode     = solution.modeList[0]
    
    # Create an empty mode
    importMode      = Mode(FEMSpaces,mesh,name=[],isStateVector=True)

    # Instantiate a reader and use to load mode
    r               = Reader(
        needInterpolation   = False,
        originalMeshFile    = None,
        isComplex           = True,
        cacheData           = True,
    )
    importMode.importData(
        r,
    )
