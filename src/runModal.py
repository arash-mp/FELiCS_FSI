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
    import  FELiCS.SpaceDisc.DefineFEMSpaces    as DefineFEMSpaces
    from    FELiCS.Fields.meanFlowClass         import meanFlowClass
    from    FELiCS.Equation.EquationCollection  import EquationCollectionClass
    from    FELiCS.Solvers.LinearSolver         import LinearSolver 
    from    FELiCS.Fields.ModeCollection        import ModeCollection

    logger.info("Running Modal analysis")
    #-----------------------------------------------------------------------
    ## INITIALIZATION
    #-----------------------------------------------------------------------
    # mesh
    mesh = param.getMesh()
    
    # FEMSpaces
    FEMSpaces = DefineFEMSpaces.FEMSpacesClass(
                param,
                mesh,
                )
     
    # read in mean flow and export to h5-file
    meanFlow = meanFlowClass(param, FEMSpaces, mesh)
    meanFlow.importDataFromFileAndExportToH5()
    
    # equation
    equation = EquationCollectionClass(
                                      param,
                                      FEMSpaces,
                                      meanFlow,
                                      mesh
                                      )

    #-----------------------------------------------------------------------
    ## MAIN PART
    #-----------------------------------------------------------------------
    # get matrices for eigenproblem
    A = equation.getLinearOperator(meanFlow)
    B = equation.getWeightMatrix  (meanFlow)

    # get parameters for eigenproblem
    guesses  = param.Numerics.EigenValueGuess
    nSol     = param.Numerics.nSolut
    adjoint  = param.Case.CalculateAdjoint

    # track time
    start= time.time()

    # solve eigenproblem for each guess
    solution = ModeCollection(FEMSpaces.VMixed, mesh)
    for guess in guesses:
        
        logger.info("Solving direct GEVP for guess: omega = " + str(guess))
        tmp     = LinearSolver.solveGeneralEigenproblem(A,
                                                        B,
                                                        guess,
                                                        nSol,
                                                        )

        solution.appendSolutionOfEigenProblem(tmp, guess)

        if adjoint:
            
            logger.info("Solving adjoint GEVP for guess: omega = " + str(guess))
            tmp = LinearSolver.solveGeneralEigenproblem(A,
                                                        B,
                                                        guess,
                                                        nSol,
                                                        adjoint=True)

            solution.appendSolutionOfEigenProblem(tmp, guess, adjoint=True)

    # end tracking time
    end = time.time() - start
    logger.info('Solving the general eigenproblem took %4g s' % end)
    residuum_max = solution.getMaximumError()
    logger.debug('Maximum residuum of all solutions:  %12g' % (residuum_max))

    #-----------------------------------------------------------------------
    ## EXPORT SOLUTION
    #-----------------------------------------------------------------------
    fluctSolutList = solution.getOldSolutionObject(meanFlow, param, FEMSpaces)
    ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)
