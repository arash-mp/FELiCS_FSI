import  time
from 	FELiCS.Misc.logging import Logger

# Get the logger
logger = Logger.get_logger("felics")

def runInputOutput(param):
    '''This function runs the calculations preset in param
    Input:
        param: Parameter objects (see parameters.py), defining the case
    '''
    # import FELiCS.IO.Import as Import
    from    FELiCS.IO.ExportSolution            import ExportFromFile 
    from    FELiCS.SpaceDisc.FEMSpaces          import FEMSpaces
    from    FELiCS.Fields.meanFlowClass         import meanFlowClass
    from    FELiCS.Equation.EquationCollection  import EquationCollectionClass
    from    FELiCS.Misc.functions               import printDebug
    from    FELiCS.Solvers.LinearSolver         import LinearSolver 
    from    FELiCS.Fields.ModeCollection        import ModeCollection
    # from   FELiCS.Fields.fluctuationClass import fluctuationSolutions
    # from   FELiCS.Fields.Mode import Mode

    logger.info("Running input/output analysis")
    #-----------------------------------------------------------------------
    ## INITIALIZATION
    #-----------------------------------------------------------------------
    # mesh
    mesh = param.getMesh()
    
    # FEMSpaces
    FEMSpaces = FEMSpaces(
        param,
        mesh,
    )
    
    # read in mean flow and export to h5-file
    meanFlow = meanFlowClass(param, FEMSpaces, mesh)
    meanFlow.importData()
    
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
    A       = equation.getLinearOperator(meanFlow)
    B       = equation.getWeightMatrix  (meanFlow)
    forcing = equation.getForcingForInputOutput(meanFlow) 

    # get parameters for input/output analysis
    omegas   = param.IOResolvent.Omegas
    names    = param.Case.StateVectorVariables

    # track time
    start= time.time()

    # solve equation system 
    solution = ModeCollection(FEMSpaces.VMixed, mesh)
    for omega in omegas:
        # define operator
        operator = A.copy()
        operator.axpy(-omega, B) #petsc command: operator = A - omega*B

        solutionVector = LinearSolver.solveEquationSystem(operator, forcing)
        
        solution.appendModeFromVector(solutionVector, frequency = omega, gain = 1) 


    # end tracking time
    end = time.time() - start
    logger.info('Solving the input/output problem took %4g s' % end)

    #-----------------------------------------------------------------------
    ## EXPORT SOLUTION
    #-----------------------------------------------------------------------
    fluctSolutList = solution.getOldSolutionObject(meanFlow, param, FEMSpaces)
    ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)
