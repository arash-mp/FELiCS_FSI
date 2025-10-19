import  time
from 	FELiCS.Misc.logging import Logger

# Get the logger
logger                      = Logger.get_logger("felics")

def runInputOutput(param):
    '''This function runs the calculations preset in param
    Input:
        param: Parameter objects (see parameters.py), defining the case
    '''
    from    FELiCS.IO.ExportSolution            import ExportFromFile 
    from    FELiCS.SpaceDisc.FEMSpaces          import FEMSpaces
    from    FELiCS.Fields.meanFlowClass         import meanFlowClass
    from    FELiCS.Equation.EquationCollection  import EquationCollectionClass
    from    FELiCS.Solvers.LinearSolver         import LinearSolver 
    from    FELiCS.Fields.ModeCollection        import ModeCollection

    logger.info("Running input/output analysis")
    
    # Get FELiCS objects required for analysis
    mesh                    = param.getMesh()
    FEMSpaces               = FEMSpaces(param, mesh)
    meanFlow                = meanFlowClass(param, FEMSpaces, mesh)
    meanFlow.importDataFromFileAndExportToH5()
    equation                = EquationCollectionClass(param,FEMSpaces,meanFlow,mesh)

    # Get operators required for input/output analysis
    A                       = equation.getLinearOperator(meanFlow)
    B                       = equation.getWeightMatrix  (meanFlow)
    forcing                 = equation.getForcingForInputOutput(meanFlow) 

    # Get parameters for input/output analysis
    omegas                  = param.IOResolvent.Omegas

    # Start tracking time
    start                   = time.time()

    # Solve equation system for each frequency 
    solution                = ModeCollection(FEMSpaces.VMixed, mesh)
    for omega in omegas:
        # define operator
        operator            = A.copy()
        operator.axpy(-omega, B) #petsc command: operator = A - omega*B
        solutionVector      = LinearSolver.solveEquationSystem(operator, forcing)
        solution.appendModeFromVector(solutionVector, frequency = omega, gain = 1) 

    # End tracking time
    logger.info(f"Solving the input/output problem took {time.time()-start:.4g} s")

    # Exporting the solution to file
    fluctSolutList          = solution.getOldSolutionObject(
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
