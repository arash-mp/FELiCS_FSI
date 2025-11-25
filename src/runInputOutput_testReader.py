import  time
from 	FELiCS.Misc.logging import Logger
from    FELiCS.IO.Reader        import Reader
from    FELiCS.Fields.Mode      import Mode

# Get the logger
logger                      = Logger.get_logger("felics")

def runInputOutput(param):
    '''This function runs the calculations preset in param
    Input:
        param: Parameter objects (see parameters.py), defining the case
    '''
    from    FELiCS.IO.Writer                    import Writer
    from    FELiCS.SpaceDisc.FEMSpaces          import FEMSpaces
    from    FELiCS.Fields.meanFlowClass         import meanFlowClass
    from    FELiCS.Equation.EquationCollection  import EquationCollectionClass
    from    FELiCS.Solvers.LinearSolver         import LinearSolver 
    from    FELiCS.Fields.ModeCollection        import ModeCollection

    logger.info("Running input/output analysis")
    
    #-----------------------------------------------------------------------
    ## INITIALIZATION
    #-----------------------------------------------------------------------
    # Get FELiCS objects required for analysis
    # Mesh
    mesh                    = param.getMesh()

    # FEMSpaces
    FEMSpaces               = FEMSpaces(
        param, 
        mesh,
    )
     
    # writer to export the results in files
    writer                  = Writer(
        mesh, 
        param.Export.ExportFolder,
    )

    # read in mean flow and export to h5-file
    meanFlow                = meanFlowClass(
        param, 
        FEMSpaces, 
        mesh,
    )
    meanFlow.importDataFromFileAndExportToH5(writer)

    # equation
    equation                = EquationCollectionClass(
        param,
        FEMSpaces,
        meanFlow,
        mesh,
    )

    #-----------------------------------------------------------------------
    ## MAIN PART
    #-----------------------------------------------------------------------
    # Get operators required for input/output analysis
    A                       = equation.getLinearOperator(meanFlow)
    B                       = equation.getWeightMatrix  (meanFlow)
    forcing                 = equation.getForcingForInputOutput(meanFlow) 

    # Get parameters for input/output analysis
    omegas                  = param.IOResolvent.Omegas

    # Start tracking time
    start                   = time.time()
    
    # Prepare solution object
    solution                = ModeCollection(
        FEMSpaces.VMixed, 
        mesh, 
        analysisType = "input_output",
    )

    # Solve equation system for each frequency 
    for omega in omegas:
        # define operator
        operator            = A.copy()
        operator.axpy(-omega, B) #petsc command: operator = A - omega*B
        solutionVector      = LinearSolver.solveEquationSystem(operator, forcing)
        
        solution.appendModeFromVector(
            solutionVector, 
            frequency = omega, 
            gain = 1
        ) 

        # Export last mode in collection
        solution.exportModes(writer, onlyNewN = 1)

    # End tracking time
    logger.info(f"Solving the input/output problem took {time.time()-start:.4g} s")
    
    #-----------------------------------------------------------------------
    ## TEST IMPORTING A MODE COLLECTION
    #-----------------------------------------------------------------------
    
    #-----------------------------------------------------------------------
    ## TEST IMPORTING A SINGLE MODE
    #-----------------------------------------------------------------------
    logger.info("Testing importing a single mode from file")
    
    # Create an empty mode to import into 
    # NOTE set isStateVector to False to test setting subnames manually
    importMode              = Mode(
        FEMSpaces.VMixed, 
        mesh,
        name                = 'q_hat_import',
        isStateVector       = False,    # So we can set subnames manually
        analysisType        = 'input_output',
    )

    # Test if the setSubFieldNames method works
    importMode.setNamesOfSubFields(solution.modeList[0].getNamesOfSubFields())
    
    # Set the frequency and mode number of the mode we want to import
    importMode.frequency    = solution.modeList[0].frequency
    
    # Instantiate a reader and use to load mode
    r                       = Reader(
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
    
    # Set the gain number of the imported mode to distinguish it
    # from others in the collection
    importMode.frequency   += 100.0
    solution.appendMode(importMode)
    
    # Export the new mode
    solution.exportModes(writer, onlyNewN = 1)

