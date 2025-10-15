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
     
    # Read in mean flow and export to h5-file
    meanFlow        = meanFlowClass(
        param, 
        FEMSpaces, 
        mesh
    )
    meanFlow.importDataFromFileAndExportToH5()
   


    # get one field from meanflow
    writer          = Writer(mesh, exportFolder="OutputTest", exportMesh = FEMSpaces.exportMesh)
    field_scalar    = meanFlow._fieldDict['nulam']
    field_scalar.setCoefficientArray(field_scalar.getCoefficientArray()*2. + 1j*field_scalar.getCoefficientArray())
    field_scalar.exportToH5(writer)
    field_vector    = meanFlow._fieldDict['u']
    field_vector.exportToH5(writer)
    exit()



    exportField     = writer.getFieldsOnExportMesh(field_vector)
    writer.writeFieldToXDMF(exportField, "test_u")
    exit()
    field_scalar    = meanFlow._fieldDict['nulam']
    exportField     = writer.getFieldsOnExportMesh(field_scalar)
    writer.writeFieldToXDMF(exportField, "test_nulam")

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
    # get matrices for eigenproblem
    A               = equation.getLinearOperator(meanFlow)
    B               = equation.getWeightMatrix  (meanFlow)

    # get parameters for eigenproblem
    guesses         = param.Numerics.EigenValueGuess
    nSol            = param.Numerics.nSolut
    adjoint         = param.Case.CalculateAdjoint
    names           = param.Case.StateVectorVariables

    # track time
    start           = time.time()

    # solve eigenproblem for each guess
    solution        = ModeCollection(
        FEMSpaces.VMixed, 
        mesh, 
        names=names
    )
    
    for guess in guesses:
        
        logger.info("Solving direct GEVP for guess: omega = " + str(guess))
        tmp         = LinearSolver.solveGeneralEigenproblem(
            A,
            B,
            guess,
            nSol,
        )

        solution.appendSolutionOfEigenProblem(tmp, guess)
        
        mode        = solution.modeList[0]
        list_modeSubFields = mode.getListOfSubFields()
        exportField = writer.getFieldsOnExportMesh(mode)
        writer.writeFieldToXDMF(exportField, "test_mode")
        mode_u      = list_modeSubFields[0]
        exportField = writer.getFieldsOnExportMesh(mode_u)
        writer.writeFieldToXDMF(exportField, "test_mode_u")

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
