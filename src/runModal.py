import pdb
import numpy as np
import copy
import time

def runModal(param, useGUI):
    '''This function runs the calculations preset in param
    Input:
        param: Parameter objects (see parameters.py), defining the case
        useGUI: Boolean, True if program is run using GUI, False if run from
        terminal directly
    '''
    import FELiCS.IO.Import as Import
    from   FELiCS.IO.ExportSolution import ExportGUI,ExportFromFile

    import FELiCS.SpaceDisc.DefineFEMSpaces as DefineFEMSpaces
    from   FELiCS.Fields.meanFlowClass import meanFlowClass
    from   FELiCS.Fields.fluctuationClass import fluctuationSolutions
    from   FELiCS.Equation.EquationCollection import EquationCollectionClass
    from   FELiCS.Misc.functions import printDebug

    from   FELiCS.Solvers.LinearSolver import LinearSolver 

    #-----------------------------------------------------------------------
    ## INITIALIZATION
    #-----------------------------------------------------------------------
    # mesh
    mesh=param.BCs.getMesh()
    # FEMSpaces
    FEMSpaces = DefineFEMSpaces.FEMSpacesClass(
                param,
                mesh,
                )
    # read in mean flow
    meanFlow = meanFlowClass(param, FEMSpaces, mesh)
    meanFlow.importDataFromFile()
    if not param.FlowInput.MeanFlowFilePath.split('.')[-1] == 'hdf5':
        meanFlow.exportBaseFlowAsHDF5()
    # export mean flow in "h5" file
    meanflowFilename = 'meanflow.h5'
    meanFlow.mapToExportMeshAndExport(FEMSpaces, meanflowFilename)
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
    #solution = ...
    for guess in guesses:
        solution_eigenProblem = LinearSolver.solveGeneralEigenproblem(A,
                                                                      B,
                                                                      guess,
                                                                      nSol,
                                                                      adjoint=False)
        #solution.append(solution_eigenProblem)....
        if adjoint==True:
            solution_eigenProblem_adjoint = LinearSolver.solveGeneralEigenproblem(A,
                                                                      B,
                                                                      guess,
                                                                      nSol,
                                                                      adjoint=True)

            #solution.append(solution_eigenProblem_adjoint)....

    # end tracking time
    end = time.time() - start
    printDebug(True, '-- Solving the general eigenproblem took %4g s' % end)
             

    LinearAlgebraObj = equation.DiscretizeFlow()
    fluctSolutList = LinearAlgebraObj.\
                     solveGEVP(param.Case.CalculateAdjoint)


    #-----------------------------------------------------------------------
    ## EXPORT SOLUTION
    #-----------------------------------------------------------------------
    if useGUI:
        ExportGUI(param, fluctSolutList, meanFlow,FEMSpaces, equation,mesh)
    else:
        ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)

