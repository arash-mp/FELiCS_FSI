import pdb
import numpy as np
import copy
import time

def runInputOutput(param, useGUI):
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
    from   FELiCS.Fields.ModeCollection import ModeCollection
    from   FELiCS.Fields.Mode import Mode


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
    # export mean flow in "h5" file
    if not param.FlowInput.MeanFlowFilePath.split('.')[-1] == 'hdf5':
        meanFlow.exportBaseFlowAsHDF5()
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
    A       = equation.getLinearOperator(meanFlow)
    B       = equation.getWeightMatrix  (meanFlow)
    forcing = equation.getForcingForInputOutput(meanFlow) 

    # get parameters for input/output analysis
    omegas   = param.IOResolvent.Omegas

    exit()

    # track time
    start= time.time()

    # solve equation system 
    solution = ModeCollection(FEMSpaces.VMixed, mesh)
    for omega in omegas:
        # define operator
        operator = A.copy()
        operator.axpy(-omega, B) #petsc command: operator = A - omega*B

        solutionVector = LinearSolver.solveEquationSystem(operator, forcing)
        
        solution.appendModeFromVector(solutionVector, gain=1) 

#                # Use PETSC to solve linear system (with LU decomposition)
#                eigenvectors_c   = self.__solveEquationSystem(OP, self.__matrix_dict_petsc['b_forcing'])
#                gains[:,i]       = 1
#                responses[:,0,i] = eigenvectors_c
#
#        for i in range(len(self.__param.IOResolvent.Omegas)):
#            # construct for each EVal and EVec a fluctuationSolution
#            for gainNumb in range(responses.shape[1]):
#                fluctSolutObjList.append(fluctuationSolutions(
#                        self.__param,
#                        self.__meanFlow,
#                        self.__FEMSpaces,
#                        self.__param.IOResolvent.Omegas[i],
#                        responses[:,gainNumb,i],
#                        True,
#                        gainNumb,
#                        gains[gainNumb,i]
#                        ))
#



    # end tracking time
    end = time.time() - start
    printDebug(True, '-- Solving the input/output problem took %4g s' % end)
    residuum_max = solution.getMaximumError()
    printDebug(True, '-- Maximum residuum of all solutions:  %12g' % (residuum_max))




    #-----------------------------------------------------------------------
    ## EXPORT SOLUTION
    #-----------------------------------------------------------------------
    fluctSolutList = solution.getOldSolutionObject(meanFlow, param, FEMSpaces)
    if useGUI:
        ExportGUI(param, fluctSolutList, MeanFlow,FEMSpaces, equationColl,mesh)
    else:
        ExportFromFile(param,FEMSpaces,fluctSolutList,MeanFlow)
