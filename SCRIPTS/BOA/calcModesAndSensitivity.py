import pdb
import numpy as np
import copy
import time
import sys
import dolfinx

from   FELiCS.Parameters.parameters import parameters

import FELiCS.IO.Import as Import
from   FELiCS.IO.ExportSolution import ExportGUI,ExportFromFile

import FELiCS.SpaceDisc.DefineFEMSpaces as DefineFEMSpaces
from   FELiCS.Fields.meanFlowClass import meanFlowClass
from   FELiCS.Fields.fluctuationClass import fluctuationSolutions
from   FELiCS.Equation.EquationCollection import EquationCollectionClass
from   FELiCS.Misc.functions import printDebug

from   FELiCS.Solvers.LinearSolver import LinearSolver 
from   FELiCS.Fields.ModeCollection import ModeCollection
from   FELiCS.Fields.Field import Field
from   FELiCS.Fields.Mode import Mode
from   FELiCS.Misc.tensorUtils import Tensor

from   CaseHandler import CaseHandler


def calculateModesAndSensitivity(settingsFileName, baseFlow_array, baseFlowSensitivity, optimizerParameters):

    #-----------------------------------------------------------------------
    ## INITIALIZATION 
    #-----------------------------------------------------------------------
    # read parameters
    param=parameters()
    param.importFromFile(settingsFileName)
    param.getOldParameters()

    # mesh
    mesh=param.BCs.getMesh()
   
    # FEMSpaces
    FEMSpaces = DefineFEMSpaces.FEMSpacesClass(
                param,
                mesh,
                )
    
    
    # initialize mean flow class
    meanFlow = meanFlowClass(param, FEMSpaces, mesh)
    
    #import mean flow data from file
    meanFlow.importDataFromFile()

    # load base flow into meanFlow object
    if np.linalg.norm(meanFlow._fieldDict['u'].x.array[:]) < 1.e-8:
        baseFlow = Field(FEMSpaces.VMixed, mesh)
        baseFlow.setCoefficientArray(baseFlow_array)
        [u,p] = baseFlow.getListOfSingleFields()
        meanFlow._fieldDict['u'] = u.function

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

    caseHandler = CaseHandler(settingsFileName, param, mesh, equation.boundaries)

    #-----------------------------------------------------------------------
    ## SOLVE EIGENPROBLEM AND SCALE LEADING MODES
    #-----------------------------------------------------------------------
    # get matrices for eigenproblem
    A  = equation.getLinearOperator(meanFlow)
    B  = equation.getWeightMatrix  (meanFlow)
    
    # get parameters for eigenproblem
    guesses  = param.Numerics.EigenValueGuess
    nSol     = param.Numerics.nSolut
    adjoint  = param.Case.CalculateAdjoint
    
    # track time
    start= time.time()
    
    # solve direct eigenproblem for each guess
    solution = ModeCollection(FEMSpaces.VMixed, mesh)
    for guess in guesses:
        tmp     = LinearSolver.solveGeneralEigenproblem(A,
                                                        B,
                                                        guess,
                                                        nSol,
                                                        )
        solution.appendSolutionOfEigenProblem(tmp, guess)
    
    
    mode_direct       = solution.getLeadingMode()
    guess             = mode_direct.getEigenValue()
    
    # solve adjoint eigenproblem only for the leading eigenvalue
    tmp = LinearSolver.solveGeneralEigenproblem(A,
                                                B,
                                                guess,
                                                1,
                                                adjoint=True)
    solution.appendSolutionOfEigenProblem(tmp, guess, adjoint=True)
    
    
    # end tracking time
    end = time.time() - start
    printDebug(True, '-- Solving the general eigenproblem took %4g s' % end)
    residuum_max = solution.getMaximumError()
    printDebug(True, '-- Maximum residuum of all solutions:  %12g' % (residuum_max))
    
    # get leading modes
    mode_direct  = solution.getLeadingMode(adjoint=False)
    mode_adjoint = solution.getLeadingMode(adjoint=True)
    eigenValue   = mode_direct.getEigenValue()
    printDebug(True, '------------------------------------------------ ')
    printDebug(True, '-- Leading eigenvalue:  ' + str(eigenValue))
    printDebug(True, '------------------------------------------------ ')
    
    
    
    # scale modes s.t. mode_adjoint^H * B * mode_direct = 1
    mode_direct_petsc  = mode_direct.getPetscVector()
    mode_adjoint_petsc = mode_adjoint.getPetscVector()
    temp               = mode_direct.getPetscVector() # gets a petsc vector "temp" of correct length
    B.mult(mode_direct_petsc, temp)                   # temp = B*mode_direct
    factor = temp.dot(mode_adjoint_petsc)             # factor = mode_adjoint ^H temp
    
    mode_direct_petsc.scale(1./np.sqrt(factor))
    mode_adjoint_petsc.scale(np.conj(1./np.sqrt(factor)))
    mode_direct.setCoefficientArray(mode_direct_petsc.getArray())
    mode_adjoint.setCoefficientArray(mode_adjoint_petsc.getArray())
   
  
    # export leading modes in standard felics format
    solution_onlyLeading = ModeCollection(FEMSpaces.VMixed, mesh)
    solution_onlyLeading.appendMode(mode_direct)
    solution_onlyLeading.appendMode(mode_adjoint)
    fluctSolutList = solution_onlyLeading.getOldSolutionObject(meanFlow, param, FEMSpaces)
    ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)



    #-----------------------------------------------------------------------
    ## CALCULATE GRADIENT WITH RESPECT TO GEOMETRY DEFORMATIONS - Part I
    #-----------------------------------------------------------------------
    A_0  = equation.getLinearOperator(meanFlow)
    
    # create and read fields
    a_i  = optimizerParameters 
    
    # deform mesh
    epsilon = 1.e-8
    
    geometryDeformer = caseHandler.getGeometryDeformer() 
    N_param          = geometryDeformer.getNumberOfParameters()
    sensitivity1     = np.zeros(N_param,dtype=complex)
    sensitivity2     = np.zeros(N_param,dtype=complex)
    baseFlowSens     = Field(FEMSpaces.VMixed, mesh)
    for i in range(0,N_param):
        a_i[i] = a_i[i] + epsilon
    
        geometryDeformer.deformMesh(a_i)
        
        A_deformed  = equation.getLinearOperator(meanFlow)
        A_deformed.axpy(-1., A_0)
        A_deformed.scale(1./epsilon)
        mode_direct_petsc  = mode_direct.getPetscVector()
        mode_adjoint_petsc = mode_adjoint.getPetscVector()
        result             = mode_direct.getPetscVector()
        A_deformed.mult(mode_direct_petsc, result)
        dolfinx.fem.petsc.set_bc(result, equation.BCs)
        sensitivity1[i]   = result.dot(mode_adjoint_petsc)
    

        baseFlowSens.setCoefficientArray(baseFlowSensitivity[i])
        [u_direct, p_direct]              = baseFlowSens.getListOfSingleFields()
        meanFlow._fieldDict['u_bilinear'] = u_direct.function
        BL                                = equation.getBilinearOperator(meanFlow)
        BL.mult(mode_direct_petsc, result)
        dolfinx.fem.petsc.set_bc(result, equation.BCs)
        sensitivity2[i]  = result.dot(mode_adjoint_petsc)


        a_i[i] = a_i[i] - epsilon
        geometryDeformer.restoreMesh()
    
    
    printDebug(True, '------------------------------------------------ ')
    printDebug(True, '-- sensitivities part 1:  ' + str((sensitivity1)))
    printDebug(True, '------------------------------------------------ ')
    printDebug(True, '------------------------------------------------ ')
    printDebug(True, '-- sensitivities part 2:  ' + str((sensitivity2)))
    printDebug(True, '------------------------------------------------ ')


    return sensitivity1, sensitivity2, mode_direct.getEigenValue()


    

