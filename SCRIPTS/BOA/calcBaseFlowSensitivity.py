import pdb
import numpy as np
import copy
import time
import sys
import dolfinx

from   FELiCS.Parameters.parameters import parameters

import FELiCS.IO.Import as Import
from   FELiCS.IO.ExportSolution import ExportGUI,ExportFromFile
import FELiCS.IO.ExportSolution as Export 

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

def calculateBaseFlowSensitivity(settingsFileName, baseFlow_array, optimizerParameters, deformed = False):

    param=parameters()
    param.importFromFile(settingsFileName)
    param.getOldParameters()
    
    
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
    
    
    # initialize mean flow class
    meanFlow = meanFlowClass(param, FEMSpaces, mesh)
    
    #import mean flow data from file
    meanFlow.importDataFromFile()


    # export mean flow in "h5" file
    if not param.FlowInput.MeanFlowFilePath.split('.')[-1] == 'hdf5':
        meanFlow.exportBaseFlowAsHDF5()
    meanflowFilename = 'meanflow.h5'
    meanFlow.mapToExportMeshAndExport(FEMSpaces, meanflowFilename)


    # load base flow into meanFlow object if the mean flow field is zero
    baseFlow = Field(FEMSpaces.VMixed, mesh)
    if np.linalg.norm(meanFlow._fieldDict['u'].x.array[:]) < 1.e-8:
        baseFlow.setCoefficientArray(baseFlow_array)
        [u,p] = baseFlow.getListOfSingleFields()
        meanFlow._fieldDict['u'] = u.function
        meanFlow._fieldDict['p'] = p.function
    elif np.linalg.norm(meanFlow._fieldDict['p'].x.array[:]) < 1.e-8:
        print('Error: for the sensitivities the pressure field is needed, please provide it alongside the velocity field.')
        exit()

    # equation
    equation = EquationCollectionClass(
                                      param,
                                      FEMSpaces,
                                      meanFlow,
                                      mesh
                                      )
    
    caseHandler = CaseHandler(settingsFileName, param, mesh, equation.boundaries)

    # deform mesh at beginning - only for testing /debugging
    if deformed == True:
        geometryDeformer = caseHandler.getGeometryDeformer()
        geometryDeformer.deformMesh(optimizerParameters)
        geometryDeformer.isDeformed = False


    #-----------------------------------------------------------------------
    ## PUT MEAN FIELD IN BASE FLOW OBJECT (WHEN CALCULATED EXTERNALLY)
    #-----------------------------------------------------------------------
    # interlopate pressure onto P1 space
    p_mean = dolfinx.fem.Function(FEMSpaces.P1)
    p_mean.interpolate(meanFlow._fieldDict['p'])
    meanFlow._fieldDict['p'] = p_mean

    # fill base flow field object
    mapping            = baseFlow.space.sub(0).sub(0).collapse()[1]
    mappingu           = meanFlow._fieldDict['u']._V.sub(0).collapse()[1]
    baseFlow.function.x.array[mapping] = meanFlow._fieldDict['u'].x.array[mappingu]
    mapping            = baseFlow.space.sub(0).sub(1).collapse()[1]
    mappingu           = meanFlow._fieldDict['u']._V.sub(1).collapse()[1]
    baseFlow.function.x.array[mapping] = meanFlow._fieldDict['u'].x.array[mappingu]
    mapping            = baseFlow.space.sub(1).collapse()[1]
    baseFlow.function.x.array[mapping] = meanFlow._fieldDict['p'].x.array[:]



    #-----------------------------------------------------------------------
    ## CALCULATE delta N / delta a_i and dq/da_i: 
    #-----------------------------------------------------------------------
    # set target function for nonlinear sponge
    targetValues = caseHandler.getTargetValuesForSponge()
    [u_t,p_t] = baseFlow.getListOfSingleFields()
    mapping_ux = u_t.space.sub(0).collapse()[1]
    mapping_uy = u_t.space.sub(1).collapse()[1]
    u_t.function.x.array[mapping_ux] = targetValues[0]
    u_t.function.x.array[mapping_uy] = targetValues[1]
    p_t.function.x.array[:]          = targetValues[2]
    # set wrong target values s.t. the sponge has no effect
    meanFlow._fieldDict['u_target'] = meanFlow._fieldDict['u']
    meanFlow._fieldDict['p_target'] = meanFlow._fieldDict['p']
    #meanFlow._fieldDict['u_target'] = u_t.function
    #meanFlow._fieldDict['p_target'] = p_t.function

    p_mean = dolfinx.fem.Function(FEMSpaces.P1)
    p_mean.interpolate(meanFlow._fieldDict['p'])
    meanFlow._fieldDict['p'] = p_mean

    # get reference nonlinear expression
    N_0  = equation.getNonlinearExpression(meanFlow)    

    ##########
    # smooth and treat with sponge => does not help for 1st order;  
    D        = equation.getFEMDiffusionMatrix(3.e-6, sponge=None)#meanFlow._fieldDict['spg'])
    solver_D = LinearSolver.createEquationSystemSolver(D)
    #residuum_vec  = LinearSolver.solveEquationSystemWithPredefinedSolver(solver_D, N_0)
    #residuum = Field(FEMSpaces.VMixed,mesh)
    #residuum.setCoefficientArray(residuum_vec)
    #N_0  = equation.getFullRHS(residuum.function)
    ##########

    # get linear solver  
    A        = equation.getLinearOperator(meanFlow)
    solver_A = LinearSolver.createEquationSystemSolver(A)

    printDebug(True, "-------------------------------------------------------------" )
    printDebug(True, "-- Calculating new base flow...")
    printDebug(True, "-------------------------------------------------------------" )

    # 1st order
    baseFlow_sens_1_array  = LinearSolver.solveEquationSystemWithPredefinedSolver(solver_A, N_0)

    baseFlow_sens_1 = Field(FEMSpaces.VMixed, mesh)
    baseFlow_sens_1.setCoefficientArray(baseFlow_sens_1_array)

    ###########
    ## smooth and treat with sponge 
    #rhs = equation.getFullRHS(baseFlow_sens_1.function)
    #baseFlow_sens_1_array  = LinearSolver.solveEquationSystemWithPredefinedSolver(solver_D, rhs)
    #baseFlow_sens_1.setCoefficientArray(baseFlow_sens_1_array)
    ###########


    # 2nd order
    [u_sens, p_sens]                  = baseFlow_sens_1.getListOfSingleFields()
    meanFlow._fieldDict['u_bilinear'] = u_sens.function
    A_q = equation.getBilinearOperator(meanFlow)
    vec1, vec2 = A_q.getVecs() 
    vec1.setValues(range(0,len(baseFlow_sens_1_array[:])), baseFlow_sens_1_array[:])
    A_q.multTranspose(vec1, vec2)

    ############
    ### smooth and treat with sponge 
    ##residuum_vec  = LinearSolver.solveEquationSystemWithPredefinedSolver(solver_D, vec2)
    ##residuum = Field(FEMSpaces.VMixed,mesh)
    ##residuum.setCoefficientArray(residuum_vec)
    ##vec2  = equation.getFullRHS(residuum.function)
    ############

    baseFlow_sens_2_array  = LinearSolver.solveEquationSystemWithPredefinedSolver(solver_A, vec2)
    baseFlow_sens_2 = Field(FEMSpaces.VMixed, mesh)
    baseFlow_sens_2.setCoefficientArray(baseFlow_sens_2_array)

    ###########
    ## smooth and treat with sponge => does not help for 1st order;  
    #rhs = equation.getFullRHS(baseFlow_sens_2.function)
    #baseFlow_sens_2_array  = LinearSolver.solveEquationSystemWithPredefinedSolver(solver_D, rhs)
    #baseFlow_sens_2.setCoefficientArray(baseFlow_sens_2_array)
    ###########


    # update base flow and mean flow
    baseFlow.setCoefficientArray(baseFlow.getCoefficientArray() + baseFlow_sens_1_array[:] + 0.5*baseFlow_sens_2_array[:])
    #baseFlow.setCoefficientArray(baseFlow.getCoefficientArray() + baseFlow_sens_1_array[:])
    [u,p] = baseFlow.getListOfSingleFields()
    meanFlow._fieldDict['u'] = u.function
    meanFlow._fieldDict['p'] = p.function
 
    # now: set correct target values for sponge
    meanFlow._fieldDict['u_target'] = u_t.function
    meanFlow._fieldDict['p_target'] = p_t.function


    #######################################
    ######## export as resolvent mode #####
    ########   (for debugging)        #####
    #######################################
    #omegas   = param.IOResolvent.Omegas
    #nSol     = 1 
    #nDofs    = A.getSizes()[0][0]
    #nOmegas  = len(omegas)
    #gains             = np.zeros((nSol,nOmegas),'complex')
    #forcings          = np.zeros((nDofs,nSol,nOmegas),'complex')
    #responses         = np.zeros((nDofs,nSol,nOmegas),'complex')
    #responses[:,0,0]  = baseFlow.getCoefficientArray()
    #forcings[:,0,0]   = baseFlow_sens_1_array[:] + 1j*baseFlow_sens_2_array[:]
    #fluctSolutList    = []
    #for i, omega in enumerate(omegas):
    #    for gainNumb in range(responses.shape[1]):
    #        fluctSolutForcing = fluctuationSolutions(
    #                            param,
    #                            meanFlow,
    #                            FEMSpaces,
    #                            omegas[i],
    #                            forcings[:,gainNumb,i],
    #                            False,
    #                            gainNumb,
    #                            gains[gainNumb, i],
    #                            )
    #        fluctSolutResponse = fluctuationSolutions(
    #                            param,
    #                            meanFlow,
    #                            FEMSpaces,
    #                            omegas[i],
    #                            responses[:,gainNumb, i],
    #                            True,
    #                            gainNumb,
    #                            gains[gainNumb, i],
    #                            )

    #        fluctSolutList.append(fluctSolutForcing)
    #        fluctSolutList.append(fluctSolutResponse)


    #Export.ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)
    #exit()
    #######################################
    #######################################
    
    # get reference nonlinear expression
    N_0  = equation.getNonlinearExpression(meanFlow)    

    # get linear solver  
    A        = equation.getLinearOperator(meanFlow)
    solver_A = LinearSolver.createEquationSystemSolver(A)

    # deform mesh
    epsilon = 1.e-8
    a_i     = optimizerParameters
    
    geometryDeformer    = caseHandler.getGeometryDeformer() 
    N_param             = geometryDeformer.getNumberOfParameters()
    baseFlowSensitivity = [None]*N_param
    for i in range(0,N_param):
        printDebug(True, "-------------------------------------------------------------" )
        printDebug(True, "-- Calculating base flow sensitivity to parameter number "+str(i+1)+"...")
        printDebug(True, "-------------------------------------------------------------" )

        a_i[i] = a_i[i] + epsilon
    
        geometryDeformer.deformMesh(a_i)
        
        N_deformed  = equation.getNonlinearExpression(meanFlow)
        N_deformed.axpy( -1., N_0) # does not need to be substracted, is zero (only for base flow...)
        N_deformed.scale(1./epsilon)

        baseFlowSensitivity[i]     = Field(FEMSpaces.VMixed, mesh)
        baseFlowSensitivity_array  = LinearSolver.solveEquationSystemWithPredefinedSolver(solver_A, N_deformed)
        baseFlowSensitivity[i].setCoefficientArray(baseFlowSensitivity_array)

        a_i[i] = a_i[i] - epsilon
        geometryDeformer.restoreMesh()
   

    ## export leading modes in standard felics format
    #baseFlow_0   = Mode(FEMSpaces.VMixed, mesh)
    #baseFlow_0.setEigenValue(0.)
    #baseFlow_0.setCoefficientArray(baseFlowSensitivity[0].getCoefficientArray())
    #baseFlow_1   = Mode(FEMSpaces.VMixed, mesh)
    #baseFlow_1.setEigenValue(0.)
    #baseFlow_1.setCoefficientArray(baseFlowSensitivity[1].getCoefficientArray())
    #baseFlow_1.isAdjoint=True
    #solution_onlyLeading = ModeCollection(FEMSpaces.VMixed, mesh)
    #solution_onlyLeading.appendMode(baseFlow_0)
    #solution_onlyLeading.appendMode(baseFlow_1)
    #fluctSolutList = solution_onlyLeading.getOldSolutionObject(meanFlow, param, FEMSpaces)
    #ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)


    baseFlowSensitivity_array=[]
    for sens in baseFlowSensitivity:
        baseFlowSensitivity_array.append(sens.getCoefficientArray())

    return baseFlowSensitivity_array, baseFlow.getCoefficientArray()




