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

def calculateNavierStokesResiduum(settingsFileName, baseFlow_array, optimizerParameters, deformed = False):

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
    ## SMOOTH MEAN FLOW 
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


    ## get smoothing matrix and rhs, and smooth mean flow
    #D        = equation.getFEMDiffusionMatrix(1.e-9)
    #solver_D = LinearSolver.createEquationSystemSolver(D)
    #rhs      = equation.getFullRHS(baseFlow.function)
 
    #baseFlow.setCoefficientArray(LinearSolver.solveEquationSystemWithPredefinedSolver(solver_D, rhs))

    [u,p] = baseFlow.getListOfSingleFields()
    meanFlow._fieldDict['u'] = u.function
    meanFlow._fieldDict['p'] = p.function
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
    meanFlow._fieldDict['u_target'] = meanFlow._fieldDict['u']#u_t.function
    meanFlow._fieldDict['p_target'] = meanFlow._fieldDict['p']#p_t.function
    #meanFlow._fieldDict['u_target'] = u_t.function
    #meanFlow._fieldDict['p_target'] = p_t.function
   # get reference nonlinear expression
    N_0  = equation.getNonlinearExpression(meanFlow)    


    # get Forcing vector 
    D        = equation.getFEMDiffusionMatrix(1.e-6, sponge=meanFlow._fieldDict['spg'])
    solver_D = LinearSolver.createEquationSystemSolver(D)

    residuum_vec  = LinearSolver.solveEquationSystemWithPredefinedSolver(solver_D, N_0)

    #####################################
    ###### export as resolvent mode #####
    ######   (for debugging)        #####
    #####################################
    #omegas   = param.IOResolvent.Omegas
    #nSol     = 1 
    #nDofs    = W.getSizes()[0][0]
    #nOmegas  = len(omegas)
    #gains             = np.zeros((nSol,nOmegas),'complex')
    #forcings          = np.zeros((nDofs,nSol,nOmegas),'complex')
    #responses         = np.zeros((nDofs,nSol,nOmegas),'complex')
    #responses[:,0,0]  = residuum_vec[:]
    #forcings[:,0,0]   = residuum_vec[:]
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

    #####################################
    #####################################

    W        = equation.getFEMWeightMatrix()

    vec1, vec2 = W.getVecs()
    vec1.setValues(range(0,len(residuum_vec[:])), residuum_vec[:])
    W.mult(vec1, vec2)

    residuum_norm = vec2.dot(vec1) 
    
   
    return residuum_norm 




