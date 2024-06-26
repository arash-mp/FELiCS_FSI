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

from   GeometryDeformer import CylinderBSpline2Pts 


def calculateBaseFlowSensitivity(settingsFileName, baseFlow_array):

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


    ## export mean flow in "h5" file
    #if not param.FlowInput.MeanFlowFilePath.split('.')[-1] == 'hdf5':
    #    meanFlow.exportBaseFlowAsHDF5()
    #meanflowFilename = 'meanflow.h5'
    #meanFlow.mapToExportMeshAndExport(FEMSpaces, meanflowFilename)


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
    
    #-----------------------------------------------------------------------
    ## CALCULATE delta N / delta a_i and dq/da_i: 
    #-----------------------------------------------------------------------
    # set target function for nonlinear sponge
    [u_t,p_t] = baseFlow.getListOfSingleFields()
    mapping_ux = u_t.space.sub(0).collapse()[1]
    mapping_uy = u_t.space.sub(1).collapse()[1]
    u_t.function.x.array[mapping_ux] = 1.
    u_t.function.x.array[mapping_uy] = 0.
    p_t.function.x.array[:]          = 0.
    meanFlow._fieldDict['u_target'] = u_t.function
    meanFlow._fieldDict['p_target'] = p_t.function
   

    # get reference nonlinear expression
    A    = equation.getLinearOperator(meanFlow)
    N_0  = equation.getNonlinearExpression(meanFlow)    

    # create and read fields
    a_i  = np.load("params.npy")
    
    # deform mesh
    epsilon = 1.e-8
    
    geometryDeformer = CylinderBSpline2Pts(mesh, equation.boundaries)
    N_param          = geometryDeformer.getNumberOfParameters()
    baseFlowSensitivity = [None]*N_param
    for i in range(0,N_param):
        a_i[i] = a_i[i] + epsilon
    
        geometryDeformer.deformMesh(a_i)
        
        N_deformed  = equation.getNonlinearExpression(meanFlow)
        #N_deformed.axpy( -1., N_0) # does not need to be substracted, is zero
        N_deformed.scale(1./epsilon)

        baseFlowSensitivity[i] = Field(FEMSpaces.VMixed, mesh)
        baseFlowSensitivity[i].setCoefficientArray(LinearSolver.solveEquationSystem(A, N_deformed))
    
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


    return baseFlowSensitivity




