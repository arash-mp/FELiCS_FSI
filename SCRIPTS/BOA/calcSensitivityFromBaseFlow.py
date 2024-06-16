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


def calculateSensitivityFromBaseFlow(settingsFileName, baseFlow_array, rhs):

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
   
    # load base flow into meanFlow object
    baseFlow = Field(FEMSpaces.VMixed, mesh)
    baseFlow.setCoefficientArray(baseFlow_array)
    [u,p] = baseFlow.getListOfSingleFields()
    meanFlow._fieldDict['u'] = u.function
    meanFlow._fieldDict['p'] = p.function

    # equation
    equation = EquationCollectionClass(
                                      param,
                                      FEMSpaces,
                                      meanFlow,
                                      mesh
                                      )
    
    
    #-----------------------------------------------------------------------
    ## CALCULATE ADJOINT BASEFLOW
    #-----------------------------------------------------------------------
    A  = equation.getLinearOperator(meanFlow)
    
    # solve the adjoint equation system
    rhs.conjugate()
    baseFlow_adjoint = Field(FEMSpaces.VMixed, mesh)
    baseFlow_adjoint.setCoefficientArray(LinearSolver.solveTransposeEquationSystem(A, rhs))
    baseFlow_adjoint.conjugate()
    baseFlow_petsc = baseFlow_adjoint.getPetscVector()
    dolfinx.fem.petsc.set_bc(baseFlow_petsc, equation.BCs)
    baseFlow_adjoint.setCoefficientArray(baseFlow_petsc.getArray())
    

    #-----------------------------------------------------------------------
    ## CALCULATE GRADIENT WITH RESPECT TO GEOMETRY DEFORMATIONS - Part II
    #-----------------------------------------------------------------------
    # set target function for nonlinear sponge
    [u_t,p_t] = baseFlow_adjoint.getListOfSingleFields()
    mapping_ux = u_t.space.sub(0).collapse()[1]
    mapping_uy = u_t.space.sub(1).collapse()[1]
    u_t.function.x.array[mapping_ux] = 1.
    u_t.function.x.array[mapping_uy] = 0.
    p_t.function.x.array[:]          = 0.
    meanFlow._fieldDict['u_target'] = u_t.function
    meanFlow._fieldDict['p_target'] = p_t.function
   

    # get reference nonlinear expression
    N_0  = equation.getNonlinearExpression(meanFlow)
    
    # create and read fields
    a_i  = np.load("params.npy")
    
    # deform mesh
    epsilon = 1.e-6
    
    geometryDeformer = CylinderBSpline2Pts(mesh, equation.boundaries)
    N_param          = geometryDeformer.getNumberOfParameters()
    sensitivity2     = np.zeros(N_param,dtype=complex)
    for i in range(0,N_param):
        a_i[i] = a_i[i] + epsilon
    
        geometryDeformer.deformMesh(a_i)
        
        N_deformed  = equation.getNonlinearExpression(meanFlow)
        N_deformed.axpy( -1., N_0)
        N_deformed.scale(1./epsilon)
        baseFlow_adjoint_petsc  = baseFlow_adjoint.getPetscVector()
        sensitivity2[i] = N_deformed.dot(baseFlow_adjoint_petsc) 
    
        a_i[i] = a_i[i] - epsilon
        geometryDeformer.restoreMesh()
    
    
    printDebug(True, '------------------------------------------------ ')
    printDebug(True, '-- sensitivities part 2:  ' + str((sensitivity2)))
    printDebug(True, '------------------------------------------------ ')
    
    
    return sensitivity2 




