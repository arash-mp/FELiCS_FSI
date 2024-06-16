import pdb
import numpy as np
import copy
import time

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



def calculateBaseFlow(settingsFileName):

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

    # equation
    equation = EquationCollectionClass(
                                      param,
                                      FEMSpaces,
                                      meanFlow,
                                      mesh
                                      )
    
    #-----------------------------------------------------------------------
    ## CALCULATE BASE FLOW 
    #-----------------------------------------------------------------------
    
    # initialize field 
    baseFlow   = Field(FEMSpaces.VMixed, mesh)
    try:
        # try reading previously saved base flow file
        array = np.load('baseFlow.npy')
        if (len(array) == len(baseFlow.getCoefficientArray())):
            baseFlow.setCoefficientArray(array)
        else:
            a = notInitializedVariable #do something to trigger "except" statement
    except:
        # create initial solution by setting velocity component ux=1
        mapping = baseFlow.space.sub(0).sub(0).collapse()[1]
        baseFlow.function.x.array[mapping] = 1.
    
    
    # set boundary conditions for the velocity components at cylinder wall, identifier "1001"
    baseFlow.function.sub(0).sub(0).x.array[dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(0), 1, equation.boundaries.find(1001))] = 0. 
    baseFlow.function.sub(0).sub(1).x.array[dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(1), 1, equation.boundaries.find(1001))] = 0. 
    
    
    # set target function for nonlinear sponge
    [u_t,p_t] = baseFlow.getListOfSingleFields()
    mapping_ux = u_t.space.sub(0).collapse()[1]
    mapping_uy = u_t.space.sub(1).collapse()[1]
    u_t.function.x.array[mapping_ux] = 1.
    u_t.function.x.array[mapping_uy] = 0.
    p_t.function.x.array[:]          = 0.
    meanFlow._fieldDict['u_target'] = u_t.function
    meanFlow._fieldDict['p_target'] = p_t.function
    
    
    ## start Newton solver
    
    # track time
    start= time.time()
    
    i=0
    residuum=1.
    while(residuum > 3.e-11 and i<20):
        i+=1
    
        [u,p] = baseFlow.getListOfSingleFields()
        meanFlow._fieldDict['u'] = u.function
        meanFlow._fieldDict['p'] = p.function
    
    
        # solve equation system
        N = equation.getNonlinearExpression(meanFlow)
        L = equation.getLinearOperator(meanFlow)
        newtonSummand_array = LinearSolver.solveEquationSystem(L,N)
    
        # update baseFlow
        baseFlow_array = baseFlow.getCoefficientArray() + newtonSummand_array
        baseFlow.setCoefficientArray(baseFlow_array)
        
        residuum = np.linalg.norm(newtonSummand_array)
        printDebug(True, "-- Base flow iteration: "+str(i)+"; Residuum: %4g " % residuum)
    
    
    
    # end tracking time
    end = time.time() - start
    printDebug(True, '-- Solving the base flow problem took %4g s' % end)
    printDebug(True, '-- Residuum:  %12g' % (residuum))
    
    # save base flow as numpy file, to accelerate future base flow calculations
    np.save("baseFlow.npy",  baseFlow.getCoefficientArray())


    return baseFlow


