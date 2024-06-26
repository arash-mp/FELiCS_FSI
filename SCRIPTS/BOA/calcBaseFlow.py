import pdb
import numpy as np
import copy
import time
import h5py

import dolfinx

from   mpi4py import MPI

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
    ## CALCULATE BASE FLOW 
    #-----------------------------------------------------------------------
    baseFlow   = Field(FEMSpaces.VMixed, mesh)

    # if the given mean flow is not zero, return to main (no base flow will be computed) 
    if np.linalg.norm(meanFlow._fieldDict['u'].x.array[:]) > 1.e-8:
        return baseFlow

    viscosityRampFactors = param.Case.MolViscRampFactors
    
    if viscosityRampFactors[0]==1.:
        #--------------------------------------------------------------------------
        # CASE CYLINDER
        #--------------------------------------------------------------------------
        # initialize field 
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
        

        # set target function for nonlinear sponge
        [u_t,p_t] = baseFlow.getListOfSingleFields()
        mapping_ux = u_t.space.sub(0).collapse()[1]
        mapping_uy = u_t.space.sub(1).collapse()[1]
        u_t.function.x.array[mapping_ux] = 1.
        u_t.function.x.array[mapping_uy] = 0.
        p_t.function.x.array[:]          = 0.
        meanFlow._fieldDict['u_target'] = u_t.function
        meanFlow._fieldDict['p_target'] = p_t.function
        
        # set boundary conditions for the velocity components at cylinder wall, identifier "1001"
        baseFlow.function.sub(0).sub(0).x.array[dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(0), 1, equation.boundaries.find(1001))] = 0. 
        baseFlow.function.sub(0).sub(1).x.array[dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(1), 1, equation.boundaries.find(1001))] = 0. 
        #--------------------------------------------------------------------------
        #--------------------------------------------------------------------------
   
    else:

        #--------------------------------------------------------------------------
        # CASE PROFILE 
        #--------------------------------------------------------------------------
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
            baseFlow.function.x.array[mapping] = 1.5
        

        # set target function for nonlinear sponge
        [u_t,p_t] = baseFlow.getListOfSingleFields()
        mapping_ux = u_t.space.sub(0).collapse()[1]
        mapping_uy = u_t.space.sub(1).collapse()[1]
        u_t.function.x.array[mapping_ux] = 1.5
        u_t.function.x.array[mapping_uy] = 0.
        p_t.function.x.array[:]          = 0.
        meanFlow._fieldDict['u_target'] = u_t.function
        meanFlow._fieldDict['p_target'] = p_t.function
        

        baseFlow.function.sub(0).sub(0).x.array[dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(0), 1, equation.boundaries.find(1))] = 0. 
        baseFlow.function.sub(0).sub(0).x.array[dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(0), 1, equation.boundaries.find(2))] = 0. 
        baseFlow.function.sub(0).sub(0).x.array[dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(0), 1, equation.boundaries.find(3))] = 0. 
        baseFlow.function.sub(0).sub(0).x.array[dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(0), 1, equation.boundaries.find(4))] = 0. 
        baseFlow.function.sub(0).sub(0).x.array[dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(0), 1, equation.boundaries.find(5))] = 0. 
        #--------------------------------------------------------------------------
        #--------------------------------------------------------------------------


    ## start Newton solver
    target_residuum = 1.e-11
    # track time
    start= time.time()
   
    nulam_target = meanFlow._fieldDict['nulam'].x.array[:]
    for viscFactor in viscosityRampFactors:
        # set the viscosity for this ramp loop
        meanFlow._fieldDict['nulam'].x.array[:] = nulam_target * viscFactor
        printDebug(True, "-- viscosity factor for base flow run:  %4g " % viscFactor)

        # calculate the base flow
        i=0
        residuum=1.
        [u,p] = baseFlow.getListOfSingleFields()
        meanFlow._fieldDict['u'] = u.function
        meanFlow._fieldDict['p'] = p.function
        N = equation.getNonlinearExpression(meanFlow)
        while(residuum > target_residuum and i<10):
            i+=1
        
            # solve equation system
            L = equation.getLinearOperator(meanFlow)
            newtonSummand_array = LinearSolver.solveEquationSystem(L,N)
        
            # update baseFlow
            baseFlow_array = baseFlow.getCoefficientArray() + newtonSummand_array
            baseFlow.setCoefficientArray(baseFlow_array)
            
            # calculate nonlinear expression & residuum
            [u,p] = baseFlow.getListOfSingleFields()
            meanFlow._fieldDict['u'] = u.function
            meanFlow._fieldDict['p'] = p.function
            N = equation.getNonlinearExpression(meanFlow)
            residuum = np.linalg.norm(N.getArray())
            printDebug(True, "-- Base flow iteration: "+str(i)+"; Residuum: %4g " % residuum)

    
    # end tracking time
    end = time.time() - start
    printDebug(True, '-- Solving the base flow problem took %4g s' % end)
    printDebug(True, '-- Residuum:  %12g' % (residuum))
    
    # save base flow as numpy file, to accelerate future base flow calculations
    np.save("baseFlow.npy",  baseFlow.getCoefficientArray())


    # export base flow in modes-format 
    baseFlow_0   = Mode(FEMSpaces.VMixed, mesh)
    baseFlow_0.setEigenValue(0.)
    baseFlow_0.setCoefficientArray(baseFlow.getCoefficientArray())
    baseFlow_1   = Mode(FEMSpaces.VMixed, mesh)
    baseFlow_1.setEigenValue(0.)
    baseFlow_1.isAdjoint = True
    baseFlow_1.setCoefficientArray(baseFlow.getCoefficientArray())
    solution = ModeCollection(FEMSpaces.VMixed, mesh)
    solution.appendMode(baseFlow_0)
    solution.appendMode(baseFlow_1)
    fluctSolutList = solution.getOldSolutionObject(meanFlow, param, FEMSpaces)
    ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)


    return baseFlow


