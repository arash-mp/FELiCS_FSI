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

from ufl import VectorElement, SpatialCoordinate, exp
from dolfinx.fem import Function, FunctionSpace, Expression
from dolfinx.mesh import  locate_entities
from   CaseHandler import CaseHandler

def saveToFelFile(mesh,FEMSpaces, u_original, sponge):
    #Define functions for source of interpolation in first order
    x_source = Function(FEMSpaces.P1,dtype=float)
    y_source = Function(FEMSpaces.P1,dtype=float)
    
    #Define functions for the target of the interpolation in P2
    x_target = Function(FEMSpaces.P2,dtype=float)
    y_target = Function(FEMSpaces.P2,dtype=float)
    ux_target = Function(FEMSpaces.P2)
    uy_target = Function(FEMSpaces.P2)
    
    # Get coordinates in P1 FEM spaces
    x_source.x.array[:] = mesh.geometry.x[:,0]
    y_source.x.array[:] = mesh.geometry.x[:,1]
    # Interpolate coordinates on P2 FEM space
    x_target.interpolate(x_source)
    y_target.interpolate(y_source)
  
    # Interpolate from P2 vector FEM space to P2 FEM space
    ux_target.interpolate(u_original.function.sub(0))
    uy_target.interpolate(u_original.function.sub(1))

    # write base flow file
    with h5py.File("base_flow_for_FELiCS.fel", 'w') as f:
        f.create_dataset('/MeanFlow/x', data=x_target.x.array)
        f.create_dataset('/MeanFlow/y', data=y_target.x.array)
        f.create_dataset('/MeanFlow/ux', data=ux_target.x.array)
        f.create_dataset('/MeanFlow/uy', data=uy_target.x.array)
        f.create_dataset('/MeanFlow/spg', data=sponge.x.array)


def saveBaseflow(mesh, u_original):

    # Interpolate on P1 elements
    P1_first = VectorElement('CG', mesh.ufl_cell(), 1)
    Sol      = FunctionSpace(mesh, P1_first)
    u        = Function(Sol)
    u.interpolate(u_original.function)

    # Write felics baseflow file
    with h5py.File("base_flow_4_plot.fel", 'w') as f:
        f.create_dataset('/MeanFlow/x', data=mesh.geometry.x[:,0])
        f.create_dataset('/MeanFlow/y', data=mesh.geometry.x[:,1])
        f.create_dataset('/MeanFlow/ux', data=u.sub(0).collapse().x.array)
        f.create_dataset('/MeanFlow/uy', data=u.sub(1).collapse().x.array)

    # Write solution
    xdmf = dolfinx.io.XDMFFile(MPI.COMM_WORLD, "base_flow.xdmf", "w")
    xdmf.write_mesh(mesh)
    xdmf.write_function(u)
    xdmf.close()


def calculateBaseFlow(settingsFileName, optimizerParameters = None, deformed = False):

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

    caseHandler = CaseHandler(settingsFileName, param, mesh, equation.boundaries)

    # deform mesh at beginning - only for testing /debugging
    if deformed == True:
        geometryDeformer = caseHandler.getGeometryDeformer()
        geometryDeformer.deformMesh(optimizerParameters)
        geometryDeformer.isDeformed = False


    #-----------------------------------------------------------------------
    ## CALCULATE BASE FLOW 
    #-----------------------------------------------------------------------
    baseFlow   = Field(FEMSpaces.VMixed, mesh)

    # if the given mean flow is not zero, return to main (no base flow will be computed) 
    if np.linalg.norm(meanFlow._fieldDict['u'].x.array[:]) > 1.e-8:
        return baseFlow.getCoefficientArray()


    # ---------------- INITIALIZE CASE --------------------------------------
    try:
        # try reading previously saved base flow file
        array = np.load('baseFlow.npy')
        if (len(array) == len(baseFlow.getCoefficientArray())):
            baseFlow.setCoefficientArray(array)
        else:
            a = notInitializedVariable #do something to trigger "except" statement
    except:
        # create initial solution 
        initialValues = caseHandler.getInitialValuesForBaseFlow()
        mapping_ux = baseFlow.space.sub(0).sub(0).collapse()[1]
        mapping_uy = baseFlow.space.sub(0).sub(1).collapse()[1]
        mapping_p  = baseFlow.space.sub(1).collapse()[1]
        baseFlow.function.x.array[mapping_ux] = initialValues[0]
        baseFlow.function.x.array[mapping_uy] = initialValues[1]
        baseFlow.function.x.array[mapping_p ] = initialValues[2]
    

    # set target function for nonlinear sponge
    targetValues = caseHandler.getTargetValuesForSponge()
    [u_t,p_t] = baseFlow.getListOfSingleFields()
    mapping_ux = u_t.space.sub(0).collapse()[1]
    mapping_uy = u_t.space.sub(1).collapse()[1]
    u_t.function.x.array[mapping_ux] = targetValues[0]
    u_t.function.x.array[mapping_uy] = targetValues[1]
    p_t.function.x.array[:]          = targetValues[2]
    meanFlow._fieldDict['u_target'] = u_t.function
    meanFlow._fieldDict['p_target'] = p_t.function
   

    # set velocity components at wall (id=1001 and id=1002, which are the default wall ids for BOA) to zero
    baseFlow.function.sub(0).sub(0).x.array[dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(0), 1, equation.boundaries.find(1001))] = 0. 
    baseFlow.function.sub(0).sub(0).x.array[dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(0), 1, equation.boundaries.find(1002))] = 0. 
    baseFlow.function.sub(0).sub(0).x.array[dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(1), 1, equation.boundaries.find(1001))] = 0. 
    baseFlow.function.sub(0).sub(0).x.array[dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(1), 1, equation.boundaries.find(1002))] = 0. 


    # ---------------- START LOOP ---------------------------------------------- 
    ## start Newton solver
    target_residuum = 1.e-11
    # track time
    start= time.time()
  
    nulam_target = Field(FEMSpaces.P2, mesh)
    nulam_target.setCoefficientArray(meanFlow._fieldDict['nulam'].x.array[:])


    for viscFactor in param.Case.MolViscRampFactors:
        # set the viscosity for this Newton-Ramp loop
        meanFlow._fieldDict['nulam'].x.array[:] = nulam_target.getCoefficientArray() * viscFactor
        printDebug(True, "-------------------------------------------------------------" )
        printDebug(True, "-- viscosity factor for base flow run:  %4g " % viscFactor)
        printDebug(True, "-------------------------------------------------------------" )

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
            L.destroy()
            N.destroy()
        
            # update baseFlow
            baseFlow_array = baseFlow.getCoefficientArray() + newtonSummand_array
            baseFlow.setCoefficientArray(baseFlow_array)
            
            # calculate nonlinear expression & residuum
            [u,p] = baseFlow.getListOfSingleFields()
            meanFlow._fieldDict['u'] = u.function
            meanFlow._fieldDict['p'] = p.function
            N = equation.getNonlinearExpression(meanFlow)
            residuum = np.linalg.norm(N.getArray())
            printDebug(True, "-------------------------------------------------------------" )
            printDebug(True, "-- Base flow iteration: "+str(i)+"; Residuum: %4g " % residuum)
            printDebug(True, "-------------------------------------------------------------" )


    
    # end tracking time
    end = time.time() - start
    printDebug(True, '-- Solving the base flow problem took %4g s' % end)
    printDebug(True, '-- Residuum:  %12g' % (residuum))
    
    # save base flow as numpy file, to accelerate future base flow calculations
    np.save("baseFlow.npy",  baseFlow.getCoefficientArray())


    # Interpolate solution on first order Lagrange Functionspace for XDMF export
    saveBaseflow(mesh.dolfinxMesh, u)
    saveToFelFile(mesh,FEMSpaces,u, meanFlow._fieldDict['spg'])
    # export base flow in modes-format 
    # baseFlow_0   = Mode(FEMSpaces.VMixed, mesh)
    # baseFlow_0.setEigenValue(0.)
    # baseFlow_0.setCoefficientArray(baseFlow.getCoefficientArray())
    # baseFlow_1   = Mode(FEMSpaces.VMixed, mesh)
    # baseFlow_1.setEigenValue(0.)
    # baseFlow_1.isAdjoint = True
    # solution = ModeCollection(FEMSpaces.VMixed, mesh)
    # solution.appendMode(baseFlow_0)
    # solution.appendMode(baseFlow_1)
    # fluctSolutList = solution.getOldSolutionObject(meanFlow, param, FEMSpaces)
    # ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)


    return baseFlow.getCoefficientArray()
