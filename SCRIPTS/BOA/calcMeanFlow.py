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

from   ufl import VectorElement, SpatialCoordinate, exp
from   dolfinx.fem import Function, FunctionSpace, Expression
from   dolfinx.mesh import  locate_entities

from   CaseHandler import CaseHandler

from   TimeStepping.CrankNicolson_inc import CrankNicolson_inc

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
    with h5py.File("mean_flow_for_FELiCS.fel", 'w') as f:
        f.create_dataset('/MeanFlow/x', data=x_target.x.array)
        f.create_dataset('/MeanFlow/y', data=y_target.x.array)
        f.create_dataset('/MeanFlow/ux', data=ux_target.x.array)
        f.create_dataset('/MeanFlow/uy', data=uy_target.x.array)
        f.create_dataset('/MeanFlow/spg', data=sponge.x.array)


def saveMeanflow(mesh, u_original):

    # Interpolate on P1 elements
    P1_first = VectorElement('CG', mesh.ufl_cell(), 1)
    Sol      = FunctionSpace(mesh, P1_first)
    u        = Function(Sol)
    u.interpolate(u_original.function)

    # Write felics baseflow file
    with h5py.File("mean_flow_4_plot.fel", 'w') as f:
        f.create_dataset('/MeanFlow/x', data=mesh.geometry.x[:,0])
        f.create_dataset('/MeanFlow/y', data=mesh.geometry.x[:,1])
        f.create_dataset('/MeanFlow/ux', data=u.sub(0).collapse().x.array)
        f.create_dataset('/MeanFlow/uy', data=u.sub(1).collapse().x.array)

    # Write solution
    xdmf = dolfinx.io.XDMFFile(MPI.COMM_WORLD, "mean_flow.xdmf", "w")
    xdmf.write_mesh(mesh)
    xdmf.write_function(u)
    xdmf.close()

def getProbeValue(f, x_0):
    tree = dolfinx.geometry.bb_tree(f.mesh._cpp_object,f.mesh.geometry.dim)
    cell_candidates = dolfinx.geometry.compute_collisions_points(tree, x_0)
    cell = dolfinx.geometry.compute_colliding_cells(f.mesh._cpp_object, cell_candidates, x_0)
    return f.function.eval(x_0,cell)


def monitor(iteration, timeDependentSolution, time):

            # write probe results in file
            probeLocation = [2.,0.,0.]
            [u,p] = timeDependentSolution.getListOfSingleFields()
            probeValue_u = getProbeValue(u, probeLocation)
            probeValue_p = getProbeValue(p, probeLocation)
            outputString = str(time) +"  "+ str(np.real(probeValue_u[0]))+"  " + str(np.real(probeValue_u[1])) +"  "+ str(np.real(probeValue_p[0]))
            print(outputString)
            outputFile = open('probe.txt', 'a')
            outputFile.write(outputString + "\n")
            outputFile.close()
    
            # write solution in file
            if np.mod(iteration,10)==1:
                # save solution as numpy file, to accelerate future calculations
                np.save("timeDependentSolution.npy",  timeDependentSolution.getCoefficientArray())

            printDebug(True, "-------------------------------------------------------------" )
            printDebug(True, "-- Iteration: "+str(iteration)+"; Time: %4g " % time)
            printDebug(True, "-------------------------------------------------------------" )


def calculateMeanFlow(settingsFileName, optimizerParameters = None, deformed = False):
    # uses implicit Crank-Nicolson-method of order 2 for time stepping:
    # (q_new-q)/dt = 0.5*(NavierStokes(q_new) + NavierStokes(q)).
    # q_new is calculated with the Newton method


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
    ## CALCULATE MEAN FLOW 
    #-----------------------------------------------------------------------
    averagedSolution       = Field(FEMSpaces.VMixed, mesh)
    q_init                 = Field(FEMSpaces.VMixed, mesh)

    # ---------------- INITIALIZE CASE --------------------------------------
    # read previously saved solution file
    array = np.load('timeDependentSolution.npy')
    q_init.setCoefficientArray(array)

    ## if the given mean flow is not zero, set it as initial value 
    #if np.linalg.norm(meanFlow._fieldDict['u'].x.array[:]) > 1.e-8:
    #    ux = Function(FEMSpaces.P2)
    #    uy = Function(FEMSpaces.P2)
    #    p  = Function(FEMSpaces.P2)
    #    mapping_ux = meanFlow._fieldDict['u'].function_space.sub(0).collapse()[1]
    #    mapping_uy = meanFlow._fieldDict['u'].function_space.sub(1).collapse()[1]
    #    ux.x.array[:] = meanFlow._fieldDict['u'].x.array[mapping_ux]
    #    uy.x.array[:] = meanFlow._fieldDict['u'].x.array[mapping_uy]
    #    p.interpolate(meanFlow._fieldDict['p'])
    #    mapping_ux = q_init.space.sub(0).sub(0).collapse()[1]
    #    mapping_uy = q_init.space.sub(0).sub(1).collapse()[1]
    #    mapping_p  = q_init.space.sub(1).collapse()[1]
    #    q_init.function.x.array[mapping_ux] = ux.x.array[:] 
    #    q_init.function.x.array[mapping_uy] = uy.x.array[:] 
    #    q_init.function.x.array[mapping_p]  = p.x.array[:] 

    #    
    #    timeDependentSolution.setCoefficientArray(q_init.getCoefficientArray())
    #    #nDofs = len(timeDependentSolution.function.x.array[:])
    #    #timeDependentSolution.setCoefficientArray(timeDependentSolution.getCoefficientArray() + (np.random.rand(nDofs)-0.5)*1.e-2)

    #else:
    #    # create initial solution
    #    initialValues = caseHandler.getInitialValuesForBaseFlow()
    #    mapping_ux = timeDependentSolution.space.sub(0).sub(0).collapse()[1]
    #    mapping_uy = timeDependentSolution.space.sub(0).sub(1).collapse()[1]
    #    mapping_p  = timeDependentSolution.space.sub(1).collapse()[1]
    #    timeDependentSolution.function.x.array[mapping_ux] = initialValues[0]
    #    timeDependentSolution.function.x.array[mapping_uy] = initialValues[1]
    #    timeDependentSolution.function.x.array[mapping_p ] = initialValues[2]
    

    # set target function for nonlinear sponge
    targetValues = caseHandler.getTargetValuesForSponge()
    [u_t,p_t] = q_init.getListOfSingleFields()
    mapping_ux = u_t.space.sub(0).collapse()[1]
    mapping_uy = u_t.space.sub(1).collapse()[1]
    u_t.function.x.array[mapping_ux] = targetValues[0]
    u_t.function.x.array[mapping_uy] = targetValues[1]
    p_t.function.x.array[:]          = targetValues[2]
    meanFlow._fieldDict['u_target'] = u_t.function
    meanFlow._fieldDict['p_target'] = p_t.function
   


    # ---------------- START TIMESTEPPING ---------------------------------------------- 
    timeStepper = CrankNicolson_inc(FEMSpaces.VMixed, mesh, meanFlow, equation)

    # track computing time
    start= time.time()

    timeStepper.setMonitoringMethod(monitor)
    solution = timeStepper.doTimeStepping(q_init = q_init, t_start = 1000., dt = 0.01, t_end = 1001.)

    # end tracking time
    end = time.time() - start
    printDebug(True, '-- calculating the mean flow took %4g s' % end)
    
    # save base flow as numpy file, to accelerate future base flow calculations
    np.save("timeDependentSolution.npy",  solution.getCoefficientArray())


    [u,p] = solution.getListOfSingleFields()
    #[u,p] = averagedSolution.getListOfSingleFields()
    # Interpolate solution on first order Lagrange Functionspace for XDMF export
    saveMeanflow(mesh.dolfinxMesh, u)
    saveToFelFile(mesh,FEMSpaces,u, meanFlow._fieldDict['spg'])
    # export base flow in modes-format 
    baseFlow_0   = Mode(FEMSpaces.VMixed, mesh)
    baseFlow_0.setEigenValue(0.)
    baseFlow_0.setCoefficientArray(solution.getCoefficientArray())
    baseFlow_1   = Mode(FEMSpaces.VMixed, mesh)
    baseFlow_1.setEigenValue(0.)
    baseFlow_1.isAdjoint = True
    solution = ModeCollection(FEMSpaces.VMixed, mesh)
    solution.appendMode(baseFlow_0)
    solution.appendMode(baseFlow_1)
    fluctSolutList = solution.getOldSolutionObject(meanFlow, param, FEMSpaces)
    ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)

    exit()

    return averagedSolution.getCoefficientArray()
