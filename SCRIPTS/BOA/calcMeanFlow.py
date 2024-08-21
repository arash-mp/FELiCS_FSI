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


def getOperatorCrankNicolson(q_new,dt,B,equation,meanFlow):
    # Calculates the operator of the Crank Nicolson method, calculated with Newton-Steps:
    # op = B + 0.5*dt*LinearOperator_new

    # linear operator with q_new
    [u_new,p_new] = q_new.getListOfSingleFields()
    meanFlow._fieldDict['u'] = u_new.function
    meanFlow._fieldDict['p'] = p_new.function
    L = equation.getLinearOperator(meanFlow)

    # L = B + L*dt/2
    L.aypx(-dt/2., B) # L is negative, due to the weak form
     
    return L 
                           




def getRightHandSideCrankNicolson(q,q_new,dt,B,equation, meanFlow):
    # Calculates the right hand side of the Crank Nicolson method, calculated with Newton-Steps:
    # rhs = (q_new-q) + 0.5*dt*(NavierStokes+NavierStokes_new)

    [u,p] = q.getListOfSingleFields()
    meanFlow._fieldDict['u'] = u.function
    meanFlow._fieldDict['p'] = p.function
    N = equation.getNonlinearExpression(meanFlow)

    [u_new,p_new] = q_new.getListOfSingleFields()
    meanFlow._fieldDict['u'] = u_new.function
    meanFlow._fieldDict['p'] = p_new.function
    N_new = equation.getNonlinearExpression(meanFlow)

    # dq = u_new - u
    q_new  = q_new.getPetscVector()
    dq     = q.getPetscVector()
    q_new.axpy(-1.,dq)

    B.mult(q_new, dq)

    # N = N + N_new
    N.axpy(1.,N_new)

    # N = dq + N*dt/2
    N.aypx(dt/2., dq) # N has the opposite sign to L 
   
    N_new.destroy()
    q_new.destroy()
    dq.destroy()

    return N
                           

def getProbeValue(function, x0, mesh):
    tree = dolfinx.geometry.bb_tree(mesh._cpp_object,mesh.geometry.dim)
    cell_candidates = dolfinx.geometry.compute_collisions_points(tree, x0)
    cell = dolfinx.geometry.compute_colliding_cells(mesh._cpp_object, cell_candidates, x0)
    #cells = dolfinx.geometry.compute_first_entity_collision(tree, mesh, x0)
    return function.eval(x0,cell)


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
    timeDependentSolution  = Field(FEMSpaces.VMixed, mesh)
    averagedSolution  = Field(FEMSpaces.VMixed, mesh)

    # if the given mean flow is not zero, return to main (no base flow will be computed) 
    if np.linalg.norm(meanFlow._fieldDict['u'].x.array[:]) > 1.e-8:
        return timeDependentSolution.getCoefficientArray()


    # ---------------- INITIALIZE CASE --------------------------------------

    # ---------------- INITIALIZE CASE --------------------------------------
    try:
        # try reading previously saved base flow file
        array = np.load('baseFlow.npy')
        if (len(array) == len(timeDependentSolution.getCoefficientArray())):
            timeDependentSolution.setCoefficientArray(array)
        else:
            a = notInitializedVariable
    except: 
        # create initial solution
        initialValues = caseHandler.getInitialValuesForBaseFlow()
        mapping_ux = timeDependentSolution.space.sub(0).sub(0).collapse()[1]
        mapping_uy = timeDependentSolution.space.sub(0).sub(1).collapse()[1]
        mapping_p  = timeDependentSolution.space.sub(1).collapse()[1]
        timeDependentSolution.function.x.array[mapping_ux] = initialValues[0]
        timeDependentSolution.function.x.array[mapping_uy] = initialValues[1]
        timeDependentSolution.function.x.array[mapping_p ] = initialValues[2]
    

    # set target function for nonlinear sponge
    targetValues = caseHandler.getTargetValuesForSponge()
    [u_t,p_t] = timeDependentSolution.getListOfSingleFields()
    mapping_ux = u_t.space.sub(0).collapse()[1]
    mapping_uy = u_t.space.sub(1).collapse()[1]
    u_t.function.x.array[mapping_ux] = targetValues[0]
    u_t.function.x.array[mapping_uy] = targetValues[1]
    p_t.function.x.array[:]          = targetValues[2]
    meanFlow._fieldDict['u_target'] = u_t.function
    meanFlow._fieldDict['p_target'] = p_t.function
   

    # set velocity components at wall (id=1001 and id=1002, which are the default wall ids for BOA) to zero
    timeDependentSolution.function.sub(0).sub(0).x.array[dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(0), 1, equation.boundaries.find(1001))] = 0. 
    timeDependentSolution.function.sub(0).sub(0).x.array[dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(0), 1, equation.boundaries.find(1002))] = 0. 
    timeDependentSolution.function.sub(0).sub(0).x.array[dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(1), 1, equation.boundaries.find(1001))] = 0. 
    timeDependentSolution.function.sub(0).sub(0).x.array[dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(1), 1, equation.boundaries.find(1002))] = 0. 



    # ---------------- START TIMESTEPPING ---------------------------------------------- 
    probeLocation = [0.2,0.,0.]

    ## solver parameiters
    numberOfNewtonIterations = 5
    timeStep     = 1.e-4
    timeStep_max = 1.e-4
    timeEnd  = 1.
    residuum_max = 1.e-10
    residuum_min = 1.e-14
  
    i =0
    t = 0.
    dq    = Field(FEMSpaces.VMixed, mesh)
    q_new = Field(FEMSpaces.VMixed, mesh)
    q_old = Field(FEMSpaces.VMixed, mesh)
    q_new.setCoefficientArray(timeDependentSolution.getCoefficientArray()) #q_new = timeDependentSolution

    # create time derivative matrix. Has to be independent of mean flow! (only for incompressible)
    B        = equation.getWeightMatrix(meanFlow)
    W        = equation.getFEMWeightMatrix()
    D        = equation.getFEMDiffusionMatrix(1.e-9)
    solver_D = LinearSolver.createEquationSystemSolver(D)
    q_petsc1 = timeDependentSolution.getPetscVector() 
    q_petsc2 = timeDependentSolution.getPetscVector()
    dofs     = range(0,q_petsc1.getSizes()[0]) 

    # track computing time
    start= time.time()
    while t < timeEnd:

        dq.setConstantValue(0.)  

        for j in range(numberOfNewtonIterations):
            q_new.setCoefficientArray(timeDependentSolution.getCoefficientArray() + dq.getCoefficientArray()) #q_new = q_new - dq_Newton

            A  = getOperatorCrankNicolson(q_new, timeStep, B, equation, meanFlow)
            b  = getRightHandSideCrankNicolson(timeDependentSolution, q_new, timeStep, B, equation, meanFlow)  
           
            dq_Newton = LinearSolver.solveEquationSystem(A,b,destroy=True)

            dq.setCoefficientArray(dq.getCoefficientArray() - dq_Newton[:]) #q_new = q_new - dq_Newton

            # calc residuum = q^T W q
            q_petsc1.setValues(dofs,dq_Newton)
            W.mult(q_petsc1, q_petsc2)
            residuum = np.real(q_petsc1.dot(q_petsc2))
            print(residuum)

            if residuum < residuum_max:
                break

        # smooth solution 
        rhs      = equation.getFullRHS(dq.function)
        dq.setCoefficientArray(LinearSolver.solveEquationSystemWithPredefinedSolver(solver_D, rhs, destroy=True))

        q_old.setCoefficientArray(timeDependentSolution.getCoefficientArray())
        timeDependentSolution.setCoefficientArray(q_old.getCoefficientArray() + dq.getCoefficientArray())

        # adaptive time step to make sure the Newton method converges (the solution has to be sufficiently near the solution of the next time step):
        if residuum > residuum_max:
            #reduce time step and begin anew
            timeStep = timeStep/2.
            timeDependentSolution.setCoefficientArray(q_old.getCoefficientArray())
            q_new.setCoefficientArray(q_old.getCoefficientArray())
            continue
        elif j==0 and residuum < residuum_min and timeStep < timeStep_max / 1.2:
            timeStep = timeStep*1.2

        # write probe results in file
        [u,p] = timeDependentSolution.getListOfSingleFields()
        probeValue_u = getProbeValue(u.function, probeLocation, mesh)
        probeValue_p = getProbeValue(p.function, probeLocation, mesh)
        outputString = str(t) +"  "+ str(np.real(probeValue_u[0]))+"  " + str(np.real(probeValue_u[1])) +"  "+ str(np.real(probeValue_p[0]))
        print(outputString)
        outputFile = open('probe.txt', 'a')
        outputFile.write(outputString + "\n")
        outputFile.close()

        # write solution in file
        if np.mod(i,10)==0:
            [u,p] = timeDependentSolution.getListOfSingleFields()
            #[u,p] = averagedSolution.getListOfSingleFields()
            # Interpolate solution on first order Lagrange Functionspace for XDMF export
            saveMeanflow(mesh.dolfinxMesh, u)
            saveToFelFile(mesh,FEMSpaces,u, meanFlow._fieldDict['spg'])
            # export base flow in modes-format 
            baseFlow_0   = Mode(FEMSpaces.VMixed, mesh)
            baseFlow_0.setEigenValue(0.)
            baseFlow_0.setCoefficientArray(timeDependentSolution.getCoefficientArray())
            baseFlow_1   = Mode(FEMSpaces.VMixed, mesh)
            baseFlow_1.setEigenValue(0.)
            baseFlow_1.isAdjoint = True
            solution = ModeCollection(FEMSpaces.VMixed, mesh)
            solution.appendMode(baseFlow_0)
            solution.appendMode(baseFlow_1)
            fluctSolutList = solution.getOldSolutionObject(meanFlow, param, FEMSpaces)
            ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)


        t = t + timeStep
        i = i + 1

        printDebug(True, "-------------------------------------------------------------" )
        printDebug(True, "-- Iteration: "+str(i)+"; Time: %4g " % t)
        printDebug(True, "-- Newton It: "+str(j+1)+"; Residuum: %4g " % residuum)
        printDebug(True, "-- Time Step: %4g " % timeStep)
        printDebug(True, "-------------------------------------------------------------" )



    # end tracking time
    end = time.time() - start
    printDebug(True, '-- calculating the mean flow took %4g s' % end)
    
    # save base flow as numpy file, to accelerate future base flow calculations
    np.save("timeDependentSolution.npy",  timeDependentSolution.getCoefficientArray())


    [u,p] = timeDependentSolution.getListOfSingleFields()
    #[u,p] = averagedSolution.getListOfSingleFields()
    # Interpolate solution on first order Lagrange Functionspace for XDMF export
    saveMeanflow(mesh.dolfinxMesh, u)
    saveToFelFile(mesh,FEMSpaces,u, meanFlow._fieldDict['spg'])
    # export base flow in modes-format 
    baseFlow_0   = Mode(FEMSpaces.VMixed, mesh)
    baseFlow_0.setEigenValue(0.)
    baseFlow_0.setCoefficientArray(timeDependentSolution.getCoefficientArray())
    baseFlow_1   = Mode(FEMSpaces.VMixed, mesh)
    baseFlow_1.setEigenValue(0.)
    baseFlow_1.isAdjoint = True
    solution = ModeCollection(FEMSpaces.VMixed, mesh)
    solution.appendMode(baseFlow_0)
    solution.appendMode(baseFlow_1)
    fluctSolutList = solution.getOldSolutionObject(meanFlow, param, FEMSpaces)
    ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)

    return averagedSolution.getCoefficientArray()
