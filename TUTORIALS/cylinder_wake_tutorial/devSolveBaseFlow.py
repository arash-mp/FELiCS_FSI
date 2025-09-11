import pdb
import numpy as np
import copy
import time
import h5py

import dolfinx

from   mpi4py import MPI

from   FELiCS.Parameters.config import config

import FELiCS.IO.Import as Import

import FELiCS.SpaceDisc.FEMSpaces as DefineFEMSpaces
from   FELiCS.Fields.meanFlowClass import meanFlowClass
from   FELiCS.Fields.fluctuationClass import fluctuationSolutions
from   FELiCS.Equation.EquationCollection import EquationCollectionClass
from   FELiCS.Misc.functions import printDebug

from   FELiCS.Solvers.LinearSolver import LinearSolver 
from   FELiCS.Fields.ModeCollection import ModeCollection
from   FELiCS.Fields.Field import Field
from   FELiCS.Fields.Mode import Mode
from   FELiCS.Misc.tensorUtils import Tensor

# from ufl import VectorElement, SpatialCoordinate, exp, FiniteElement
from dolfinx.fem import Function, FunctionSpace, Expression
from dolfinx.mesh import  locate_entities
from scipy.interpolate import griddata
# from   CaseHandler import CaseHandler

def interpolateFieldsFromMeanfield(field, meanField, fieldsList = ['u', 'p']):
    # Get dimension of case from meanField
    dim = meanField._param.Case.nDim
    for i, fieldName in enumerate(fieldsList):
        # Extract field from meanField
        fieldMean = meanField._fieldDict[fieldName]
        if fieldName == 'u':
            for iu in range(dim):
                # Get subspace and mappings for both meanField and Field (field)
                space, mapping = field.space.sub(i).sub(iu).collapse()
                spaceMean, mappingMean = fieldMean.space.sub(iu).collapse()
                # Get coordinates from meanField
                initx = spaceMean.tabulate_dof_coordinates()[:,0]
                inity = spaceMean.tabulate_dof_coordinates()[:,1]
                # Get coordinates from Field
                tgtx = space.tabulate_dof_coordinates()[:,0]
                tgty = space.tabulate_dof_coordinates()[:,1]
                # Get value from meanField
                initValue = fieldMean.function.x.array[mappingMean]
                # Interpolate value to target mesh
                tgtValue = griddata((initx, inity), initValue, (tgtx, tgty), method = 'nearest')
                # Give interpolated value to Field
                field.function.x.array[mapping] = tgtValue
        else:
            # Get subspace and mappings for Field (field)
            space, mapping = field.space.sub(i).collapse()
            # Get coordinates from meanField
            initx = fieldMean.space.tabulate_dof_coordinates()[:,0]
            inity = fieldMean.space.tabulate_dof_coordinates()[:,1]
            # Get coordinates from Field
            tgtx = space.tabulate_dof_coordinates()[:,0]
            tgty = space.tabulate_dof_coordinates()[:,1]
            # Get value from meanField
            initValue = fieldMean.function.x.array[:]
            # Interpolate value to target mesh
            tgtValue = griddata((initx, inity), initValue, (tgtx, tgty), method = 'nearest')
            # Give interpolated value to Field
            field.function.x.array[mapping] = tgtValue
    pass

# def saveScalarField(mesh, p_original, paraname, fieldname):
#     # Interpolate on P1 elements
#     P1_first = FiniteElement('CG', mesh.ufl_cell(), 1)
#     Sol      = FunctionSpace(mesh, P1_first)
#     p        = Function(Sol)
#     try:
#         p.interpolate(p_original.function)
#     except AttributeError:
#         p.interpolate(p_original)
#     p.name = paraname

#     # Write solution
#     xdmf = dolfinx.io.XDMFFile(MPI.COMM_WORLD, fieldname, "w")
#     xdmf.write_mesh(mesh)
#     xdmf.write_function(p)
#     xdmf.close()

# def saveToFelFile(mesh,FEMSpaces, u_original, p_original, nu_original, sponge):
#     #Define functions for source of interpolation in first order
#     x_source = Function(FEMSpaces.P1,dtype=float)
#     y_source = Function(FEMSpaces.P1,dtype=float)
    
#     #Define functions for the target of the interpolation in P2
#     x_target = Function(FEMSpaces.P2,dtype=float)
#     y_target = Function(FEMSpaces.P2,dtype=float)
#     ux_target = Function(FEMSpaces.P2)
#     uy_target = Function(FEMSpaces.P2)
#     p_target = Function(FEMSpaces.P2)
#     nu_target = Function(FEMSpaces.P2)
    
#     # Get coordinates in P1 FEM spaces
#     x_source.x.array[:] = mesh.geometry.x[:,0]
#     y_source.x.array[:] = mesh.geometry.x[:,1]
#     # Interpolate coordinates on P2 FEM space
#     x_target.interpolate(x_source)
#     y_target.interpolate(y_source)
  
#     # Interpolate from P2 vector FEM space to P2 FEM space
#     ux_target.interpolate(u_original.function.sub(0))
#     uy_target.interpolate(u_original.function.sub(1))
#     p_target.interpolate(p_original.function)
#     try:
#         nu_target.interpolate(nu_original.function)
#     except(AttributeError):
#         nu_target.interpolate(nu_original)
        
#     print(f'Maximum visicosity is {max(nu_target.x.array)}')

#     # write base flow file
#     with h5py.File("base_flow_for_FELiCS.fel", 'w') as f:
#         f.create_dataset('/MeanFlow/x', data=x_target.x.array)
#         f.create_dataset('/MeanFlow/y', data=y_target.x.array)
#         f.create_dataset('/MeanFlow/ux', data=ux_target.x.array)
#         f.create_dataset('/MeanFlow/uy', data=uy_target.x.array)
#         f.create_dataset('/MeanFlow/p', data=p_target.x.array)
#         f.create_dataset('/MeanFlow/nulam', data=nu_target.x.array)
#         # f.create_dataset('/MeanFlow/spg', data=sponge.x.array)


# def saveBaseflow(mesh, u_original,param):
#     # Read coordinate system and define dimension
#     if param.Case.CoordinateSystem == 'Cartesian':
#         dimSet = 2
#     elif param.Case.CoordinateSystem == 'Cylindrical':
#         dimSet = 3
#     # Interpolate on P1 elements
#     P1_first = VectorElement('CG', mesh.ufl_cell(), 1, dim=dimSet)
#     Sol      = FunctionSpace(mesh, P1_first)
#     u        = Function(Sol)
#     u.interpolate(u_original.function)

#     # Write felics baseflow file
#     with h5py.File("base_flow_4_plot.fel", 'w') as f:
#         f.create_dataset('/MeanFlow/x', data=mesh.geometry.x[:,0])
#         f.create_dataset('/MeanFlow/y', data=mesh.geometry.x[:,1])
#         f.create_dataset('/MeanFlow/ux', data=u.sub(0).collapse().x.array)
#         f.create_dataset('/MeanFlow/uy', data=u.sub(1).collapse().x.array)

#     # Write solution
#     xdmf = dolfinx.io.XDMFFile(MPI.COMM_WORLD, "base_flow.xdmf", "w")
#     xdmf.write_mesh(mesh)
#     xdmf.write_function(u)
#     xdmf.close()

# def save_as_a_h5File(param, meanFlow):
#         # export mean flow in "h5" file
#     if not param.FlowInput.MeanFlowFilePath.split('.')[-1] == 'hdf5':
#         meanFlow.exportBaseFlowAsHDF5()
#     meanflowFilename = 'meanflow.h5'

def calculateBaseFlow(settingsFileName, optimizerParameters = None, deformed = False):

    #-----------------------------------------------------------------------
    ## INITIALIZATION 
    #-----------------------------------------------------------------------
    # read parameters
    param=config()
    param.importFromFile(settingsFileName)
    # param.getOldParameters()

    # mesh
    mesh=param.__mesh__
    
    # FEMSpaces
    FEMSpaces = DefineFEMSpaces.FEMSpaces(
                param,
                mesh,
                )
    
    
    # initialize mean flow class
    meanFlow = meanFlowClass(param, FEMSpaces, mesh)
    
    #import mean flow data from file
    meanFlow.importDataFromFileAndExportToH5()

    # export mean flow in "h5" file
    # if not param.FlowInput.MeanFlowFilePath.split('.')[-1] == 'hdf5':
    #     meanFlow.exportBaseFlowAsHDF5()
    # meanflowFilename = 'meanflow.h5'
    # meanFlow.mapToExportMeshAndExport(FEMSpaces, meanflowFilename)
    
    # Read mean flow field from file if doing optimization
    if not optimizerParameters == None:
        print(f'Base flow solver: Optimization mode. Perform on {optimizerParameters}')
        tempField = Field(FEMSpaces.VMixed, mesh).getListOfSingleFields()[0].getListOfSingleFields()[0]
        nulamField = Field(tempField.function.function_space, mesh)
        nulamField.function.x.array[:] = np.load(optimizerParameters + '_field.npy')
        meanFlow._fieldDict[optimizerParameters].interpolate(nulamField.function)

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
    
    # Define a Field class to save calculating results
    baseFlow   = Field(FEMSpaces.VMixed, mesh)
    # Initialize fields in baseFlow with meanFlow fields
    interpolateFieldsFromMeanfield(baseFlow, meanFlow)

    # ---------------- START LOOP ---------------------------------------------- 
    ## start Newton solver
    target_residuum = 1.e-11
    # track time
    start= time.time()
    # calculate the base flow
    i=0
    residuum=1.
    [u,p] = baseFlow.getListOfSingleFields()
    meanFlow._fieldDict['u'] = u.function
    meanFlow._fieldDict['p'] = p.function
    N = equation.getNonlinearExpression(meanFlow)
    while(residuum > target_residuum and i<0):
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
        N = equation.getNonlinearExpression(baseFlow)
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
    # saveBaseflow(mesh.dolfinxMesh, u, param)
    # saveToFelFile(mesh,FEMSpaces,u, p, meanFlow._fieldDict['nulam'], meanFlow._fieldDict['spg'])
    
    # # Save scalarFields as xdmf file
    # saveScalarField(mesh.dolfinxMesh, p, 'p', 'BaseFlow_p.xdmf')
    # saveScalarField(mesh.dolfinxMesh, meanFlow._fieldDict['nulam'], 'nu', 'BaseFlow_nu.xdmf')

    return baseFlow

# if __name__ == "__main__":
configFilename       = 'Re50.json'
calculateBaseFlow(configFilename)