#  ___________________________________   _______________________________________________________
# /-----------------------------------\ /-------------------------------------------------------\
# |   (         (                (     |  This source code is part of FELiCS                     |
# |   )\ )     ) )        (      )\ )  |  (F)inite (E)lement (Li)nearized (C)ombustion (S)olver  |  
# |  (()/(  (  (()/( (    )\   (()/(   |                                                         |  
# |  /(_)) )\  /(_)))\  (((_)  /(_))   |  Licensed under the GNU GPLv3                           |
# |  (_)_)((_) (_)) ((_) )\___ (_))    |                                                         |
# |  | __|| __|| |   (_)((/ __|/ __|   |  (C) 2018-2025: The FELiCS Developers (www.felics.eu)   |
# |  | _| | _| | |__ | | | (__ \__ \   |  Visit          www.felics.eu                           |
# |  |_|  |___||____||_|  \___||___/   |  Contact        info@felics.eu                          |
# \___________________________________/ \_______________________________________________________/
#
# TODO: 
# - Create first test for Field class => Xiuyang
# - do all other functions  => Anant

import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from FELiCS.Fields.Field import Field
from FELiCS.Misc.tensorUtils import *
from FELiCS.Equation.UflDecorator import UflDecorator

import numpy as np
from mpi4py              import MPI
from petsc4py.PETSc      import ScalarType
from ufl                 import (
    dx,TestFunction,TrialFunction,conj,SpatialCoordinate, Dx,dot, inner, grad,
    TestFunctions, TrialFunctions
)
import dolfinx
from dolfinx             import mesh
from FELiCS.Fields.Field import Field
from FELiCS.Fields.Mode  import Mode
from tests.RandomCaseHandler import RandomCaseHandler
from tests.UnitTestHelper    import UnitTestHelper
from FELiCS.Misc.logging import Logger

# initialize the logger
logger = Logger(logger_name="felics_unit_test")
logger = Logger.get_logger("felics_unit_test")
#################################################################
##### define necessary functions ################################
#################################################################
class FieldTestHandler(RandomCaseHandler):
    def __init__(self, dim_vector):
        super().__init__(dim_vector)
        self.scalar_field = self.createFELiCSField("scalar")
        self.vector_field = self.createFELiCSField("vector")
        self.mixed_field  = self.createFELiCSField("mixed")
        pass
    
    

#################################################################
##### initialize tests ##########################################
#################################################################

# 0. define parameters
dim_vector           = 3 #dimension of vector function space, 2 or 3 (at the moment: is fixed to 3; that should test the 2, too?) 

# 1. create Test class
randomField = FieldTestHandler(dim_vector)
# TODO: So far the tests can affect the whole class, which is not expected,
# This need to be solved in the future
# --------------------------------------------------------------
# Tests for Field methods
# --------------------------------------------------------------

def test_calculateL2Norm():
    logger.info("Testing L2 norm calculation")
    # 1. define expressions
    for field in [randomField.scalar_field, randomField.vector_field, randomField.mixed_field]:
        logger.info(
        " - Testing field of type: ",
        field._name,
        )
        list = field.get_list_of_sub_fields()
        J_hat = field.mesh.coordinate_system.J_hat
        for i, subfield in enumerate(list):
            fieldTens = subfield.get_tensor()
            if i == 0:
                validExpr = (i_dot(
                fieldTens,
                i_conj(fieldTens),
                )).ufl_tens*J_hat*dx
            else:
                validExpr += (i_dot(
                fieldTens,
                i_conj(fieldTens),
                )).ufl_tens*J_hat*dx 
        # 2. check alignment
        computedL2Norm = field.calculate_l2_norm()
        validL2Norm = np.sqrt(dolfinx.fem.assemble_scalar(dolfinx.fem.form(validExpr)))
        assert np.abs(computedL2Norm - validL2Norm) < 1.e-14
    logger.info("... passed.")

def test_getVorticityField():
    # test for 2-D vector
    dim_vector = 2
    logger.info("Testing computation of the vorticity field in 2D")
    randomField     = FieldTestHandler(dim_vector)
    field           = randomField.vector_field
    components      = field.get_list_of_sub_fields()
    componentList   = []  
    for subfield in components:
        componentList.append(subfield.get_gradient_field())
    dvdx = componentList[1].get_list_of_sub_fields()[0]    
    dudy = componentList[0].get_list_of_sub_fields()[1] 
    vorticity_field_valid = dvdx - dudy
    vorticity_field_test = field.get_vorticity_field()
    assert np.linalg.norm(vorticity_field_valid.function.x.array[:] - vorticity_field_test.function.x.array[:]) < 1e-13
    # Test for 3D vector
    # NOTE: we do this in the last step
    dim_vector = 3
    logger.info("Testing computation of the vorticity field in 3D")
    randomField     = FieldTestHandler(dim_vector)
    field           = randomField.vector_field
    components      = field.get_list_of_sub_fields()
    componentList   = []  
    for subfield in components:
        componentList.append(subfield.get_gradient_field())
    dwdy = componentList[2].get_list_of_sub_fields()[1]
    dvdz = componentList[1].get_list_of_sub_fields()[2]

    dudz = componentList[0].get_list_of_sub_fields()[2]
    dwdx = componentList[2].get_list_of_sub_fields()[0]

    dvdx = componentList[1].get_list_of_sub_fields()[0]
    dudy = componentList[0].get_list_of_sub_fields()[1]
    
    vorticity_x = dwdy - dvdz
    vorticity_y = dudz - dwdx
    vorticity_z = dvdx - dudy
    
    vorticity_field_valid = randomField.createFELiCSField("vector")
    vorticity_field_valid.set_list_of_sub_fields([vorticity_x, vorticity_y, vorticity_z])
    vorticity_field_test = field.get_vorticity_field()
    assert np.linalg.norm(vorticity_field_valid.function.x.array[:] - vorticity_field_test.function.x.array[:]) < 1e-12

def test_evaluateUflExpression():
    logger.info("testing evaluate Ufl expression")
    from FELiCS.Solvers.LinearSolver import LinearSolver
    # NOTE: So far only scalar field is tested
    # 0 Define a simple ufl expression
    random_scalar_function = randomField.createDolfinxFunction(dim=1)
    expr = random_scalar_function * conj(randomField.test_scalar) * dx
    temp_scalar_field = randomField.createFELiCSField("scalar")
    # 1 Assemble validation array
    testFunc_field = TestFunctions(temp_scalar_field.space)
    trialFunc_field = TrialFunctions(temp_scalar_field.space)
    matrix_ufl = UflDecorator()
    i=0
    for test in testFunc_field:
        try:
            j=0
            for subTest in test:
                matrix_ufl.add(conj(subTest)*trialFunc_field[i][j]*dx)
                j+=1
        except:
            matrix_ufl.add(conj(test)*trialFunc_field[i]*dx)
        i+=1
    matrix = matrix_ufl.get_assembled_matrix(
    temp_scalar_field.mesh,
    bcs =[],
    )
    temp_scalar_field.space.FEMWeightSolver = LinearSolver.create_equation_system_solver(matrix)
    expr_ufl = UflDecorator(expr)
    petscVec = expr_ufl.get_assembled_vector(
    temp_scalar_field.mesh,
    bcs=[],
    )
    validation_array = LinearSolver.solve_equation_system_with_predefined_solver(
    temp_scalar_field.space.FEMWeightSolver,
    petscVec,
    )
    # 2 Check alignment
    randomField.scalar_field.evaluate_ufl_expression(expr)
    test_array = randomField.scalar_field.function.x.array.copy()
    # NOTE: This need a higher tolerance than 1e-14
    assert np.linalg.norm(validation_array-test_array) < 1e-13
    logger.info("... passed.")
    
def test_setBoundaryConditions():
    # NOTE: Only scalar field
    # 0 Create a bc list
    randomField.mesh.topology.create_connectivity(
    randomField.mesh.topology.dim - 1,
    randomField.mesh.topology.dim,
    )
    boundary_facets = dolfinx.mesh.exterior_facet_indices(randomField.mesh.topology)
    boundary_dofs = dolfinx.fem.locate_dofs_topological(
    randomField.space_scalar,
    randomField.mesh.topology.dim - 1,
    boundary_facets,
    
    )
    u_bc = randomField.createDolfinxFunction(dim=1)
    u_bc.x.array[:]=1.0+1.0j
    bc = dolfinx.fem.dirichletbc(
    u_bc,
    boundary_dofs,
    )
    bcs = [bc]
    
    # 1 Assemble validation array
    petscArray = randomField.scalar_field.get_petsc_vector().copy()
    dolfinx.fem.petsc.set_bc(
    petscArray,
    bcs,
    )
    validation_array = petscArray.getArray()
    # 2 Check alignment
    randomField.scalar_field.set_boundary_conditions(bcs)
    test_array = randomField.scalar_field.function.x.array[:]
    assert np.linalg.norm(validation_array-test_array) < 1e-14
    logger.info("... passed.")
            
def test_smoothUflTensorExpression():
    logger.info("testing evaluate Ufl expression with smooth")
    from FELiCS.Solvers.LinearSolver import LinearSolver
    # NOTE: So far only scalar field is tested
    # 0 Define a simple ufl expression and smooth factor
    random_scalar_function = randomField.createDolfinxFunction(dim=1)
    expr = random_scalar_function * conj(randomField.test_scalar) * dx
    temp_scalar_field = randomField.createFELiCSField("scalar")
    smoothFactor = 1e-2
    # 1 Assemble validation array
    testFunc_field = TestFunctions(temp_scalar_field.space)
    trialFunc_field = TrialFunctions(temp_scalar_field.space)
    matrix_ufl = UflDecorator()
    coordinate_system = temp_scalar_field.mesh.coordinate_system
    J_hat = coordinate_system.J_hat
    i=0
    for test in testFunc_field:
        iTest = Tensor(
        test,
        coordinate_system,
        mayHaveSpectralDimension=True,
        )
        iFluc = Tensor(
        trialFunc_field[i],
        coordinate_system,
        mayHaveSpectralDimension=True,
        )
        if iTest.order  == 1:
            matrix_ufl.add( ( i_dot(
            iFluc,
            i_conj(iTest),
            ) ).ufl_tens*J_hat*dx)
            matrix_ufl.add((smoothFactor* i_inner(
            i_grad(iFluc),
            i_grad(iTest),
            )).ufl_tens*J_hat*dx)
        elif iTest.order == 0:
            matrix_ufl.add( ( iFluc * i_conj(iTest)).ufl_tens*J_hat*dx)
            matrix_ufl.add((smoothFactor* i_dot(
            i_grad(iFluc),
            i_grad(i_conj(iTest)),
             )).ufl_tens*J_hat*dx)
        else:
            logger.info("Wrong order for test function")
        i+=1
    try:
        matrix_ufl.setCorrectMeshObject(temp_scalar_field.mesh)
    except:
        pass
    matrix = dolfinx.fem.petsc.assemble_matrix(
    dolfinx.fem.form(matrix_ufl.lhs),
    bcs = [],
    )
    matrix.assemble()
    temp_scalar_field.space.FEMSmoothSolver = LinearSolver.create_equation_system_solver(matrix)
    expr_ufl = UflDecorator(expr)
    petscVec = expr_ufl.get_assembled_vector(
    temp_scalar_field.mesh,
    bcs=[],
    )
    validation_array = LinearSolver.solve_equation_system_with_predefined_solver(
    temp_scalar_field.space.FEMSmoothSolver,
    petscVec,
    )
    # 2 Check alignment
    randomField.scalar_field.smooth_ufl_tensor_expression(
    expr,
    smoothFactor,
    )
    test_array = randomField.scalar_field.function.x.array.copy()
    assert np.linalg.norm(validation_array-test_array) < 1e-14
    logger.info("... passed")
    
def test_smooth():
    logger.info("testing smoothing functionality")
    from FELiCS.Solvers.LinearSolver import LinearSolver
    # NOTE: So far only scalar field is tested
    # 0 Define a simple ufl expression and smooth factor
    temp_scalar_field = randomField.createFELiCSField("scalar")
    temp_scalar_field.function.x.array[:] = randomField.scalar_field.function.x.array.copy()
    smoothFactor = 1e-2
    # 1 Assemble validation array
    testFunc_field = TestFunctions(temp_scalar_field.space)
    trialFunc_field = TrialFunctions(temp_scalar_field.space)
    matrix_ufl = UflDecorator()
    coordinate_system = temp_scalar_field.mesh.coordinate_system
    J_hat = coordinate_system.J_hat
    i=0
    for test in testFunc_field:
        iTest = Tensor(
        test,
        coordinate_system,
        mayHaveSpectralDimension=True,
        )
        iFluc = Tensor(
        trialFunc_field[i],
        coordinate_system,
        mayHaveSpectralDimension=True,
        )
        if iTest.order  == 1:
            matrix_ufl.add( ( i_dot(
            iFluc,
            i_conj(iTest),
            ) ).ufl_tens*J_hat*dx)
            matrix_ufl.add((smoothFactor* i_inner(
            i_grad(iFluc),
            i_grad(iTest),
            )).ufl_tens*J_hat*dx)
        elif iTest.order == 0:
            matrix_ufl.add( ( iFluc * i_conj(iTest)).ufl_tens*J_hat*dx)
            matrix_ufl.add((smoothFactor* i_dot(
            i_grad(iFluc),
            i_grad(i_conj(iTest)),
             )).ufl_tens*J_hat*dx)
        else:
            logger.info("Wrong order for test function")
        i+=1
    matrix = matrix_ufl.get_assembled_matrix(
    temp_scalar_field.mesh,
    bcs=[],
    )
    temp_scalar_field.space.FEMSmoothSolver = LinearSolver.create_equation_system_solver(matrix)
    expr_ufl = UflDecorator()
    listOfFields = temp_scalar_field.get_list_of_sub_fields()
    test_FEM     = TestFunctions(temp_scalar_field.space)
    coordinate_system = temp_scalar_field.mesh.coordinate_system
    J_hat = coordinate_system.J_hat
    for i in range(len(listOfFields)):
        iTest = Tensor(
        test_FEM[i],
        coordinate_system,
        mayHaveSpectralDimension=True,
        )
        field     = listOfFields[i]
        expr_ufl += (i_dot(
        field.get_tensor(),
        i_conj(iTest),
        )).ufl_tens*J_hat*dx
    petscVec = expr_ufl.get_assembled_vector(
    temp_scalar_field.mesh,
    bcs=[],
    )
    validation_array = LinearSolver.solve_equation_system_with_predefined_solver(
    temp_scalar_field.space.FEMSmoothSolver,
    petscVec,
    )
    # 2 Check alignment
    randomField.scalar_field.smooth(smoothFactor)
    test_array = randomField.scalar_field.function.x.array.copy()
    assert np.linalg.norm(validation_array-test_array) < 1e-14
    logger.info("... passed")
    
#################################################################
##### moving log files to TESTS folder ##########################
#################################################################
Logger.change_log_location("./TESTS/Fields/logs")
