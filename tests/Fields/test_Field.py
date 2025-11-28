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
# --------------------------------------------------------------
# Tests for Field methods
# --------------------------------------------------------------

def test_calculateL2Norm():
    print("Testing L2 norm calculation")
    # 1. define expressions
    for field in [randomField.scalar_field, randomField.vector_field, randomField.mixed_field]:
        print(" - Testing field of type: ", field._name)
        list = field.getListOfSubFields()
        J_hat = field.mesh.coordinateSystem.J_hat
        for i, subfield in enumerate(list):
            fieldTens = subfield.getTensor()
            if i == 0:
                validExpr = (iDot(fieldTens, iConj(fieldTens))).ufl_tens*J_hat*dx
            else:
                validExpr += (iDot(fieldTens, iConj(fieldTens))).ufl_tens*J_hat*dx 
        # 2. check alignment
        computedL2Norm = field.calculateL2Norm()
        validL2Norm = np.sqrt(dolfinx.fem.assemble_scalar(dolfinx.fem.form(validExpr)))
        assert np.abs(computedL2Norm - validL2Norm) < 1.e-14
    print("... passed.")

def test_getVorticityField():
    # test for 2-D vector
    dim_vector = 2
    print("Testing computation of the vorticity field in 2D")
    randomField     = FieldTestHandler(dim_vector)
    field           = randomField.vector_field
    components      = field.getListOfSubFields()
    componentList   = []  
    for subfield in components:
        componentList.append(subfield.getGradientField())
    dvdx = componentList[1].getListOfSubFields()[0]    
    dudy = componentList[0].getListOfSubFields()[1] 
    vorticity_field_valid = dvdx - dudy
    vorticity_field_test = field.getVorticityField()
    assert np.linalg.norm(vorticity_field_valid.function.x.array[:] - vorticity_field_test.function.x.array[:]) < 1e-14
    # Test for 3D vector
    # NOTE: we do this in the last step
    # dim_vector = 3
    # print("Testing computation of the vorticity field in 2D")
    # randomField     = FieldTestHandler(dim_vector)
    # field           = randomField.vector_field
    # components      = field.getListOfSubFields()
    # componentList   = []  
    # for subfield in components:
    #     componentList.append(subfield.getGradientField())
    # dwdy = componentList[2].getListOfSubFields()[1]
    # dvdz = componentList[1].getListOfSubFields()[2]

    # dudz = componentList[0].getListOfSubFields()[2]
    # dwdx = componentList[2].getListOfSubFields()[0]

    # dvdx = componentList[1].getListOfSubFields()[0]
    # dudy = componentList[0].getListOfSubFields()[1]
    
    # vorticity_x = dwdy - dvdz
    # vorticity_y = dudz - dwdx
    # vorticity_z = dvdx - dudy
    
    # vorticity_field_valid = randomField.createFELiCSField("vector")
    # vorticity_field_valid.setListOfSubFields([vorticity_x, vorticity_y, vorticity_z])
    # vorticity_field_test = field.getVorticityField()
    # assert np.linalg.norm(vorticity_field_valid.function.x.array[:] - vorticity_field_test.function.x.array[:]) < 1e-14

def test_evaluateUflExpression():
    print("testing evaluate Ufl expression")
    from FELiCS.Solvers.LinearSolver import LinearSolver
    # NOTE: So far only scalar field is tested
    # 0 Define a simple ufl expression
    random_scalar_function = randomField.createDolfinxFunction(dim=1)
    expr = random_scalar_function * conj(randomField.test_scalar) * dx
    temp_scalar_field = randomField.createFELiCSField("scalar")
    # 1 Assemble validation array
    if not hasattr(temp_scalar_field.space, 'FEMWeightSolver'):
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
        matrix = matrix_ufl.getAssembledMatrix(temp_scalar_field.mesh,bcs =[])
        temp_scalar_field.space.FEMWeightSolver = LinearSolver.createEquationSystemSolver(matrix)
    expr_ufl = UflDecorator(expr)
    petscVec = expr_ufl.getAssembledVector(temp_scalar_field.mesh, bcs=[])
    validation_array = LinearSolver.solveEquationSystemWithPredefinedSolver(temp_scalar_field.space.FEMWeightSolver, petscVec)
    # 2 Check alignment
    randomField.scalar_field.evaluateUflExpression(expr)
    test_array = randomField.scalar_field.function.x.array.copy()
    # NOTE: This need a higher tolerance than 1e-14
    assert np.linalg.norm(validation_array-test_array) < 1e-13
    print("... passed")
    
def test_setBoundaryConditions():
    # NOTE: Only scalar field
    # 0 Create a bc list
    randomField.mesh.topology.create_connectivity(randomField.mesh.topology.dim - 1, randomField.mesh.topology.dim)
    boundary_facets = dolfinx.mesh.exterior_facet_indices(randomField.mesh.topology)
    boundary_dofs = dolfinx.fem.locate_dofs_topological(
        randomField.space_scalar, randomField.mesh.topology.dim - 1, boundary_facets
    )
    u_bc = randomField.createDolfinxFunction(dim=1)
    u_bc.x.array[:]=1.0+1.0j
    bc = dolfinx.fem.dirichletbc(u_bc, boundary_dofs)
    bcs = [bc]
    
    # 1 Assemble validation array
    petscArray = randomField.scalar_field.getPetscVector().copy()
    dolfinx.fem.petsc.set_bc(petscArray,bcs)
    validation_array = petscArray.getArray()
    # 2 Check alignment
    randomField.scalar_field.setBoundaryConditions(bcs)
    test_array = randomField.scalar_field.function.x.array[:]
    assert np.linalg.norm(validation_array-test_array) < 1e-14

            
