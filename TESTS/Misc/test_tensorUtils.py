# TODO: 
# - do all other functions (only iGrad here) => Xiuyang and Kai?
# - implement logging and replace print statements => ? 
# - create bigger test architecture, at best object oriented, to conduct all tests most efficiently => Sophie


from FELiCS.Misc.tensorUtils import*

import numpy as np
from mpi4py              import MPI
from petsc4py.PETSc      import ScalarType
from ufl                 import (
    dx,TestFunction,TrialFunction,conj,SpatialCoordinate, Dx,dot, inner, grad
)
from dolfinx.fem         import (
    Function,
    functionspace,
    dirichletbc,
    Constant,
    form,
    locate_dofs_topological,
    petsc
)
from dolfinx             import mesh
from FELiCS.Fields.Field import Field
from FELiCS.Fields.Mode  import Mode

#################################################################
##### initialize tests ##########################################
#################################################################

# 0. define parameters
dim_vector           = 3 #dimension of vector function space, 2 or 3 (at the moment: is fixed to 3; that should test the 2, too?) 

# 1. create test mesh and coordinate system
testMesh = mesh.create_rectangle(comm=MPI.COMM_WORLD,
                            points=((0.0, 0.0), (1.0, 1.0)), n=(10, 10),
                            cell_type=mesh.CellType.triangle,
                            ghost_mode=mesh.GhostMode.none)

# 2. create function spaces and their test functions
space_scalar = functionspace(testMesh, ("CG", 2))
space_vector = functionspace(testMesh, ("CG", 2,(dim_vector,)))

test_scalar  = TestFunction(space_scalar)
test_vector  = TestFunction(space_vector)
trial_scalar = TrialFunction(space_scalar)
trial_vector = TrialFunction(space_vector)


# 3. create functions with random coefficient arrays and smooth them
smoothFactor           = 1.e-8
func_scalar            = Function(space_scalar)
func_vector1           = Function(space_vector)
func_vector2           = Function(space_vector)
func_scalar.x.array[:]  = np.random.rand(len(func_scalar.x.array[:])) + 1j*np.random.rand(len(func_scalar.x.array[:]))
func_vector1.x.array[:] = np.random.rand(len(func_vector1.x.array[:]))+ 1j*np.random.rand(len(func_vector1.x.array[:]))
func_vector2.x.array[:] = np.random.rand(len(func_vector2.x.array[:]))+ 1j*np.random.rand(len(func_vector1.x.array[:]))
# scalar smoothing
rhs = func_scalar*conj(test_scalar)*dx
lhs = conj(test_scalar)*trial_scalar*dx + smoothFactor*inner(grad(trial_scalar),grad(test_scalar))*dx
problem = petsc.LinearProblem(lhs, rhs, bcs=[],
                                  petsc_options={"ksp_type": "preonly",
                                                 "pc_type": "lu"})
func_scalar = problem.solve()
# vector smoothing
rhs1 = inner(func_vector1,  test_vector)*dx
rhs2 = inner(func_vector2,  test_vector)*dx
lhs = inner(trial_vector, test_vector)*dx + smoothFactor*inner(grad(trial_vector),grad(test_vector))*dx
problem = petsc.LinearProblem(lhs, rhs1, bcs=[], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
func_vector1 = problem.solve()
problem = petsc.LinearProblem(lhs, rhs2, bcs=[], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
func_vector2 = problem.solve()


def initializeTensorUtils():
    # create tensorUtils specific things
    global J_hat, itest_scalar, itest_vector, tens_vector1, tens_vector2
    testCoordinateSystem = CoordinateSystem(
                             SpatialCoordinate(testMesh),
                             coordinateSystemName,
                             m=m,
                             mesh_dims=(1,1,0),
                             )
    J_hat        = testCoordinateSystem.J_hat
    itest_scalar = Tensor(test_scalar, testCoordinateSystem,  hasSpectralDimension =True) 
    itest_vector = Tensor(test_vector, testCoordinateSystem,  hasSpectralDimension =True) 
    tens_vector1 = Tensor(func_vector1, testCoordinateSystem, hasSpectralDimension = True)
    tens_vector2 = Tensor(func_vector2, testCoordinateSystem, hasSpectralDimension = True)


def getValidGrad():
    t11 = Dx(func_vector1[0],0)
    t12 = Dx(func_vector1[0],1)
    t21 = Dx(func_vector1[1],0)
    t22 = Dx(func_vector1[1],1)
    t31 = Dx(func_vector1[2],0)
    t32 = Dx(func_vector1[2],1)
    if coordinateSystemName == "cartesian": 
        t13 = 1j*m*func_vector1[0]
        t23 = 1j*m*func_vector1[1]
        t33 = 1j*m*func_vector1[2]
    elif coordinateSystemName == "cylindricalfelics": 
        t13 = 1j*m*func_vector1[0]/r
        t23 = 1j*m*func_vector1[1]/r - func_vector1[2]/r
        t33 = 1j*m*func_vector1[2]/r + func_vector1[1]/r
    return t11, t12, t13, t21, t22, t23, t31, t32, t33


def checkIfVectorsAlign(petscVec1, petscVec2):
    assert np.linalg.norm(petscVec1.getArray())> 0
    assert np.linalg.norm(petscVec2.getArray())> 0
    assert np.linalg.norm(petscVec1.getArray()-petscVec2.getArray()) < 1.e-14


def validateGrad():
    t11, t12, t13, t21, t22, t23, t31, t32, t33 = getValidGrad()
    tensor_expr = (iDot((iDot(iGrad(tens_vector1),tens_vector2)),iConj(itest_vector))).ufl_tens*J_hat*dx
    valid_expr  = (t11*func_vector2[0] + t12*func_vector2[1] + t13*func_vector2[2])*conj(test_vector[0])*r*dx
    valid_expr += (t21*func_vector2[0] + t22*func_vector2[1] + t23*func_vector2[2])*conj(test_vector[1])*r*dx
    valid_expr += (t31*func_vector2[0] + t32*func_vector2[1] + t33*func_vector2[2])*conj(test_vector[2])*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2


def test_iGrad():
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validateGrad()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validateGrad()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validateGrad()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validateGrad()





  
    
