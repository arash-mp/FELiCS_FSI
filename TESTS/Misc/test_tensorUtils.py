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
space_dyade = functionspace(testMesh, ("CG", 2,(dim_vector,dim_vector)))

test_scalar  = TestFunction(space_scalar)
test_vector  = TestFunction(space_vector)
test_dyade   = TestFunction(space_dyade)

trial_scalar = TrialFunction(space_scalar)
trial_vector = TrialFunction(space_vector)
trial_dyade   = TrialFunction(space_dyade)


# 3. create functions with random coefficient arrays and smooth them
smoothFactor           = 1.e-8
random_float          = np.random.random() + 1j*np.random.random()

func_scalar1            = Function(space_scalar)
func_scalar2            = Function(space_scalar)
func_vector1           = Function(space_vector)
func_vector2           = Function(space_vector)
func_dyade1            = Function(space_dyade)  
func_dyade2            = Function(space_dyade)

func_scalar1.x.array[:]  = np.random.rand(len(func_scalar1.x.array[:])) + 1j*np.random.rand(len(func_scalar1.x.array[:]))
func_scalar2.x.array[:]  = np.random.rand(len(func_scalar2.x.array[:])) + 1j*np.random.rand(len(func_scalar1.x.array[:]))
func_vector1.x.array[:] = np.random.rand(len(func_vector1.x.array[:]))+ 1j*np.random.rand(len(func_vector1.x.array[:]))
func_vector2.x.array[:] = np.random.rand(len(func_vector2.x.array[:]))+ 1j*np.random.rand(len(func_vector1.x.array[:]))
func_dyade1.x.array[:]  = np.random.rand(len(func_dyade1.x.array[:])) + 1j*np.random.rand(len(func_dyade1.x.array[:]))
func_dyade2.x.array[:]  = np.random.rand(len(func_dyade2.x.array[:])) + 1j*np.random.rand(len(func_dyade2.x.array[:]))

# scalar smoothing
rhs1 = func_scalar1*conj(test_scalar)*dx
rhs2 = func_scalar2*conj(test_scalar)*dx
lhs = conj(test_scalar)*trial_scalar*dx + smoothFactor*inner(grad(trial_scalar),grad(test_scalar))*dx
problem = petsc.LinearProblem(lhs, rhs1, bcs=[],
                                  petsc_options={"ksp_type": "preonly",
                                                 "pc_type": "lu"})
func_scalar1 = problem.solve()
problem = petsc.LinearProblem(lhs, rhs2, bcs=[],
                                  petsc_options={"ksp_type": "preonly",
                                                 "pc_type": "lu"})
func_scalar2 = problem.solve()
# vector smoothing
rhs1 = inner(func_vector1,  test_vector)*dx
rhs2 = inner(func_vector2,  test_vector)*dx
lhs = inner(trial_vector, test_vector)*dx + smoothFactor*inner(grad(trial_vector),grad(test_vector))*dx
problem = petsc.LinearProblem(lhs, rhs1, bcs=[], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
func_vector1 = problem.solve()
problem = petsc.LinearProblem(lhs, rhs2, bcs=[], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
func_vector2 = problem.solve()
# dyade smoothing
rhs1 = inner(func_dyade1,  test_dyade)*dx
rhs2 = inner(func_dyade2,  test_dyade)*dx
lhs = inner(trial_dyade, test_dyade)*dx + smoothFactor*inner(grad(trial_dyade),grad(test_dyade))*dx
problem = petsc.LinearProblem(lhs, rhs1, bcs=[], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
func_dyade1 = problem.solve()
problem = petsc.LinearProblem(lhs, rhs2, bcs=[], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
func_dyade2 = problem.solve()


def initializeTensorUtils():
    # NOTE: Failed on defining dyade tensors.
    
    # create tensorUtils specific things
    # global J_hat, itest_scalar, itest_vector, itest_dyade, tens_vector1, tens_vector2, tens_scalar1, tens_scalar2, tens_dyade1, tens_dyade2
    global J_hat, itest_scalar, itest_vector, tens_vector1, tens_vector2, tens_scalar1, tens_scalar2
    testCoordinateSystem = CoordinateSystem(
                             SpatialCoordinate(testMesh),
                             coordinateSystemName,
                             m=m,
                             mesh_dims=(1,1,0),
                             )
    J_hat        = testCoordinateSystem.J_hat
    itest_scalar = Tensor(test_scalar, testCoordinateSystem,  hasSpectralDimension =True) 
    itest_vector = Tensor(test_vector, testCoordinateSystem,  hasSpectralDimension =True) 
    # itest_dyade   = Tensor(test_dyade, testCoordinateSystem,  hasSpectralDimension =True)
    
    tens_vector1 = Tensor(func_vector1, testCoordinateSystem, hasSpectralDimension = True)
    tens_vector2 = Tensor(func_vector2, testCoordinateSystem, hasSpectralDimension = True)
    tens_scalar1 = Tensor(func_scalar1, testCoordinateSystem, hasSpectralDimension = True)
    tens_scalar2 = Tensor(func_scalar2, testCoordinateSystem, hasSpectralDimension = True)
    # tens_dyade1  = Tensor(func_dyade1, testCoordinateSystem, hasSpectralDimension = True)
    # tens_dyade2  = Tensor(func_dyade2, testCoordinateSystem, hasSpectralDimension = True)


def getValidGrad(func_vector):
    t11 = Dx(func_vector[0],0)
    t12 = Dx(func_vector[0],1)
    t21 = Dx(func_vector[1],0)
    t22 = Dx(func_vector[1],1)
    t31 = Dx(func_vector[2],0)
    t32 = Dx(func_vector[2],1)
    if coordinateSystemName == "cartesian": 
        t13 = 1j*m*func_vector[0]
        t23 = 1j*m*func_vector[1]
        t33 = 1j*m*func_vector[2]
    elif coordinateSystemName == "cylindricalfelics": 
        t13 = 1j*m*func_vector[0]/r
        t23 = 1j*m*func_vector[1]/r - func_vector[2]/r
        t33 = 1j*m*func_vector[2]/r + func_vector[1]/r
    return t11, t12, t13, t21, t22, t23, t31, t32, t33

def getValidScalarGrad(func_scalar):
    t1 = Dx(func_scalar,0)
    t2 = Dx(func_scalar,1)
    if coordinateSystemName == "cartesian": 
        t3 = 1j*m*func_scalar
    elif coordinateSystemName == "cylindricalfelics": 
        t3 = 1j*m*func_scalar/r
    return t1, t2, t3


def checkIfVectorsAlign(petscVec1, petscVec2):
    assert np.linalg.norm(petscVec1.getArray())> 0
    assert np.linalg.norm(petscVec2.getArray())> 0
    assert np.linalg.norm(petscVec1.getArray()-petscVec2.getArray()) < 1.e-14

# --------------------------------------------------------------
# Tests for tensor algebra
# --------------------------------------------------------------
def validate_scalar_minus_float():
    tensor_expr = ((tens_scalar1 - random_float) * iConj(itest_scalar)).ufl_tens*J_hat*dx
    valid_expr  = (func_scalar1 - random_float) * conj(test_scalar)*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_scalar_minus_float():
    print("Testing scalar minus float")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_minus_float()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_minus_float()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_minus_float()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_minus_float()

def validate_float_minus_scalar():
    tensor_expr = ((-1 * tens_scalar1 + random_float) * iConj(itest_scalar)).ufl_tens*J_hat*dx
    valid_expr  = (-1 * func_scalar1 + random_float) * conj(test_scalar)*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_float_minus_scalar():
    print("Testing float minus scalar")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_float_minus_scalar()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_float_minus_scalar()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_float_minus_scalar()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_float_minus_scalar()
    
def validate_scalar_minus_scalar():
    tensor_expr = ((tens_scalar1 - tens_scalar2) * iConj(itest_scalar)).ufl_tens*J_hat*dx
    valid_expr  = (func_scalar1 - func_scalar2) * conj(test_scalar)*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_scalar_minus_scalar():
    print("Testing scalar minus scalar")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_minus_scalar()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_minus_scalar()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_minus_scalar()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_minus_scalar()
    
def validate_scalar_plus_scalar():
    tensor_expr = ((tens_scalar1 + tens_scalar2) * iConj(itest_scalar)).ufl_tens*J_hat*dx
    valid_expr  = (func_scalar1 + func_scalar2) * conj(test_scalar)*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_scalar_plus_scalar():
    print("Testing scalar plus scalar")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_plus_scalar()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_plus_scalar()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_plus_scalar()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_plus_scalar()
    
def validate_float_times_scalar():
    tensor_expr = ((random_float * tens_scalar1) * iConj(itest_scalar)).ufl_tens*J_hat*dx
    valid_expr  = (random_float * func_scalar1) * conj(test_scalar)*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_float_times_scalar():
    print("Testing float times scalar")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_float_times_scalar()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_float_times_scalar()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_float_times_scalar()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_float_times_scalar()
    
def validate_scalar_times_float():
    tensor_expr = ((tens_scalar1 * random_float) * iConj(itest_scalar)).ufl_tens*J_hat*dx
    valid_expr  = (func_scalar1 * random_float) * conj(test_scalar)*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_scalar_times_float():
    print("Testing scalar times float.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_times_float()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_times_float()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_times_float()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_times_float()
    
def validate_scalar_times_scalar():
    tensor_expr = ((tens_scalar1 * tens_scalar2) * iConj(itest_scalar)).ufl_tens*J_hat*dx
    valid_expr  = (func_scalar1 * func_scalar2) * conj(test_scalar)*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_scalar_times_scalar():
    print("Testing scalar times scalar.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_times_scalar()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_times_scalar()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_times_scalar()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_times_scalar()
    
def validate_scalar_divide_scalar():
    tensor_expr = ((tens_scalar1 / tens_scalar2) * iConj(itest_scalar)).ufl_tens*J_hat*dx
    valid_expr  = (func_scalar1 / func_scalar2) * conj(test_scalar)*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_scalar_divide_scalar():
    # NOTE: The unit test didn't pass through in the cylindrical coordinate system
    print("Testing scalar divide scalar.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_divide_scalar()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_divide_scalar()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_divide_scalar()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_divide_scalar()

def validate_scalar_exp_to_positive_float():
    tensor_expr = ((tens_scalar1 ** random_float) * iConj(itest_scalar)).ufl_tens*J_hat*dx
    valid_expr  = (func_scalar1 ** random_float) * conj(test_scalar)*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_scalar_exp_to_positive_float():
    # NOTE: The unit test didn't pass through in the cylindrical coordinate system
    print("Testing scalar exponent to positive float.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_exp_to_positive_float()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_exp_to_positive_float()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_exp_to_positive_float()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_exp_to_positive_float()

def validate_scalar_exp_to_negetive_float():
    tensor_expr = ((tens_scalar1 ** (-1 * random_float)) * iConj(itest_scalar)).ufl_tens*J_hat*dx
    valid_expr  = (func_scalar1 ** (-1 * random_float)) * conj(test_scalar)*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_scalar_exp_to_negetive_float():
    # NOTE: The unit test didn't pass through in the cylindrical coordinate system
    print("Testing scalar exponent to negetive float.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_exp_to_negetive_float()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_exp_to_negetive_float()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_exp_to_negetive_float()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_exp_to_negetive_float()

def validate_vector_dot_vector():
    tensor_expr = (iDot(tens_vector1, tens_vector2) * iConj(itest_scalar)).ufl_tens*J_hat*dx
    valid_expr  = (dot(func_vector1, func_vector2)) * conj(test_scalar)*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_vector_dot_vector():
    print("Testing vector dot vector.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_dot_vector()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_dot_vector()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_dot_vector()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_dot_vector()
    
def validate_scalar_times_vector():
    tensor_expr = (iDot(tens_scalar1 * tens_vector1, iConj(itest_vector))).ufl_tens*J_hat*dx
    valid_expr  = (dot(func_scalar1 * func_vector1, conj(test_vector)))*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_scalar_times_vector():
    print("Testing scalar times vector.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_times_vector()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_times_vector()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_times_vector()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_times_vector()
    
def validate_vector_times_scalar():
    tensor_expr = (iDot(tens_vector1 * tens_scalar1, iConj(itest_vector))).ufl_tens*J_hat*dx
    valid_expr  = (dot(func_vector1 * func_scalar1, conj(test_vector)))*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_vector_times_scalar():
    print("Testing vector times scalar.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_times_scalar()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_times_scalar()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_times_scalar()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_times_scalar()

def validate_one_divide_scalar_times_vector():
    # NOTE: Don't know how to define 1 / scalar
    tensor_expr = (iDot((1 + 1j) / tens_scalar1 * tens_vector1, iConj(itest_vector))).ufl_tens*J_hat*dx
    valid_expr  = (dot((1 + 1j) / func_scalar1 * func_vector1, conj(test_vector)))*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_one_divide_scalar_times_vector():
    print("Testing one divide by scalar times vector.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_one_divide_scalar_times_vector()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_one_divide_scalar_times_vector()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_one_divide_scalar_times_vector()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_one_divide_scalar_times_vector()
    
def validate_vector_divided_by_scalar():
    # NOTE: Don't know how to define 1 / scalar
    tensor_expr = (iDot(tens_vector1 / tens_scalar1 , iConj(itest_vector))).ufl_tens*J_hat*dx
    valid_expr  = (dot(func_vector1 / func_scalar1 , conj(test_vector)))*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_vector_divided_by_scalar():
    print("Testing vector divided by scalar.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_divided_by_scalar()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_divided_by_scalar()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_divided_by_scalar()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_divided_by_scalar()

def validate_vector_dot_dyade():
    t11, t12, t13, t21, t22, t23, t31, t32, t33 = getValidGrad(func_vector2)
    tensor_expr = (iDot(iDot(tens_vector1, iGrad(tens_vector2)) , iConj(itest_vector))).ufl_tens*J_hat*dx
    valid_expr  = (func_vector1[0] * t11 + func_vector1[1] * t21 + func_vector1[2] * t31)*conj(test_vector[0])*r*dx
    valid_expr += (func_vector1[0] * t12 + func_vector1[1] * t22 + func_vector1[2] * t32)*conj(test_vector[1])*r*dx
    if not (coordinateSystemName == "cartesian" and m == 0):
        valid_expr += (func_vector1[0] * t13 + func_vector1[1] * t23 + func_vector1[2] * t33)*conj(test_vector[2])*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_vector_dot_dyade():
    print("Testing vector dot dyade.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_dot_dyade()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_dot_dyade()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_dot_dyade()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_dot_dyade()

def validate_dyade_dot_vector():
    t11, t12, t13, t21, t22, t23, t31, t32, t33 = getValidGrad(func_vector1)
    tensor_expr = (iDot(iDot(iGrad(tens_vector1), tens_vector2) , iConj(itest_vector))).ufl_tens*J_hat*dx
    valid_expr  = (t11 * func_vector2[0]+ t12 * func_vector2[1] + t13 * func_vector2[2])*conj(test_vector[0])*r*dx
    valid_expr += (t21 * func_vector2[0]+ t22 * func_vector2[1] + t23 * func_vector2[2])*conj(test_vector[1])*r*dx
    valid_expr += (t31 * func_vector2[0]+ t32 * func_vector2[1] + t33 * func_vector2[2])*conj(test_vector[2])*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_dyade_dot_vector():
    print("Testing dyade dot vector.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_dot_vector()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_dot_vector()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_dot_vector()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_dot_vector()

def validate_vector_dot_dyade_dot_vector():
    t11, t12, t13, t21, t22, t23, t31, t32, t33 = getValidGrad(func_vector1)
    tensor_expr = ((iDot(iDot(tens_vector1,iGrad(tens_vector1)), tens_vector2)) * iConj(itest_scalar)).ufl_tens*J_hat*dx
    valid_expr  = (func_vector1[0] * t11 + func_vector1[1] * t21 + func_vector1[2] * t31) * func_vector2[0] * conj(test_scalar)*r*dx
    valid_expr += (func_vector1[0] * t12 + func_vector1[1] * t22 + func_vector1[2] * t32) * func_vector2[1] * conj(test_scalar)*r*dx
    if not (coordinateSystemName == "cartesian" and m == 0):
        valid_expr += (func_vector1[0] * t13 + func_vector1[1] * t23 + func_vector1[2] * t33) * func_vector2[2] * conj(test_scalar)*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_vector_dot_dyade_dot_vector():
    print("Testing vector dot dyade dot vector.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_dot_dyade_dot_vector()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_dot_dyade_dot_vector()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_dot_dyade_dot_vector()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_dot_dyade_dot_vector()
    
def validate_dyade_dot_dyade():
    # NOTE: Not implemented yet
    pass

def test_dyade_dot_dyade():
    print("Testing dyade dot dyade.")
    raise NotImplementedError("Not implemented yet.")
    
def validate_dyade_inner_dyade():
    t11_1, t12_1, t13_1, t21_1, t22_1, t23_1, t31_1, t32_1, t33_1 = getValidGrad(func_vector1)
    t11_2, t12_2, t13_2, t21_2, t22_2, t23_2, t31_2, t32_2, t33_2 = getValidGrad(func_vector2)
    tensor_expr = (iInner(iGrad(tens_vector1),iGrad(tens_vector2)) * iConj(itest_scalar)).ufl_tens*J_hat*dx
    valid_expr = ((t11_1 * t11_2 + t21_1 * t21_2 + t31_1 * t31_2) + (t12_1 * t12_2 + t22_1 * t22_2 + t32_1 * t32_2) + (t13_1 * t13_2 + t23_1 * t23_2 + t33_1 * t33_2)) * conj(test_scalar)*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_dyade_inner_dyade():
    print("Testing dyade inner dyade.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_inner_dyade()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_inner_dyade()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_inner_dyade()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_inner_dyade()
    
def validate_dyade_transpose():
    t11, t12, t13, t21, t22, t23, t31, t32, t33 = getValidGrad(func_vector1)
    tensor_expr = (iDot(iDot(iT(iGrad(tens_vector1)),tens_vector2), iConj(itest_vector))).ufl_tens*J_hat*dx
    valid_expr  = (t11 * func_vector2[0]+ t21 * func_vector2[1] + t31 * func_vector2[2])*conj(test_vector[0])*r*dx
    valid_expr += (t12 * func_vector2[0]+ t22 * func_vector2[1] + t32 * func_vector2[2])*conj(test_vector[1])*r*dx
    if not (coordinateSystemName == "cartesian" and m == 0):
        valid_expr += (t13 * func_vector2[0]+ t23 * func_vector2[1] + t33 * func_vector2[2])*conj(test_vector[2])*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_dyade_transpose():
    print("Testing dyade transpose.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_transpose()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_transpose()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_transpose()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_transpose()
    
def validate_dyade_trace():
    t11, t12, t13, t21, t22, t23, t31, t32, t33 = getValidGrad(func_vector1)
    tensor_expr = (iTr(iGrad(tens_vector1)) * iConj(itest_scalar)).ufl_tens*J_hat*dx
    valid_expr  = (t11 * t11+ t22 * t22 + t33 * t33)*conj(test_scalar)*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_dyade_trace():
    print("Testing dyade trace.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_trace()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_trace()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_trace()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_trace()
    
def validate_scalar_conj():
    tensor_expr = (iConj(tens_scalar1) * iConj(itest_scalar)).ufl_tens*J_hat*dx
    valid_expr  = (conj(func_scalar1))*conj(test_scalar)*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_scalar_conj():
    print("Testing conjugate of scalar.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_conj()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_conj()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_conj()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_conj()
    
def validate_vector_conj():
    tensor_expr = (iDot(iConj(tens_vector1), iConj(itest_vector))).ufl_tens*J_hat*dx
    valid_expr  = dot(conj(func_vector1),conj(test_vector))*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_vector_conj():
    print("Testing conjugate of vector.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_conj()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_conj()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_conj()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_vector_conj()
    
def validate_dyade_hermitian():
    t11, t12, t13, t21, t22, t23, t31, t32, t33 = getValidGrad(func_vector1)
    tensor_expr = (iDot(iDot(iT(iConj(iGrad(tens_vector1))), tens_vector2), iConj(itest_vector))).ufl_tens*J_hat*dx
    valid_expr = (conj(t11) * func_vector2[0]+ conj(t21) * func_vector2[1] + conj(t31) * func_vector2[2])*conj(test_vector[0])*r*dx
    valid_expr += (conj(t12) * func_vector2[0]+ conj(t22) * func_vector2[1] + conj(t32) * func_vector2[2])*conj(test_vector[1])*r*dx
    if not (coordinateSystemName == "cartesian" and m == 0):
        valid_expr += (conj(t13) * func_vector2[0]+ conj(t23) * func_vector2[1] + conj(t33) * func_vector2[2])*conj(test_vector[2])*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_dyade_hermitian():
    print("Testing hermitian transpose of dyade.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_hermitian()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_hermitian()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_hermitian()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_hermitian()
    
def validate_identity():
    # NOTE: Not implemented yet
    pass

def test_identity():
    print("Testing identity functionality.")
    raise NotImplementedError("Not implemented yet.")
    
def validate_grad_scalar_dot_grad_scalar():
    t1_1, t2_1, t3_1 = getValidScalarGrad(func_scalar1)
    t1_2, t2_2, t3_2 = getValidScalarGrad(func_scalar2)
    tensor_expr = (iDot(iGrad(tens_scalar1), iGrad(tens_scalar2)) * iConj(itest_scalar)).ufl_tens*J_hat*dx
    valid_expr  = (t1_1*t1_2 + t2_1*t2_2 + t3_1*t3_2) * conj(test_scalar)*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_grad_scalar_dot_grad_scalar():
    print("Testing grad of scalar dot grad of scalar.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_grad_scalar_dot_grad_scalar()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_grad_scalar_dot_grad_scalar()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_grad_scalar_dot_grad_scalar()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_grad_scalar_dot_grad_scalar()
    
def validate_scalar_times_grad_scalar():
    t1_2, t2_2, t3_2 = getValidScalarGrad(func_scalar2)
    tensor_expr = (iDot(tens_scalar1 * iGrad(tens_scalar2) , iConj(itest_vector))).ufl_tens*J_hat*dx
    valid_expr  = (func_scalar1 * t1_2) * conj(test_vector[0])*r*dx
    valid_expr  += (func_scalar1 * t2_2) * conj(test_vector[1])*r*dx
    if not (m == 0):
        valid_expr  += (func_scalar1 * t3_2) * conj(test_vector[2])*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_scalar_times_grad_scalar():
    print("Testing scalar times grad of scalar.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_times_grad_scalar()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_times_grad_scalar()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_times_grad_scalar()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_scalar_times_grad_scalar()
    
def validate_grad_scalar_times_scalar():
    t1_1, t2_1, t3_1 = getValidScalarGrad(func_scalar1)    
    tensor_expr = (iDot(iGrad(tens_scalar1) * tens_scalar2 , iConj(itest_vector))).ufl_tens*J_hat*dx
    valid_expr  = (t1_1 * func_scalar2) * conj(test_vector[0])*r*dx
    valid_expr  += (t2_1 * func_scalar2) * conj(test_vector[1])*r*dx
    if not (m == 0):
        valid_expr  += (t3_1 * func_scalar2) * conj(test_vector[2])*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_grad_scalar_times_scalar():
    print("Testing grad of scalar times scalar.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_grad_scalar_times_scalar()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_grad_scalar_times_scalar()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_grad_scalar_times_scalar()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_grad_scalar_times_scalar()
    
def validate_one_divide_scalar_times_grad_scalar():
    # NOTE: Not implemented yet
    pass

def test_one_divide_scalar_times_grad_scalar():
    print("Testing one divide scalar times grad of scalar.")
    raise NotImplementedError("Not implemented yet.")

def validate_grad_scalar_divided_by_scalar():
    # NOTE: Not implemented yet
    pass

def test_grad_scalar_divided_by_scalar():
    print("Testing grad of scalar divided by scalar.")
    raise NotImplementedError("Not implemented yet.")

def validate_grad_scalar_dot_dyade():
    t11, t12, t13, t21, t22, t23, t31, t32, t33 = getValidGrad(func_vector1)
    t1, t2, t3 = getValidScalarGrad(func_scalar1)
    tensor_expr = (iDot(iDot(iGrad(tens_scalar1), iGrad(tens_vector1)) , iConj(itest_vector))).ufl_tens*J_hat*dx
    valid_expr  = (t11 * t1 + t21 * t2 + t31 * t3)*conj(test_vector[0])*r*dx
    valid_expr += (t12 * t1 + t22 * t2 + t32 * t3)*conj(test_vector[1])*r*dx
    if not (coordinateSystemName == "cartesian" and m == 0):
        valid_expr += (t13 * t1 + t23 * t2 + t33 * t3)*conj(test_vector[2])*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_grad_scalar_dot_dyade():
    print("Testing grad of scalar dot dyade")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_grad_scalar_dot_dyade()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_grad_scalar_dot_dyade()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_grad_scalar_dot_dyade()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_grad_scalar_dot_dyade()

def validate_dyade_dot_grad_scalar():
    t11, t12, t13, t21, t22, t23, t31, t32, t33 = getValidGrad(func_vector1)
    t1, t2, t3 = getValidScalarGrad(func_scalar2)
    tensor_expr = (iDot(iDot(iGrad(tens_vector1), iGrad(tens_scalar2)) , iConj(itest_vector))).ufl_tens*J_hat*dx
    valid_expr  = (t11 * t1 + t12 * t2 + t13 * t3)*conj(test_vector[0])*r*dx
    valid_expr += (t21 * t1 + t22 * t2 + t23 * t3)*conj(test_vector[1])*r*dx
    valid_expr += (t31 * t1+ t32 * t2 + t33 * t3)*conj(test_vector[2])*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_dyade_dot_grad_scalar():
    print("Testing dyade dot grad of scalar.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_dot_grad_scalar()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_dot_grad_scalar()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_dot_grad_scalar()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_dyade_dot_grad_scalar()

def validate_grad_scalar_dot_dyade_dot_grad_scalar():
    # NOTE: Not implemented yet
    pass

def test_grad_scalar_dot_dyade_dot_grad_scalar():
    print("Testing grad of scalar dot dyade dot grad of scalar.")
    raise NotImplementedError("Not implemented yet.")

def validate_grad_scalar_conj():
    t1, t2, t3 = getValidScalarGrad(func_scalar1)
    tensor_expr = (iDot(iConj(iGrad(tens_scalar1)), iConj(itest_vector))).ufl_tens*J_hat*dx
    valid_expr  = (conj(t1) * conj(test_vector[0]))*r*dx
    valid_expr  += (conj(t2) * conj(test_vector[1]))*r*dx
    if not (m == 0):
        valid_expr  += (conj(t3) * conj(test_vector[2]))*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2

def test_grad_scalar_conj():
    print("Testing conjugate of grad of scalar.")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_grad_scalar_conj()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_grad_scalar_conj()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_grad_scalar_conj()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_grad_scalar_conj()

# ---------------------------------------------------------------
# Tests for tensor analysis
# ---------------------------------------------------------------

def validate_div_vector():
    t11, t12, t13, t21, t22, t23, t31, t32, t33 = getValidGrad(func_vector1)
    tensor_expr = (iDiv(tens_vector1) * iConj(itest_scalar)).ufl_tens*J_hat*dx
    valid_expr  = (t11 + t22 + t33)*conj(test_scalar)*r*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()

    checkIfVectorsAlign(res1, res2)
    return res1, res2


def test_iDiv():
    print("Testing iDiv")
    global coordinateSystemName, m, r
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0; r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_div_vector()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20); r = 1.
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_div_vector()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0; r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_div_vector()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20); r   = SpatialCoordinate(testMesh)[1]
    initializeTensorUtils()
    print("\n\n  #####, ", coordinateSystemName, m)
    validate_div_vector()





  
    
