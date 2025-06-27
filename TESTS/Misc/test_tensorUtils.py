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
x   = SpatialCoordinate(testMesh)
r   = x[1]

# 2. create function spaces and their test functions
space_scalar = functionspace(testMesh, ("CG", 2))
space_vector = functionspace(testMesh, ("CG", 2,(dim_vector,)))

test_scalar  = TestFunction(space_scalar)
test_vector  = TestFunction(space_vector)
trial_scalar = TrialFunction(space_scalar)
trial_vector = TrialFunction(space_vector)


# 3. create functions with random coefficient arrays and smooth them
smoothFactor           = 1.
func_scalar            = Function(space_scalar)
func_vector1           = Function(space_vector)
func_vector2           = Function(space_vector)
func_scalar.x.array[:]  = np.random.rand(len(func_scalar.x.array[:]))
func_vector1.x.array[:] = np.random.rand(len(func_vector1.x.array[:]))
func_vector2.x.array[:] = np.random.rand(len(func_vector2.x.array[:]))
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
    assert np.linalg.norm(petscVec1.getArray()-petscVec2.getArray()) == 0


def validateGrad():
    t11, t12, t13, t21, t22, t23, t31, t32, t33 = getValidGrad()
    tensor_expr = (iDot((iDot(iGrad(tens_vector1),tens_vector2)),iConj(itest_vector))).ufl_tens*J_hat*dx
    valid_expr  = (t11*func_vector2[0] + t12*func_vector2[1] + t13*func_vector2[2])*conj(test_vector[0])*dx
    valid_expr += (t21*func_vector2[0] + t22*func_vector2[1] + t23*func_vector2[2])*conj(test_vector[1])*dx
    valid_expr += (t31*func_vector2[0] + t32*func_vector2[1] + t33*func_vector2[2])*conj(test_vector[2])*dx
    res1 = petsc.assemble_vector(form(tensor_expr))
    res2 = petsc.assemble_vector(form(valid_expr))
    res1.assemble()
    res2.assemble()
    checkIfVectorsAlign(res1, res2)
    return res1, res2


def test_iGrad():
    global coordinateSystemName, m
    # 1. Cartesian Cordinates, m=0
    coordinateSystemName = "cartesian"; m = 0
    initializeTensorUtils()
    #validateGrad()
    # 2. Cartesian Cordinates, m= random integer
    coordinateSystemName = "cartesian"; m = np.random.randint(20)
    initializeTensorUtils()
    validateGrad()
    # 3. Cylindrical Cordinates, m= 0
    coordinateSystemName = "cylindricalfelics"; m = 0
    initializeTensorUtils()
    validateGrad()
    # 4. Cylindrical Cordinates, m= random integer
    coordinateSystemName = "cylindricalfelics"; m = np.random.randint(20)
    initializeTensorUtils()
    validateGrad()



#    space  = u1.space
#    space2 = u.space
#    test   = TestFunction(space)
#    test2  = TestFunction(space2)
#    itest  = Tensor(test, equation._coordinateSystem, containsTestFunction=True)
#    itest2 = Tensor(test2, equation._coordinateSystem, containsTestFunction=True)
#
#    u_mean = Tensor(u_old, equation._coordinateSystem, containsFluctuation = True)
#
#
#    q_old_field = Field(FEMSpaces.VMixed,mesh)
#    [u_,p_] = q_old_field.getListOfSingleFields()
#    mapping_u1 = q_old_field.space.sub(0).sub(0).collapse()[1]
#    mapping_u2 = q_old_field.space.sub(0).sub(1).collapse()[1]
#    mapping_u3 = q_old_field.space.sub(0).sub(2).collapse()[1]
#    mapping_u1_ = meanFlow._fieldDict['u'].function_space.sub(0).collapse()[1]
#    mapping_u2_ = meanFlow._fieldDict['u'].function_space.sub(1).collapse()[1]
#    mapping_u3_ = meanFlow._fieldDict['u'].function_space.sub(2).collapse()[1]
#    q_old_field.function.x.array[mapping_u1] = 1.e4*0.2903858 * meanFlow._fieldDict['u'].x.array[mapping_u2_]
#    q_old_field.function.x.array[mapping_u2] = 1.e4*meanFlow._fieldDict['u'].x.array[mapping_u2_]
#    q_old_field.function.x.array[mapping_u3] = 1.e4*0.59012834*meanFlow._fieldDict['u'].x.array[mapping_u2_]
#
#    [u_field, p_old_field] = q_old_field.getListOfSingleFields()
#    [u1_field, u2_field, u3_field]  = u_field.getListOfSingleFields()
#
#    u_mean = Tensor(u_field.function, equation._coordinateSystem, containsFluctuation = True)
#
#    expression = (iDot(u_mean, iConj(itest2))).ufl_tens*J_hat*dx
#    u_field.smoothUflTensorExpression(expression, 100.)
#
#    u_mean = Tensor(u_field.function, equation._coordinateSystem, containsFluctuation = True)
#
#
#    m   = equation.m
#    print("#### m: ", m)
#    x   = SpatialCoordinate(mesh)
#    r   = x[1]
#    ##############################################################################################
#    print("")
#    print("#### iGrad ####")
#    tensor = (iDot(iGrad(u_mean),iConj(itest2)))
#
#    # values calculated without tensor utils
#    t11 = Dx(u_field.function[0],0)
#    t12 = Dx(u_field.function[0],1)
#    t13 = 1j*m*u_field.function[0]/r
#    t21 = Dx(u_field.function[1],0)
#    t22 = Dx(u_field.function[1],1)
#    t23 = 1j*m*u_field.function[1]/r - u_field.function[2]/r
#    t31 = Dx(u_field.function[2],0)
#    t32 = Dx(u_field.function[2],1)
#    t33 = (1j*m*u_field.function[2]/r + u_field.function[1]/r)
#
#    expression   = tensor.ufl_tens[0]*J_hat*dx
#    expression_T = (t11*conj(test2[0]) + t12*conj(test2[1]) + t13*conj(test2[2]))*dx
#    u_big1.evaluateUflTensorExpression(expression)
#    u_big2.evaluateUflExpression(expression_T)
#    print("1:   ", np.linalg.norm(u_big1.getCoefficientArray() - u_big2.getCoefficientArray()), np.linalg.norm(u_big1.getCoefficientArray()), np.linalg.norm(u_big2.getCoefficientArray()))
#
#    expression   = tensor.ufl_tens[1]*J_hat*r*dx
#    expression_T = (t21*conj(test2[0]) + t22*conj(test2[1]) + t23*conj(test2[2]))*r*dx
#    u_big1.evaluateUflTensorExpression(expression)
#    u_big2.evaluateUflExpression(expression_T)
#    print("2:   ", np.linalg.norm(u_big1.getCoefficientArray() - u_big2.getCoefficientArray()), np.linalg.norm(u_big1.getCoefficientArray()), np.linalg.norm(u_big2.getCoefficientArray()))
#
#    expression   = tensor.ufl_tens[2]*J_hat*r*r*dx
#    expression_T = (t31*conj(test2[0]) + t32*conj(test2[1]) + t33*conj(test2[2]))*r*dx
#    u_big1.evaluateUflTensorExpression(expression)
#    u_big2.evaluateUflExpression(expression_T)
#    print("3:   ", np.linalg.norm(u_big1.getCoefficientArray() - u_big2.getCoefficientArray()), np.linalg.norm(u_big1.getCoefficientArray()), np.linalg.norm(u_big2.getCoefficientArray()))
#
#
#    tensor = (iInner(iGrad(u_mean), iGrad(u_mean)))
#    expression   = (tensor*iConj(itest)).ufl_tens*J_hat*r*r*dx
#    expression_T = (t11*t11 + t12*t12 + t13*t13 + t21*t21 + t22*t22 + t23*t23 + t31*t31 + t32*t32 + t33*t33)*conj(test)*r*r*dx
#    u1.evaluateUflTensorExpression(expression)
#    u2.evaluateUflExpression(expression_T)
#    print("all: ", np.linalg.norm(u1.getCoefficientArray() - u2.getCoefficientArray()), np.linalg.norm(u1.getCoefficientArray()), np.linalg.norm(u2.getCoefficientArray()))
#
#
#    tensor = iDot((iDot(iGrad(u_mean),u_mean)),iConj(itest2))
#    expression   = tensor.ufl_tens*J_hat*r*dx
#    expression_T = (t11*u_field.function[0] + t12*u_field.function[1] + t13*u_field.function[2])*conj(test2[0])
#    expression_T+= (t21*u_field.function[0] + t22*u_field.function[1] + t23*u_field.function[2])*conj(test2[1])
#    expression_T+= (t31*u_field.function[0] + t32*u_field.function[1] + t33*u_field.function[2])*conj(test2[2])
#    u_big1.evaluateUflTensorExpression(expression)
#    u_big2.evaluateUflExpression(expression_T*r*dx)
#    print("all2:", np.linalg.norm(u_big1.getCoefficientArray() - u_big2.getCoefficientArray()), np.linalg.norm(u_big1.getCoefficientArray()), np.linalg.norm(u_big2.getCoefficientArray()))
#
#    # test function values calculated without tensor utils
#    test2_ = conj(test2)
#    t11_ = Dx(test2_[0],0)
#    t12_ = Dx(test2_[0],1)
#    t13_ = -1j*m*test2_[0]/r
#    t21_ = Dx(test2_[1],0)
#    t22_ = Dx(test2_[1],1)
#    t23_ = -1j*m*test2_[1]/r - test2_[2]/r
#    t31_ = Dx(test2_[2],0)
#    t32_ = Dx(test2_[2],1)
#    t33_ = (-1j*m*test2_[2]/r + test2_[1]/r)
#
#    #t11  = Dx(u_field.function[0],0)
#    #t12  = Dx(u_field.function[0],1)
#    #t13  = 1j*m*u_field.function[0]/r
#    #t21  = Dx(u_field.function[1],0)
#    #t22  = Dx(u_field.function[1],1)
#    #t23  = 1j*m*u_field.function[1]/r - u_field.function[2]/r
#    #t31  = Dx(u_field.function[2],0)
#    #t32  = Dx(u_field.function[2],1)
#    #t33  = (1j*m*u_field.function[2]/r + u_field.function[1]/r)
#
#
#
#
#    tensor = iInner(iGrad(u_mean), iGrad(iConj(itest2)))
#    expression   = tensor.ufl_tens*J_hat*r*r*dx
#    expression_T = (t11*t11_ + t12*t12_ + t13*t13_)
#    expression_T+= (t21*t21_ + t22*t22_ + t23*t23_)
#    expression_T+= (t31*t31_ + t32*t32_ + t33*t33_)
#    u_big1.evaluateUflTensorExpression(expression)
#    u_big2.evaluateUflExpression(expression_T*r*r*dx)
#    print("all3:", np.linalg.norm(u_big1.getCoefficientArray() - u_big2.getCoefficientArray()), np.linalg.norm(u_big1.getCoefficientArray()), np.linalg.norm(u_big2.getCoefficientArray()))
#
#    # => OPEN QUESTION: why is the 3rd row devided by r? Does that make any sense?
#    # => repeat with finer grid: why is the error so high?
#
#    ##############################################################################################
#    print("")
#    print("#### iT(iGrad) ####")
#
#    tensor = (iInner((iT(iGrad(u_mean))), iGrad(u_mean)))
#
#    expression   = (tensor*iConj(itest)).ufl_tens*J_hat*r*r*dx
#    expression_T = (t11*t11 + t21*t12 + t31*t13 + t12*t21 + t22*t22 + t32*t23 + t13*t31 + t23*t32 + t33*t33)*conj(test)*r*r*dx
#    u1.evaluateUflTensorExpression(expression)
#    u2.evaluateUflExpression(expression_T)
#    print("all: ", np.linalg.norm(u1.getCoefficientArray() - u2.getCoefficientArray()), np.linalg.norm(u1.getCoefficientArray()), np.linalg.norm(u2.getCoefficientArray()))
#
#    tensor = iDot((iDot(iT(iGrad(u_mean)),u_mean)),iConj(itest2))
#    expression   = tensor.ufl_tens*J_hat*dx
#    expression_T = (t11*u_field.function[0] + t21*u_field.function[1] + t31*u_field.function[2])*conj(test2[0])
#    expression_T+= (t12*u_field.function[0] + t22*u_field.function[1] + t32*u_field.function[2])*conj(test2[1])
#    expression_T+= (t13*u_field.function[0] + t23*u_field.function[1] + t33*u_field.function[2])*conj(test2[2])
#    u_big1.evaluateUflTensorExpression(expression)
#    u_big2.evaluateUflExpression(expression_T*dx)
#    print("all2:", np.linalg.norm(u_big1.getCoefficientArray() - u_big2.getCoefficientArray()), np.linalg.norm(u_big1.getCoefficientArray()), np.linalg.norm(u_big2.getCoefficientArray()))
#
#    #tensor = (iInner(iT(iGrad(u_mean)) + iGrad(u_mean), iGrad(u_mean)))
#    tensor = (iInner((iGrad(u_mean)) + iT(iGrad(u_mean)), iGrad(u_mean)))
#    #tensor = (iInner(iT(iGrad(u_mean)), iGrad(u_mean))) + iInner(iGrad(u_mean),iGrad(u_mean))
#
#    expression   = (tensor*iConj(itest)).ufl_tens*J_hat*r*r*dx
#    expression_T = (2.*t11*t11 + (t12+t21)*t12 + (t13+t31)*t13 + (t21+t12)*t21 + 2.*t22*t22 + (t23+t32)*t23 + (t31+t13)*t31 + (t32+t23)*t32 + 2.*t33*t33)*conj(test)*r*r*dx
#    u1.evaluateUflTensorExpression(expression)
#    u2.evaluateUflExpression(expression_T)
#    print("all3:", np.linalg.norm(u1.getCoefficientArray() - u2.getCoefficientArray()), np.linalg.norm(u1.getCoefficientArray()), np.linalg.norm(u2.getCoefficientArray()))
#
#
#    tensor = iDot((iDot(iT(iGrad(u_mean))+ iGrad(u_mean),u_mean)),iConj(itest2))
#    #tensor  = iDot((iDot(iT(iGrad(u_mean)),u_mean)),iConj(itest2))
#    #tensor += iDot((iDot(iGrad(u_mean),u_mean)),iConj(itest2))
#    expression   = tensor.ufl_tens*J_hat*dx
#    expression_T = (2.*t11   *u_field.function[0] + (t21+t12)*u_field.function[1] + (t31+t13)*u_field.function[2])*conj(test2[0])
#    expression_T+= ((t12+t21)*u_field.function[0] + 2.*t22   *u_field.function[1] + (t32+t23)*u_field.function[2])*conj(test2[1])
#    expression_T+= ((t13+t31)*u_field.function[0] + (t23+t32)*u_field.function[1] + 2.*t33   *u_field.function[2])*conj(test2[2])
#    u_big1.evaluateUflTensorExpression(expression)
#    u_big2.evaluateUflExpression(expression_T*dx)
#    print("all4:", np.linalg.norm(u_big1.getCoefficientArray() - u_big2.getCoefficientArray()), np.linalg.norm(u_big1.getCoefficientArray()), np.linalg.norm(u_big2.getCoefficientArray()))
#
#
#
#    #### substraction!
#    #tensor = (iInner(iT(iGrad(u_mean)) + iGrad(u_mean), iGrad(u_mean)))
#    tensor = (iInner((iGrad(u_mean)) - iT(iGrad(u_mean)), iGrad(u_mean)))
#    #tensor = (iInner(iGrad(u_mean), iGrad(u_mean))) - iInner(iT(iGrad(u_mean)),iGrad(u_mean))
#
#    expression   = (tensor*iConj(itest)).ufl_tens*J_hat*r*r*dx
#    expression_T = ((t12-t21)*t12 + (t13-t31)*t13 + (t21-t12)*t21 + (t23-t32)*t23 + (t31-t13)*t31 + (t32-t23)*t32 )*conj(test)*r*r*dx
#    u1.evaluateUflTensorExpression(expression)
#    u2.evaluateUflExpression(expression_T)
#    print("all5:", np.linalg.norm(u1.getCoefficientArray() - u2.getCoefficientArray()), np.linalg.norm(u1.getCoefficientArray()), np.linalg.norm(u2.getCoefficientArray()))
#
#    #tensor = (iInner(iT(iGrad(u_mean)) + iGrad(u_mean), iGrad(u_mean)))
#    tensor = (iInner(iT(iGrad(u_mean)) - (iGrad(u_mean)), iGrad(u_mean)))
#    #tensor = (iInner(iGrad(u_mean), iGrad(u_mean))) - iInner(iT(iGrad(u_mean)),iGrad(u_mean))
#
#    expression   = (tensor*iConj(itest)).ufl_tens*J_hat*r*r*dx
#    expression_T = ((t12-t21)*t12 + (t13-t31)*t13 + (t21-t12)*t21 + (t23-t32)*t23 + (t31-t13)*t31 + (t32-t23)*t32 )*conj(test)*r*r*dx
#    u1.evaluateUflTensorExpression(expression)
#    u2.evaluateUflExpression(-expression_T)
#    print("all6:", np.linalg.norm(u1.getCoefficientArray() - u2.getCoefficientArray()), np.linalg.norm(u1.getCoefficientArray()), np.linalg.norm(u2.getCoefficientArray()))
#
#   ##############################################################################################
#    print("")
#    print("#### power ####")
#
#    exponent = meanFlow.nulam
#
#    tensor = meanFlow.nulam ** exponent
#
#    expression   = (tensor*iConj(itest)).ufl_tens*J_hat*dx
#    u1.evaluateUflTensorExpression(expression)
#    meanFlow._fieldDict["nulam"].x.array[:] = meanFlow._fieldDict["nulam"].x.array[:]**exponent
#    expression_T = (meanFlow._fieldDict["nulam"]*conj(test)*dx)
#    u2.evaluateUflExpression(expression_T)
#    print("all: ", np.linalg.norm(u1.getCoefficientArray() - u2.getCoefficientArray()), np.linalg.norm(u1.getCoefficientArray()), np.linalg.norm(u2.getCoefficientArray()))



  
    
