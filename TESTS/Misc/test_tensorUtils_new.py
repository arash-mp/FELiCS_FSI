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
from random_fields_handle import RandomFieldsHandle
from unit_test_handle import UnitTestHandle

#################################################################
##### define necessary functions ################################
#################################################################

class TensorUtilsTestHandle(RandomFieldsHandle):
    def __init__(self,dim_vector):
        super().__init__(dim_vector)
        self.func_scalar1 = self.createDolfinxFunction(dim=1)
        self.func_scalar2 = self.createDolfinxFunction(dim=1)
        self.func_vector1 = self.createDolfinxFunction(dim=dim_vector)
        self.func_vector2 = self.createDolfinxFunction(dim=dim_vector)
        pass
    
    def updateTensorCoordinate(self, coordinateSystemName, m):
        # 1. create CoordinateSystem
        testCoordinateSystem = CoordinateSystem(
                             SpatialCoordinate(self.mesh),
                             coordinateSystemName,
                             m=m,
                             gdim = 2,
                             trueDim = 3
                             )
        self.J_hat = testCoordinateSystem.J_hat
        if coordinateSystemName == "cartesian":
            self.r = 1.0
        elif coordinateSystemName == "cylindricalfelics":
            self.r   = SpatialCoordinate(self.mesh)[1]
        # 2. create TensorUtils for test functions
        self.test_scalar_T = Tensor(self.test_scalar, testCoordinateSystem,  mayHaveSpectralDimension =True) 
        self.test_vector_T = Tensor(self.test_vector, testCoordinateSystem,  mayHaveSpectralDimension =True) 
        # 3. create TensorUtils for random functions
        self.tens_vector1 = Tensor(self.func_vector1, testCoordinateSystem, mayHaveSpectralDimension = True)
        self.tens_vector2 = Tensor(self.func_vector2, testCoordinateSystem, mayHaveSpectralDimension = True)
        self.tens_scalar1 = Tensor(self.func_scalar1, testCoordinateSystem, mayHaveSpectralDimension = True)
        self.tens_scalar2 = Tensor(self.func_scalar2, testCoordinateSystem, mayHaveSpectralDimension = True)
        # 4. create spatial gradients
        self.grad_vector1 = UnitTestHandle.getVecGrad(self.func_vector1, coordinateSystemName, m, self.r)
        self.grad_vector2 = UnitTestHandle.getVecGrad(self.func_vector2, coordinateSystemName, m, self.r)
        self.grad_scalar1 = UnitTestHandle.getScalarGrad(self.func_scalar1, coordinateSystemName, m, self.r)
        self.grad_scalar2 = UnitTestHandle.getScalarGrad(self.func_scalar2, coordinateSystemName, m, self.r)
        pass
    
    def checkExpressionInAllCoordinateSystems(self, expressionFunction):
        coordinateSystemList = ["cartesian", "cylindricalfelics"]
        m_list = [0,np.random.randint(20)]
        for coordinateSystemName in coordinateSystemList:
            for m in m_list:
                print(f"\n\n Checking alignment for coordinate system {coordinateSystemName} and m={m}")
                self.updateTensorCoordinate(coordinateSystemName, m)
                tensor_expr, valid_expr = expressionFunction(self)
                UnitTestHandle.checkVectorExpressionAlignment(tensor_expr, valid_expr)
        pass

#################################################################
##### initialize tests ##########################################
#################################################################

# 0. define parameters
dim_vector           = 3 #dimension of vector function space, 2 or 3 (at the moment: is fixed to 3; that should test the 2, too?) 

# 1. create Test class
randomField = TensorUtilsTestHandle(dim_vector)

# --------------------------------------------------------------
# Tests for tensor algebra
# --------------------------------------------------------------
def test_scalar_minus_float():
    print("Testing scalar minus float")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = ((fd_field.tens_scalar1 - fd_field.random_float) * iConj(fd_field.test_scalar_T)).ufl_tens* fd_field.J_hat*dx
        valid_expr  = (fd_field.func_scalar1 - fd_field.random_float) * conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    print("... passed")
    
def test_float_minus_scalar():
    print("Testing float minus scalar")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = ((-1 * fd_field.tens_scalar1 + fd_field.random_float) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (-1 * fd_field.func_scalar1 + fd_field.random_float) * conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    print("... passed")
    

def test_iDiv():
    print("Testing iDiv")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iDiv(fd_field.tens_vector1) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.grad_vector1[0] + fd_field.grad_vector1[4] + fd_field.grad_vector1[8])*conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    print("... passed")