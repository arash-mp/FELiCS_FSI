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
# - do all other functions (only iGrad here) => Xiuyang and Kai?
# - implement logging and replace print statements => ? 
# - create bigger test architecture, at best object oriented, to conduct all tests most efficiently => Sophie

import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from FELiCS.Misc.tensorUtils import*

import numpy as np
from mpi4py              import MPI
from petsc4py.PETSc      import ScalarType
from ufl                 import (
    dx,TestFunction,TrialFunction,conj,SpatialCoordinate, Dx,dot, inner, grad, outer
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
from tests.RandomCaseHandler import RandomCaseHandler
from tests.UnitTestHelper    import UnitTestHelper
from FELiCS.Misc.logging import Logger

# initialize the logger
logger = Logger(logger_name="felics_unit_test")
logger = Logger.get_logger("felics_unit_test")

#################################################################
##### define necessary functions ################################
#################################################################

class TensorUtilsTestHandle(RandomCaseHandler):
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
        self.m = m
        if coordinateSystemName == "cartesian":
            self.r = 1.0
            self.coordinateSystemName = "cartesian"
        elif coordinateSystemName == "cylindricalfelics":
            self.r   = SpatialCoordinate(self.mesh)[1]
            self.coordinateSystemName = "cylindricalfelics"
        # 2. create TensorUtils for test functions
        self.test_scalar_T = Tensor(self.test_scalar, testCoordinateSystem,  mayHaveSpectralDimension =True) 
        self.test_vector_T = Tensor(self.test_vector, testCoordinateSystem,  mayHaveSpectralDimension =True) 
        # 3. create TensorUtils for random functions
        self.tens_vector1 = Tensor(self.func_vector1, testCoordinateSystem, mayHaveSpectralDimension = True)
        self.tens_vector2 = Tensor(self.func_vector2, testCoordinateSystem, mayHaveSpectralDimension = True)
        self.tens_scalar1 = Tensor(self.func_scalar1, testCoordinateSystem, mayHaveSpectralDimension = True)
        self.tens_scalar2 = Tensor(self.func_scalar2, testCoordinateSystem, mayHaveSpectralDimension = True)
        # 4. create spatial gradients
        self.grad_vector1 = UnitTestHelper.getVecGrad(self.func_vector1, coordinateSystemName, m, self.r)
        self.grad_vector2 = UnitTestHelper.getVecGrad(self.func_vector2, coordinateSystemName, m, self.r)
        self.grad_scalar1 = UnitTestHelper.getScalarGrad(self.func_scalar1, coordinateSystemName, m, self.r)
        self.grad_scalar2 = UnitTestHelper.getScalarGrad(self.func_scalar2, coordinateSystemName, m, self.r)
        pass
    
    def checkExpressionInAllCoordinateSystems(self, expressionFunction, tol=1.e-14):
        coordinateSystemList = ["cartesian", "cylindricalfelics"]
        m_list = [0, np.random.randint(20)+1]
        for coordinateSystemName in coordinateSystemList:
            for m in m_list:
                logger.info(f"Checking alignment for coordinate system {coordinateSystemName} and m={m}")
                self.updateTensorCoordinate(coordinateSystemName, m)
                tensor_expr, valid_expr = expressionFunction(self)
                UnitTestHelper.checkVectorExpressionAlignment(tensor_expr, valid_expr, tol)
                logger.info("... passed")
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
    logger.info("Testing scalar minus float")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = ((fd_field.tens_scalar1 - fd_field.random_float) * iConj(fd_field.test_scalar_T)).ufl_tens* fd_field.J_hat*dx
        valid_expr  = (fd_field.func_scalar1 - fd_field.random_float) * conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_float_minus_scalar():
    logger.info("Testing float minus scalar")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = ((-1 * fd_field.tens_scalar1 + fd_field.random_float) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (-1 * fd_field.func_scalar1 + fd_field.random_float) * conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_scalar_minus_scalar():
    logger.info("Testing scalar minus scalar")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = ((fd_field.tens_scalar1 - fd_field.tens_scalar2) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.func_scalar1 - fd_field.func_scalar2) * conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_scalar_plus_scalar():
    logger.info("Testing scalar plus scalar")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = ((fd_field.tens_scalar1 + fd_field.tens_scalar2) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.func_scalar1 + fd_field.func_scalar2) * conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    

def test_float_times_scalar():
    logger.info("Testing float times scalar")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = ((fd_field.random_float * fd_field.tens_scalar1) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.random_float * fd_field.func_scalar1) * conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_scalar_times_float():
    logger.info("Testing scalar times float.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = ((fd_field.tens_scalar1 * fd_field.random_float) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.func_scalar1 * fd_field.random_float) * conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_scalar_times_scalar():
    logger.info("Testing scalar times scalar.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = ((fd_field.tens_scalar1 * fd_field.tens_scalar2) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.func_scalar1 * fd_field.func_scalar2) * conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_scalar_divide_scalar():
    logger.info("Testing scalar divide scalar.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = ((fd_field.tens_scalar1 / fd_field.tens_scalar2) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.func_scalar1 / fd_field.func_scalar2) * conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_scalar_exp_to_positive_float():
    logger.info("Testing scalar exponent to positive float.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = ((fd_field.tens_scalar1 ** fd_field.random_float) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.func_scalar1 ** fd_field.random_float) * conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_scalar_exp_to_negetive_float():
    logger.info("Testing scalar exponent to negetive float.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = ((fd_field.tens_scalar1 ** (-1 * fd_field.random_float)) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.func_scalar1 ** (-1 * fd_field.random_float)) * conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_vector_dot_vector():
    logger.info("Testing vector dot vector.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iDot(fd_field.tens_vector1, fd_field.tens_vector2) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (dot(fd_field.func_vector1, fd_field.func_vector2)) * conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_scalar_times_vector():
    logger.info("Testing scalar times vector.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iDot(fd_field.tens_scalar1 * fd_field.tens_vector1, iConj(fd_field.test_vector_T))).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (dot(fd_field.func_scalar1 * fd_field.func_vector1, conj(fd_field.test_vector)))*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_vector_times_scalar():
    logger.info("Testing vector times scalar.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iDot(fd_field.tens_vector1 * fd_field.tens_scalar1, iConj(fd_field.test_vector_T))).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (dot(fd_field.func_vector1 * fd_field.func_scalar1, conj(fd_field.test_vector)))*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_float_divide_scalar_times_vector():
    logger.info("Testing one divide by scalar times vector.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iDot((fd_field.random_float) / fd_field.tens_scalar1 * fd_field.tens_vector1, iConj(fd_field.test_vector_T))).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (dot((fd_field.random_float) / fd_field.func_scalar1 * fd_field.func_vector1, conj(fd_field.test_vector)))*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_vector_divided_by_scalar():
    logger.info("Testing vector divided by scalar.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iDot(fd_field.tens_vector1 / fd_field.tens_scalar1 , iConj(fd_field.test_vector_T))).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (dot(fd_field.func_vector1 / fd_field.func_scalar1 , conj(fd_field.test_vector)))*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_vector_dot_dyade():
    logger.info("Testing vector dot dyade.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iDot(iDot(fd_field.tens_vector1, iGrad(fd_field.tens_vector2)) , iConj(fd_field.test_vector_T))).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.func_vector1[0] * fd_field.grad_vector2[0][0] + fd_field.func_vector1[1] * fd_field.grad_vector2[1][0] + fd_field.func_vector1[2] * fd_field.grad_vector2[2][0])*conj(fd_field.test_vector[0])*fd_field.r*dx
        valid_expr += (fd_field.func_vector1[0] * fd_field.grad_vector2[0][1] + fd_field.func_vector1[1] * fd_field.grad_vector2[1][1] + fd_field.func_vector1[2] * fd_field.grad_vector2[2][1])*conj(fd_field.test_vector[1])*fd_field.r*dx
        if not (fd_field.coordinateSystemName == "cartesian" and fd_field.m == 0):
            valid_expr += (fd_field.func_vector1[0] * fd_field.grad_vector2[0][2] + fd_field.func_vector1[1] * fd_field.grad_vector2[1][2] + fd_field.func_vector1[2] * fd_field.grad_vector2[2][2])*conj(fd_field.test_vector[2])*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_dyade_dot_vector():
    logger.info("Testing dyade dot vector.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iDot(iDot(iGrad(fd_field.tens_vector1), fd_field.tens_vector2) , iConj(fd_field.test_vector_T))).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.grad_vector1[0][0] * fd_field.func_vector2[0]+ fd_field.grad_vector1[0][1] * fd_field.func_vector2[1] + fd_field.grad_vector1[0][2] * fd_field.func_vector2[2])*conj(fd_field.test_vector[0])*fd_field.r*dx
        valid_expr += (fd_field.grad_vector1[1][0] * fd_field.func_vector2[0]+ fd_field.grad_vector1[1][1] * fd_field.func_vector2[1] + fd_field.grad_vector1[1][2] * fd_field.func_vector2[2])*conj(fd_field.test_vector[1])*fd_field.r*dx
        valid_expr += (fd_field.grad_vector1[2][0] * fd_field.func_vector2[0]+ fd_field.grad_vector1[2][1] * fd_field.func_vector2[1] + fd_field.grad_vector1[2][2] * fd_field.func_vector2[2])*conj(fd_field.test_vector[2])*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_vector_dot_dyade_dot_vector():
    logger.info("Testing vector dot dyade dot vector.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = ((iDot(iDot(fd_field.tens_vector1,iGrad(fd_field.tens_vector1)), fd_field.tens_vector2)) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.func_vector1[0] * fd_field.grad_vector1[0][0] + fd_field.func_vector1[1] * fd_field.grad_vector1[1][0] + fd_field.func_vector1[2] * fd_field.grad_vector1[2][0]) * fd_field.func_vector2[0] * conj(fd_field.test_scalar)*fd_field.r*dx
        valid_expr += (fd_field.func_vector1[0] * fd_field.grad_vector1[0][1] + fd_field.func_vector1[1] * fd_field.grad_vector1[1][1] + fd_field.func_vector1[2] * fd_field.grad_vector1[2][1]) * fd_field.func_vector2[1] * conj(fd_field.test_scalar)*fd_field.r*dx
        if not (fd_field.coordinateSystemName == "cartesian" and fd_field.m == 0):
            valid_expr += (fd_field.func_vector1[0] * fd_field.grad_vector1[0][2] + fd_field.func_vector1[1] * fd_field.grad_vector1[1][2] + fd_field.func_vector1[2] * fd_field.grad_vector1[2][2]) * fd_field.func_vector2[2] * conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_dyade_dot_dyade():
    logger.info("Testing dyade dot dyade.")
    # 1. define expressions
    def expression(fd_field):
        [[t11_1, t12_1, t13_1], [t21_1, t22_1, t23_1], [t31_1, t32_1, t33_1]] = fd_field.grad_vector1
        [[t11_2, t12_2, t13_2], [t21_2, t22_2, t23_2], [t31_2, t32_2, t33_2]] = fd_field.grad_vector2
        tensor_expr = (iDot(iDot(iDot(iGrad(fd_field.tens_vector1), iGrad(fd_field.tens_vector2)),fd_field.tens_vector2 ), iConj(fd_field.test_vector_T))).ufl_tens*fd_field.J_hat*dx
        valid_expr  = ((t11_1 * t11_2 + t12_1 * t21_2 + t13_1 * t31_2) * fd_field.func_vector2[0]+ \
            (t11_1 * t12_2 + t12_1 * t22_2 + t13_1 * t32_2) * fd_field.func_vector2[1] + \
            (t11_1 * t13_2 + t12_1 * t23_2 + t13_1 * t33_2) * fd_field.func_vector2[2])*conj(fd_field.test_vector[0])*fd_field.r*dx
        valid_expr += ((t21_1 * t11_2 + t22_1 * t21_2 + t23_1 * t31_2) * fd_field.func_vector2[0]+ \
            (t21_1 * t12_2 + t22_1 * t22_2 + t23_1 * t32_2) * fd_field.func_vector2[1] + \
            (t21_1 * t13_2 + t22_1 * t23_2 + t23_1 * t33_2) * fd_field.func_vector2[2])*conj(fd_field.test_vector[1])*fd_field.r*dx
        valid_expr += ((t31_1 * t11_2 + t32_1 * t21_2 + t33_1 * t31_2) * fd_field.func_vector2[0]+ \
            (t31_1 * t12_2 + t32_1 * t22_2 + t33_1 * t32_2) * fd_field.func_vector2[1] + \
            (t31_1 * t13_2 + t32_1 * t23_2 + t33_1 * t33_2) * fd_field.func_vector2[2])*conj(fd_field.test_vector[2])*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_dyade_inner_dyade():
    logger.info("Testing dyade inner dyade.")
    # 1. define expressions
    def expression(fd_field):
        [[t11_1, t12_1, t13_1], [t21_1, t22_1, t23_1], [t31_1, t32_1, t33_1]] = fd_field.grad_vector1
        [[t11_2, t12_2, t13_2], [t21_2, t22_2, t23_2], [t31_2, t32_2, t33_2]] = fd_field.grad_vector2
        tensor_expr = (iInner(iGrad(fd_field.tens_vector1),iGrad(fd_field.tens_vector2)) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr = ((t11_1 * conj(t11_2) + t21_1 * conj(t21_2) + t31_1 * conj(t31_2)) + (t12_1 * conj(t12_2) + t22_1 * conj(t22_2) + t32_1 * conj(t32_2)) + (t13_1 * conj(t13_2) + t23_1 * conj(t23_2) + t33_1 * conj(t33_2))) * conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_dyade_transpose():
    logger.info("Testing dyade transpose.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iDot(iDot(iT(iGrad(fd_field.tens_vector1)),fd_field.tens_vector2), iConj(fd_field.test_vector_T))).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.grad_vector1[0][0] * fd_field.func_vector2[0]+ fd_field.grad_vector1[1][0] * fd_field.func_vector2[1] + fd_field.grad_vector1[2][0] * fd_field.func_vector2[2])*conj(fd_field.test_vector[0])*fd_field.r*dx
        valid_expr += (fd_field.grad_vector1[0][1] * fd_field.func_vector2[0]+ fd_field.grad_vector1[1][1] * fd_field.func_vector2[1] + fd_field.grad_vector1[2][1] * fd_field.func_vector2[2])*conj(fd_field.test_vector[1])*fd_field.r*dx
        if not (fd_field.coordinateSystemName == "cartesian" and fd_field.m == 0):
            valid_expr += (fd_field.grad_vector1[0][2] * fd_field.func_vector2[0]+ fd_field.grad_vector1[1][2] * fd_field.func_vector2[1] + fd_field.grad_vector1[2][2] * fd_field.func_vector2[2])*conj(fd_field.test_vector[2])*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_dyade_trace():
    logger.info("Testing dyade trace.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iTr(iGrad(fd_field.tens_vector1)) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.grad_vector1[0][0] + fd_field.grad_vector1[1][1] + fd_field.grad_vector1[2][2])*conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_scalar_conj():
    logger.info("Testing conjugate of scalar.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iConj(fd_field.tens_scalar1) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (conj(fd_field.func_scalar1))*conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_vector_conj():
    logger.info("Testing conjugate of vector.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iDot(iConj(fd_field.tens_vector1), iConj(fd_field.test_vector_T))).ufl_tens*fd_field.J_hat*dx
        valid_expr  = dot(conj(fd_field.func_vector1),conj(fd_field.test_vector))*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_dyade_hermitian():
    logger.info("Testing hermitian transpose of dyade.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iDot(iDot(iT(iConj(iGrad(fd_field.tens_vector1))), fd_field.tens_vector2), iConj(fd_field.test_vector_T))).ufl_tens*fd_field.J_hat*dx
        valid_expr = (conj(fd_field.grad_vector1[0][0]) * fd_field.func_vector2[0]+ conj(fd_field.grad_vector1[1][0]) * fd_field.func_vector2[1] + conj(fd_field.grad_vector1[2][0]) * fd_field.func_vector2[2])*conj(fd_field.test_vector[0])*fd_field.r*dx
        valid_expr += (conj(fd_field.grad_vector1[0][1]) * fd_field.func_vector2[0]+ conj(fd_field.grad_vector1[1][1]) * fd_field.func_vector2[1] + conj(fd_field.grad_vector1[2][1]) * fd_field.func_vector2[2])*conj(fd_field.test_vector[1])*fd_field.r*dx
        if not (fd_field.coordinateSystemName == "cartesian" and fd_field.m == 0):
            valid_expr += (conj(fd_field.grad_vector1[0][2]) * fd_field.func_vector2[0]+ conj(fd_field.grad_vector1[1][2]) * fd_field.func_vector2[1] + conj(fd_field.grad_vector1[2][2]) * fd_field.func_vector2[2])*conj(fd_field.test_vector[2])*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_identity():
    logger.info("Testing identity functionality.")
    # 1. define expressions
    def expression(fd_field):
        test_tensor  = iGrad(fd_field.tens_vector1)
        valid_tensor = Identity(3)
        test_tensor  = iIdentity(test_tensor)
        
        tensor_expr  = (iDot(iDot(test_tensor, iConj(fd_field.test_vector_T)),fd_field.tens_vector1)).ufl_tens*fd_field.J_hat*dx
        valid_expr = dot(dot(valid_tensor, conj(fd_field.test_vector)),fd_field.func_vector1)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_grad_scalar_dot_grad_scalar():
    logger.info("Testing grad of scalar dot grad of scalar.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iDot(iGrad(fd_field.tens_scalar1), iGrad(fd_field.tens_scalar2)) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.grad_scalar1[0]*fd_field.grad_scalar2[0] + fd_field.grad_scalar1[1]*fd_field.grad_scalar2[1] + fd_field.grad_scalar1[2]*fd_field.grad_scalar2[2]) * conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_scalar_times_grad_scalar():
    logger.info("Testing scalar times grad of scalar.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iDot(fd_field.tens_scalar1 * iGrad(fd_field.tens_scalar2) , iConj(fd_field.test_vector_T))).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.func_scalar1 * fd_field.grad_scalar2[0]) * conj(fd_field.test_vector[0])*fd_field.r*dx
        valid_expr  += (fd_field.func_scalar1 * fd_field.grad_scalar2[1]) * conj(fd_field.test_vector[1])*fd_field.r*dx
        if not (fd_field.m == 0):
            valid_expr  += (fd_field.func_scalar1 * fd_field.grad_scalar2[2]) * conj(fd_field.test_vector[2])*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_grad_scalar_times_scalar():
    logger.info("Testing grad of scalar times scalar.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iDot(iGrad(fd_field.tens_scalar1) * fd_field.tens_scalar2 , iConj(fd_field.test_vector_T))).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.grad_scalar1[0] * fd_field.func_scalar2) * conj(fd_field.test_vector[0])*fd_field.r*dx
        valid_expr  += (fd_field.grad_scalar1[1] * fd_field.func_scalar2) * conj(fd_field.test_vector[1])*fd_field.r*dx
        if not (fd_field.m == 0):
            valid_expr  += (fd_field.grad_scalar1[2] * fd_field.func_scalar2) * conj(fd_field.test_vector[2])*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_float_divide_scalar_times_grad_scalar():
    logger.info("Testing one divide scalar times grad of scalar.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iDot((fd_field.random_float) / fd_field.tens_scalar1 * iGrad(fd_field.tens_scalar2), iConj(fd_field.test_vector_T))).ufl_tens*fd_field.J_hat*dx
        valid_expr  = ((fd_field.random_float) / fd_field.func_scalar1 * fd_field.grad_scalar2[0]) * conj(fd_field.test_vector[0])*fd_field.r*dx
        valid_expr  += ((fd_field.random_float) / fd_field.func_scalar1 * fd_field.grad_scalar2[1]) * conj(fd_field.test_vector[1])*fd_field.r*dx
        if not (fd_field.m == 0):
            valid_expr  += ((fd_field.random_float) / fd_field.func_scalar1 * fd_field.grad_scalar2[2]) * conj(fd_field.test_vector[2])*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_grad_scalar_divided_by_scalar():
    logger.info("Testing grad of scalar divided by scalar.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iDot(iGrad(fd_field.tens_scalar1) / fd_field.tens_scalar2 , iConj(fd_field.test_vector_T))).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.grad_scalar1[0] / fd_field.func_scalar2) * conj(fd_field.test_vector[0])*fd_field.r*dx
        valid_expr  += (fd_field.grad_scalar1[1] / fd_field.func_scalar2) * conj(fd_field.test_vector[1])*fd_field.r*dx
        if not (fd_field.m == 0):
            valid_expr  += (fd_field.grad_scalar1[2] / fd_field.func_scalar2) * conj(fd_field.test_vector[2])*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_grad_scalar_dot_dyade():
    logger.info("Testing grad of scalar dot dyade")
    # 1. define expressions
    def expression(fd_field):
        [[t11, t12, t13], [t21, t22, t23], [t31, t32, t33]] = fd_field.grad_vector1
        t1, t2, t3 = fd_field.grad_scalar1
        tensor_expr = (iDot(iDot(iGrad(fd_field.tens_scalar1), iGrad(fd_field.tens_vector1)) , iConj(fd_field.test_vector_T))).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (t11 * t1 + t21 * t2 + t31 * t3)*conj(fd_field.test_vector[0])*fd_field.r*dx
        valid_expr += (t12 * t1 + t22 * t2 + t32 * t3)*conj(fd_field.test_vector[1])*fd_field.r*dx
        if not (fd_field.coordinateSystemName == "cartesian" and fd_field.m == 0):
            valid_expr += (t13 * t1 + t23 * t2 + t33 * t3)*conj(fd_field.test_vector[2])*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_dyade_dot_grad_scalar():
    logger.info("Testing dyade dot grad of scalar.")
    # 1. define expressions
    def expression(fd_field):
        [[t11, t12, t13], [t21, t22, t23], [t31, t32, t33]] = fd_field.grad_vector1
        t1, t2, t3 = fd_field.grad_scalar2
        tensor_expr = (iDot(iDot(iGrad(fd_field.tens_vector1), iGrad(fd_field.tens_scalar2)) , iConj(fd_field.test_vector_T))).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (t11 * t1 + t12 * t2 + t13 * t3)*conj(fd_field.test_vector[0])*fd_field.r*dx
        valid_expr += (t21 * t1 + t22 * t2 + t23 * t3)*conj(fd_field.test_vector[1])*fd_field.r*dx
        valid_expr += (t31 * t1+ t32 * t2 + t33 * t3)*conj(fd_field.test_vector[2])*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_grad_scalar_dot_dyade_dot_grad_scalar():
    logger.info("Testing grad of scalar dot dyade dot grad of scalar.")
    # 1. define expressions
    def expression(fd_field):
        [[t11, t12, t13], [t21, t22, t23], [t31, t32, t33]] = fd_field.grad_vector1
        t1_1, t2_1, t3_1 = fd_field.grad_scalar1
        t1_2, t2_2, t3_2 = fd_field.grad_scalar2
        tensor_expr = ((iDot(iDot(iGrad(fd_field.tens_scalar1),iGrad(fd_field.tens_vector1)), iGrad(fd_field.tens_scalar2))) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (t1_1 * t11 + t2_1 * t21 + t3_1 * t31) * t1_2 * conj(fd_field.test_scalar)*fd_field.r*dx
        valid_expr += (t1_1 * t12 + t2_1 * t22 + t3_1 * t32) * t2_2 * conj(fd_field.test_scalar)*fd_field.r*dx
        # if not (coordinateSystemName == "cartesian" and m == 0):
        if not (fd_field.m == 0):
            logger.info("Including z-component in validation.")
            valid_expr += (t1_1 * t13 + t2_1 * t23 + t3_1* t33) * t3_2 * conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    # randomField.checkExpressionInAllCoordinateSystems(expression, tol=1e-10)
    randomField.checkExpressionInAllCoordinateSystems(expression)

    # NOTE: This case needs a higher tolerance, possibly due to the multiple gradients involved and the functions are not smoothed.
    
    
def test_grad_scalar_conj():
    logger.info("Testing conjugate of grad of scalar.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iDot(iConj(iGrad(fd_field.tens_scalar1)), iConj(fd_field.test_vector_T))).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (conj(fd_field.grad_scalar1[0]) * conj(fd_field.test_vector[0]))*fd_field.r*dx
        valid_expr  += (conj(fd_field.grad_scalar1[1]) * conj(fd_field.test_vector[1]))*fd_field.r*dx
        if not (fd_field.m == 0):
            valid_expr  += (conj(fd_field.grad_scalar1[2]) * conj(fd_field.test_vector[2]))*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
    
def test_div_vector():
    logger.info("Testing divergence of vector.")
    # 1. define expressions
    def expression(fd_field):
        tensor_expr = (iDiv(fd_field.tens_vector1) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (fd_field.grad_vector1[0][0] + fd_field.grad_vector1[1][1] + fd_field.grad_vector1[2][2])*conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
def test_vector_outer_vector():
    logger.info("Testing outer product of vectors.")
    # 1. define expressions
    def expression(fd_field):
        t1_1, t2_1, t3_1 = fd_field.func_vector1
        t1_2, t2_2, t3_2 = fd_field.func_vector2
        tensor_expr = (iInner(iOuter(fd_field.tens_vector1, fd_field.tens_vector2), iOuter(fd_field.tens_vector1, fd_field.tens_vector2)) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (inner(outer(fd_field.func_vector1, fd_field.func_vector2), outer(fd_field.func_vector1, fd_field.func_vector2))) * conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    randomField.checkExpressionInAllCoordinateSystems(expression)
    
def test_grad_scalar_outer_grad_scalar():
    logger.info("Testing outer product of grad of two scalars.")
    # 1. define expressions
    def expression(fd_field):
        t1_1, t2_1, t3_1 = fd_field.grad_scalar1
        t1_2, t2_2, t3_2 = fd_field.grad_scalar2
        tensor_expr = (iInner(iOuter(iGrad(fd_field.tens_scalar1), iGrad(fd_field.tens_scalar2)), iOuter(iGrad(fd_field.tens_scalar1), iGrad(fd_field.tens_scalar2))) * iConj(fd_field.test_scalar_T)).ufl_tens*fd_field.J_hat*dx
        valid_expr  = (t1_1 * conj(t1_1) * t1_2 * conj(t1_2) + \
                        t1_1 * conj(t1_1) * t2_2 * conj(t2_2) + \
                        t1_1 * conj(t1_1) * t3_2 * conj(t3_2) + \
                        t2_1 * conj(t2_1) * t1_2 * conj(t1_2) + \
                        t2_1 * conj(t2_1) * t2_2 * conj(t2_2) + \
                        t2_1 * conj(t2_1) * t3_2 * conj(t3_2) + \
                        t3_1 * conj(t3_1) * t1_2 * conj(t1_2) + \
                        t3_1 * conj(t3_1) * t2_2 * conj(t2_2) + \
                        t3_1 * conj(t3_1) * t3_2 * conj(t3_2) \
                        )*conj(fd_field.test_scalar)*fd_field.r*dx
        return tensor_expr, valid_expr
    # 2. check alignment
    # randomField.checkExpressionInAllCoordinateSystems(expression, tol=1e-11)
    randomField.checkExpressionInAllCoordinateSystems(expression)


#################################################################
##### moving log files to TESTS folder ##########################
#################################################################
Logger.change_log_location("./TESTS/Misc/logs")