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

import io
import contextlib
from FELiCS.Misc.tensorUtils import *
from FELiCS.Equation.UflDecorator import UflDecorator

import numpy as np
from mpi4py              import MPI
from petsc4py.PETSc      import ScalarType
from ufl                 import (
    dx,TestFunction,TrialFunction,conj,SpatialCoordinate, Dx,dot, inner, grad
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
class UflDecoratorTestHandler(RandomCaseHandler):
    def __init__(self, dim_vector):
        super().__init__(dim_vector)
        # NOTE: So far only test expression using scalar field
        self.scalar_field = self.createDolfinxFunction(dim=1)
        pass
    
    

#################################################################
##### initialize tests ##########################################
#################################################################

# 0. define parameters
dim_vector           = 3 #dimension of vector function space, 2 or 3 (at the moment: is fixed to 3; that should test the 2, too?) 

# 1. create Test class
randomField = UflDecoratorTestHandler(dim_vector)

# 2. create Ufl expression used for testing
# NOTE: Test with Poisson problem
lhs_expr = inner(grad(randomField.trial_scalar), grad(randomField.test_scalar)) * dx
rhs_expr = randomField.scalar_field * randomField.test_scalar * dx
full_expr = lhs_expr - rhs_expr
# --------------------------------------------------------------
# Tests for UflDocorator methods
# --------------------------------------------------------------
def test_printExpression():
    for expr in [None, lhs_expr, rhs_expr, full_expr]:
        # 1. define validation (no need here)
        
        # 2. check alignment
        # create UflDecorator object with the expression
        UflDeco = UflDecorator(expr)
        # get output string
        string_buffer = io.StringIO()
        with contextlib.redirect_stdout(string_buffer):
            UflDeco.printExpression()
        testString = string_buffer.getvalue()
        print(testString)
        
        if expr is None:
            assert "Ufl expression is zero." in testString
        else:
            print(f"validating expression: {expr}")
            assert f"{expr}" in testString
            print(f"validating arguments: {expr.arguments()}")
            assert f"{expr.arguments()}" in testString
            print(f"validating number of arguments: {len(expr.arguments())}")
            assert f"{len(expr.arguments())}" in testString
        
    print("... passed.")

def test_lhsIsZero():
    print("Validating if the left hand side of the ufl expression equals zero")
    from ufl import lhs, rhs
    for expr in [None, lhs_expr, rhs_expr, full_expr]:
        UflDeco = UflDecorator(expr)
        lhs_is_zero_test    = UflDeco.lhsIsZero()
        if UflDeco.isZero():
            lhs_is_zero_valid = True
        else: 
            temp_expression     = lhs(UflDeco._expression)
            if len(temp_expression.arguments())<2:
                lhs_is_zero_valid = True
            else:
                lhs_is_zero_valid = False
        print(f"Validating the expression: {expr}")
        assert lhs_is_zero_test == lhs_is_zero_valid
    print("... passed.")

def test_rhsIsZero():
    print("Validating if the right hand side of the ufl expression equals zero")
    from ufl import lhs, rhs
    for expr in [None, lhs_expr, rhs_expr, full_expr]:
        UflDeco = UflDecorator(expr)
        rhs_is_zero_test    = UflDeco.rhsIsZero()
        if UflDeco.isZero():
            rhs_is_zero_valid = True
        else: 
            temp_expression     = rhs(UflDeco._expression)
            if len(temp_expression.arguments())<1:
                rhs_is_zero_valid = True
            else:
                rhs_is_zero_valid = False
        print(f"Validating the expression: {expr}")
        assert rhs_is_zero_test == rhs_is_zero_valid
    print("... passed.")
