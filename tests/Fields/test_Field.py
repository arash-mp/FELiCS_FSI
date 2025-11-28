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



            
