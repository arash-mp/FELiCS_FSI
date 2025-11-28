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
from FELiCS.Fields.ModeCollection import ModeCollection

from tests.RandomCaseHandler import RandomCaseHandler
from tests.UnitTestHelper    import UnitTestHelper

#################################################################
##### define necessary functions ################################
#################################################################
class ModeCollectionTestHandler(RandomCaseHandler):
    def __init__(self, dim_vector):
        super().__init__(dim_vector)
        self.mode_collection = ModeCollection(self.space_scalar, self.felics_mesh,)
        
    def appendModesInCollection(self, analysisType, num_modes=5):
        for i in range(num_modes):
            if analysisType == "Modal":
                mode1 = self.createFELiCSMode(name=f"fake_mode_direct_{i}", analysisType=analysisType, modeType="Direct")
                mode1.isAdjoint = False
                mode1.eigenValue = np.random.randint(9999) + 1j*np.random.randint(9999)
                mode2 = self.createFELiCSMode(name=f"fake_mode_adjoint_{i}", analysisType=analysisType, modeType="Adjoint")
                mode2.function.x.array[:] = mode1.function.x.array[:].conj()
                mode2.isAdjoint = True
                mode2.eigenValue = mode1.eigenValue
                self.mode_collection.appendMode(mode1)
                self.mode_collection.appendMode(mode2)
            elif analysisType == "Resolvent" or analysisType == "Input_Output":
                mode1 = self.createFELiCSMode(name=f"fake_mode_forcing_{i}", analysisType=analysisType, modeType="Forcing")
                mode2 = self.createFELiCSMode(name=f"fake_mode_response{i}", analysisType=analysisType, modeType="Response")
                mode1.isAdjoint = False
                mode2.isAdjoint = False
                mode1.eigenValue = np.random.randint(9999) + 1j*np.random.randint(9999)
                mode2.eigenValue = np.random.randint(9999) + 1j*np.random.randint(9999)
                self.mode_collection.appendMode(mode1)
                self.mode_collection.appendMode(mode2)            
        
    
    

#################################################################
##### initialize tests ##########################################
#################################################################

# 0. define parameters
dim_vector           = 3 #dimension of vector function space, 2 or 3 (at the moment: is fixed to 3; that should test the 2, too?) 

# 1. create Test class
randomModeCollector = ModeCollectionTestHandler(dim_vector)
# --------------------------------------------------------------
# Tests for Field methods
# --------------------------------------------------------------
def test_getDirectEigenValueSpectrum():
    print("Testing get direct eigenvalue spectrum")
    # NOTE: This test only makes sense for Modal analysis
    analysisType = "Modal"
    randomModeCollector.appendModesInCollection(analysisType, num_modes=5)
    # 1. define validation 
    validation_eigenvalues = []
    for mode in randomModeCollector.mode_collection.modeList:
        if mode.isAdjoint == False:
            validation_eigenvalues.append(mode.eigenValue)
    # 2. compare eigenvalues from the method
    test_eigenvalues = randomModeCollector.mode_collection.getDirectEigenValueSpectrum()
    
    for i, eigenvalue in enumerate(test_eigenvalues):
        print(" - Testing eigenvalue ", i, ": ", eigenvalue, " vs ", validation_eigenvalues[i])
        assert np.abs(eigenvalue - validation_eigenvalues[i]) < 1e-14
    pass