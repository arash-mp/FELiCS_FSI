import pdb
import numpy as np
import copy
import time
import sys

from   ufl import (TrialFunctions,
                   TestFunctions,
                   )

from   FELiCS.Parameters.parameters import parameters

import FELiCS.IO.Import as Import
from   FELiCS.IO.ExportSolution import ExportGUI,ExportFromFile

import FELiCS.SpaceDisc.DefineFEMSpaces as DefineFEMSpaces
from   FELiCS.Fields.meanFlowClass import meanFlowClass
from   FELiCS.Fields.fluctuationClass import fluctuationSolutions
from   FELiCS.Equation.EquationCollection import EquationCollectionClass
from   FELiCS.Misc.functions import printDebug

from   FELiCS.Solvers.LinearSolver import LinearSolver 
from   FELiCS.Fields.ModeCollection import ModeCollection
from   FELiCS.Fields.Field import Field
from   FELiCS.Fields.Mode import Mode
from   FELiCS.Misc.tensorUtils import Tensor





#-----------------------------------------------------------------------
## READ PARAMETER FILE 
#-----------------------------------------------------------------------
# Replace spurious double quotation marks which occur in MobaXterm when executing shell scripts
for i_argument,argument in enumerate(sys.argv):
        sys.argv[i_argument] = argument.replace(chr(8221), '')

if len(sys.argv) == 1:
    useGUI= True
elif sys.argv[1] == '-file': useGUI=False
else: useGUI=True

SettingsFileName = sys.argv[2]
param=parameters()
param.importFromFile(SettingsFileName)
param.getOldParameters()



#-----------------------------------------------------------------------
## INITIALIZATION 
#-----------------------------------------------------------------------
# mesh
mesh=param.BCs.getMesh()

# FEMSpaces
FEMSpaces = DefineFEMSpaces.FEMSpacesClass(
            param,
            mesh,
            )


# initialize mean flow class
meanFlow = meanFlowClass(param, FEMSpaces, mesh)

#import mean flow data from file
meanFlow.importDataFromFile()

# export mean flow in "h5" file
if not param.FlowInput.MeanFlowFilePath.split('.')[-1] == 'hdf5':
    meanFlow.exportBaseFlowAsHDF5()
meanflowFilename = 'meanflow.h5'
meanFlow.mapToExportMeshAndExport(FEMSpaces, meanflowFilename)

# equation
equation = EquationCollectionClass(
                                  param,
                                  FEMSpaces,
                                  meanFlow,
                                  mesh
                                  )


#-----------------------------------------------------------------------
## CALCULATE ADJOINT BASEFLOW
#-----------------------------------------------------------------------
A  = equation.getLinearOperator(meanFlow)

# create and read fields
mode_direct  = Mode(FEMSpaces.VMixed, mesh) 
mode_adjoint = Mode(FEMSpaces.VMixed, mesh) 
rhs          = Field(FEMSpaces.VMixed, mesh) 

mode_direct.setCoefficientArray(np.load("mode_direct.npy"))
mode_adjoint.setCoefficientArray(np.load("mode_adjoint.npy"))
rhs.setCoefficientArray(np.load("rhs.npy"))

# solve the adjoint equation system
rhs.conjugate()
baseFlow_adjoint = Field(FEMSpaces.VMixed, mesh)
baseFlow_adjoint.setCoefficientArray(LinearSolver.solveTransposeEquationSystem(A, rhs.getPetscVector()))
baseFlow_adjoint.conjugate()


#-----------------------------------------------------------------------
## EXPORT 
#-----------------------------------------------------------------------
np.save("baseFlow_adjoint.npy",  baseFlow_adjoint.getCoefficientArray())





