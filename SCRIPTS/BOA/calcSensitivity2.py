import pdb
import copy
import time
import sys
import numpy as np


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


from   GeometryDeformer import CylinderBSpline2Pts 




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
## CALCULATE GRADIENT WITH RESPECT TO GEOMETRY DEFORMATIONS 
#-----------------------------------------------------------------------
N_0  = equation.getNonlinearExpression(meanFlow)


epsilon = 1.e-4 #epsilon for finite difference mesh deformation


# create and read fields
a_i  = np.load("params.npy")
baseFlow_adjoint  = Field(FEMSpaces.VMixed, mesh) 
baseFlow_adjoint.setCoefficientArray(np.load("baseFlow_adjoint.npy"))

factor = np.load("factor.npy")

# deform mesh
epsilon = 1.e-6

geometryDeformer = CylinderBSpline2Pts(mesh, equation.boundaries)
N = geometryDeformer.getNumberOfParameters()
sensitivity2 = np.zeros(N,dtype=complex)
for i in range(0,N):
    a_i[i] = a_i[i] + epsilon

    geometryDeformer.deformMesh(a_i)
    
    N_deformed  = equation.getNonlinearExpression(meanFlow)
    N_deformed.axpy( -1., N_0)
    N_deformed.scale(1./epsilon)
    baseFlow_adjoint_petsc  = baseFlow_adjoint.getPetscVector()
    sensitivity2[i] = N_deformed.dot(baseFlow_adjoint_petsc) 

    a_i[i] = a_i[i] - epsilon
    geometryDeformer.restoreMesh()


sensitivity2 = sensitivity2 / factor

printDebug(True, '------------------------------------------------ ')
printDebug(True, '-- sensitivities part 2:  ' + str((sensitivity2)))
printDebug(True, '------------------------------------------------ ')


np.save("sensitivity2.npy", sensitivity2)


