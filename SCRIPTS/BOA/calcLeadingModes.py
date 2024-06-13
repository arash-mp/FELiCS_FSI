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
## SOLVE EIGENPROBLEM AND CALCULATE RHS OF ADJOINT BASEFLOW EQUATION
#-----------------------------------------------------------------------
# get matrices for eigenproblem
A  = equation.getLinearOperator(meanFlow)
B  = equation.getWeightMatrix  (meanFlow)

# get parameters for eigenproblem
guesses  = param.Numerics.EigenValueGuess
nSol     = param.Numerics.nSolut
adjoint  = param.Case.CalculateAdjoint

# track time
start= time.time()

# solve direct eigenproblem for each guess
solution = ModeCollection(FEMSpaces.VMixed, mesh)
for guess in guesses:
    tmp     = LinearSolver.solveGeneralEigenproblem(A,
                                                    B,
                                                    guess,
                                                    nSol,
                                                    )
    solution.appendSolutionOfEigenProblem(tmp, guess)


mode_direct       = solution.getLeadingMode()
guess             = mode_direct.getEigenValue()

# solve adjoint eigenproblem only for the leading eigenvalue
tmp = LinearSolver.solveGeneralEigenproblem(A,
                                            B,
                                            guess,
                                            1,
                                            adjoint=True)
solution.appendSolutionOfEigenProblem(tmp, guess, adjoint=True)


# end tracking time
end = time.time() - start
printDebug(True, '-- Solving the general eigenproblem took %4g s' % end)
residuum_max = solution.getMaximumError()
printDebug(True, '-- Maximum residuum of all solutions:  %12g' % (residuum_max))

# get leading modes
mode_direct  = solution.getLeadingMode(adjoint=False)
mode_adjoint = solution.getLeadingMode(adjoint=True)
eigenValue   = mode_direct.getEigenValue()
printDebug(True, '------------------------------------------------ ')
printDebug(True, '-- Leading eigenvalue:  ' + str(eigenValue))
printDebug(True, '------------------------------------------------ ')

# calculate scaling factor mode_adjoint^H * B * mode_direct



# scale modes s.t. mode_adjoint^H * B * mode_direct = 1
mode_direct_petsc  = mode_direct.getPetscVector()
mode_adjoint_petsc = mode_adjoint.getPetscVector()
temp               = mode_direct.getPetscVector() # gets a petsc vector "temp" of correct length
B.mult(mode_direct_petsc, temp)                   # temp = B*mode_direct
factor = temp.dot(mode_adjoint_petsc)             # factor = mode_adjoint ^H temp

#mode_direct_petsc.scale(1./np.sqrt(factor))
#mode_adjoint_petsc.scale(np.conj(1./np.sqrt(factor)))
#mode_direct.setCoefficientArray(mode_direct_petsc.getArray())
#mode_adjoint.setCoefficientArray(mode_adjoint_petsc.getArray())


# calculate rhs of adjoint base flow equation system (because this set of BCs is needed)
[u_direct, p_direct]              = mode_direct.getListOfSingleFields()
meanFlow._fieldDict['u_bilinear'] = u_direct.function
BL                                = equation.getBilinearOperator(meanFlow)

mode_adjoint_petsc = mode_adjoint.getPetscVector()
rhs                = mode_adjoint.getPetscVector() # get petsc vector of correct size
BL.multHermitian(mode_adjoint_petsc, rhs)          # rhs = BL^H * mode_adjoint


#-----------------------------------------------------------------------
## EXPORT SOLUTION
#-----------------------------------------------------------------------
# f is the growth rate (imaginary part) of the leading eigenvalue
f =  np.imag(eigenValue)
np.save("f.npy", f)

np.save("factor.npy", factor)

np.save("mode_direct.npy",  mode_direct.getCoefficientArray())
np.save("mode_adjoint.npy", mode_adjoint.getCoefficientArray())
np.save("rhs.npy", rhs.getArray())

## export modes in standard felics format
#solution_onlyLeading = ModeCollection(FEMSpaces.VMixed, mesh)
#solution_onlyLeading.appendMode(mode_direct)
#solution_onlyLeading.appendMode(mode_adjoint)
#fluctSolutList = solution_onlyLeading.getOldSolutionObject(meanFlow, param, FEMSpaces)
#ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)


