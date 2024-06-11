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
## MAIN PART
#-----------------------------------------------------------------------
# get matrices for eigenproblem
A = equation.getLinearOperator(meanFlow)
B = equation.getWeightMatrix  (meanFlow)

# get parameters for eigenproblem
guesses  = param.Numerics.EigenValueGuess
nSol     = param.Numerics.nSolut
adjoint  = param.Case.CalculateAdjoint

# track time
start= time.time()

# solve eigenproblem for each guess
solution = ModeCollection(FEMSpaces.VMixed, mesh)
for guess in guesses:
    tmp     = LinearSolver.solveGeneralEigenproblem(A,
                                                    B,
                                                    guess,
                                                    nSol,
                                                    )

    solution.appendSolutionOfEigenProblem(tmp, guess)

    if adjoint==True:
        tmp = LinearSolver.solveGeneralEigenproblem(A,
                                                    B,
                                                    guess,
                                                    nSol,
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

field_list = mode_direct.getListOfSingleFields()

print(len(field_list))



# solve adjoint equation system to get adjoint baseFlow
# 0.: field: getListOfFields, getPetscVector, getCoefficientArray
# 1. get u_dir from mode (method in field?)
# 2. u_dir = eps*u_dir + u_mean
# 3. save u_mean, put u_dir in dict
#A_1 = equation.getLinearOperator(meanFlow)
#BL  = (A_1 - L) / eps
# BL^T * conj(mode_adjoint )
#solution = LinearSolver.solveTransposeEquationSystem(A, b)
# base flow = conj(solution)




f = 1.
df = [0.3,0.4]
np.save("f.npy", f)
np.save("df.npy", df)


#-----------------------------------------------------------------------
## EXPORT SOLUTION
#-----------------------------------------------------------------------
fluctSolutList = solution.getOldSolutionObject(meanFlow, param, FEMSpaces)
ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)

