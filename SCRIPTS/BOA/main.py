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
## SOLVE EIGENPROBLEM
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

#-----------------------------------------------------------------------
## SOLVE ADJOINT BASEFLOW EQUATION 
#-----------------------------------------------------------------------
### manipulate velocity and get right-hand-side of base flow equation from a finite difference (epsilon = 1.e-8)
# TODO: Sophie: this is a quick (and dirty) implementation. Re-write once the restructuring of FELiCS has progressed sufficiently
epsilon        = 1.e-4
u_mean         = meanFlow._fieldDict['u']
[u_dir,p_dir]  = mode_direct.getListOfSingleFields()
coeff          = u_dir.getCoefficientArray()
coeff          = u_mean.x.array[:] + epsilon * coeff
u_dir.setCoefficientArray(coeff)
meanFlow._fieldDict['u'] = u_dir.function

# calculate disturbed operator and return "u" to its old value
A_1 = equation.getLinearOperator(meanFlow)
meanFlow._fieldDict['u'] = u_mean

# calculate finite difference of operator and multiplicate its transpose with the adjoint eigenvector
A_1.axpy(-1., A)      #A_1 = A_1 - A
A_1.scale(1./epsilon) #A_1 = A_1 / epsilon
mode_adjoint_petsc = mode_adjoint.getPetscVector()
rhs                = mode_adjoint.getPetscVector() #gets a petsc vector of correct length
mode_adjoint_petsc.conjugate()
A_1.multTranspose(mode_adjoint_petsc, rhs)

baseFlow_adjoint_array = LinearSolver.solveTransposeEquationSystem(A, rhs)
baseFlow_adjoint       = Field(FEMSpaces.VMixed, mesh)
baseFlow_adjoint.setCoefficientArray(baseFlow_adjoint_array)
baseFlow_adjoint.conjugate()


solution.popList()
newMode = Mode(FEMSpaces.VMixed, mesh)
newMode.function.x.array[:] = baseFlow_adjoint.function.x.array[:]
#newMode.function.x.array[:] = mode_adjoint.function.x.array[:]
newMode.setEigenValue(0.0)
newMode.isAdjoint = True
solution.appendMode(newMode)

#-----------------------------------------------------------------------
## CALCULATE GRADIENT WITH RESPECT TO GEOMETRY DEFORMATIONS 
#-----------------------------------------------------------------------
#geometryDeformer = CylinderBSpline()
#
#N = geometryDeformer.getNumberOfParameters()
#for i in range(N):
#    geometryDeformer.changeMesh(parameterIndex = i)
#    # change of nonlinear operator
#
#    # change of linear operator
#
#    geometryDeformer.changeMeshBackToOriginalState(parameterIndex = i)



#-----------------------------------------------------------------------
## WRITE FUNCITON VALUE AND FUNCTION GRADIENT INTO FILES 
#-----------------------------------------------------------------------

# f is the growth rate (imaginary part) of the leading eigenvalue
f =  np.imag(mode_direct.getEigenValue())
df = [0.3,0.4]
np.save("f.npy", f)
np.save("df.npy", df)


#-----------------------------------------------------------------------
## EXPORT SOLUTION
#-----------------------------------------------------------------------
fluctSolutList = solution.getOldSolutionObject(meanFlow, param, FEMSpaces)
ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)

