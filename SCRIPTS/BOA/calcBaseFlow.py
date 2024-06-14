import pdb
import numpy as np
import copy
import time
import sys

import dolfinx

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
## CALCULATE BASE FLOW 
#-----------------------------------------------------------------------
#create nonlinear boundary conditions by hand (hard-coded):
#TODO: read nonlinear boundary file
from petsc4py.PETSc import ScalarType
bcs_nonlin = []
ft = equation.boundaries
#ux
#bcs_nonlin.append( dolfinx.fem.dirichletbc(ScalarType(1.), dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(0), 1, ft.find(1)),    FEMSpaces.VMixed.sub(0).sub(0)) ) # inlet
##bcs_nonlin.append( dolfinx.fem.dirichletbc(ScalarType(1.), dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(0), 1, ft.find(3)),    FEMSpaces.VMixed.sub(0).sub(0)) ) # outlet
##bcs_nonlin.append( dolfinx.fem.dirichletbc(ScalarType(1.), dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(0), 1, ft.find(4)),    FEMSpaces.VMixed.sub(0).sub(0)) ) # lateral
bcs_nonlin.append( dolfinx.fem.dirichletbc(ScalarType(0.), dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(0), 1, ft.find(1001)), FEMSpaces.VMixed.sub(0).sub(0)) )
#uy
#bcs_nonlin.append( dolfinx.fem.dirichletbc(ScalarType(0.), dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(1), 1, ft.find(1)),    FEMSpaces.VMixed.sub(0).sub(1)) ) # inlet
#bcs_nonlin.append( dolfinx.fem.dirichletbc(ScalarType(0.), dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(1), 1, ft.find(2)),    FEMSpaces.VMixed.sub(0).sub(1)) ) # symm
##bcs_nonlin.append( dolfinx.fem.dirichletbc(ScalarType(0.), dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(1), 1, ft.find(3)),    FEMSpaces.VMixed.sub(0).sub(1)) ) # outlet
#bcs_nonlin.append( dolfinx.fem.dirichletbc(ScalarType(0.), dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(1), 1, ft.find(4)),    FEMSpaces.VMixed.sub(0).sub(1)) ) # lateral
bcs_nonlin.append( dolfinx.fem.dirichletbc(ScalarType(0.), dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(0).sub(1), 1, ft.find(1001)), FEMSpaces.VMixed.sub(0).sub(1)) )
#p
#bcs_nonlin.append( dolfinx.fem.dirichletbc(ScalarType(0.), dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(1), 1, ft.find(1)),    FEMSpaces.VMixed.sub(1)) )
##bcs_nonlin.append( dolfinx.fem.dirichletbc(ScalarType(0.), dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(1), 1, ft.find(3)),    FEMSpaces.VMixed.sub(1)) )
##bcs_nonlin.append( dolfinx.fem.dirichletbc(ScalarType(0.), dolfinx.fem.locate_dofs_topological(FEMSpaces.VMixed.sub(1), 1, ft.find(4)),    FEMSpaces.VMixed.sub(1)) )


# initialize field with "1"
#TODO: alternative: read mean flow and take that as start field
baseFlow   = Field(FEMSpaces.VMixed, mesh)
baseFlow.setConstantValue(1.)
baseFlow.setBoundaryConditions(bcs_nonlin)


## start Newton solver

# track time
start= time.time()

i=0
residuum=1.
while(residuum > 3.e-11 and i<2):
    i+=1

    # put base flow arrays in meanFlow-obejct
    [u,p] = baseFlow.getListOfSingleFields()
    meanFlow._fieldDict['u'].x.array[:] = u.getCoefficientArray()
    meanFlow._fieldDict['p'].x.array[:] = p.getCoefficientArray()

    # solve equation system
    N = equation.getNonlinearExpression(meanFlow)
    L = equation.getLinearOperator(meanFlow)
    newtonSummand_array = LinearSolver.solveEquationSystem(L,N)

    # update baseFlow
    baseFlow_array = baseFlow.getCoefficientArray() - newtonSummand_array
    baseFlow.setCoefficientArray(baseFlow_array)
    baseFlow.setBoundaryConditions(bcs_nonlin)
    
    residuum = np.linalg.norm(newtonSummand_array)
    printDebug(True, "-- Base flow iteration: "+str(i)+"; Residuum: %4g " % residuum)


# end tracking time
end = time.time() - start
printDebug(True, '-- Solving the base flow problem took %4g s' % end)
printDebug(True, '-- Residuum:  %12g' % (residuum))

np.save("baseFlow.npy",  baseFlow.getCoefficientArray())

mode_direct = Mode(FEMSpaces.VMixed, mesh)
mode_adjoint = Mode(FEMSpaces.VMixed, mesh)
mode_adjoint.isAdjoint=True
mode_direct.setCoefficientArray(baseFlow.getCoefficientArray())
mode_adjoint.setCoefficientArray(baseFlow.getCoefficientArray())

mode_direct.setEigenValue(0.)
mode_adjoint.setEigenValue(0.)

# export modes in standard felics format
solution_onlyLeading = ModeCollection(FEMSpaces.VMixed, mesh)
solution_onlyLeading.appendMode(mode_direct)
solution_onlyLeading.appendMode(mode_adjoint)
fluctSolutList = solution_onlyLeading.getOldSolutionObject(meanFlow, param, FEMSpaces)
ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)


