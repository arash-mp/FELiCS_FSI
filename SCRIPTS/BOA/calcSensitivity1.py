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
## CALCULATE GRADIENT WITH RESPECT TO GEOMETRY DEFORMATIONS 
#-----------------------------------------------------------------------
A_0  = equation.getLinearOperator(meanFlow)

# deform mesh
epsilon = 1.e-4
x = mesh.geometry.x
bcDict        = param.BCs.getBCsDict()
bcs           = equation.boundaries
bc_facet_dofs = bcs.indices[bcs.values==1001]

import dolfinx

mesh.topology.create_connectivity(mesh.topology.dim-1, mesh.topology.dim)
#boundary_facets    = dolfinx.mesh.exterior_facet_indices(mesh.topology)
bc_dofs  = dolfinx.mesh.compute_incident_entities(mesh.topology,bc_facet_dofs, mesh.topology.dim-1, 0)
#vertex_to_geometry = dolfinx.cpp.mesh.entities_to_geometry(mesh._cpp_object, 0, boundary_vertices, False)

print(x[bc_dofs])




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



