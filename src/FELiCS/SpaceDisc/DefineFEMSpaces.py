#/* Copyright (C) 2019 FLOW group TU Berlin - All Rights Reserved
# * You may NOT use, distribute or modify this code without explicit
# * opermission of the copyright owner, the FLOW group at TU Berlin
# * However, permissions to use and modify the code are generally
# * granted when asked for.
# * To ask for permission please contact t.kaiser@tu-berlin.de
# */
'''
# **********************************************************************
# * This file provides a class to which all FiniteElement Spaces are
# * inherent. This should be changed to a Dict and this File should be
# * then deleted
# * This file was created by Thomas L. Kaiser. Significant contributions
# * were made by
# * -
# *
# *
# ********************
'''
from dolfinx.fem import (   FunctionSpace,
                            VectorFunctionSpace)

from dolfinx.mesh import (refine)

from ufl import (FiniteElement,
                 MixedElement,
                 triangle,
                 VectorElement,
                 tetrahedron)

from os import (listdir,
                remove,
                mkdir)

from FELiCS.Misc.functions import (printDebug,
                                   printError)
from FELiCS.SpaceDisc.FELiCSMesh import FELiCSMesh

from FELiCS.IO.Mapping import Mapping
import pdb

class FEMSpacesClass():
    """Class containing the FE-spaces.
    When declared the object obtains the attributes
    \t -P1: Continuous Galerkin FEM-Space of first order test functions
    \t -P2: Continuous Galerkin FEM-Space of second order test functions
    \t -VMixed: More dimensional Continuous Galerkin FEM-Space containing all the unknowns
    Necessary function arguments:
    \t -parameter object
    \t -mesh object"""
    def __init__(self,param,mesh, degree=2):

        printDebug(True,'--------------------------------')
        printDebug(True,'-- Defining FEMSpaces...')

        if param.Case.nDim==2:
            element_shape = triangle
        elif param.Case.nDim==3:
            element_shape = tetrahedron
        elementTypeStr='CG'
        self._nVelocityComponents = param.BoundaryCondition.nVelocityComponents

        # refine the mesh to get the exportMesh:
        exportMesh_dolfinx = refine(mesh.dolfinxMesh)
        exportMesh         = FELiCSMesh(param.Case.CoordinateSystem,inputMesh=exportMesh_dolfinx)
        exportMesh.gdim    = param.Case.nDim
        meshfileName       = f'{param.Case.AnalysisMode}_mesh.h5'

        exportMesh.saveInFELiCSFormat(f'{param.Export.ExportFolder}/{meshfileName}')

        self.exportMesh = exportMesh
        
        if 'u' in param.getTransportedQuantityList():
            velocityOrder = param.Numerics.PolynomialOrder['u']
        else:
            velocityOrder = 2
        self.FunctionSpaceVectorVelocity = VectorFunctionSpace(
            mesh.dolfinxMesh,
            (self.elementTypeStr,
             velocityOrder),
            dim = self._nVelocityComponents,
            )
        self.FunctionSpaceVectorVelocityExport = VectorFunctionSpace(
            exportMesh.dolfinxMesh,
            (self.elementTypeStr,
             1),
            dim = self._nVelocityComponents
            )

        # Collect all of them in a list
        self.MixedList=[]

        # Then a scalar for all the remaining quantities
        self.FunctionSpaceList=[]
        self.FunctionSpaceListExport=[]
        for name in param.getTransportedQuantityList():
            printDebug(param.debug,"-- Adding finite element space of order "+
              f"{param.Numerics.PolynomialOrder[name]} for "+
              f"{name}")
            if name=='u':
                FE = VectorElement(
                                self.elementTypeStr,
                                self.element_shape,
                                param.Numerics.PolynomialOrder['u'],
                        dim = self._nVelocityComponents,
                                )   
                    
            else:
                FE=FiniteElement(
                                 self.elementTypeStr,
                                 self.element_shape,
                                 param.Numerics.PolynomialOrder[name],
                                 )
                    
            self.MixedList.append(FE)        

        printDebug(param.debug,'-- The mixed finite element list is: ' + str(self.MixedList))

        # Create a element of the mixed function space
        MixedFE=MixedElement(self.MixedList)

        # Create a function space containing of mixed elements on the given mesh
        self.VMixed = FunctionSpace(mesh.dolfinxMesh,MixedFE)

        self.FunctionSpaceVectorVelocityP1 = VectorFunctionSpace(
            mesh.dolfinxMesh,
            (self.elementTypeStr,
             1),
            dim=self._nVelocityComponents,
        )
        self.FunctionSpaceListExport=[]

#       # Get function spaces for first order and second order elements...
        self.P1=FunctionSpace(mesh.dolfinxMesh,FiniteElement(self.elementTypeStr, self.element_shape, 1))
        self.P2=FunctionSpace(mesh.dolfinxMesh,FiniteElement(self.elementTypeStr, self.element_shape, 2))

#       # Create a function space containing of mixed elements on the given mesh
        self.VMixedExport = FunctionSpace(exportMesh.dolfinxMesh,MixedFE)

        # Get a function space for the velocity components
        # Get function spaces for first order and second order elements...
        self.P1Export = FunctionSpace(exportMesh.dolfinxMesh,FiniteElement(self.elementTypeStr, self.element_shape, 1))
        self.P2Export = FunctionSpace(exportMesh.dolfinxMesh,FiniteElement(self.elementTypeStr, self.element_shape, 2))

        self.mappingObj = Mapping(self)

    def addCustomScalarSpaceToMixedSpace(self,mesh,order):
        FE=FiniteElement(
                         self.elementTypeStr,
                         self.element_shape,
                         order,
                         )
        self.MixedList.append(FE)

        # Create a element of the mixed function space
        MixedFE=MixedElement(self.MixedList)
        # Create a function space containing of mixed elements on the given mesh
        self.VMixed = FunctionSpace(mesh.dolfinxMesh,MixedFE)
#       # Create a function space containing of mixed elements on the export mesh
        self.VMixedExport = FunctionSpace(self.exportMesh.dolfinxMesh,MixedFE)

        #self.mappingObj = Mapping(self)

    def _projectField2allFEMSpaces(self, field, nfluctvar, nDim):
        """
        This function is used to project the field of a FEM space to all the
        FEM spaces used for fluctuations.
        For example, this is useful to project the forcing/response limiter
        in resolvent analyses to all the fluctuations fields.

        INPUTS:
            field: should be a FEM function, to be projected
            nfluctvar: total number of fluctuation variables
            nDim: number of spatial dimension = number of velocity components
   
        (ToDo: We need a P1 version of the field for cases where some of the fluctuations are P1)
        Update Sophie: Completed. Any polynomial order for any field can be used now (in this method)
        """

        from dolfinx.fem import Function
        import sys
        import numpy as np

        fieldVMixed = Function(self.VMixed)
        for i in range(self.VMixed.num_sub_spaces):
            
            if self.VMixed.sub(i).num_sub_spaces > 0:   # Vector field, we go through sub-space

                for j in range(self.VMixed.sub(i).num_sub_spaces):
                    space_ii, map_ii = self.VMixed.sub(i).sub(j).collapse()
                    targetFunction = Function(space_ii)
                    targetFunction.interpolate(field)
                    fieldVMixed.x.array[map_ii] = targetFunction.x.array
                    
            else:
                space_i, map_i = self.VMixed.sub(i).collapse() 
                targetFunction = Function(space_i)
                targetFunction.interpolate(field)
                fieldVMixed.x.array[map_i] = targetFunction.x.array



        # First we project to all velocity components
        #fieldvelocity = Function(self.FunctionSpaceVectorVelocity)
        #spacevel = []
        #fieldvel = []
        #for i in range(nDim):
        #   spacevel.append(self.P2) # Here we assume P2 by default!
    #       fieldvel.append(field)
            # fieldvelocity.sub(i).collapse().x.array[:] = field.x.array[:]
        #   fieldvelocity.sub(i).x.array[:] = field.x.array[:]

        # Then we project to the whole mixed space, including remaining scalar fields
        #nscalar = nfluctvar - nDim
        #fieldVMixed = Function(self.VMixed)
        #spacemixed = [self.FunctionSpaceVectorVelocity]
        #fieldmixed = [fieldvelocity]
        #for i in range(nscalar):
        #   spacemixed.append(self.P2) # Here we assume P2 by default!
    #       fieldmixed.append(field)
            # fieldVMixed.sub(i+1).collapse().x.array[:] = field.x.array[:]
    #       fieldVMixed.sub(i+1).x.array[:] = field.x.array[:]

        return fieldVMixed


