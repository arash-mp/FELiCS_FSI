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
from dolfinx.fem import (
							FunctionSpace,
							VectorFunctionSpace,
							)

from dolfinx.mesh import (
							refine,
						)

from ufl import (
					FiniteElement,
					MixedElement,
					triangle,
					VectorElement,
					tetrahedron,
					)
from os import (
	listdir,
	remove,
	mkdir,
	)

from functions import (
						printDebug,
						printError,
					)
from GUI.BCsSettingsClass import (
	FELiCSMesh,
	)

from Mapping import Mapping
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

		if param.BCs.dim==2:
			element_shape=triangle
		elif param.BCs.dim==3:
			element_shape = tetrahedron
		elementTypeStr='CG'
		self._degree = 2
		self._nVelocityComponents = len(param.Case.getVelocityComponents())
		# Create a finite element for a scalar
		FE_scalar=FiniteElement(
			elementTypeStr,
			element_shape,
			self._degree,
			)
		FE_scalar_p=FiniteElement(
			elementTypeStr,
			element_shape,
			self._degree,
			)
		# Create a finite element for a vector (velocity in x-y-plane)
		FE_vector=VectorElement(
			elementTypeStr,
			element_shape,
			self._degree,
			dim = self._nVelocityComponents,
			)

		self.FunctionSpaceVectorVelocity = VectorFunctionSpace(
			mesh,
			(elementTypeStr,
			 self._degree),
			dim = self._nVelocityComponents,
			)

		# refine the mesh to get the exportMesh:
		exportMesh = refine(mesh)
		exportMesh = FELiCSMesh(inputMesh=exportMesh)
		exportMesh.gdim = param.Case.nDim
		meshfileName = f'{param.Case.AnalysisMode}_mesh.h5'

		exportMesh.saveInFELiCSFormat(f'{param.Export.ExportFolder}/{meshfileName}')

		self.FunctionSpaceVectorVelocityExport = VectorFunctionSpace(
			exportMesh,
			(elementTypeStr,
			 1),
			dim = self._nVelocityComponents
			)



		# Collect all of them in a list
		MixedList=[]
		# Then a scalar for all the remaining quantities
		self.FunctionSpaceList=[]
		print(FE_vector)
		for name in param.SolutionList:
			printDebug(param.debug,'Adding finite element space for '+name+'...')
			if name=='u':
				#In order to stabilize the equations, the shape functions for pressure must be of one order lower than those of the velocity
				MixedList.append(FE_vector)
				self.FunctionSpaceList.append(FunctionSpace(mesh,FE_vector))
			else:
				MixedList.append(FE_scalar)
				self.FunctionSpaceList.append(FunctionSpace(mesh,FE_scalar))
		printDebug(param.debug,'The mixed finite element list is: ' + str(MixedList))
		# Create a element of the mixed function space
		element=MixedElement(MixedList)
		# Create a function space containing of mixed elements on the given mesh
		self.VMixed = FunctionSpace(mesh,element)

		self.FunctionSpaceVectorVelocityP1 = VectorFunctionSpace(
			mesh,
			(elementTypeStr,
			 1),
			dim=self._nVelocityComponents,
		)
		# Get a function space for the velocity components
		# Get function spaces for first order and second order elements...
		self.P1=FunctionSpace(mesh,FiniteElement(elementTypeStr, element_shape, 1))
		self.P2=FunctionSpace(mesh,FiniteElement(elementTypeStr, element_shape, 2))

		# Collect all of them in a list
		MixedList=[]
		# Then a scalar for all the remaining quantities
		self.FunctionSpaceListExport=[]
		for name in param.SolutionList:
			printDebug(param.debug,'Adding finite element space for '+name+'...')
			if name=='u':
				#In order to stabilize the equations, the shape functions for pressure must be of one order lower than those of the velocity
				MixedList.append(FE_vector)
				self.FunctionSpaceListExport.append(FunctionSpace(exportMesh,FE_vector))
			else:
				MixedList.append(FE_scalar)
				self.FunctionSpaceListExport.append(FunctionSpace(exportMesh,FE_scalar))
		printDebug(param.debug,'The mixed finite element list is: ' + str(MixedList))
		# Create a element of the mixed function space
		element=MixedElement(MixedList)
		# Create a function space containing of mixed elements on the given mesh
		self.VMixedExport = FunctionSpace(exportMesh,element)

		# Get a function space for the velocity components
		# Get function spaces for first order and second order elements...
		self.P1Export = FunctionSpace(exportMesh,FiniteElement(elementTypeStr, element_shape, 1))
		self.P2Export = FunctionSpace(exportMesh,FiniteElement(elementTypeStr, element_shape, 2))

		self.mappingObj = Mapping(self)

	def _projectField2allFEMSpaces(self, field, nfluctvar, nDim):
		"""
		This function is used to project the field on a FEM space to all the
		FEM spaces used for fluctuations.
		For example, this is useful to project the forcing/response limiter
		in resolvent analyses to all the fluctuations fields.

		NOTE: So far this function uses P2 for all variables!

		INPUTS:
			field: should be a FEM function, to be projected
			nfluctvar: total number of fluctuation variables
			nDim: number of spatial dimension = number of velocity components
		"""

		from dolfinx.fem import Function
		import sys
		import numpy as np

		fieldVMixed = Function(self.VMixed)
		for i in range(self.VMixed.num_sub_spaces):
			
			if self.VMixed.sub(i).num_sub_spaces > 0:

				for j in range(self.VMixed.sub(i).num_sub_spaces):
					space_ii, map_ii = self.VMixed.sub(i).sub(j).collapse()
					fieldVMixed.x.array[map_ii] = field.x.array
					
			else:
				space_i, map_i = self.VMixed.sub(i).collapse()
				fieldVMixed.x.array[map_i] = field.x.array



		# First we project to all velocity components
		#fieldvelocity = Function(self.FunctionSpaceVectorVelocity)
		#spacevel = []
		#fieldvel = []
		#for i in range(nDim):
		#	spacevel.append(self.P2) # Here we assume P2 by default!
	#		fieldvel.append(field)
			# fieldvelocity.sub(i).collapse().x.array[:] = field.x.array[:]
		#	fieldvelocity.sub(i).x.array[:] = field.x.array[:]

		# Then we project to the whole mixed space, including remaining scalar fields
		#nscalar = nfluctvar - nDim
		#fieldVMixed = Function(self.VMixed)
		#spacemixed = [self.FunctionSpaceVectorVelocity]
		#fieldmixed = [fieldvelocity]
		#for i in range(nscalar):
		#	spacemixed.append(self.P2) # Here we assume P2 by default!
	#		fieldmixed.append(field)
			# fieldVMixed.sub(i+1).collapse().x.array[:] = field.x.array[:]
	#		fieldVMixed.sub(i+1).x.array[:] = field.x.array[:]

		return fieldVMixed

