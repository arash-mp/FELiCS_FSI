#/* Copyright (C) 2019 FLOW group TU Berlin - All Rights Reserved
# * You may NOT use, distribute or modify this code without explicit
# * permission of the copyright owner, the FLOW group at TU Berlin!
# * However, permissions to use and modify the code are generally
# * granted when asked for.
# * To ask for permission please contact t.kaiser@tu-berlin.de.
# */
from h5py import File, string_dtype, AttributeManager
import pdb
import numpy as np
from functions import getLastGitCommit

class parameters():
	def __init__(self):
		from GUI.CaseSettingsClass import CaseSettingsClass
		self.Case=CaseSettingsClass()
		from GUI.FlowInputSettingsClass import FlowInputSettingsClass
		self.FlowInput=FlowInputSettingsClass()
		from GUI.BCsSettingsClass import BCsSettingsClass
		self.BCs = BCsSettingsClass(self.Case.MeshFilePath)
		from GUI.IOResolventSettingsClass import IOResolventSettingsClass
		self.IOResolvent = IOResolventSettingsClass()
		from GUI.NumericsSettingsClass import NumericsSettingsClass
		self.Numerics = NumericsSettingsClass()
		from GUI.ExportSettingsClass import ExportSettingsClass
		self.Export = ExportSettingsClass()
		pass
	def export(self, filestring):
		'''function exporting the parameters to a file
		\t Input:
		\t -filestring: path of parameter file'''
		from inspect import isclass

		# Writing the parameters to a file
		if '.h5' in filestring:
			file = File(filestring, 'a')

			if 'parameters' not in file.keys():
				paramGroup = file.create_group('parameters')
			else:
				paramGroup = file['parameters']
			#dt = string_dtype()
			# add the fenics-version:
			if 'FELiCSVersion' not in paramGroup.keys():
				felicsVersion = paramGroup.create_group('FELiCSVersion')
			else:
				felicsVersion = paramGroup['FELiCSVersion']
			felicsVersion.attrs.create('FELiCSVersion', data=getLastGitCommit())
		else:
			file = open(filestring,'w')

		# Iterate over all attributes of the object
		for group in dir(self):
			# Only not internal (marked by _) and not callable attributes should be exported
			if not group.startswith('_') and not callable(eval('self.'+group)):
				# Iterate over all attributes of the object
				for parameter in dir(eval('self.'+group)):
					# Only not internal (marked by _) and not callable or class attributes should be exported
					if not parameter.startswith('_') and not callable(eval('self.'+group+'.'+parameter)) and not isclass(parameter):

						# Check if the attributes are string, float, int, bool or list and export them
						if type(eval('self.'+group+'.'+parameter)) == str:
							if isinstance(file, File):
								if group not in list(paramGroup.keys()):
									currentGroup = paramGroup.create_group(group)
									# stringInArray = np.array([eval(f'self.{group}.{parameter}')])
									#

									currentGroup.attrs.create(parameter, data=eval(f'self.{group}.{parameter}'))
								else:

									#tringInArray = np.array(str(eval(f'self.{group}.{parameter}')))
									file[f'parameters/{group}'].attrs.create(parameter, data=eval(f'self.{group}.{parameter}'))
							else:
								file.writelines(parameter+'='+'\''+str(eval('self.'+group+'.'+parameter))+'\''+'\n')
						elif type(eval('self.'+group+'.'+parameter)) in [float,int,bool,list,dict]:
							if isinstance(file, File):
								if group not in list(paramGroup.keys()):
									currentGroup = paramGroup.create_group(group)
									#try:
									currentGroup.attrs.create(parameter, data=str(eval(f'self.{group}.{parameter}')))
									# except:
									# 	# stringInArray = np.array(str(eval(f'self.{group}.{parameter}')))
									#
									# 	currentGroup.attrs.create(parameter, data=str(eval(f'self.{group}.{parameter}')))
								else:
									#try:
									file[f'parameters/{group}'].attrs.create(parameter, data=str(eval(f'self.{group}.{parameter}')))
									# except:
									# 	stringInArray = np.array(str(eval(f'self.{group}.{parameter}')))
									# 	file[f'param/{group}'].create_dataset(parameter, data=stringInArray)
							else:
								file.writelines(parameter+'='+str(eval('self.'+group+'.'+parameter))+'\n')


		file.close()
	def importFromFile(self,filestring):
		''' Loading parameters from file
		\t Input:
		\t -filestring: path of parameter file
		'''
		from MixtureClass import MixtureClass
		if '.h5' in filestring:
			self.Case.importFromH5File(filestring)
			#The Mixture is not loaded but constructed from the inputs
			self.Case.Mixture=MixtureClass(self.Case.MixtureFilePath,self.Case.SpeciesFilePath)
			self.FlowInput.importFromH5File(filestring)
			self.BCs.importFromH5File(filestring)
			self.BCs.readDomainData(self.Case.MeshFilePath,
				self.Case.getExtendedTransportedQuantityList())
			self.IOResolvent.importFromH5File(filestring)
			self.Numerics.importFromH5File(filestring)
			self.Export.importFromH5File(filestring)
		else:

			self.Case.importSettings(filestring)
			#The Mixture is not loaded but constructed from the inputs
			self.Case.Mixture=MixtureClass(self.Case.MixtureFilePath,self.Case.SpeciesFilePath)
			self.FlowInput.importSettings(filestring)
			self.BCs.importSettings(filestring)
			self.BCs.readDomainData(self.Case.MeshFilePath,
				self.Case.nDim,
				self.Case.getExtendedTransportedQuantityList())
			self.IOResolvent.importSettings(filestring)
			self.Numerics.importSettings(filestring)
			self.Export.importSettings(filestring)
			self.getOldParameters()
	def complete(self):
		'''Check all parameters for completeness and consistency'''
		#Check all sub parameter groups for completeness.
		CompleteBool=False
		if self.Case.complete():
			CompleteBool = self.FlowInput.complete() and\
			   self.BCs.complete(self.Case.getExtendedTransportedQuantityList) and\
			   self.Numerics.complete(self.Case.getTransportedQuantityList())
		# IOResolvent is only checked if Analysis mode is wither of IO or Resolvent
			if self.Case.AnalysisMode in ['Input-Output','Resolvent']:
				CompleteBool= CompleteBool and self.IOResolvent.complete()
		return CompleteBool

	def getOldParameters(self):
		'''This function provides the parameters in the \'old\' fashion for compatibility with the rest of the code. This is redundant information and as a consequence the with the adaptations in the rest of the code the \'old\' parameters as set here in the function will disappear. In the end this function will be obsolete and must be deleted. '''
		# AdditionalVelocityComponents
		from copy import copy
		from functions import printDebug, printError
		#from fenics import DirichletBC,MeshFunction
		from dolfinx.fem import dirichletbc as DirichletBC
		from DefineFEMSpaces import FEMSpacesClass

		#nVelocityComponents
		self.nVelocityComponents=len(self.Case.getInternalVelocityComponents())
		#SolutionList
		self.SolutionList=self.Case.getTransportedQuantityList()
		#
		self.VelocityComponents = self.Case.getInternalVelocityComponents()
		#SpeciesList
		self.SpeciesList=self.Case.Mixture.getSpeciesList('transported')

		#ChemistryModel
		if len(self.Case.Mixture.ReactionMechanism.keys()):
			self.ChemistryModel = self.Case.Mixture.ReactionMechanism['type']
		else:
			self.ChemistryModel = ''

		self.ExportMode='both'

		#debug
		self.debug=True
		#NumericalScheme
		self.NumericalScheme = 'Continuous Galerkin'
		#Additional Species
		self.additionalSpecies=[]
