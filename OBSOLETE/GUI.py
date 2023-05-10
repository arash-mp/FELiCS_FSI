#/* Copyright (C) 2019 FLOW group TU Berlin - All Rights Reserved
# * You may NOT use, distribute or modify this code without explicit
# * permission of the copyright owner, the FLOW group at TU Berlin!
# * However, permissions to use and modify the code are generally 
# * granted when asked for.
# * To ask for permission please contact t.kaiser@tu-berlin.de.
# */
'''
# **********************************************************************
# * This file provides a GUI for using the Felics code. All necessary
# * settings needed for the calculation can be chosen. Furthermore, 
# * setting files (*.set) can be exported and imported, making a rund 
# * reproducable and increasing usability for the end user
# * This file was created by Mario Casel. Significant contributions 
# * were done by...
# **********************************************************************
'''
  
import tkinter as tk   # gui module
from tkinter import filedialog
#import matplotlib
import os
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import matplotlib.pyplot as plt
plt.close('all')
#%%
import numpy as np
from functools import partial # partial is used to pass arguments to command-based subfunctions within a tkinter-gui ('partial(sub_func, arg1,arg2,...)')
from functions import *

class parameters():
	def __init__(self, paramsTMP):
		# Write entries from paramsTMP (actually the GUI) into the class of parameters
		self.MeshPath   = paramsTMP.meshFileLstr.get()
		self.BaseFlowPath = paramsTMP.baseFlowLstr.get()
		self.BCsFile = paramsTMP.BCsFilestr.get()
		self.AveragingMode = paramsTMP.avgModestr.get()
		self.AveragingFromAxis = paramsTMP.avgAxFromstr.get()
#		self.AveragingToAxis = paramsTMP.avgAxTostr.get()
		if self.BaseFlowPath.split('.')[-1] == 'h5':
			self.BaseFlowMesh = paramsTMP.baseFlowMeshLstr.get()
        
		self.FlowMode = paramsTMP.FlowModestr.get()
		self.SolutionDirectory = paramsTMP.outDirLstr.get()
           
		self.molVisc                 = float(paramsTMP.molViscStr.get())
		self.m                  = float(paramsTMP.mstr.get())
		self.nSolut             = int(paramsTMP.nSolstr.get())
		self.EigenvalueGuess    = float(paramsTMP.eigGuessstr.get())
		self.nCPU               = int(paramsTMP.ncpustr.get())
		self.LimitRe                = paramsTMP.LimitReVar.get()
		if self.LimitRe:
			self.ReLimit 				= float(paramsTMP.ReLimitVar.get())
           
		self.NumericalScheme         = paramsTMP.NumericalSchemestr.get()
		self.CoordinateSystem   = paramsTMP.coordSysstr.get()
		self.ReynoldsNumberSource           = paramsTMP.ReSrcstr.get()
		if self.ReynoldsNumberSource == 'TKE-based':
			self.TKEnutConst = float(paramsTMP.TKEnutConstVar.get())
		#self.BoundaryDescription            = paramsTMP.BCdescrstr.get()
		self.FlowModeM          = paramsTMP.FlowModestr.get()
		self.AnalysisMode       = paramsTMP.anaModestr.get()
		self.ExportMode         = paramsTMP.videoExpstr.get()
		self.VideoExport            = paramsTMP.videoVar.get() 
		self.forcingMode        = paramsTMP.fModestr.get()
		#self.BCfile             = paramsTMP.BCsaveVar.get()
		#self.additionalSpecies  = list(paramsTMP.additionalSpeciesEntry.get().split(","))
		self.additionalSpecies  = paramsTMP.additionalSpeciesEntry.get().replace(']','').replace('[','').replace('\'','').split(',')
		#if '' in self.additionalSpecies:
		#	self.additionalSpecies.remove('')
		
		self.nonDim 			= paramsTMP.nonDimVar.get()
		if self.nonDim:
			self.Uinf 				= float(paramsTMP.UinfVar.get())
			self.charL 				= float(paramsTMP.charLVar.get())
		
		#Remove empty strings from list
		if '' in self.additionalSpecies:
			self.additionalSpecies.remove('')
		if self.AnalysisMode in ['Resolvent', 'Input-Output'] :
			self.omegas = list(set(paramsTMP.omegas + paramsTMP.omegas1 + paramsTMP.omegas2))
			self.omegas = sorted(paramsTMP.omegas)
			self.omegas = paramsTMP.omegas
			self.response_coeff = paramsTMP.response_coeff
			self.forcing_coeff = paramsTMP.forcing_coeff
            
		self.ChemistryModel = ' '
		if self.FlowMode in ['Reacting','LowMachReactingEnthalpy']:
			self.ChemistryModel = paramsTMP.ChemModelstr.get()

		## Debug Switch should be given to user
		self.debug=True
		## The following should be either asked from the user, or better, read out frmo the mesh data. 
		printWarning("Two dimensional mesh is assumed. If this is not the case, then at this position in the code the number of mesh dimensions must be determined...")
		self.nMeshDim=2
		self.getVelocityComponents()
		#self.getSpeciesListSolution()
		#self.getPreSpeciesListMean()
		#self.getSpeciesListMean()
		#printDebug(self.debug,"Mean species list is " +str(self.SpeciesList))
		self.getMeanFlowFieldNames()
		self.getSolutionInfo()
		if self.forcingMode == 'Boundary':
			self.BFidList = getListFromString(paramsTMP.bndFliststr.get())
			self.BFidList = [int(i) for i in self.BFidList]
	def getVelocityComponents(self):
		''' This function provides the velocity components that need to be solved for '''
		self.VelocityComponents=[]
		self.AdditionalVelocityComponents=[]
		self.VelocityComponents.append('x')
		if self.CoordinateSystem in ['Cartesian']:
			self.VelocityComponents.append('y')
		elif self.CoordinateSystem in ['Cylindrical']:
			self.VelocityComponents.append('r')
			self.AdditionalVelocityComponents.append('t')
		self.nVelocityComponents= len(self.VelocityComponents)
		printDebug(self.debug,"The velocity components are " +str(self.VelocityComponents))
		printDebug(self.debug,"The number of velocity components is  " + str(self.nVelocityComponents))

	def getMeanFlowFieldNames(self):
		''' This function provides the mean fields which mus be read in '''
		#import parameters as param
		self.MeanList=[]
		# Add velocity components
		self.MeanList.append('u')
		# If the number of dimensions of the mesh is smaller, than the number of velocities (e.g. homogeneous velocity in third dimension for 2D grids),
		# the additional momentum must be solved in an additional momentum equation. Then the respective velocity is not part of the velocity vector...
		# to be consistent not only for the solution but also for the mean field.
		for component in self.AdditionalVelocityComponents:
			self.MeanList.append('u'+component)
		#if self.ReynoldsNumberSource in ['Boussinesq', 'TKE-based', 'Boussinesq(xr)'] and self.CoordinateSystem == 'Cylindrical': 
		#	self.MeanList.extend(['rstxx', 'rstrr', 'rsttt', 'rstxr', 'rstxt', 'rstrt'])
		#elif self.ReynoldsNumberSource in ['Boussinesq', 'TKE-based', 'Boussinesq(xr)'] and self.CoordinateSystem == 'Cartesian': 
		#	self.MeanList.extend(['rstxx', 'rstyy', 'rstzz', 'rstxy', 'rstxz', 'rstyz'])
			
		# If necessary, add density
		if self.FlowMode in ['LowMach','Reacting','LowMachEnthalpy','LowMachReactingEnthalpy']:
			self.MeanList.append('rho')
		# Add species
		for specie in self.getMixtureDict()['species']:
			self.MeanList.append(specie)

		self.MeanList.extend(self.getMixtureDict()['species'])
		# If Input-Output analysis is used, the forcing must be read in (at least curently) for
		# every conservative variable ()...
		if self.AnalysisMode in ['Input-Output']:
			CurrentList=self.MeanList.copy()
			for entry in CurrentList:
				# If the field is the velocity vector, all components must be considered
				if entry == 'u':
					#for component in self.VelocityComponents:
					self.MeanList.append(entry+'_forcing_r')#jvs
					self.MeanList.append(entry+'_forcing_i')
				else:
					self.MeanList.append(entry+'_forcing_r')
					self.MeanList.append(entry+'_forcing_i')
		self.MeanList.append('nuSGS')
		## If necessary, add eddy viscosity Should not be used, since user can calculate Re number on his own!
		#if self.ReynoldsNumberSource=='EddyViscosity':
		#	self.MeanList.append('nuturb')	
		# If necessary, add temperature
		if self.FlowMode in ['Reacting','LowMachEnthalpy','LowMachReactingEnthalpy']:
			self.MeanList.append('T')
			
		if self.ReynoldsNumberSource in ['File']:
			
			self.MeanList.append('nuturb')
		if self.ReynoldsNumberSource in ['Boussinesq', 'TKE-based', 'Boussinesq(xr)'] and self.CoordinateSystem == 'Cylindrical': 
			self.MeanList.extend(['rstxx', 'rstrr', 'rsttt', 'rstxr', 'rstxt', 'rstrt','rstyy', 'rstzz', 'rstxy', 'rstxz', 'rstyz'])
		elif self.ReynoldsNumberSource in ['Boussinesq', 'TKE-based', 'Boussinesq(xr)'] and self.CoordinateSystem == 'Cartesian': 
			self.MeanList.extend(['rstxx', 'rstyy', 'rstzz', 'rstxy', 'rstxz', 'rstyz'])
		if self.FlowMode in ['LowMachEnthalpy','LowMachReactingEnthalpy']:
			self.MeanList.append('cp')
			self.MeanList.append('alpha')
			self.MeanList.append('he')
			self.MeanList.append('T')
		if self.FlowMode in ['LowMachReactingEnthalpy']:
			self.MeanList.append('dQ')
		printDebug(self.debug,"Mean flow fields to be read are " +str(self.MeanList))
	def getMixtureDict(self):
		
		#self.speciesList=[]
		try:
			Mixture='Mixture'
			fileMixture = open(Mixture,'r')
			MixtureDict = eval(fileMixture.read())
		#for species in list(SpeciesDict.keys()):
		#	self.speciesList.append(species)
		except: 
			MixtureDict={}
			MixtureDict['species']=[]
			MixtureDict['species_calc']=[]
		return(MixtureDict)
		#if self.ChemistryModel=='Global':
		#	#self.speciesList.append(param.ReactionProgressMarker )
		#	self.speciesList.append('CH4')
		#	self.speciesList.append('O2')
		#elif self.ChemistryModel=='BFER':
		#	self.speciesList.append('CH4' )
		#	self.speciesList.append('O2' )
		#	self.speciesList.append('CO2' )
		#	self.speciesList.append('CO')
		#	self.speciesList.append('H2O')
		#	self.speciesList.append('N2')
		## Add user specified additional species
		#for specie in self.additionalSpecies:
		#	self.speciesList.append(specie)
	def getSolutionInfo(self):
		''' This function provides the solution fields which are solved for '''
		from itertools import compress
		self.SolutionList=[]
		self.AllVelocityList=[]
		# Add velocity components
		self.SolutionList.append('u')
		for component in self.VelocityComponents:
			self.AllVelocityList.append('u'+component)
		# If the number of dimensions of the mesh is smaller, than the number of velocities (e.g. homogeneous velocity in third dimension for 2D grids),
		# the additional momentum must be solved in an additional momentum equation. Then the respective velocity is not part of the velocity vector...
		for component in self.AdditionalVelocityComponents:
			self.SolutionList.append('u'+component)
			self.AllVelocityList.append('u'+component)
		self.ExtendedSolutionList=self.AllVelocityList.copy()
		# Add pressure
		self.SolutionList.append('p')
		self.ExtendedSolutionList.append('p')
		# If necessary, add density
		if self.FlowMode in ['LowMach', 'Reacting','LowMachEnthalpy','LowMachReactingEnthalpy']:
			self.SolutionList.append('rho')
			self.ExtendedSolutionList.append('rho')
		# Add species
		mixtureDict=self.getMixtureDict()
		transportedSpecies=list(compress(mixtureDict['species'], [i == 'transported' for i in mixtureDict['species_calc']]))
		self.SolutionList.extend(transportedSpecies)
		self.ExtendedSolutionList.extend(list(compress(mixtureDict['species'], mixtureDict['species_calc'])))
		#for species in list(compress(mixtureDict['species'], mixtureDict['species_calc'])):
		#	if speciesDict[species]['calc']=='transported':
		#		self.SolutionList.append(species)
		#		self.ExtendedSolutionList.append(species)
		self.nSolutions= len(self.SolutionList)
		self.nExtendedSolutions= len(self.ExtendedSolutionList)
		printDebug(self.debug,"All velocity fields are " +str(self.SolutionList))
		printDebug(self.debug,"Solution fields are " +str(self.SolutionList))
		printDebug(self.debug,"Extended solution fields are " +str(self.ExtendedSolutionList))
	#def getSpeciesListSolution(self):
	#	''' This function provides a list for the species that are solved for '''
	#	self.SpeciesList=[]
	#	Species='Species'
	#	fileSpecies = open(Species,'r')
	#	SpeciesDict = eval(fileSpecies.read())
	#	for species in list(SpeciesDict.keys()):
	#		self.SpeciesList.append(species)
	#	##import parameters as param
	#	#self.SpeciesList=[]
	#	#if self.ChemistryModel=='Global':
	#	#	self.SpeciesList.append('CH4' )
	#	#elif self.ChemistryModel=='BFER':
	#	#	self.SpeciesList.append('CH4' )
	#	#	self.SpeciesList.append('O2' )
	#	#	self.SpeciesList.append('CO2' )
	#	#	self.SpeciesList.append('CO')
	#	#	self.SpeciesList.append('H2O')
	#	## Add user specified additional species
	#	#for specie in self.additionalSpecies:
	#	#	self.SpeciesList.append(specie)
	#	## Finally get the number of species
	#	#self.nSpecies=len(self.SpeciesList)	
	#	#printDebug(self.debug,"Solution species list is " +str(self.SpeciesList))
			
def getListFromString(string):
	tmpVar = ''
	List = []
	for k in range(0,len(string)):
			if not string[k] == '[' and not string[k] == ']' and not string[k] == ',':
				tmpVar += string[k]
				
			if string[k] == ',' or string[k] == ']':
				List += [float(tmpVar)]
				tmpVar = ''
	
	List = [round(k,3) for k in List]			
	return List
#%%
#def initBCdict(param, lineIDs, SolutionDict):
#	
#	bcDict = {}
#	for k in param.SolutionList:
#		if k == 'u':
#			for kk in param.VelocityComponents:
#				bcDict[k+kk] = []
#				for kkk in lineIDs:
#					bcDict[k+kk].append({'ID':kkk, 'type': 'Dirichlet', 'value':0})	
#		else:
#			bcDict[k] = []
#			for kkk in lineIDs:
#				bcDict[k].append({'ID':kkk, 'type': 'Dirichlet', 'value':0})	
#				
#	return bcDict

# Main GUI function 	
   
headerFont  = ("calibri", 16, "bold")
labelFont   = ('calibri', 10, 'bold') 
 
#%%
#Inside the Setup class all the gui related settings are implemented, such as the
#definition of input, label and interactive objects. Some of these objects are connected
#to certain events. E.g. "press button" -> "do something"
#
#GUI.py is structured as follows:
#   - Block1: initialization of the main window and all objects required
#   - Block2: Variable initialization, gridding, parameter initialization
#   - Block3: importing old settings from external file
#   - Block4: sub-routines called by events or other sub-routines  
#   - Block5: main routine - former main.py
#   - BOTTOM: setup class is called
#%% Block1: initialization of the main window and all objects required
class Setup:
	def __init__(self):
		# Get working directory
		self.workDir = os.getcwd()
			#Initialize Default values
		self.win = tk.Tk()      # creates the outer main frame of GUI
		self.win.resizable(True, True)
		self.win.title('Felics Setup')
		self.initVars()
		self.win.winfo_x
		## Definition of Input, Labels, Boxes, and Switches (Objects in GUI)
		self.importSetB     = tk.Button(self.win, text='Import Settings', \
					command=self.askForSettings)
		self.importSetLstr  = tk.StringVar()
		self.importSetL         = tk.Label(self.win, textvariable=self.importSetLstr)
		#self.refreshSetB    = tk.Button(self.win, text='Refresh',\
		#                        command=self.refreshSettings)
		######################## Path Settings #####################################
		self.pathSetL   = tk.Label(self.win, text='Path Settings', font = headerFont)
		   
		self.baseFlowB  = tk.Button(self.win, text='Baseflow File', \
		                      command=self.askForBaseFlow)
		self.baseFlowL   = tk.Label(self.win, textvariable=self.baseFlowLstr)
		    
		 
		self.BaseFlowMeshB = tk.Button(self.win, text='Import Base Flow Mesh', \
		                   command=self.askForBaseFlowMesh)  
		self.baseFlowMeshL   = tk.Label(self.win, textvariable=self.baseFlowMeshLstr)
		self.BCsFileB = tk.Button(self.win, text='Import BCs', \
		                   command=self.askForBCsFile)  
		self.BCsFileL   = tk.Label(self.win, textvariable=self.BCsFilestr)
		
		self.assBcB     = tk.Button(self.win, text='Assign BC\'s', \
                              command=self.assignBC)
		
#		self.impBcB     = tk.Button(self.win, text='Import BC\'s', \
#							  command=self.importBC)
		
		self.meshFileB  = tk.Button(self.win, text='Mesh File', \
		                      command=self.askForMeshFile)
		self.meshFileL    = tk.Label(self.win, textvariable=self.meshFileLstr)
			
		avgCH 	 = {'None','Azimuthal','Straight'}
		self.avgM 	= tk.OptionMenu(self.win, self.avgModestr, *avgCH)
			
		axisFromCH 	 = {'x','y','z'}
		self.axisFromCH 	= tk.OptionMenu(self.win, self.avgAxFromstr, *axisFromCH)
		                         
		self.outDirB    =  tk.Button(self.win, text='Output Directory', \
		                      command=self.createOutDir)
		self.outDirL    = tk.Label(self.win, textvariable=self.outDirLstr)
		####################### Base Flow Settings ################################
		self.bFlowSetL  = tk.Label(self.win, text='Baseflow Settings', font = headerFont)
		   
		self.NumericalSchemeL    = tk.Label(self.win, text='Numerical Scheme', font=labelFont)
		NumericalSchemestCH = {'Continuous Galerkin', 'Discontinuous Galerkin'}
		self.NumericalSchemeM    = tk.OptionMenu(self.win, self.NumericalSchemestr, *NumericalSchemestCH)
		   
		self.coordSysL  = tk.Label(self.win, text='Coordinate System', font=labelFont)
		coordSysCH = {'Cartesian', 'Cylindrical'}
		self.coordSysM  = tk.OptionMenu(self.win, self.coordSysstr, *coordSysCH,
		                                command=self.checkForFRcoeffs)
		   
		self.ReSrcL     = tk.Label(self.win, text='Local Re - calculation', font=labelFont)
		ReSrcCH = {'Constant', 'File', 'Sutherland', 'PowerLaw_AVBP', 'EddyViscosity', 'Boussinesq', 'Boussinesq(xr)' ,'TKE-based'}
		self.ReSrcM     = tk.OptionMenu(self.win, self.ReSrcstr, *ReSrcCH,\
		                       command=self.getReynoldsNumberParams)
		   
		self.ReL    = tk.Label(self.win, text='mol. visc.[m2/s]', font=labelFont)
		self.ReE        = tk.Entry(self.win, textvariable = self.molViscStr, width = 5)
		   
		self.mL     = tk.Label(self.win, text='Azim. wave number', font=labelFont)
		self.mE         = tk.Entry(self.win, textvariable = self.mstr, width = 5)
		   
		#self.BCdescrL   = tk.Label(self.win, text='Boundary Description', font=labelFont)
		#BCdescrCH = {'Boundary', 'Body', 'None', 'None'}
		#self.BCdescrM   = tk.OptionMenu(self.win, self.BCdescrstr, *BCdescrCH)
		   
		self.FlowModeL  = tk.Label(self.win, text='Flow Mode', font=labelFont)
		FlowModeCH = {'ColdFlow', 'ColdFlow Variable Density', 'Multi-Species', 'Reacting', 'LowMach','LowMachEnthalpy','LowMachReactingEnthalpy'}
		self.FlowModeM  = tk.OptionMenu(self.win, self.FlowModestr, *FlowModeCH,\
		                           command=self.checkForMode)
			
		self.LimitReL     = tk.Label(self.win, text='Limit Re', font=labelFont)
		self.LimitReCB = tk.Checkbutton(self.win, text='', variable = self.LimitReVar,\
		                       command=self.getReLimit)
		
		self.LimitReE     = tk.Entry(self.win, textvariable = self.ReLimitVar, width = 5)
		   
		self.anaModeL   = tk.Label(self.win, text='Analysis Mode', font=labelFont)
		anaModeCH = {'Modal', 'Resolvent','Input-Output'}
		self.anaModeM   = tk.OptionMenu(self.win, self.anaModestr, *anaModeCH,\
		                           command=self.checkForAnaMode)
	
		# Add label and entry for additional species
		self.additionalSpeciesLabel = tk.Label(self.win, text='AdditionalSpecies', font=labelFont)
		self.additionalSpeciesEntry = tk.Entry(self.win,textvariable=self.additionalSpeciesstr, width=10)
		
		self.nonDimBox = tk.Checkbutton(self.win, text='non-dimensionalized', variable = self.nonDimVar,\
		                       command=self.getnonDimParam)
		self.UinfL     = tk.Label(self.win, text='U_inf', font=labelFont)
		self.charLL     = tk.Label(self.win, text='char. Length', font=labelFont)
		self.charLE         = tk.Entry(self.win, textvariable = self.charLVar, width = 5)
		self.UinfE         = tk.Entry(self.win, textvariable = self.UinfVar, width = 5)
		
		self.TKEnutConstL     = tk.Label(self.win, text='TKEnut c', font=labelFont)
		self.TKEnutConstE     = tk.Entry(self.win, textvariable = self.TKEnutConstVar, width = 5)
		
		if self.nonDimVar.get() == True:
			self.ReL['text'] = 'mol. visc.[m2/s]'
		
		############ General Settings for computation and output ####################
		self.CompSetL   = tk.Label(self.win, text='Computation Settings', font = headerFont)
		self.nSolL  = tk.Label(self.win, text='Number of Solutions', font=labelFont)
		self.nSolE      = tk.Entry(self.win, textvariable = self.nSolstr, width = 5)
		   
		self.eigGuessL  = tk.Label(self.win, text='Eigenvalue Guess', font=labelFont)
		self.eigGuessE      = tk.Entry(self.win, textvariable = self.eigGuessstr, width = 5)
		   
		self.ncpuL  = tk.Label(self.win, text='Number of CPUs', font=labelFont)
		self.ncpuE      = tk.Entry(self.win, textvariable = self.ncpustr, width = 5)
		   
		self.videoL     = tk.Label(self.win, text='Export Video', font=labelFont)
		self.videoCB = tk.Checkbutton(self.win, text='', variable = self.videoVar,\
		                       command=self.getVideoState)
		   
		self.videoExpL  = tk.Label(self.win, text='Export Mode', font=labelFont)
		videoExpCH = {'vtk', 'mat','hdf5','both'}
		self.videoExpM  = tk.OptionMenu(self.win, self.videoExpstr, *videoExpCH)
		   
		self.saveSetB = tk.Button(self.win, text='Save Settings', \
                                command=self.saveSettingsFile)
		########################################################################
		self.runB   = tk.Button(self.win, text="RUN", font = ('calibri', 12, 'bold', 'underline'), command=self.runMain)
		self.quitB   = tk.Button(self.win, text="QUIT", font = ('calibri', 12, 'bold', 'underline'), command=self.win.destroy)
		self.pltFlow   = tk.Button(self.win, text="Plot Mean Flow", font = ('calibri', 12, 'bold', 'underline'), command=self.plotMeanFlows)
		 
		 
		self.progessStat = tk.Label(self.win, textvariable = self.progessText, font = ('calibri', 12, 'bold', 'underline'))
		 
		 
		self.gridding(key='basic')
		self.win.mainloop()
#%% Block2: Variable initialization, gridding, parameter initialization
#		- variables for the tkinter format are defined. these are called within the further workflow
#		- gridding is a subfunction which calls the tkinter's grid function or grid_forget function
#			- only after placing an object on the main windows's grid, the object appeas visible
#			- grid_forget let the object disappear
#		- paramter initialization is important for assigning the parameters set up inside the gui to 
#		  the variables names inside the main routine
		  
	def initVars(self):
		self.progessText = tk.StringVar()
	     
		self.baseFlowLstr     = tk.StringVar()
		self.BCsFilestr        = tk.StringVar()
		self.baseFlowMeshLstr = tk.StringVar()
		self.meshFileLstr     = tk.StringVar()
		self.avgModestr     = tk.StringVar()
		self.avgAxFromstr   = tk.StringVar()
		self.avgAxTostr     = tk.StringVar()
		self.outDirLstr     = tk.StringVar()
		self.NumericalSchemestr  = tk.StringVar()
		self.NumericalSchemestr.set('Continuous Galerkin')
		self.coordSysstr    = tk.StringVar()
		self.coordSysstr.set('Cartesian')
		self.ReSrcstr     	= tk.StringVar()
		self.ReSrcstr.set('Constant')
		self.ReLimitVar = tk.StringVar()
		self.molViscStr     		= tk.StringVar()
		self.molViscStr.set('NaN')
		self.mstr     		= tk.StringVar()
		self.mstr.set('NaN')
		self.bndFliststr   	= tk.StringVar()
		#self.BCdescrstr     = tk.StringVar()
		self.FlowModestr    = tk.StringVar()
		self.ChemModelstr   = tk.StringVar()
		self.anaModestr     = tk.StringVar()
		self.nSolstr     	= tk.StringVar()
		self.nSolstr.set('10')
		self.eigGuessstr    = tk.StringVar()
		self.eigGuessstr.set('NaN')
		self.ncpustr     	= tk.StringVar()
		self.ncpustr.set('1')
		self.plotVar     	= tk.IntVar()
		self.videoVar     	= tk.IntVar()
		self.videoExpstr    = tk.StringVar() 
		self.LimitReVar     	= tk.BooleanVar()
		self.omega1str  		= tk.StringVar()
		self.omega2strstart 	= tk.StringVar()
		self.omega2strstep 	= tk.StringVar()
		self.omega2strend 	= tk.StringVar()
		self.omegas2str 		= tk.StringVar()
		self.fModestr 		= tk.StringVar()
		self.fCoeffstr  		= tk.StringVar()
		self.rCoeffstr  		= tk.StringVar()
		self.omegasstr 		= tk.StringVar()
		self.chBuxVar 		= tk.IntVar()
		self.chBuyVar 		= tk.IntVar()
		self.chBurVar 		= tk.IntVar()
		self.chButVar 		= tk.IntVar()
		self.chBrhoVar 		= tk.IntVar()
		self.chBpVar 		= tk.IntVar()
		self.chBuxVarR 		= tk.IntVar()
		self.chBuyVarR 		= tk.IntVar()
		self.chBurVarR 		= tk.IntVar()
		self.chButVarR 		= tk.IntVar()
		self.chBrhoVarR 		= tk.IntVar()
		self.chBpVarR 		= tk.IntVar()
		self.saveSetdir     = tk.StringVar()
		self.BCsaveVar 		= tk.StringVar()
		self.additionalSpeciesstr  		= tk.StringVar()
	       
		self.plusBcnt = 1
		self.addcnt     = 1
		self.plotFlag = 0
		  
		self.linStepListV   = []
		self.linEndListV    = []
		self.linStartListV  = []
		  
		self.linStepListE = []
		self.linEndListE = []
		self.linStartListE = []
		  
		self.omegas2plus = []
		self.omegas2Eplus = []
		self.omPlusStr = tk.StringVar()
		
		self.nonDimVar = tk.BooleanVar()
		self.UinfVar = tk.StringVar()
		self.charLVar = tk.StringVar()
		
		self.TKEnutConstVar = tk.StringVar()
##########################################################################################
#The gridding function sub-divides in smaller sub-routines called by key word. At the 
#begging the constantly appearing objects are gridded. However, some objects' appearance
#rely on the user's input. Hence, variables objects are connected to events, letting them appear 
#or disappear
	def gridding(self, key):
		if key == 'basic':
			self.importSetB.grid        (row = 1, column = 2, rowspan = 1, columnspan = 1, sticky="E")
			self.importSetL.grid        (row = 1, column = 3, rowspan = 1, columnspan = 5, sticky="W")
			#self.refreshSetB.grid       (row = 1, column = 8, rowspan = 1, columnspan = 1)
			############################################################################
			self.pathSetL.grid          (row = 2, column = 1, rowspan = 1, columnspan = 3, pady=(20,0), sticky="W")
			self.baseFlowB.grid         (row = 3, column = 1, rowspan = 1, columnspan = 1, sticky="E")
			self.baseFlowL.grid         (row = 3, column = 2, rowspan = 1, columnspan = 10, sticky="W")
			self.meshFileB.grid         (row = 4, column = 1, rowspan = 1, columnspan = 1, sticky="E")
			self.meshFileL.grid         (row = 4, column = 2, rowspan = 1, columnspan = 10, sticky="W")
			self.BCsFileB.grid          (row = 5, column = 1, rowspan = 1, columnspan = 1, sticky="E")
			self.BCsFileL.grid          (row = 5, column = 2, rowspan = 1, columnspan = 10, sticky="W")
			self.assBcB.grid            (row = 4, column = 8, rowspan = 1, columnspan = 1, sticky="")
#			self.impBcB.grid            (row = 4, column = 9, rowspan = 1, columnspan = 1, sticky="W")
			self.avgM.grid              (row = 4, column = 6, rowspan = 1, columnspan = 1, sticky="E")
			self.axisFromCH.grid        (row = 4, column = 7, rowspan = 1, columnspan = 1, sticky="W")
			#            self.axisToCH.grid             (row = 4, column = 8, rowspan = 1, columnspan = 1, sticky="W")
			self.outDirB.grid           (row = 6, column = 1, rowspan = 1, columnspan = 1, sticky="E")
			self.outDirL.grid           (row = 6, column = 2, rowspan = 1, columnspan = 10, sticky="W")
			############################################################################
			self.bFlowSetL.grid         (row = 7, column = 1, rowspan = 1, columnspan = 3, pady=(20,0), sticky="W")
			self.NumericalSchemeL.grid       (row = 8, column = 1, rowspan = 1, columnspan = 1)
			self.NumericalSchemeM.grid       (row = 9, column = 1, rowspan = 1, columnspan = 1)
			self.coordSysL.grid         (row = 8, column = 2, rowspan = 1, columnspan = 1)
			self.coordSysM.grid         (row = 9, column = 2, rowspan = 1, columnspan = 1)
			self.ReSrcL.grid            (row = 8, column = 3, rowspan = 1, columnspan = 1)
			self.ReSrcM.grid            (row = 9, column = 3, rowspan = 1, columnspan = 1)
			self.LimitReL.grid    (row = 8, column = 4, rowspan = 1, columnspan = 1)
			self.LimitReCB.grid   (row = 9, column = 4, rowspan = 1, columnspan = 1)
			self.ReL.grid               (row = 8, column = 5, rowspan = 1, columnspan = 1, padx=(5,5))
			self.ReE.grid               (row = 9, column = 5, rowspan = 1, columnspan = 1)
			self.mL.grid                (row = 8, column = 6, rowspan = 1, columnspan = 1, padx=(5,5))
			self.mE.grid                (row = 9, column = 6, rowspan = 1, columnspan = 1)
			#self.BCdescrL.grid          (row = 8, column = 7, rowspan = 1, columnspan = 1)
			#self.BCdescrM.grid          (row = 9, column = 7, rowspan = 1, columnspan = 1)
			self.FlowModeL.grid         (row = 8, column = 7, rowspan = 1, columnspan = 1)
			self.FlowModeM.grid         (row = 9, column = 7, rowspan = 1, columnspan = 1)  
			self.anaModeL.grid          (row = 8, column = 8, rowspan = 1, columnspan = 1)
			self.anaModeM.grid          (row = 9, column = 8, rowspan = 1, columnspan = 1)
			self.additionalSpeciesLabel.grid          (row = 8, column = 9, rowspan = 1, columnspan = 1)
			self.additionalSpeciesEntry.grid          (row = 9, column = 9, rowspan = 1, columnspan = 1)
			############################################################################
			self.CompSetL.grid          (row = 35, column = 1, rowspan = 1, columnspan = 3, pady=(20,0), sticky="W")
			self.nSolL.grid             (row = 36, column = 1, rowspan = 1, columnspan = 1)
			self.nSolE.grid             (row = 37, column = 1, rowspan = 1, columnspan = 1)
			self.eigGuessL.grid         (row = 36, column = 2, rowspan = 1, columnspan = 1)
			self.eigGuessE.grid         (row = 37, column = 2, rowspan = 1, columnspan = 1)
			self.ncpuL.grid             (row = 36, column = 3, rowspan = 1, columnspan = 1)
			self.ncpuE.grid             (row = 37, column = 3, rowspan = 1, columnspan = 1)
			self.videoL.grid            (row = 36, column = 5, rowspan = 1, columnspan = 1)
			self.videoCB.grid           (row = 37, column = 5, rowspan = 1, columnspan = 1)
			self.videoL.grid            (row = 36, column = 6, rowspan = 1, columnspan = 1)
			self.videoCB.grid           (row = 37, column = 6, rowspan = 1, columnspan = 1)
			self.videoExpL.grid         (row = 36, column = 7, rowspan = 1, columnspan = 1)
			self.videoExpM.grid         (row = 37, column = 7, rowspan = 1, columnspan = 1)
			self.saveSetB.grid          (row = 1, column = 1, rowspan = 1, columnspan = 1)
			#            self.saveSetCH.grid         (row = 49, column = 1, rowspan = 1, columnspan = 1)
			self.nonDimBox.grid         (row = 7, column = 5, rowspan = 1, columnspan = 1, pady=(0,0), sticky="W")
			  
			self.runB.grid              (row = 50, column = 8, rowspan = 1, columnspan = 1, sticky="W", pady = (30,0))
			self.pltFlow.grid           (row = 50, column = 7, rowspan = 1, columnspan = 1, sticky="E", pady = (30,0))
			
		elif key == 'chemModel':
			if self.FlowMode == 'Reacting':
				self.ChemModelL.grid    (row = 10, column = 7, rowspan = 1, columnspan = 1)
				self.ChemModelM.grid    (row = 11,column = 7, rowspan = 1, columnspan = 1)
            
			if not self.FlowMode == 'Reacting' and hasattr(self, 'ChemModelL'):
				self.ChemModelL.grid_forget()
				self.ChemModelM.grid_forget()
				
		elif key=='nonDim':
			self.charLL.grid        (row = 7, column = 6, rowspan = 1, columnspan = 1, pady=(0,0), sticky="")
			self.charLE.grid        (row = 7, column = 7, rowspan = 1, columnspan = 1, pady=(0,0), sticky="")
			self.UinfL.grid        (row = 7, column = 8, rowspan = 1, columnspan = 1, pady=(0,0), sticky="")
			self.UinfE.grid        (row = 7, column = 9, rowspan = 1, columnspan = 1, pady=(0,0), sticky="")
		elif key=='nonDim_forget':
			self.charLL.grid_forget()
			self.UinfL.grid_forget()
			self.charLE.grid_forget()
			self.UinfE.grid_forget()
			
		elif key=='getReLimit':
			self.LimitReE.grid(row = 10, column = 4, rowspan = 1, columnspan = 1)
		elif key=='getReLimit_forget':
			self.LimitReE.grid_forget()
			
		elif key=='ReynoldsNumberSource':
			self.TKEnutConstL.grid(row = 7, column = 3, rowspan = 1, columnspan = 1, pady=(0,0), sticky="")
			self.TKEnutConstE.grid(row = 7, column = 4, rowspan = 1, columnspan = 1, pady=(0,0), sticky="")
		elif key=='ReynoldsNumberSource_forget':
			self.TKEnutConstL.grid_forget()
			self.TKEnutConstE.grid_forget()
				
		elif key == 'resolvent':
			self.omegaL.grid        (row = 12, column = 2, rowspan = 1, columnspan = 1)
			self.omega1L.grid       (row = 13, column = 2, rowspan = 1, columnspan = 1, sticky="E")
			self.omega1E.grid       (row = 13, column = 3, rowspan = 1, columnspan = 5, sticky="W")
			self.omega2Lstart.grid  (row = 15, column = 2, rowspan = 1, columnspan = 1, sticky="W")
			self.omega2Lstep.grid   (row = 15, column = 2, rowspan = 1, columnspan = 1)
			self.omega2Lend.grid    (row = 15, column = 2, rowspan = 1, columnspan = 1, sticky="E")
			self.omega2L.grid       (row = 15, column = 3, rowspan = 1, columnspan = 4, sticky="W")
			self.omLinspaceL.grid   (row = 16, column = 1, rowspan = 1, columnspan = 1, sticky="E")
			self.omega2Estart.grid  (row = 16, column = 2, rowspan = 1, columnspan = 1, sticky="W")
			self.omega2Estep.grid   (row = 16, column = 2, rowspan = 1, columnspan = 1)
			self.omega2Eend.grid    (row = 16, column = 2, rowspan = 1, columnspan = 1, sticky="E")
			self.omegasL.grid       (row = 26, column = 2, rowspan = 1, columnspan = 1, pady=(0,15), sticky="E")
			self.refreshOmegaB.grid (row = 12, column = 4, rowspan = 2, columnspan = 3)
			self.clearOmegasB.grid  (row = 12, column = 5, rowspan = 2, columnspan = 3)
			self.fModeL.grid        (row = 31, column = 1, rowspan = 1, columnspan = 1)
			self.fModeM.grid        (row = 32, column = 1, rowspan = 1, columnspan = 1)
			self.fCoeffL.grid       (row = 31, column = 2, rowspan = 1, columnspan = 1)
			self.rCoeffL.grid       (row = 31, column = 3, rowspan = 1, columnspan = 1)
			
		elif key == 'carthesianCoeffs':
			self.chBux.grid     (row = 32, column = 2, rowspan = 1, columnspan = 1, sticky="W")
			self.chBuy.grid     (row = 33, column = 2, rowspan = 1, columnspan = 1, sticky="W")
			self.chBuxR.grid    (row = 32, column = 3, rowspan = 1, columnspan = 1, sticky="W")
			self.chBuyR.grid    (row = 33, column = 3, rowspan = 1, columnspan = 1, sticky="W")
               
			if self.FlowModestr.get() == 'Reacting' or self.FlowModestr.get() == 'LowMach' or self.FlowModestr.get() == 'LowMachEnthalpy' or self.FlowModestr.get() == 'LowMachReactingEnthalpy':
				self.chBrho.grid    (row = 33, column = 2, rowspan = 1, columnspan = 1, sticky="")
				self.chBrhoR.grid   (row = 33, column = 3, rowspan = 1, columnspan = 1, sticky="")
				
			self.chBp.grid      (row = 32, column = 2, rowspan = 1, columnspan = 1, sticky="")  
			self.chBpR.grid     (row = 32, column = 3, rowspan = 1, columnspan = 1, sticky="")
			
		elif key == 'cylinderCoeffs':
			self.chBux.grid     (row = 32, column = 2, rowspan = 1, columnspan = 1, sticky="W")
			self.chBur.grid     (row = 33, column = 2, rowspan = 1, columnspan = 1, sticky="W")
			self.chBut.grid     (row = 34, column = 2, rowspan = 1, columnspan = 1, sticky="W")
			self.chBuxR.grid    (row = 32, column = 3, rowspan = 1, columnspan = 1, sticky="W")
			self.chBurR.grid    (row = 33, column = 3, rowspan = 1, columnspan = 1, sticky="W")
			self.chButR.grid    (row = 34, column = 3, rowspan = 1, columnspan = 1, sticky="W")

			if self.FlowModestr.get() == 'Reacting' or self.FlowModestr.get() == 'LowMach' or self.FlowModestr.get() == 'LowMachEnthalpy' or self.FlowModestr.get() == 'LowMachReactingEnthalpy':
				self.chBrho.grid    (row = 33, column = 2, rowspan = 1, columnspan = 1, sticky="")
				self.chBrhoR.grid   (row = 33, column = 3, rowspan = 1, columnspan = 1, sticky="")
				
			self.chBp.grid      (row = 32, column = 2, rowspan = 1, columnspan = 1, sticky="")  
			self.chBpR.grid     (row = 32, column = 3, rowspan = 1, columnspan = 1, sticky="")
			
		elif key == 'resolvent_forget':
			self.checkForFRcoeffs()
			try:
				self.omLinspaceL.grid_forget()
				self.omegaL.grid_forget()
				self.omega1L.grid_forget()
				self.omega2Lstart.grid_forget()
				self.omega2Lstep.grid_forget()
				self.omega2Lend.grid_forget()
				self.omega2L.grid_forget()
				self.omega2Estart.grid_forget()
				self.omega2Estep.grid_forget()
				self.omega2Eend.grid_forget()
				self.omegasL.grid_forget()
				self.refreshOmegaB.grid_forget()
				self.clearOmegasB.grid_forget()
				self.fModeL.grid_forget()
				self.fModeM.grid_forget()
				self.fCoeffL.grid_forget()
				self.rCoeffL.grid_forget()
				self.omega1E.grid_forget()
				self.omegas2E.grid_forget()
				self.omegasE.grid_forget()
				self.plusB.grid_forget()
				
			except: AttributeError
			
		elif key == 'Omegas' and self.plusBcnt <2:
			self.omegas2E.grid      (row = 16, column = 3, rowspan = 1, columnspan = 6, sticky="W")
			self.omegasE.grid       (row = 26, column = 3, rowspan = 1, columnspan = 6, sticky="W")
       
 
#%% Block3: importing old settings from external file
#The sub-routines in this section make the user put in information of the desired file to save the settings
#and ask for a path to store the file in.
#For importing old settings from a file the user had to choose file right at the beginning. After pressing
#REFRESH the function refreshSettings is called, copening the file, reading line after line. The lines
#inside the .set files are handled as a command in the terminal being executed.
	def saveSettingsFile(self):
		# Ask the user in which file to store the settings and write it directly to the GUI
		self.importSetL=filedialog.asksaveasfilename(initialdir = "./",title = "Select file",filetypes = (("settings files","*.set *.py"),("all files","*.*")))
		self.importSetLstr.set(self.importSetL)
		param = parameters(self) 
		paramNames = dir(param)
		file = open(self.importSetL,'w')
		for k in range(0, len(paramNames)):
			if paramNames[k].startswith('_') == False:
				if type(eval('param.'+paramNames[k])) == str:
					file.writelines(paramNames[k]+'='+'\''+eval('param.'+paramNames[k])+'\''+'\n')
				elif type(eval('param.'+paramNames[k])) == float\
				or type(eval('param.'+paramNames[k])) == int\
				or type(eval('param.'+paramNames[k])) == bool\
				or type(eval('param.'+paramNames[k])) == list:
					file.writelines(paramNames[k]+'='+str(eval('param.'+paramNames[k]))+'\n')
		file.close()
		
	def refreshSettings(self,filestring):
		nan='nan'
		#name = self.importSetLstr.get()
		file = open(filestring,'r')
		readFlag = 1
		while readFlag == 1:
			line = file.readline()
			try: 
				exec(line)
			except SyntaxError:
				pass
			if line == '':
				readFlag = 0
		       
		file.close()
		try: exec('self.baseFlowLstr.set(str(BaseFlowPath))')
		except NameError : exec('self.baseFlowLstr.set(str(None))')
		try: exec('self.BCsFilestr.set(str(BCsFile))')
		except NameError : exec('self.BCsFilestr.set(str(None))')
		try: exec('self.baseFlowMeshLstr.set(str(BaseFlowMesh))')
		except NameError : exec('self.baseFlowMeshLstr.set(str(None))')
		try: exec('self.meshFileLstr.set(str(MeshPath))')
		except NameError : exec('self.meshFileLstr.set(str(None))')
		#try: exec('self.BCsaveVar.set(str(BCfile))')
		#except NameError : exec('self.BCsaveVar.set(str(None))')
		try: exec('self.avgModestr.set(str(AveragingMode))')
		except NameError : exec('self.avgModestr.set(str(None))')
		try: exec('self.avgAxFromstr.set(str(AveragingFromAxis))')
		except NameError : exec('self.avgAxFromstr.set(str(None))')
		try: exec('self.outDirLstr.set(str(SolutionDirectory))')
		except NameError : exec('self.outDirLstr.set(str(None))')
		try: exec('self.NumericalSchemestr.set(str(NumericalScheme))')
		except NameError : exec('self.NumericalSchemestr.set(str(None))')
		try: exec('self.coordSysstr.set(str(CoordinateSystem))')
		except NameError : exec('self.coordSysstr.set(str(None))')
		try: exec('self.ReSrcstr.set(str(ReynoldsNumberSource))')
		except NameError : exec('self.ReSrcstr.set(str(None))')
		try: exec('self.molViscStr.set(str(molVisc))')
		except NameError : exec('self.molViscStr.set(str(None))')
		try: exec('self.additionalSpeciesstr.set(str(additionalSpecies))')
		except NameError : exec('self.additionalSpeciesstr.set(\'[]\')')
		try: exec('self.mstr.set(str(m))')
		except NameError : exec('self.mstr.set(str(None))')
		#try: exec('self.BCdescrstr.set(str(BoundaryDescription))')
		#except NameError : exec('self.BCdescrstr.set(str(None))')
		try: exec('self.FlowModestr.set(str(FlowMode))')
		except NameError : exec('self.FlowModestr.set(str(None))')
		try: exec('self.anaModestr.set(str(AnalysisMode))')
		except NameError : exec('self.anaModestr.set(str(None))')
		try: exec('self.nSolstr.set(str(nSolut))')
		except NameError : exec('self.nSolstr.set(str(None))')
		try: exec('self.eigGuessstr.set(str(EigenvalueGuess))')
		except NameError : exec('self.eigGuessstr.set(str(None))')
		try: exec('self.ncpustr.set(str(nCPU))')
		except NameError : exec('self.ncpustr.set(str(None))')
		try: exec('self.videoVar.set(int(VideoExport))')
		except NameError : exec('self.videoVar.set(int(0))')
		try: exec('self.videoExpstr.set(str(ExportMode))')
		except NameError : exec('self.videoExpstr.set(str(None))')
		try: exec('self.ChemModelstr.set(str(ChemistryModel))')
		except NameError : exec('self.ChemModelstr.set(str(None))')
		try: exec('self.fModestr.set(str(forcingMode))')
		except NameError : exec('self.fModestr.set(str(None))')
		try: exec('self.LimitReVar.set(bool(LimitRe))')
		except NameError : exec('self.LimitReVar.set(bool(False))')
		  
		try: exec('self.omega1str.set(str(omegas1))')
		except NameError : exec('self.omega1str.set(str(None))')
		try: exec('self.omega2strstart.set(str(omega2start))')
		except NameError : exec('self.omega2strstart.set(str(None))')
		try: exec('self.omega2strstep.set(str(omega2step))')
		except NameError : exec('self.omega2strstep.set(str(None))')
		try: exec('self.omega2strend.set(str(omega2end))')
		except NameError : exec('self.omega2strend.set(str(None))')
		try: exec('self.omegas2str.set(str(omegas2))')
		except NameError : exec('self.omegas2str.set(str(None))')
		try: exec('self.omegasstr.set(str(omegas))')
		except NameError : exec('self.omegasstr = tk.StringVar()')
		 
		try: exec('self.nonDimVar.set(bool(nonDim))')
		except NameError : exec('self.nonDimVar.set(bool(False))')
		try: exec('self.UinfVar.set(str(Uinf))')
		except NameError : exec('self.UinfVar.set(str(1))')
		try: exec('self.charLVar.set(str(charL))')
		except NameError : exec('self.charLVar.set(str(1))')
		
		try: exec('self.TKEnutConstVar.set(str(TKEnutConst))')
		except NameError : exec('self.TKEnutConstVar.set(str(0.2))')
		
		try: exec('self.ReLimitVar.set(str(ReLimit))')
		except NameError : exec('self.ReLimitVar.set(str())')
		
		try: exec('self.bndFliststr.set(str(BFidList))')
		except NameError : exec('self.bndFliststr = tk.StringVar()')
		 
		try: exec('self.response_coeff = response_coeff')
		except NameError : exec('self.response_coeff = [0,0,0,0]')
		try: exec('self.forcing_coeff = forcing_coeff')
		except NameError : exec('self.forcing_coeff = [0,0,0,0]')
		
		if self.baseFlowLstr.get().split('.')[-1] == 'h5':
			self.BaseFlowMeshB.grid(row = 3, column = 4, rowspan = 1, columnspan = 2, sticky="")
			self.baseFlowMeshL.grid(row = 3, column = 5, rowspan = 1, columnspan = 3, sticky="E")
			
			
		self.getnonDimParam()     
		self.getReynoldsNumberParams() 
		self.getReLimit()
		if self.anaModestr.get() in ['Resolvent','Input-Output']:
			self.coeffCheckBoxes()
			self.checkForAnaMode()
			self.addOmegas()
			self.checkForMode()


#%% Block4: sub-routines called by events or other sub-routines	
	def getReLimit(self,*args):
		if self.LimitReVar.get():
			self.gridding(key='getReLimit')
		else:
			self.gridding(key='getReLimit_forget')
			
	def getReynoldsNumberParams(self,*args):
		if self.ReSrcstr.get() == 'TKE-based':
			self.gridding(key='ReynoldsNumberSource')
		else:
			self.gridding(key='ReynoldsNumberSource_forget')
			
	def getnonDimParam(self):
		if self.nonDimVar.get() == True:
			self.gridding(key='nonDim')
			self.ReL['text'] = 'mol. visc.[m2/s]'
		else:
			self.gridding(key='nonDim_forget')
			self.ReL['text'] = 'mol. visc.[m2/s]'
	    
	def getVideoState(self):
		if self.videoVar.get() == 1:
			self.VideoExport = True
		else:
			self.VideoExport = False
	       
	def checkForAnaMode(self, *args):
		self.AnalysisMode       = self.anaModestr.get()
		self.checkForFRcoeffs()
		self.setFRcoeffs()
		if self.anaModestr.get() in ['Resolvent','Input-Output']:
	    ## Definition of Input, Labels, Boxes, and Switches (Objects in GUI)
			self.omegaL         	= tk.Label(self.win, text='Omegas', font = labelFont)
			self.omega1L    		= tk.Label(self.win, text='Omegas(1)', font = labelFont)
			self.omega1E        = tk.Entry(self.win, textvariable = self.omega1str, width=35)
			self.omLinspaceL    = tk.Label(self.win, text='linspace(start,step,stop)', font=labelFont)
			self.omega2L 		= tk.Label(self.win, text='Omegas(2)', font = labelFont)
			self.omega2Lstart 	= tk.Label(self.win, text='start', font = labelFont)
			self.omega2Lstep 	= tk.Label(self.win, text='step', font = labelFont)
			self.omega2Lend 		= tk.Label(self.win, text='end', font = labelFont)
			self.omega2Estart 	= tk.Entry(self.win, textvariable = self.omega2strstart, width=5)
			self.omega2Estep 	= tk.Entry(self.win, textvariable = self.omega2strstep, width=5)
			self.omega2Eend 		= tk.Entry(self.win, textvariable = self.omega2strend, width=5)
			self.omegas2E 		= tk.Entry(self.win, width=100)
			self.omegasE 		= tk.Entry(self.win, width=100)
			   
			self.refreshOmegaB 	= tk.Button(self.win, text='Add Omegas',\
			                       command=self.addOmegas)
			self.clearOmegasB 	= tk.Button(self.win, text='Clear Omegas',\
			                       command=self.clearOmegas)
			self.omegasL 		= tk.Label(self.win, text='Omegas(total)=(1)+(2)', font=labelFont)
			   
			self.fModeL     		= tk.Label(self.win, text='Forcing Mode', font=labelFont)
			fModeCH = {'Body', 'Boundary'}
			self.fModeM     		= tk.OptionMenu(self.win, self.fModestr, *fModeCH,\
			                    command=self.setfMode)
			self.forcingMode 	= self.fModestr.get()
			self.fCoeffL    		= tk.Label(self.win, text='Forcing Coefficients', font=labelFont)
			   
			   
			self.rCoeffL    		= tk.Label(self.win, text='Response Coefficients', font=labelFont)
			 
			self.gridding(key='resolvent')
	          
		if not self.AnalysisMode in ['Resolvent','Input-Output']:
			self.gridding(key='resolvent_forget')
	        
	def clearOmegas(self):
		self.omega1str.set('')
		self.omega2strstart.set('')
		self.omega2strstep.set('')
		self.omega2strend.set('')
		self.omegas2str.set('')
		self.omegasstr.set('')
		self.omegas1 = []
		self.omegas = []
	        
	def addOmegas(self):
		omega1TMP = self.omega1str.get()
		tmpVar = ''
		self.omegas1 = []
		for k in range(0,len(omega1TMP)):
			if not omega1TMP[k] == '[' and not omega1TMP[k] == ']' and not omega1TMP[k] == ',':
				tmpVar += omega1TMP[k]
				
			if omega1TMP[k] == ',' or omega1TMP[k] == ']':
				self.omegas1 += [float(tmpVar)]
				tmpVar = ''
				
			self.omegas1 = [round(k,3) for k in self.omegas1]
		try:        
			omega2TMPstart = float(self.omega2strstart.get())
			omega2TMPstep = float(self.omega2strstep.get())
			omega2TMPend = float(self.omega2strend.get())
			tmpVar2 = np.arange(omega2TMPstart, omega2TMPend, omega2TMPstep)
		except ValueError: 
			tmpVar2 = []
		if self.addcnt == 1 and type(self.omegasstr.get()) == str:
			self.omegas = list(self.omegasstr.get())
			tmpVar = ''
			omegasTMP = []
			for k in range(0,len(self.omegas)):
				if not self.omegas[k] == '[' and not self.omegas[k] == ']' and not self.omegas[k] == ',':
					tmpVar += self.omegas[k]
            
				if self.omegas[k] == ',' or self.omegas[k] == ']':
					omegasTMP += [float(tmpVar)]
					tmpVar = ''
            
			self.omegas = [round(k,3) for k in omegasTMP]

	        
		elif self.addcnt == 1 and self.omegasstr.get() == None:
			self.omegas = []
	    
		self.omegas2 = [round(k,3) for k in tmpVar2]
		self.omegas2str.set(str(self.omegas2))
		self.omegas2E = tk.Entry(self.win, textvariable = self.omegas2str, width=120)
		self.omegas = list(set(self.omegas + self.omegas1 + self.omegas2))
		self.omegas = sorted(self.omegas)
		self.omegasstr.set(str(self.omegas))
		self.addcnt+=1
	    
		self.omegasE = tk.Entry(self.win, textvariable = self.omegasstr, width=120)
		self.gridding(key='Omegas')
	            
	def setfMode(self, *args):
			self.forcingMode = self.fModestr.get()
			if self.fModestr.get() == 'Boundary':
				self.bndFidEntry = tk.Entry(self.win, textvariable = self.bndFliststr, width = 10)
				self.bndFidEntry.grid        (row = 33, column = 1, rowspan = 1, columnspan = 1)
	           
	def checkForMode(self, *args):
		self.checkForFRcoeffs()
		self.setfMode()
		self.FlowMode       = self.FlowModestr.get()
		if self.FlowMode in ['Reacting','LowMachReactingEnthalpy']:
			self.ChemModelL     = tk.Label(self.win, text='Model', font = labelFont)
			ChemModelCH = {'Global', 'BFER'}
			self.ChemModelM     = tk.OptionMenu(self.win, self.ChemModelstr, *ChemModelCH)
			self.ChemistryModel         = self.ChemModelstr.get()
			self.gridding(key = 'chemModel')
		elif not self.FlowMode in ['Reacting','LowMachReactingEnthalpy']:
			self.gridding(key = 'chemModel')
	          
	def askForMeshFile(self):
		self.MeshPath = os.path.relpath(filedialog.askopenfilename(filetypes = (("Mesh files","*.msh"),("all files","*.*"))),self.workDir)
		self.meshFileLstr.set(self.MeshPath)
		#self.meshFileLstr.set(filedialog.askopenfilename(\
		#				   filetypes = (("xml files","*.xml *.msh"),("all files","*.*")))\
	        #                      )
		#self.MeshPath = self.meshFileLstr.get()
		#self.ChooseTheFilePath(self.MeshPath)
	       
	def readBCFile(self, bcDict, flag):
		# MARIO: Was ist flag? Bitte so programmieren, dass man auf anhieb weiss, was welche Groesse ist...
		import DefineFEMSpaces
		from fenics import MeshFunction, DirichletBC
		# if boundaries were created by CREATE!-button in gui window -> user has possiblity to save settings
#		def saveBC(bcDict, mshPath):    # settings will be stored in bcDict and saved as such under a user-defined file (ending: .npy)
#			self.BCsaveVar.set(filedialog.asksaveasfilename(initialdir = "./",title = "Select file",filetypes = (("BC-Files","*.npy"),("all files","*.*"))))
#			np.save(self.BCsaveVar.get(), bcDict) 
			# read again via: bcDict = np.load('saveTest.npy').item()
#			self.saveSetdir.set(filedialog.asksaveasfilename(initialdir = "./",title = "Select file",filetypes = (("BC files","*.bc"),("all files","*.*"))))
#			name = self.saveSetdir.get()
#			file = open(name,'w')
#			file.write(str(bcDict))
#			file.close()
			
		param = parameters(self)
		# Initialize FEMSpaces
		self.FEMSpaces=DefineFEMSpaces.FEMSpacesClass(param,self.mesh)
         
		# when creating the .geo file for gmsh properly, a *_face_region.xml-file is created when converting the .msh-file with dolfin-convert
		# this _facet_region-file includes the regions which were labeled in .geo
		# Meshfunction uses the these regions and creates fenics-subdomains -> used to define Boundary Conditions
		#  subdomains = MeshFunction('size_t', self.mesh, param.MeshPath.split('.')[0] + '_physical_region.xml')
		boundaries = MeshFunction('size_t', self.mesh, param.MeshPath.split('.')[0] + '_facet_region.xml')
         
        
		# Open BCs file
		file = open(param.BCsFile,'r')
		
		bcDict = eval(file.read())
		file.close()
		   
		printDebug(param.debug,"Creating list with boundary conditions...")
		self.BClist = []
		# Get the intersection of the bcDict loaded from the file and the solutions needed by the solver
		#intersectVar=list(set(list(bcDict.keys())).intersection(param.ExtendedSolutionList))
		# Iterate through all the variables in the just read BCs file
		for varname in list(bcDict.keys()):
			# only write the read in information if varname is in the extended solution list (and is needed by the solver)
			if varname in param.ExtendedSolutionList:
				# Get the index of the varible in the Extended solution list
				varindex=param.ExtendedSolutionList.index(varname)
				# Get index of equation/variable i_eqn and if needed the index of the velocity component
				if varname[0]=='u' and varname[1] in param.VelocityComponents:
					i_eqn=0
					i_component=param.VelocityComponents.index(varname[1])
				else:
					i_eqn = param.SolutionList.index(varname)
				# Iterate through all the Boundaries (BCcount)
				for BCcount in range(0,len(bcDict[varname])):
					# Write the respective values to the GUI
					self.optMenueVar[BCcount][varindex].set(bcDict[varname][BCcount]['type'])
					self.EntryVar[BCcount][varindex].set(bcDict[varname][BCcount]['value'])
						
#					if bcDict[k][mm]['type']=='Dirichlet':
#						self.BClist.append( DirichletBC(self.FEMSpaces.VMixed.sub(m), bcDict[k][mm]['value'], boundaries, bcDict[k][mm]['ID']) )     
					#if bcDict[k][mm]['type']=='Dirichlet':
					#	if m<param.nVelocityComponents:
					#		self.BClist.append( DirichletBC(self.FEMSpaces.VMixed.sub(0).sub(m), bcDict[k][mm]['value'], boundaries, bcDict[k][mm]['ID']) )
					#	else:
					#		self.BClist.append( DirichletBC(self.FEMSpaces.VMixed.sub(m-param.nVelocityComponents+1), bcDict[k][mm]['value'], boundaries, bcDict[k][mm]['ID']) )     
					# Create a Boundary condition list (is this really needed?)
					if bcDict[varname][BCcount]['type']=='Dirichlet':
						printDebug(param.debug,"Adding Dirichlet BC for "+str(varname)+ " in equation "+str(i_eqn)+" with value "+str(bcDict[varname][BCcount]['value'])+" on boundary with index "+str(bcDict[varname][BCcount]['ID']))
						if varname in ['u'+ component for component in param.VelocityComponents]:
							self.BClist.append( DirichletBC(self.FEMSpaces.VMixed.sub(i_eqn).sub(i_component), bcDict[varname][BCcount]['value'], boundaries, bcDict[varname][BCcount]['ID']) )

							
						else:
							self.BClist.append( DirichletBC(self.FEMSpaces.VMixed.sub(i_eqn), bcDict[varname][BCcount]['value'], boundaries, bcDict[varname][BCcount]['ID']) )     
					
			
		self.bcDict = bcDict
	def createBC(self, bcDict,flag,*args):
		# MARIO: Was ist flag? Bitte so programmieren, dass man auf anhieb weiss, was welche Groesse ist...
		import DefineFEMSpaces
		from fenics import MeshFunction, DirichletBC
		# if boundaries were created by CREATE!-button in gui window -> user has possiblity to save settings
#		def saveBC(bcDict, mshPath):    # settings will be stored in bcDict and saved as such under a user-defined file (ending: .npy)
#			self.BCsaveVar.set(filedialog.asksaveasfilename(initialdir = "./",title = "Select file",filetypes = (("BC-Files","*.npy"),("all files","*.*"))))
#			np.save(self.BCsaveVar.get(), bcDict) 
			# read again via: bcDict = np.load('saveTest.npy').item()
#			self.saveSetdir.set(filedialog.asksaveasfilename(initialdir = "./",title = "Select file",filetypes = (("BC files","*.bc"),("all files","*.*"))))
#			name = self.saveSetdir.get()
#			file = open(name,'w')
#			file.write(str(bcDict))
#			file.close()
			
		param = parameters(self)
		# Initialize FEMSpaces
		if not hasattr(self,'FEMSpaces'): self.FEMSpaces=DefineFEMSpaces.FEMSpacesClass(param,self.mesh)
         
		# when creating the .geo file for gmsh properly, a *_face_region.xml-file is created when converting the .msh-file with dolfin-convert
		# this _facet_region-file includes the regions which were labeled in .geo
		# Meshfunction uses the these regions and creates fenics-subdomains -> used to define Boundary Conditions
		#  subdomains = MeshFunction('size_t', self.mesh, param.MeshPath.split('.')[0] + '_physical_region.xml')
		boundaries = MeshFunction('size_t', self.mesh, param.MeshPath.split('.')[0] + '_facet_region.xml')
         
		#if flag == 1:
		#	bcLoadVar = tk.StringVar()
		#	bcLoadVar.set(filedialog.askopenfilename(\
		#					filetypes = (("BC-Files","*.npy, *.bc"),("all files","*.*")))\
                #                  )
#		#	bcDict = np.load(bcLoadVar.get()).item()    # BC's are loaded from previously saved boundary settings (*.npy-fike)
		#	file = open(bcLoadVar.get(),'r')
		#	bcDict = eval(file.read())
		#	file.close()
        
		#if flag == 2:
		# Open BCs file
		#file = open(param.BCsFile,'r')
		
		#bcDict = eval(file.read())
		#file.close()
		   
		printDebug(param.debug,"Creating list with boundary conditions...")
		self.BClist = []
		for varname,varindex in zip(list(bcDict.keys()),range(0,len(bcDict.keys()))):
			# Get index of equation/variable i_eqn and if needed the index of the velocity component
			if varname[0]=='u' and varname[1] in param.VelocityComponents:
				i_eqn=0
				i_component=param.VelocityComponents.index(varname[1])
			else:
				i_eqn = param.SolutionList.index(varname)
			
			#for BCcount in range(0,len(bcDict[varname])):
			for BCcount in range(0,len(self.optMenueVar)):
				if flag == 0:
					bcDict[varname][BCcount]['type'] = self.optMenueVar[BCcount][varindex].get()
					bcDict[varname][BCcount]['value'] = self.EntryVar[BCcount][varindex].get()
				elif flag == 1 or flag == 2: # mode 1: gui-elements are set by previously created .npy-file
					self.optMenueVar[BCcount][varindex].set(bcDict[varindex][BCcount]['type'])
					self.EntryVar[BCcount][varindex].set(bcDict[varindex][BCcount]['value'])
					
#				if bcDict[k][mm]['type']=='Dirichlet':
#					self.BClist.append( DirichletBC(self.FEMSpaces.VMixed.sub(m), bcDict[k][mm]['value'], boundaries, bcDict[k][mm]['ID']) )     
				#if bcDict[k][mm]['type']=='Dirichlet':
				#	if m<param.nVelocityComponents:
				#		self.BClist.append( DirichletBC(self.FEMSpaces.VMixed.sub(0).sub(m), bcDict[k][mm]['value'], boundaries, bcDict[k][mm]['ID']) )
				#	else:
				#		self.BClist.append( DirichletBC(self.FEMSpaces.VMixed.sub(m-param.nVelocityComponents+1), bcDict[k][mm]['value'], boundaries, bcDict[k][mm]['ID']) )     
				if bcDict[varname][BCcount]['type']=='Dirichlet':
					printDebug(param.debug,"Adding Dirichlet BC for "+str(varname)+ " in equation "+str(i_eqn)+" with value "+str(bcDict[varname][BCcount]['value'])+" on boundary with index "+str(bcDict[varname][BCcount]['ID']))
					if varname in ['u'+ component for component in param.VelocityComponents]:
						self.BClist.append( DirichletBC(self.FEMSpaces.VMixed.sub(i_eqn).sub(i_component), bcDict[varname][BCcount]['value'], boundaries, bcDict[varname][BCcount]['ID']) )

						
					else:
						self.BClist.append( DirichletBC(self.FEMSpaces.VMixed.sub(i_eqn), bcDict[varname][BCcount]['value'], boundaries, bcDict[varname][BCcount]['ID']) )     
					
			
		self.bcDict = bcDict
		
		
		#if flag == 0:
		#	self.BCsaveVar.set(filedialog.asksaveasfilename(initialdir = "./",title = "Select file",filetypes = (("BC files","*.bc"),("all files","*.*"))))
		#	name = self.BCsaveVar.get()
		#	file = open(name,'w')
		#	file.write(str(bcDict))
		#	file.close()
		#	flag = 1
		file = open(param.BCsFile,'w')
		file.write(str(bcDict))
		file.close()
#		tk.Button(self.bcWin, text='Save!', \
#				command=partial(saveBC, bcDict, param.MeshPath.split('.')[0]) ).grid(row = 12, column = m, rowspan = 1, columnspan = 3, pady=(10,0))
      
	# Function being started from GUI when pressing button 'Assign BC'
	# for the new boundary description all the files from the mesh creation with gmsh are required as follows:
	#       - fName.geo <- assigning physical surfaces(interior) & physical curves(boundaries) is required + V2 ASCII fName.msh -export DONT save all elemens
	#       - fName.msh <- this has to be called within the GUI from now on, instead of an xml-file. The dolfin-convert is called within the Felics algorithm and creates the files necessary
	# ....... dolfin-convert fName.msh fName.xml will den create:
	#               1) fName.xml    <- main mesh file as previously known, for the interior
	#               2) fName_facet_region.xml   <- has the information for creating subdomains for the boundary conditions
	#               3) fName_physical_region.xml    <- not in use yet
     
	def assignBC(self):  

		param = parameters(self)
		import Mesh
		from fenics import plot
		import meshio
		from functions import getSolutionInfo
		from readBCFile2 import initBCdict
		SolutionDict=getSolutionInfo(param)
		plt.ioff() # do not show figures after matplotlib-plot command
         
		self.bcWin = tk.Tk()    # BC-window
		self.bcWin.resizable(True, True)
		self.bcWin.title('Boundary Conditions')
		#if not hasattr(self, 'mesh'): self.mesh   = Mesh.ReadMesh(self.meshFileLstr.get())     # read the *.msh file + dolfin-convert to xml-fileS (see above)
		self.mesh   = Mesh.ReadMesh(self.meshFileLstr.get())     # read the *.msh file + dolfin-convert to xml-fileS (see above)
		gmsh        = meshio.read(self.meshFileLstr.get())   # gmsh is required to visualize the lines one the boundary
         
        # find the lines assigned within .geo-file for physical boundaries
		print(type(gmsh.cells))
		print(gmsh.cells)
		borderIDX = gmsh.cells['line']
#		borderIDX = gmsh.cells[0]
		
		print(type(gmsh.cell_data))
		print(gmsh.cell_data)
		lineIDs = np.unique(gmsh.cell_data['line']['gmsh:physical'])    # lineIDs are essential for finding sub domains in *_facet_region.xmlg
		#lineIDs = np.unique(gmsh.cell_data['gmsh:physical'])
		lineIDX = [[]]
		# lines-list includes [ [[line1-xCoords],[line1-yCoords]], [[line2-xCoords],[line2-yCoords]], ... , [[lineN-xCoords],[lineN-yCoords]] ]
		# line 1 , x -coordinates = lines[0][0] <- length = amount of points on that certain line
		# line 1 , y -coordinates = lines[0][1] <- length = amount of points on that certain line
		lines = [[[],[]]]
		lineIDX[0] = list(borderIDX[gmsh.cell_data['line']['gmsh:physical']==lineIDs[0]][:,1])
#		lineIDX[0] = list(borderIDX[gmsh.cell_data['gmsh:physical']==lineIDs[0]][:,1])
		lines[0] = [list(self.mesh.coordinates()[lineIDX[0]][:,0]),list(self.mesh.coordinates()[lineIDX[0]][:,1])]
		m = 1
                 
		fig1 = plt.figure()     # fig handle to pass ot canvas 
		fig1.add_subplot(1,1,1)
		plot(self.mesh)     # plot of xml-mesh
		# If no BCs file chosen yet, create one
		if self.BCsFilestr.get() == 'None' or self.BCsFilestr.get() == '':
			bcDict = initBCdict(param, lineIDs, SolutionDict)			
			
			print("Creating new BCs file...")
			self.BCsFile = os.path.relpath(filedialog.asksaveasfilename(initialdir = "./",
				title = "Create boundary condition file",
				filetypes = (("BC files","*.bc"),
				("all files","*.*"))),
				self.workDir)
			self.BCsFilestr.set(self.BCsFile)
			file = open(self.BCsFile,'w')
			file.write(str(bcDict))
			file.close()
			#flag = 1
			param=parameters(self)
		plt.scatter(lines[0][0],lines[0][1], label='line-ID = '+str(lineIDs[0])) # visualize boundaries as scatter points
         
		for k in lineIDs[1:]: # first line in lines was created above -> following lines are automatically created in this loop
			lineIDX.append(list(borderIDX[gmsh.cell_data['line']['gmsh:physical']==lineIDs[m]][:,1]))
			lines.append([list(self.mesh.coordinates()[lineIDX[m]][:,0]),list(self.mesh.coordinates()[lineIDX[m]][:,1])])
                 
			plt.scatter(lines[m][0],lines[m][1], label='line-ID = '+str(lineIDs[m]))
			m+=1
                 
		plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",\
							 mode="expand", borderaxespad=0, ncol=3)
          
		canvas = FigureCanvasTkAgg(fig1, master = self.bcWin)
		canvas.draw()
		toolbar_frame = tk.Frame(self.bcWin)
		toolbar_frame.grid(row=1,column=1,columnspan=10)
		tBar = NavigationToolbar2Tk( canvas, toolbar_frame )
		tBar.update()
		canvas.get_tk_widget().grid             (row = 2, column = 1, rowspan = 1, columnspan = 10, pady=(20,20), sticky="N")
		# init boundary conditions dictionary: each dictionary entry will be another dictionary
		# first level: string indicating the line-ID/id pointing on the boundary assigned by gmsh: bcDict['lineID_#'] (replace # with integer)
		# second level: for each flow variable (created by SolutionDict) there is a float entry (initially set to zero)
		# entries on second level are assigned to DirichletBC's later in subfunction 'createBC' 
		# entries are set in GUI; NaN means that no boundary will be set for this values (variable 'choice')
		bcDict = initBCdict(param, lineIDs, SolutionDict)    
		tk.Label(self.bcWin, text='Line-ID', font=labelFont).grid    (row = 3, column = 1, rowspan = 1, columnspan = 1, pady=(10,0), sticky="W")
		for k,m in zip(bcDict.keys(),range(0,len(bcDict.keys()))):
			tk.Label(self.bcWin, text=k, font=labelFont).grid(row = 5+m, column = 1, rowspan = 1, columnspan = 1, pady=(10,0), sticky="W")
		choice = {'Dirichlet','Neumann'}
		self.optMenueVar = []
		self.EntryVar = []
		from functools import partial
		for k , m in zip( lineIDs, range(0,len(lineIDs)) ):
			tk.Label        (self.bcWin, text=str(k), font=labelFont).grid    (row = 3, column = 2+2*(m+1), rowspan = 1, columnspan = 2, pady=(10,0))
			tk.Label 		(self.bcWin, text='BC-type', font=labelFont).grid (row = 4, column = 2+2*(m+1), rowspan = 1, columnspan = 1, pady=(10,0))
			tk.Label 		(self.bcWin, text='value', font=labelFont).grid (row = 4, column = 2+2*(m+1)+1, rowspan = 1, columnspan = 1, pady=(10,0))
			self.optMenueVar.append([])
			self.EntryVar.append([])
			for kk,mm in zip(range(0, len(bcDict.keys())),list(bcDict.keys())):
				self.optMenueVar[m].append(tk.StringVar(self.bcWin))
				tk.OptionMenu(self.bcWin, self.optMenueVar[m][kk], *choice).grid(row = 5+kk, column = 2+2*(m+1), rowspan = 1, columnspan = 1, pady=(0,0))
				self.optMenueVar[m][kk].trace("w",partial(self.createBC,bcDict,0))
				
				self.EntryVar[m].append(tk.DoubleVar(self.bcWin))
				tk.Entry(self.bcWin, textvariable =  self.EntryVar[m][kk], width = 5).grid(row = 5+kk, column = 2+2*(m+1)+1, rowspan = 1, columnspan = 1, pady=(0,0))
				
				bcDict[mm][m]['ID'] = k
				bcDict[mm][m]['type'] = self.optMenueVar[m][kk].get()
				bcDict[mm][m]['value'] = self.EntryVar[m][kk].get()
            
		# here, third input into partial is a flag deciding over the mode of createBC (0=Dictionary is filled gui-values, 1=import bcSettings, gui values are set by imported *.npy-file)
		#tk.Button(self.bcWin, text='Create&Save!', \
		#				command=partial(self.createBC, bcDict, 0) ).grid          (row = 12, column = 2+m, rowspan = 1, columnspan = 3, pady=(10,0))
		#tk.Button(self.bcWin, text='Import!', \
		#				command=partial(self.createBC, bcDict, 1) ).grid          (row = 12, column = 1, rowspan = 1, columnspan = 3, pady=(10,0))
		#if hasattr(param, 'BCsFile'): self.createBC(bcDict, 2)
		self.readBCFile(bcDict, 2)
		
		self.bcWin.mainloop()
     
	def testTrace(self,test,*args):
		print("Variable changed!         "+test)
	def importBC(self,mesh): # sub function for import of previously assigned boundary conditinos -> Button: 'Import BC's'
		import imp # module for loading .py-file (possiblity to load oformer param_BCs.py-file)
		import DefineFEMSpaces
		from fenics import DirichletBC, MeshFunction
		import Mesh
		# if not already imported by 'Refresh Settings', the user is asked for choosing a file to open (either .py or .npy-file)
		if self.BCsaveVar.get() == 'None':
			self.BCsaveVar.set(filedialog.askopenfilename(initialdir = "./",title = "Select file",filetypes = (("BC-Files","*.npy *.py"),("all files","*.*"))))
            
		param = parameters(self)
		self.FEMSpaces=DefineFEMSpaces.FEMSpacesClass(param,mesh)
		boundaries = MeshFunction('size_t', mesh, param.MeshPath.split('.')[0] + '_facet_region.xml')
        
		
		if param.BCsFile.split('.')[-1] == 'bc':   # similar procedure as above (see createBC)
			file = open(param.BCsFile,'r')
			bcDict = eval(file.read())
			self.bcDict = bcDict
			file.close()
			self.optMenueVar = []
			self.EntryVar = []
			self.BClist = []
			for k,m in zip(list(bcDict.keys()),range(0,len(bcDict.keys()))):
				# Get index of equation/variable i_eqn and if needed the index of the velocity component
				if k[0]=='u' and k[1] in param.VelocityComponents:
					i_eqn=0
					i_component=param.VelocityComponents.index(k[1])
				elif k in param.SolutionList:
					i_eqn = param.SolutionList.index(k)
				else:
					continue
				
				for mm in range(0,len(bcDict[k])) :
					self.optMenueVar.append([])
					self.optMenueVar[mm].append(tk.StringVar().set(bcDict[k][mm]['type']))
#					self.optMenueVar[mm][m].set(bcDict[k][mm]['type'])
					self.EntryVar.append([])
					self.EntryVar[mm].append(tk.DoubleVar().set(bcDict[k][mm]['value']))
#					self.EntryVar[mm][m].set(bcDict[k][mm]['value'])
					if bcDict[k][mm]['type']=='Dirichlet':	

				
						printDebug(param.debug,"Adding Dirichlet BC for "+str(k)+ " in equation "+str(i_eqn)+" with value "+str(bcDict[k][mm]['value'])+" on boundary with index "+str(bcDict[k][mm]['ID']))
						if k in ['u'+ component for component in param.VelocityComponents]:
							self.BClist.append( DirichletBC(self.FEMSpaces.VMixed.sub(i_eqn).sub(i_component), bcDict[k][mm]['value'], boundaries, bcDict[k][mm]['ID']) )

							
						else:
							self.BClist.append( DirichletBC(self.FEMSpaces.VMixed.sub(i_eqn), bcDict[k][mm]['value'], boundaries, bcDict[k][mm]['ID']) )     
			# in case the user wrote a .py module (formerly known as param_BCs.py) this is imported vie imp-module 
		#elif param.BCsFile.split('.')[-1] == 'py':   
		#	param_BCs = imp.load_source('module.name', param.BCsFile)
		#	self.BClist = param_BCs.getBCs(self.FEMSpaces)
  

            
	def askForBCsFile(self):
		# Get the file path from the user. os.path.relpath provides the relative path to the file
		self.BCsFile = os.path.relpath(filedialog.askopenfilename(filetypes = (("BCs files","*bc"),("all files","*.*")),initialdir="./"),self.workDir)
		# Set the file path to the GUI
		self.BCsFilestr.set(self.BCsFile)
	
	def askForBaseFlowMesh(self):
		# Get the file path from the user. os.path.relpath provides the relative path to the file
		self.BaseFlowMesh = os.path.relpath(filedialog.askopenfilename(filetypes = (("Base Flow Mesh Files","*.mesh.h5 "),("all files","*.*")),initialdir="./"),self.workDir)
		# Set the file path to the GUI
		self.baseFlowMeshLstr.set(self.BaseFlowMesh)
	    
	def askForBaseFlow(self):    
		self.BaseFlowPath=os.path.relpath(filedialog.askopenfilename(filetypes = (("Base Flow Files","*.hdf5 *.f00001 *.h5 *.cgns *.mat *.vtk *.fel"),("all files","*.*")),initialdir="./"))
		self.baseFlowLstr.set(self.BaseFlowPath)
	    
		if self.BaseFlowPath.split('.')[-1] == 'h5':
			self.BaseFlowMeshB.grid(row = 3, column = 4, rowspan = 1, columnspan = 2, sticky="")
			self.baseFlowMeshL.grid(row = 3, column = 5, rowspan = 1, columnspan = 3, sticky="E")
		else:
			try: 
				self.BaseFlowMeshB.grid_forget()
				self.baseFlowMeshL.grid_forget()
			except: 
				pass
			
	def askForSettings(self):
		# Get the path for the settings file from the user
		self.importSetPath=os.path.relpath(filedialog.askopenfilename(\
				filetypes = (("set files","*.set"),("all files","*.*"))\
				,initialdir="./"),self.workDir)
		# Write the settings file to the GUI
		self.importSetLstr.set(self.importSetPath)	
		# Read the settings from the settingsfile
		self.refreshSettings(self.importSetPath)
	       
	def createOutDir(self):
		# Get output directory from user. os.path.relpath provides the relative path to the file
		self.SolutionDirectory = os.path.relpath(filedialog.askdirectory(initialdir='./'),self.workDir)
		# Write the output directory to the GUI
		self.outDirLstr.set(self.SolutionDirectory)
		# If the directory does not exist, create it
		if not os.path.exists(self.SolutionDirectory):
			os.makedirs(self.SolutionDirectory)
	       
	def browse_dir(self):
		self.folder_path.set(filedialog.askdirectory())
	    
	def coeffCheckBoxes(self):
		if self.coordSysstr.get() == 'Cartesian':
			for k in self.forcing_coeff:
				if k  == 0: self.chBuxVar.set(1)
				elif k == 1: self.chBuyVar.set(1)
				elif k == 2: self.chBpVar.set(1)
				elif k == 3: self.chBrhoVar.set(1)
	                
			for k in self.response_coeff:
				if k  == 0: self.chBuxVarR.set(1)
				elif k == 1: self.chBuyVarR.set(1)
				elif k == 2: self.chBpVarR.set(1)
				elif k == 3: self.chBrhoVarR.set(1)
		elif self.coordSysstr.get() == 'Cylindrical':
			for k in self.forcing_coeff:
				if k  == 0: self.chBuxVar.set(1)
				elif k == 1: self.chBurVar.set(1)
				elif k == 2: self.chButVar.set(1)
				elif k == 3: self.chBpVar.set(1)
				elif k == 4: self.chBrhoVar.set(1)
	                
			for k in self.response_coeff:
				if k  == 0: self.chBuxVarR.set(1)
				elif k == 1: self.chBurVarR.set(1)
				elif k == 2: self.chButVarR.set(1)
				elif k == 3: self.chBpVarR.set(1)
				elif k == 4: self.chBrhoVarR.set(1)
				
	def setFRcoeffs(self):
		if self.coordSysstr.get() == 'Cartesian':
			tmpList = [self.chBuxVar.get(), self.chBuyVar.get(),\
				  self.chBpVar.get(), self.chBrhoVar.get()]
	          
			tmpListR = [self.chBuxVarR.get(), self.chBuyVarR.get(),\
				   self.chBpVarR.get(), self.chBrhoVarR.get()]
	          
		elif self.coordSysstr.get() == 'Cylindrical':
			tmpList = [self.chBuxVar.get(), self.chBurVar.get(), self.chButVar.get(),\
				  self.chBpVar.get(), self.chBrhoVar.get()]
	          
			tmpListR = [self.chBuxVarR.get(), self.chBurVarR.get(), self.chButVarR.get(),\
				   self.chBpVarR.get(), self.chBrhoVarR.get()]
			
		self.forcing_coeff = []
		for (k,kk) in zip(tmpList, range(0,len(tmpList))):
			if k == 1:
				self.forcing_coeff.append(kk)
				
		self.response_coeff = []
		for (k,kk) in zip(tmpListR, range(0,len(tmpListR))):
			if k == 1:
				self.response_coeff.append(kk)
	          
	def checkForFRcoeffs(self, *args):
	      
		try: self.chBux.grid_forget() 
		except AttributeError: pass
		try: self.chBuy.grid_forget() 
		except AttributeError: pass
		try: self.chBur.grid_forget() 
		except AttributeError: pass
		try: self.chBut.grid_forget() 
		except AttributeError: pass
		try: self.chBrho.grid_forget() 
		except AttributeError: pass
		try: self.chBp.grid_forget() 
		except AttributeError: pass
		  
		try: self.chBuxR.grid_forget() 
		except AttributeError: pass
		try: self.chBuyR.grid_forget() 
		except AttributeError: pass
		try: self.chBurR.grid_forget() 
		except AttributeError: pass
		try: self.chButR.grid_forget() 
		except AttributeError: pass
		try: self.chBrhoR.grid_forget() 
		except AttributeError: pass
		try: self.chBpR.grid_forget() 
		except AttributeError: pass
	      
		if self.anaModestr.get() in ['Resolvent','Input-Output']:
			if self.FlowModestr.get() == 'Reacting' or self.FlowModestr.get() == 'LowMach' or self.FlowModestr.get() == 'LowMachEnthalpy'  or self.FlowModestr.get() == 'LowMachReactingEnthalpy':
				self.chBrho  = tk.Checkbutton(self.win, text='rho', variable = self.chBrhoVar,\
								  command=self.setFRcoeffs)
				self.chBrhoR  = tk.Checkbutton(self.win, text='rho', variable = self.chBrhoVarR,\
								   command=self.setFRcoeffs)
				
			self.chBp  = tk.Checkbutton(self.win, text='p', variable = self.chBpVar,\
							   command=self.setFRcoeffs)
			self.chBpR  = tk.Checkbutton(self.win, text='p', variable = self.chBpVarR,\
								command=self.setFRcoeffs)
			if self.coordSysstr.get() == 'Cartesian':
				self.chBux  = tk.Checkbutton(self.win, text='ux', variable = self.chBuxVar,\
								 command=self.setFRcoeffs)
				self.chBuy  = tk.Checkbutton(self.win, text='uy', variable = self.chBuyVar,\
								 command=self.setFRcoeffs)
	              
				self.chBuxR  = tk.Checkbutton(self.win, text='ux', variable = self.chBuxVarR,\
								  command=self.setFRcoeffs)
				self.chBuyR  = tk.Checkbutton(self.win, text='uy', variable = self.chBuyVarR,\
								  command=self.setFRcoeffs)
	              
				self.gridding(key = 'carthesianCoeffs')
				
			elif self.coordSysstr.get() == 'Cylindrical':
				self.chBux  = tk.Checkbutton(self.win, text='ux', variable = self.chBuxVar,\
								 command=self.setFRcoeffs)
				self.chBur  = tk.Checkbutton(self.win, text='ur', variable = self.chBurVar,\
								 command=self.setFRcoeffs)
				self.chBut  = tk.Checkbutton(self.win, text='ut', variable = self.chButVar,\
								 command=self.setFRcoeffs)
	              
				self.chBuxR  = tk.Checkbutton(self.win, text='ux', variable = self.chBuxVarR,\
								  command=self.setFRcoeffs)
				self.chBurR  = tk.Checkbutton(self.win, text='ur', variable = self.chBurVarR,\
								  command=self.setFRcoeffs)
				self.chButR  = tk.Checkbutton(self.win, text='ut', variable = self.chButVarR,\
								  command=self.setFRcoeffs)
	              
				self.gridding(key = 'cylinderCoeffs')
	          
	
	def plotMeanFlows(self):
		#self.createParam()
		param = parameters(self)
		import Import
		import Mesh 
		#import functions
		import DefineFEMSpaces
		from dolfin import plot
        
		mesh=Mesh.ReadMesh(param.MeshPath)
		FEMSpaces=DefineFEMSpaces.FEMSpacesClass(param,mesh)
		MeanFlowDict=Import.ReadMeanFlow(param,FEMSpaces)
        
        	 
		plt.ioff()
		fig = plt.figure()
		plt.set_cmap('coolwarm')
		i=0
		# Determine the number of plots needed
		nPlots=len(MeanFlowDict)+param.nVelocityComponents-1
		if param.AnalysisMode in ['Input-Output']:
			nPlots += (param.nVelocityComponents-1)*2
		for name in MeanFlowDict.keys():
			print(name)
			if name[0] == 'u' and name not in ['ut','ut_forcing_r','ut_forcing_i']:
				for idx,component in enumerate(param.VelocityComponents):
					nameComponent=name[:1]+component+name[1:]
					print(nameComponent)
					dofIDX = MeanFlowDict[name].function_space().sub(idx).dofmap().dofs()
					ax = fig.add_subplot(nPlots,1,i+1)
					cs=plot(MeanFlowDict[name][idx],\
						vmin=np.min(MeanFlowDict[name].split()[idx].vector()[dofIDX]),\
						vmax=np.max(MeanFlowDict[name].split()[idx].vector()[dofIDX]))
					cbar = fig.colorbar(cs,ticks=[np.min(MeanFlowDict[name].split()[idx].vector()[dofIDX]),\
						np.max(MeanFlowDict[name].vector()[:])])
					plt.title(nameComponent)
					i+=1
			else:
				ax = fig.add_subplot(nPlots,1,i+1)
				print(name)
				cs=plot(MeanFlowDict[name])
				cbar = fig.colorbar(cs,ticks=[np.min(MeanFlowDict[name].vector()[:]),\
									  np.max(MeanFlowDict[name].vector()[:])])
				plt.title(name)
				i+=1
		
		pltWin = tk.Tk()
		pltWin.title('Mean Flow')
		canvas = FigureCanvasTkAgg(fig, master = pltWin)
		canvas.draw()
		tBar = NavigationToolbar2Tk( canvas, pltWin )
		tBar.update()
		canvas._tkcanvas.pack(side = tk.TOP, fill = tk.BOTH, expand = 1)
		pltWin.mainloop()

	def runMain(self):
		''' Function initiating the main part of the program '''
		import time
		from runCase import runCase
		
		param = parameters(self) 
		useGUI=True
		runCase(param,useGUI)	
		#import Mesh 
		#import Import
		#import DefineFEMSpaces
		#from ExportSolution import ExportGUI
		#import global_variables as glob
		#import FEM
		#from readBCFile2 import readBCFile2
		## get Reynoldsnumber (factor in front of diffusive terms which is called Re) from molecular viscoisty
		#if param.nonDim:
		#	param.Re = param.Uinf*param.charL/param.molVisc
		#else:
		#	param.Re = 1.0/param.molVisc
		#
		#
		## if not hasattr() ... and elif ... ensure that that required objects are present and that they are initialized only once
		## Read mesh
		## Here the mesh should not be a part of params!!! This must be changed here and in all functions.
		#mesh = param.mesh=Mesh.ReadMesh(param.MeshPath)
		## Get FEM Spaces
		#FEMSpaces=DefineFEMSpaces.FEMSpacesClass(param,mesh)
		## Read BCs
		#self.importBC(mesh)
		#param.BClist, param.bcDict = readBCFile2(param, FEMSpaces)
		##param.BClist = self.BClist
		##param.bcDict = self.bcDict
		## Read mean flow
		#MeanFlowDict=Import.ReadMeanFlow(param,FEMSpaces)
         
		## Discretize the linearized equations
		#glob.matrix_dict=FEM.DiscretizeFlow(param,MeanFlowDict,FEMSpaces,mesh)
		###Solve GEVP or Resolvent problem
		#start1 = time.time()
		#if param.AnalysisMode=='Modal':
		#	results=FEM.solveGEVP(param,FEMSpaces)
		#elif param.AnalysisMode=='Resolvent':
		#	results=FEM.solveResolvent(param,FEMSpaces)
		#elif param.AnalysisMode=='Input-Output':
		#	results=FEM.solveInputOutput(param,FEMSpaces)
		###Export Solution
		#print('time for solveResolvent was: ' +str(time.time()-start1)+ 's. (Excluding the Export time)')
		#
		#ExportGUI(param,results,MeanFlowDict,FEMSpaces)
		#self.runB.grid_forget()

