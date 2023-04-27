import tkinter as tk
from GUI.ToolTip import CreateToolTip
from GUI.GUISettings import labelFont
class CaseSettingsGUI():
	def __init__(self,mainGUI):
		'''Initializing the caseSettings GUI
		\t Input:
		\t -mainGUI: The FELiCS main GUI object
		'''
		self.Case=mainGUI.param.Case
		self.workDir=mainGUI.workDir
		self.BCs=mainGUI.param.BCs
		self.window =tk.Toplevel(mainGUI.window)
		self.window.title("Case Settings")
		AllSettings=self.Case.getAllSettingsDict()
		# If mesh file was changed, it needs to be reread:
		self.__NewMeshFileBool__=False

		# For choosing the number of dimensions
		nDimColumn=1
		nDimRow=1
		self.window.nDimL  = tk.Label(self.window, text='No. dimensions', font=labelFont())
		self.window.nDimL.grid         (row =nDimRow , column = nDimColumn, rowspan = 1, columnspan = 1)
		self.nDimStr    = tk.StringVar()
		nDimCH = {'1', '2','3'}
		self.window.nDimM  = tk.OptionMenu(self.window, self.nDimStr, *nDimCH,command=self.refresh)
		self.window.nDimM.grid         (row =nDimRow+1 , column = nDimColumn, rowspan = 1, columnspan = 1)
		if hasattr(self.Case,'nDim'):
			self.nDimStr.set(str(self.Case.nDim))
		else:
			self.nDimStr.set(str(AllSettings['nDim']['default']))

		# For choosing the coordiate system
		coordSysColumn=2
		coordSysRow=1
		self.window.coordSysL  = tk.Label(self.window, text='Coordinate System', font=labelFont())
		self.window.coordSysL.grid         (row =coordSysRow , column = coordSysColumn, rowspan = 1, columnspan = 1)
		self.coordSysStr    = tk.StringVar()
		coordSysCH = {'Cartesian', 'Cylindrical'}
		self.window.coordSysM  = tk.OptionMenu(self.window, self.coordSysStr, *coordSysCH,command=self.refresh)
		self.window.coordSysM.grid         (row =coordSysRow+1 , column = coordSysColumn, rowspan = 1, columnspan = 1)
		if hasattr(self.Case,'CoordinateSystem'):
			self.coordSysStr.set(self.Case.CoordinateSystem)
		else:
			self.coordSysStr.set(AllSettings['CoordinateSystem']['default'])

		# For choosing the analysis mode (Modal analysis, RA or Input-Output)
		analysisModeColumn=3
		analysisModeRow=1
		self.window.analysisModeL  = tk.Label(self.window, text='Analysis Type', font=labelFont())
		self.window.analysisModeL.grid         (row =analysisModeRow , column = analysisModeColumn, rowspan = 1, columnspan = 1)
		self.analysisModeStr    = tk.StringVar()
		analysisModeCH = {'Modal', 'Input-Output', 'Resolvent'}
		self.window.analysisModeM  = tk.OptionMenu(self.window, self.analysisModeStr, *analysisModeCH,command=self.doNothing)
		self.window.analysisModeM.grid         (row =analysisModeRow+1 , column = analysisModeColumn, rowspan = 1, columnspan = 1)
		if hasattr(self.Case,'AnalysisMode'):
			self.analysisModeStr.set(self.Case.AnalysisMode)
		else:
			self.analysisModeStr.set(AllSettings['AnalysisMode']['default'])


		# Model for Molecular Viscosity
		molViscModelColumn=4
		molViscModelRow=1
		self.window.molViscModelL  = tk.Label(self.window, text='Mol. Visc.', font=labelFont())
		self.window.molViscModelL.grid         (row =molViscModelRow , column = molViscModelColumn, rowspan = 1, columnspan = 1)
		self.molViscModelStr    = tk.StringVar()
		molViscModelCH = {'File','Constant','Sutherland'}
		self.window.molViscModelM  = tk.OptionMenu(self.window, self.molViscModelStr, *molViscModelCH,command=self.refresh)
		self.window.molViscModelM.grid         (row =molViscModelRow+1 , column = molViscModelColumn, rowspan = 1, columnspan = 1)
		self.molViscModelStr.set(self.Case.MolViscModel)
		self.molViscStr    = tk.StringVar()
		self.molViscE         = tk.Entry(self.window, textvariable = self.molViscStr, width = 5)
		self.molViscE.grid         (row = molViscModelRow+2 , column = molViscModelColumn, rowspan = 1, columnspan = 1)
		self.molViscStr.set(self.Case.MolVisc)

		# Model for Molecular Viscosity Perturbation
		molViscPerturbModelColumn=5
		molViscPerturbModelRow=1
		self.window.molViscPerturbModelL  = tk.Label(self.window, text='Mol. Visc. Pert.', font=labelFont())
		self.window.molViscPerturbModelL.grid         (row =molViscPerturbModelRow , column = molViscPerturbModelColumn, rowspan = 1, columnspan = 1)
		self.molViscPerturbModelStr    = tk.StringVar()
		molViscPerturbModelCH = {'None','Sutherland','Sutherland mean'}
		self.window.molViscPerturbModelM  = tk.OptionMenu(self.window, self.molViscPerturbModelStr, *molViscPerturbModelCH,command=self.refresh)
		self.window.molViscPerturbModelM.grid         (row =molViscPerturbModelRow+1 , column = molViscPerturbModelColumn, rowspan = 1, columnspan = 1)
		self.molViscPerturbModelStr.set(self.Case.MolViscPerturbModel)
		CreateToolTip([self.window.molViscPerturbModelL,self.window.molViscPerturbModelM],
			text= 'Define how the fluctuation in diffusion (momentum, enthalpy, species) is taken into account.\n'
			'\'None\': No fluctuation in diffusion.\n'
			'\'Sutherland\': The diffusion fluctuation is based on the Sutherland law. The base state is obtained by an evaluation of the Sutherland law at the temporally averaged temperature.\n'
			'\'Sutherland mean\': The diffusion fluctuation is based on the Sutherland law. The base state is the temporally averaged diffusion which is read from file.\n')

		# For choosing the Turbulence Model
		turbulenceModelColumn=6
		turbulenceModelRow=1
		self.window.turbulenceModelL  = tk.Label(self.window, text='Turbulence Model', font=labelFont())
		self.window.turbulenceModelL.grid         (row =turbulenceModelRow , column = turbulenceModelColumn, rowspan = 1, columnspan = 1)
		self.turbulenceModelStr    = tk.StringVar()
		turbulenceModelCH = {'None', 'File', 'Sutherland', 'PowerLaw_AVBP', 'EddyViscosity', 'Boussinesq', 'Boussinesq(xr)' ,'TKE-based'}
		self.window.turbulenceModelM  = tk.OptionMenu(self.window, self.turbulenceModelStr, *turbulenceModelCH,command=self.doNothing)
		self.window.turbulenceModelM.grid         (row =turbulenceModelRow+1 , column = turbulenceModelColumn, rowspan = 1, columnspan = 1)
		if hasattr(self.Case,'TurbulenceModel'):
			self.turbulenceModelStr.set(self.Case.TurbulenceModel)
		else:
			self.turbulenceModelStr.set(AllSettings['TurbulenceModel']['default'])

		# Azimuthal/Crossstreamwise wave number
		mColumn=7
		mRow=1
		self.mL     = tk.Label(self.window, text='Transv. wave no.', font=labelFont)
		self.mL.grid         (row =mRow , column = mColumn, rowspan = 1, columnspan = 1)
		self.mStr    = tk.StringVar()
		self.mE         = tk.Entry(self.window, textvariable = self.mStr, width = 5)
		self.mE.grid         (row =mRow+1 , column = mColumn, rowspan = 1, columnspan = 1)
		if hasattr(self.Case,'m'):
			self.mStr.set(str(self.Case.m))
		else:
			self.mStr.set(AllSettings['m']['default'])

		## For choosing the FlowMode
		#flowModeColumn=7
		#flowModeRow=1
		#self.window.flowModeL  = tk.Label(self.window, text='Flow Mode', font=labelFont())
		#self.window.flowModeL.grid         (row =flowModeRow , column = flowModeColumn, rowspan = 1, columnspan = 1)
		#self.flowModeStr    = tk.StringVar()
		#flowModeCH = {'ColdFlow', 'ColdFlow Variable Density', 'Multi-Species', 'Reacting', 'LowMach','LowMachEnthalpy'}
		#self.window.flowModeM  = tk.OptionMenu(self.window, self.flowModeStr, *flowModeCH,command=self.doNothing)
		#self.window.flowModeM.grid         (row =flowModeRow+1 , column = flowModeColumn, rowspan = 1, columnspan = 1)
		#if hasattr(self.Case,'FlowMode'):
		#	self.flowModeStr.set(self.Case.FlowMode)
		#else:
		#	self.flowModeStr.set(AllSettings['FlowMode']['default'])

		# velocity fluctuations check box
		velFlucColumn=1
		velFlucRow=4
		self.velFlucL     = tk.Label(self.window, text='Vel. Fluc.', font=labelFont)
		self.velFlucL.grid         (row = velFlucRow , column = velFlucColumn, rowspan = 1, columnspan = 1)
		self.velFlucBool = tk.BooleanVar()
		self.velFlucCB = tk.Checkbutton(self.window, text='', variable = self.velFlucBool,\
				command=self.doNothing)
		self.velFlucCB.grid         (row =velFlucRow+1 , column = velFlucColumn, rowspan = 1, columnspan = 1)
		if hasattr(self.Case,'VelFluc'):
			self.velFlucBool.set(self.Case.VelFluc)
		else:
			self.velFlucBool.set(AllSettings['VelFluc']['default'])

		# transverse velocity fluctuations check box
		TransVelFlucColumn=2
		TransVelFlucRow=4
		self.TransVelFlucL     = tk.Label(self.window, text='Trans. Vel. Fluc.', font=labelFont)
		self.TransVelFlucL.grid         (row = TransVelFlucRow , column = TransVelFlucColumn, rowspan = 1, columnspan = 1)
		self.TransVelFlucBool = tk.BooleanVar()
		self.TransVelFlucCB = tk.Checkbutton(self.window, text='', variable = self.TransVelFlucBool,\
				command=self.doNothing)
		self.TransVelFlucCB.grid         (row =TransVelFlucRow+1 , column = TransVelFlucColumn, rowspan = 1, columnspan = 1)
		#if hasattr(self.Case,'TransVelFluc'):
		#	self.TransVelFlucBool.set(self.Case.TransVelFluc)
		#else:
		#	self.TransVelFlucBool.set(AllSettings['TransVelFluc']['default'])
		self.TransVelFlucCB.config(state='disabled')
		if self.coordSysStr.get() == 'Cylindrical':
			self.TransVelFlucBool.set(True)
		else:
			self.TransVelFlucBool.set(False)
		# heat transfer check box
		heatTransferColumn=3
		heatTransferRow=4
		self.heatTransferL     = tk.Label(self.window, text='Heat Transfer', font=labelFont)
		self.heatTransferL.grid         (row = heatTransferRow , column = heatTransferColumn, rowspan = 1, columnspan = 1)
		self.heatTransferBool = tk.BooleanVar()
		self.heatTransferCB = tk.Checkbutton(self.window, text='', variable = self.heatTransferBool,\
				command=self.doNothing)
		self.heatTransferCB.grid         (row =heatTransferRow+1 , column = heatTransferColumn, rowspan = 1, columnspan = 1)
		if hasattr(self.Case,'HeatTransfer'):
			self.heatTransferBool.set(self.Case.HeatTransfer)
		else:
			self.heatTransferBool.set(AllSettings['HeatTransfer']['default'])

		# compressible check box
		compressibleColumn=4
		compressibleRow=4
		self.compressibleL     = tk.Label(self.window, text='Compressible', font=labelFont)
		self.compressibleL.grid         (row = compressibleRow , column = compressibleColumn, rowspan = 1, columnspan = 1)
		self.compressibleBool = tk.BooleanVar()
		self.compressibleCB = tk.Checkbutton(self.window, text='', variable = self.compressibleBool,\
				command=self.doNothing)
		self.compressibleCB.grid         (row =compressibleRow+1 , column = compressibleColumn, rowspan = 1, columnspan = 1)
		if hasattr(self.Case,'Compressible'):
			self.compressibleBool.set(self.Case.Compressible)
		else:
			self.compressibleBool.set(AllSettings['Compressible']['default'])

		# Reaction
		ReactionColumn=5
		ReactionRow=4
		self.ReactionL     = tk.Label(self.window, text='Reaction', font=labelFont)
		self.ReactionL.grid         (row = ReactionRow , column = ReactionColumn, rowspan = 1, columnspan = 1)
		self.ReactionBool = tk.BooleanVar()
		self.ReactionCB = tk.Checkbutton(self.window, text='', variable = self.ReactionBool,\
				command=self.doNothing)
		self.ReactionCB.grid         (row =ReactionRow+1 , column = ReactionColumn, rowspan = 1, columnspan = 1)
		if hasattr(self.Case,'Reaction'):
			self.ReactionBool.set(self.Case.Reaction)
		else:
			self.ReactionBool.set(AllSettings['Reaction']['default'])

		# Mesh File
		MeshFilePathColumn=1
		MeshFilePathRow=6
		self.MeshFilePathB = tk.Button(self.window, text='Load Mesh File', command=self.askForMeshFilePath)
		self.MeshFilePathB.grid         (row = MeshFilePathRow , column = MeshFilePathColumn, rowspan = 1, columnspan = 1)
		self.MeshFilePathStr    = tk.StringVar()
		self.MeshFilePathL   = tk.Label(self.window, textvariable=self.MeshFilePathStr)
		self.MeshFilePathL.grid         (row = MeshFilePathRow , column = MeshFilePathColumn+1, rowspan = 1, columnspan = 1)
		if hasattr(self.Case,'MeshFilePath'):
			self.MeshFilePathStr.set(self.Case.MeshFilePath)
		else:
			self.MeshFilePathStr.set(AllSettings['MeshFilePath']['default'])

		# Mixture File
		MixtureFilePathColumn=1
		MixtureFilePathRow=7
		self.MixtureFilePathB = tk.Button(self.window, text='Load Mixture File', command=self.askForMixtureFilePath)
		self.MixtureFilePathB.grid         (row = MixtureFilePathRow , column = MixtureFilePathColumn, rowspan = 1, columnspan = 1)
		self.MixtureFilePathStr    = tk.StringVar()
		self.MixtureFilePathL   = tk.Label(self.window, textvariable=self.MixtureFilePathStr)
		self.MixtureFilePathL.grid         (row = MixtureFilePathRow , column = MixtureFilePathColumn+1, rowspan = 1, columnspan = 1)
		if hasattr(self.Case,'MixtureFilePath'):
			self.MixtureFilePathStr.set(self.Case.MixtureFilePath)
		else:
			self.MixtureFilePathStr.set(AllSettings['MixtureFilePath']['default'])

		# Species File
		SpeciesFilePathColumn=1
		SpeciesFilePathRow=8
		self.SpeciesFilePathB = tk.Button(self.window, text='Load Species File', command=self.askForSpeciesFilePath)
		self.SpeciesFilePathB.grid         (row = SpeciesFilePathRow , column = SpeciesFilePathColumn, rowspan = 1, columnspan = 1)
		self.SpeciesFilePathStr    = tk.StringVar()
		self.SpeciesFilePathL   = tk.Label(self.window, textvariable=self.SpeciesFilePathStr)
		self.SpeciesFilePathL.grid         (row = SpeciesFilePathRow , column = SpeciesFilePathColumn+1, rowspan = 1, columnspan = 1)
		if hasattr(self.Case,'SpeciesFilePath'):
			self.SpeciesFilePathStr.set(self.Case.SpeciesFilePath)
		else:
			self.SpeciesFilePathStr.set(AllSettings['SpeciesFilePath']['default'])


		# Buttons to close and save and close
		self.window.CancelB     = tk.Button(self.window,text='Cancel', width=10, command=lambda: self.Cancel())
		self.window.CancelB.grid        (row = 10, column = 1, rowspan = 1, columnspan = 1, sticky="E")
		self.window.SaveNCloseB     = tk.Button(self.window,text='Save&Close', width=10, command=lambda: self.SaveNClose(mainGUI))
		self.window.SaveNCloseB.grid        (row = 10, column = 2, rowspan = 1, columnspan = 1, sticky="E")
		self.refresh()

	def ReturnParametersToMain(self,mainGUI):
		'''Passing the settings to the parameters object of the mainGUI
		\t Input:
		\t -mainGUI: Needed to pass the changes in the parameter file'''
		from MixtureClass import MixtureClass
		#from loadSpecies import loadSpecies
		mainGUI.param.Case.nDim=int(self.nDimStr.get())
		mainGUI.param.Case.CoordinateSystem=self.coordSysStr.get()
		#mainGUI.param.Case.FlowMode=self.flowModeStr.get()
		mainGUI.param.Case.TurbulenceModel=self.turbulenceModelStr.get()
		mainGUI.param.Case.m=float(self.mStr.get())
		mainGUI.param.Case.AnalysisMode=self.analysisModeStr.get()
		mainGUI.param.Case.MolViscModel=self.molViscModelStr.get()
		mainGUI.param.Case.MolViscPerturbModel=self.molViscPerturbModelStr.get()
		mainGUI.param.Case.MolVisc=float(self.molViscStr.get())
		mainGUI.param.Case.VelFluc=self.velFlucBool.get()
		mainGUI.param.Case.TransVelFluc=self.TransVelFlucBool.get()
		mainGUI.param.Case.HeatTransfer=self.heatTransferBool.get()
		mainGUI.param.Case.Compressible=self.compressibleBool.get()
		mainGUI.param.Case.Reaction=self.ReactionBool.get()
		mainGUI.param.Case.MeshFilePath=self.MeshFilePathStr.get()
		mainGUI.param.Case.MixtureFilePath=self.MixtureFilePathStr.get()
		mainGUI.param.Case.Mixture=MixtureClass(self.MixtureFilePathStr.get(),self.SpeciesFilePathStr.get())
		mainGUI.param.Case.SpeciesFilePath=self.SpeciesFilePathStr.get()

		if self.__NewMeshFileBool__:
			mainGUI.param.BCs.readDomainData(mainGUI.param.Case.MeshFilePath,
				mainGUI.param.Case.getExtendedTransportedQuantityList())

	def Cancel(self):
		'''Function closing the current window'''
		self.window.destroy()

	def SaveNClose(self,mainGUI):
		''' Function when user pushes Save and Close
		\tInput:
		\t-mainGUI: mainGUI object needed to refresh'''
		self.ReturnParametersToMain(mainGUI)
		self.window.destroy()
		mainGUI.refresh()

	def refresh(self,*args):
		'''Function refreshing the GUI, depending on the current inputs'''
		if self.molViscModelStr.get()=='Constant':
			self.molViscE.config(state='normal')
		else:
			self.molViscE.config(state='disabled')
		if self.coordSysStr.get() == 'Cylindrical':
			self.TransVelFlucBool.set(True)
		else:
			self.TransVelFlucBool.set(False)

	def doNothing(self,*args):
		pass

	def askForMeshFilePath(self):
		'''Function opening a dialog for choosing a mesh file'''
		import os
		tempstring = os.path.relpath(tk.filedialog.askopenfilename(filetypes = (("Mesh files","*.msh"),("all files","*.*"))),self.workDir)
		self.MeshFilePathStr.set(tempstring)
		self.__NewMeshFileBool__=True

	def askForMixtureFilePath(self):
		'''Function opening a dialog for choosing a mixture file'''
		import os
		tempstring = os.path.relpath(tk.filedialog.askopenfilename(filetypes = (("Mixture files","*.mix"),("all files","*.*"))),self.workDir)
		self.MixtureFilePathStr.set(tempstring)

	def askForSpeciesFilePath(self):
		import os
		tempstring = os.path.relpath(tk.filedialog.askopenfilename(filetypes = (("Species files","*.spe"),("all files","*.*"))),self.workDir)
		self.SpeciesFilePathStr.set(tempstring)
