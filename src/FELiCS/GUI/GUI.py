

class FELiCS_GUI():
	def __init__(self):
		'''Initialiizng the main GUI object'''
		import tkinter as tk
		from os import getcwd
		from FELiCS.parameters import parameters
		# Get the working directory
		self.workDir = getcwd()	
		# Initialize the parameter object
		self.param = parameters()
		# Let's create the Tkinter window.
		self.window = tk.Tk()
		self.window.title("FELiCS")


		self.importSetB     = tk.Button(text='Import Settings from File', width=20, command=self.loadParameters)
		self.importSetB.grid        (row = 1, column = 1, rowspan = 1, columnspan = 1, sticky="E")		
	
		self.exportSetB     = tk.Button(text='Create/Export Settings to File', width=20, command=self.exportParameters)
		self.exportSetB.grid        (row = 1, column = 2, rowspan = 1, columnspan = 1, sticky="E")

		self.caseSetB     = tk.Button(text='Case Settings', width=44, command=self.openCaseSettingsGUI,state=tk.DISABLED)
		self.caseSetB.grid        (row = 2, column = 1, rowspan = 1, columnspan = 2, sticky="E")
                                                                        
		self.baseFlowB     = tk.Button(text='Base Flow Data', width=44, command=self.openBaseFlowGUI,state=tk.DISABLED)
		self.baseFlowB.grid        (row = 3, column = 1, rowspan = 1, columnspan = 2, sticky="E")

		self.BCsB     = tk.Button(text='Boundary Conditions', width=44, command=self.openBCGUI,state=tk.DISABLED)
		self.BCsB.grid        (row = 4, column = 1, rowspan = 1, columnspan = 2, sticky="E")
                                                                        
		self.IOResolventB     = tk.Button(text='Input-Output/Resolvent Settings', width=44, command=self.openIOResolventSettingsGUI,state=tk.DISABLED)
		self.IOResolventB.grid        (row = 5, column = 1, rowspan = 1, columnspan = 2, sticky="E")

		self.NumericsB     = tk.Button(text='Numerics Settings', width=44, command=self.openNumericsSettingsGUI,state=tk.DISABLED)
		self.NumericsB.grid        (row = 6, column = 1, rowspan = 1, columnspan = 2, sticky="E")

		self.OutputB     = tk.Button(text='Output Settings', width=44, command=self.openExportSettingsGUI,state=tk.DISABLED)
		self.OutputB.grid        (row = 7, column = 1, rowspan = 1, columnspan = 2, sticky="E")

		self.RunB     = tk.Button(text='RUN', width=44, command=self.runMain,state=tk.DISABLED)
		self.RunB.grid        (row = 8, column = 1, rowspan = 1, columnspan = 2, sticky="E")
		
		self.param.string='origial'

		self.window.mainloop()
		
		
	def openCaseSettingsGUI(self):
		''' Function opening the case settings GUI'''
		from GUI.CaseSettingsGUI import CaseSettingsGUI
		CaseSettingsGUI(self)

	def openBaseFlowGUI(self):
		''' Function opening the input flow settings GUI'''
		from GUI.FlowInputSettingsGUI import FlowInputSettingsGUI
		FlowInputSettingsGUI(self)

	def openBCGUI(self):
		''' Function opening the boundary condition settings GUI'''
		from GUI.BCsSettingsGUI import BCsSettingsGUI
		BCsSettingsGUI(self)
	
	def openIOResolventSettingsGUI(self):
		''' Function opening the resolvent settings GUI'''
		from GUI.IOResolventSettingsGUI import IOResolventSettingsGUI
		IOResolventSettingsGUI(self)

	def openNumericsSettingsGUI(self):
		''' Function opening the numeric settings GUI'''
		from GUI.NumericsSettingsGUI import NumericsSettingsGUI
		NumericsSettingsGUI(self)
		
	def openExportSettingsGUI(self):
		''' Function opening the export settings GUI'''
		from GUI.ExportSettingsGUI import ExportSettingsGUI
		ExportSettingsGUI(self)
	
	def loadParameters(self):
		''' Function opening a dialog to load a settings file'''
		from os import path
		from tkinter import filedialog

		                # Get the path for the settings file from the user
		self.caseFile=path.relpath(filedialog.askopenfilename(\
			filetypes = (("set files","*.set"),("all files","*.*")),\
			initialdir="./"),\
			self.workDir)
		self.param.importFromFile(self.caseFile)
		self.refresh()
	def exportParameters(self):
		''' Function opening a dialog to choose a file to write'''
		from tkinter import filedialog
		self.caseFile=filedialog.asksaveasfilename(initialdir = "./",title = "Select file",filetypes = (("settings files","*.set"),("all files","*.*")))
		self.param.export(self.caseFile)
		# After parameters is successfully loaded, Case Settings are enabled 
		self.refresh()
	def refresh(self):
		''' Refresh the settings in the main GUI '''
		# First disable all Buttons
		self.caseSetB["state"]="disabled"
		self.baseFlowB["state"] = "disabled"
		self.BCsB["state"] = "disabled"
		self.IOResolventB["state"] = "disabled"
		self.NumericsB["state"] = "disabled"
		self.OutputB["state"]="disabled"
		
		if hasattr(self,'caseFile'):
			self.caseSetB["state"] = "normal"

		# If the case description is complete the remaining setting buttons can be activated.
		if self.param.Case.complete():
			self.baseFlowB["state"] = "normal"
			self.BCsB["state"] = "normal"
			if self.param.Case.AnalysisMode in ['Input-Output','Resolvent']:
				self.IOResolventB["state"] = "normal"
			self.NumericsB["state"] = "normal"
			self.OutputB["state"]="normal"
			# If all settings are complete, then the code is ready to run (IOResolvent is only necessary if Input-Output is complete)
			#if self.param.FlowInput.complete() and\
			#   self.param.BCs.complete([]) and\
			#   (self.param.IOResolvent.complete() or not self.param.Case.AnalysisMode in ['Input-Output','Resolvent']) and\
			#   self.param.Numerics.complete(self.param.Case.getTransportedQuantityList()):
			if self.param.complete():
				self.RunB["state"]= "normal"


#test=FELiCS_GUI()
	
	def runMain(self):
		''' Function initiating the main part of the program '''
		import time
		from runCase import runCase
		
		self.param.getOldParameters()
		useGUI=True
		runCase(self.param,useGUI)	
