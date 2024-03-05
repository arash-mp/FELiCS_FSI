import tkinter as tk
from FELiCS.GUI.GUISettings import labelFont
class FlowInputSettingsGUI():
	def __init__(self,mainGUI):
		'''Initializing the FlowInput GUI
		\t Input: 
		\t -mainGUI: The FELiCS main GUI object
		'''
		import numpy as np
		self.workDir=mainGUI.workDir
		self.window =tk.Toplevel(mainGUI.window)
		self.FlowInput=mainGUI.param.FlowInput
		self.window.title("Mean Flow Settings")

		# MeanFlow File
		MeanFlowFilePathColumn=1
		MeanFlowFilePathRow=1
		self.MeanFlowFilePathB = tk.Button(self.window, text='Load MeanFlow File', command=self.askForMeanFlowFilePath)
		self.MeanFlowFilePathB.grid         (row = MeanFlowFilePathRow , column = MeanFlowFilePathColumn, rowspan = 1, columnspan = 1)
		self.MeanFlowFilePathStr    = tk.StringVar()
		self.MeanFlowFilePathL   = tk.Label(self.window, textvariable=self.MeanFlowFilePathStr)
		self.MeanFlowFilePathL.grid         (row = MeanFlowFilePathRow , column = MeanFlowFilePathColumn+1, rowspan = 1, columnspan = 1)
		self.MeanFlowFilePathStr.set(self.FlowInput.MeanFlowFilePath)

		# Averaging Direction
		AveragingDirectionColumn=1
		AveragingDirectionRow=2
		self.AveragingDirectionL = tk.Label(self.window,
							   text='Averaging Mode', 
							   font = labelFont())
		self.AveragingDirectionL.grid(row = AveragingDirectionRow,
						     column = AveragingDirectionColumn,
						     rowspan = 1, 
						     columnspan = 1)
		self.AveragingDirectionStr = tk.StringVar()
		AveragingDirectionChoices = {'None','Azimuthal'}
		self.AveragingDirectionMenu = tk.OptionMenu(self.window,
								   self.AveragingDirectionStr,
								   *AveragingDirectionChoices,
							           command=self.doNothing)
		self.AveragingDirectionMenu.grid(row = AveragingDirectionRow+1,
							column = AveragingDirectionColumn,
							rowspan = 1,
							columnspan = 1)
		self.AveragingDirectionStr.set(self.FlowInput.AveragingDirection)
		# Averaging Axis
		AveragingAxisColumn=2
		AveragingAxisRow=2
		self.AveragingAxisL = tk.Label( self.window,
						text='Averaging Axis', 
						font = labelFont() )
		self.AveragingAxisL.grid(row = AveragingAxisRow,
						 column = AveragingAxisColumn,
						 rowspan = 1,
						 columnspan = 1)
		self.AveragingAxisStr = tk.StringVar()
		AveragingAxisChoices = {'x','y','z'}
		self.AveragingAxisMenu = tk.OptionMenu(self.window,
								   self.AveragingAxisStr,
								   *AveragingAxisChoices,
							           command=self.doNothing)
		self.AveragingAxisMenu.grid(row = AveragingAxisRow+1,
						 column = AveragingAxisColumn,
						 rowspan = 1,
						 columnspan = 1)
		self.AveragingAxisStr.set(self.FlowInput.AveragingAxis)

		# Buttons to close and save and close
		self.CancelB     = tk.Button(self.window,text='Cancel', width=10, command=lambda: self.Cancel())
		self.CancelB.grid        (row = 99, column = 1, rowspan = 1, columnspan = 1, sticky="E")		
		self.SaveNCloseB     = tk.Button(self.window,text='Save&Close', width=10, command=lambda: self.SaveNClose(mainGUI))
		self.SaveNCloseB.grid        (row = 99, column = 2, rowspan = 1, columnspan = 1, sticky="E")		

	def ReturnParametersToMain(self,mainGUI):
		'''Passing the settings to the parameters object of the mainGUI
		\t Input: 
		\t -mainGUI: Needed to pass the changes in the parameter file'''
		mainGUI.param.FlowInput.MeanFlowFilePath=self.MeanFlowFilePathStr.get()
		mainGUI.param.FlowInput.AveragingDirection=self.AveragingDirectionStr.get()
		mainGUI.param.FlowInput.AveragingAxis=self.AveragingAxisStr.get()
		
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
	def doNothing(self,*args):
		pass

	def askForMeanFlowFilePath(self):
		import os
		tempstring = os.path.relpath(tk.filedialog.askopenfilename(filetypes = (("FELiCS files","*.fel *.hdf5 *.mat"),("all files","*    .*"))),self.workDir)
		self.MeanFlowFilePathStr.set(tempstring)

