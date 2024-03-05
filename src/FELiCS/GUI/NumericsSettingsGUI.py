import tkinter as tk
from tkinter import scrolledtext
from FELiCS.functions import printWarning
from FELiCS.GUI.GUISettings import labelFont
class NumericsSettingsGUI():
	def __init__(self,mainGUI):
		'''Initializing the Numerics GUI
		\t Input: 
		\t -mainGUI: The FELiCS main GUI object
		'''
		self.window =tk.Toplevel(mainGUI.window)
		self.window.Numerics=mainGUI.param.Numerics
		self.window.Case=mainGUI.param.Case
		self.window.title("Numerics Settings")
		AllSettings=self.window.Numerics.getAllSettingsDict()
		
		# Number of Solutions
		nSolutColumn=1
		nSolutRow=1
		self.window.nSolutFrame = tk.LabelFrame(self.window, text="Number of solutions",font=labelFont())
		self.window.nSolutFrame.grid        (row = nSolutRow, column = nSolutColumn, rowspan = 1, columnspan = 1, sticky="")
		self.window.nSolutStr    = tk.StringVar()
		self.window.nSolutE  = tk.Entry(self.window.nSolutFrame, textvariable = self.window.nSolutStr, width = 5)
		self.window.nSolutE.pack()
		self.window.nSolutStr.set(self.window.Numerics.nSolut)
		#Eigen Value Guess
		EigenValueColumn=2
		EigenValueRow=1
		self.window.EigenValueGuessFrame = tk.LabelFrame(self.window,
			text="Eigen value guess",
			font=labelFont())
		self.window.EigenValueGuessFrame.grid        (row = EigenValueRow, column = EigenValueColumn, rowspan = 1, columnspan = 2, sticky="")
		self.window.EigenValueGuessStr    = tk.StringVar()
		self.window.EigenValueGuessE  = tk.Entry(self.window.EigenValueGuessFrame, textvariable = self.window.EigenValueGuessStr, width = 5)
		self.window.EigenValueGuessE.pack()
		self.window.EigenValueGuessStr.set(
			self.__listToStr(self.window.Numerics.EigenValueGuess)	
			)
		#nCPU
		nCPUColumn=4
		nCPURow=1
		self.window.nCPUFrame = tk.LabelFrame(self.window,
			text="Number of Threads",
			font=labelFont())
		self.window.nCPUFrame.grid        (row = nCPURow, column = nCPUColumn, rowspan = 1, columnspan = 2, sticky="")
		self.window.nCPUStr    = tk.StringVar()
		self.window.nCPUE  = tk.Entry(self.window.nCPUFrame, textvariable = self.window.nCPUStr, width = 5)
		self.window.nCPUE.pack()
		self.window.nCPUStr.set(self.window.Numerics.nCPU)

		# Set frame for the discretization
		DiscretizationColumn=1
		DiscretizationRow=2
		self.window.DiscretizationFrame = tk.LabelFrame(self.window, text="Spatial Discretization",font=labelFont())
		self.window.DiscretizationFrame.grid        (row = DiscretizationRow, column = DiscretizationColumn, rowspan = 2, columnspan = 2, sticky="")
	
		# set labels for variables
		variableLableDict={}
		variableList = mainGUI.param.Case.getTransportedQuantityList()
		for variable in variableList:
			variableLableDict[variable] = tk.Label(self.window.DiscretizationFrame, text=variable, font=labelFont())
			variableLableDict[variable].grid(row = 2+variableList.index(variable) , column = 1, rowspan = 1, columnspan = 1)
		# set label for scheme in the Discretization frame
		SchemeLabel = tk.Label(self.window.DiscretizationFrame, text='Scheme', font=labelFont())
		SchemeLabel.grid(row = 1,
			column = 2,
			rowspan = 1,
			columnspan = 1)
		# set label for the polynomial order in the Discretization frame
		PolynomialOrderLabel = tk.Label(self.window.DiscretizationFrame, text='Polynomial\nOrder', font=labelFont())
		PolynomialOrderLabel.grid(row = 1,
			column = 3,
			rowspan = 1,
			columnspan = 1)
		# Set menue and entry widgets for the scheme and the polynomial order
		self.window.SchemeMenuDict={}
		self.window.SchemeStrDict={}
		self.window.PolynomialOrderEntryDict={}
		self.window.PolynomialOrderStrDict={}
		# Currently the only choice for the scheme is continuous Galerkin
		SchemeChoices = {'CG'}
		for variable in variableList:
			self.window.SchemeStrDict[variable] = tk.StringVar()
			self.window.SchemeMenuDict[variable] = tk.OptionMenu(self.window.DiscretizationFrame,
				self.window.SchemeStrDict[variable],
				*SchemeChoices,
				command=self.doNothing)
			self.window.SchemeMenuDict[variable].grid(row =2+variableList.index(variable),
				column = 2,
				rowspan = 1,
				columnspan = 1)
			if variable in list(self.window.Numerics.Schemes.keys()):
				self.window.SchemeStrDict[variable].set(self.window.Numerics.Schemes[variable])
			else:
				self.window.SchemeStrDict[variable].set('CG')
			
			self.window.PolynomialOrderStrDict[variable] = tk.StringVar()
			self.window.PolynomialOrderEntryDict[variable] = tk.Entry(self.window.DiscretizationFrame,
				textvariable=self.window.PolynomialOrderStrDict[variable],
				width=2)
			self.window.PolynomialOrderEntryDict[variable].grid(row =2+variableList.index(variable),
				 column = 3,
				 rowspan = 1,
				 columnspan = 1)
			if variable in list(self.window.Numerics.PolynomialOrder.keys()):
				self.window.PolynomialOrderStrDict[variable].set(str(self.window.Numerics.PolynomialOrder[variable]))
			else:
				self.window.PolynomialOrderStrDict[variable].set('2')
		
		# Set frame for the Linear Algebra solver to use
		LinearAlgebraSolverColumn=3
		LinearAlgebraSolverRow=2
		self.window.LinearAlgebraSolverFrame = tk.LabelFrame(self.window, text="Lin. Algebra Solver",font=labelFont())
		self.window.LinearAlgebraSolverFrame.grid        (row = LinearAlgebraSolverRow, column = LinearAlgebraSolverColumn, rowspan = 1, columnspan = 2, sticky="")
		self.window.LinearAlgebraSolverStr = tk.StringVar()
		LinearAlgebraSolverChoices={'python','matlab','matrix export'}
		self.window.LinearAlgebraSolverMenu = tk.OptionMenu(self.window.LinearAlgebraSolverFrame,
			self.window.LinearAlgebraSolverStr,
			*LinearAlgebraSolverChoices,
			command=self.doNothing)
		self.window.LinearAlgebraSolverMenu.pack()
		self.window.LinearAlgebraSolverStr.set(self.window.Numerics.LinearAlgebraSolver)
		
		# Set frame for the preconditioning method
		PreconditionerColumn=3
		PreconditionerRow=3
		self.window.PreconditionerFrame = tk.LabelFrame(self.window, text="Preconditioner",font=labelFont())
		self.window.PreconditionerFrame.grid        (row = PreconditionerRow, column = PreconditionerColumn, rowspan = 1, columnspan = 2, sticky="")
		self.window.PreconditionerStr = tk.StringVar()
		PreconditionerChoices={'None','Sum of row'}
		self.window.PreconditionerMenu = tk.OptionMenu(self.window.PreconditionerFrame,
			self.window.PreconditionerStr,
			*PreconditionerChoices,
			command=self.doNothing)
		self.window.PreconditionerMenu.pack()
		self.window.PreconditionerStr.set(self.window.Numerics.Preconditioner)


		# Buttons to close and save and close
		self.window.CancelB     = tk.Button(self.window,text='Cancel', width=10, command=lambda: self.Cancel())
		self.window.CancelB.grid        (row = 99, column = 1, rowspan = 1, columnspan = 2, sticky="")		
		self.window.SaveNCloseB     = tk.Button(self.window,text='Save&Close', width=10, command=lambda: self.SaveNClose(mainGUI))
		self.window.SaveNCloseB.grid        (row = 99, column = 3, rowspan = 1, columnspan = 2, sticky="")		
		self.refresh()	

	def ReturnParametersToMain(self,mainGUI):
		'''Passing the settings to the parameters object of the mainGUI
		\t Input: 
		\t -mainGUI: Needed to pass the changes in the parameter file'''
		mainGUI.param.Numerics.nSolut=int(self.window.nSolutStr.get())
		mainGUI.param.Numerics.EigenValueGuess = self.__strToListOfFloat(self.window.EigenValueGuessStr.get())
		mainGUI.param.Numerics.nCPU=int(self.window.nCPUStr.get())
		Schemes={}
		for key in list(self.window.SchemeStrDict.keys()):
			Schemes[key]=self.window.SchemeStrDict[key].get()
		mainGUI.param.Numerics.Schemes=Schemes
		PolynomialOrder={}
		for key in list(self.window.PolynomialOrderStrDict.keys()):
			PolynomialOrder[key]=int(self.window.PolynomialOrderStrDict[key].get())
		mainGUI.param.Numerics.PolynomialOrder=PolynomialOrder
		mainGUI.param.Numerics.LinearAlgebraSolver=self.window.LinearAlgebraSolverStr.get()
		mainGUI.param.Numerics.Preconditioner=self.window.PreconditionerStr.get()

	def __listToStr(self,inputList):
		'''Function converting a list's entries to a comma seperated string '''
		tempString = ''
		if isinstance(inputList, list):
			if inputList:
				for item in inputList:
					tempString+=str(item)+','
				tempString=tempString[:-1]
		elif isinstance(inputList, float):
			tempString = str(inputList)
		else:
			printWarning('type of inputList undefined...')
		return tempString

	def __strToListOfFloat(self,tempString):
		''' Separates a string by comma and makes it a list of floats'''
		tempList=tempString.split(',')
		outputList=[]
		for entry in tempList:
			if not entry=='':
				outputList.append(float(entry))
		return outputList
		
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
		pass

	def doNothing(self,*args):
		pass

