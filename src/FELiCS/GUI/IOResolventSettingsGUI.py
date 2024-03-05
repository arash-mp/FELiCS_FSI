import tkinter as tk
from tkinter import scrolledtext

from FELiCS.GUI.GUISettings import labelFont
class IOResolventSettingsGUI():
	def __init__(self,mainGUI):
		'''Initializing the Input-Output/Resolvent settings GUI
		\t Input: 
		\t -mainGUI: The FELiCS main GUI object
		'''
		self.IOResolvent=mainGUI.param.IOResolvent
		self.Case=mainGUI.param.Case
		self.workDir=mainGUI.workDir
		self.window =tk.Toplevel(mainGUI.window)
		self.window.title("Input-Output/Resolvent Settings")
		AllSettings=self.IOResolvent.getAllSettingsDict()
		
		# Listing all omegas
		omegasColumn=1
		omegasRow=1
		self.window.nDimL  = tk.Label(self.window, text='Omegas', font=labelFont())
		self.window.nDimL.grid         (row =omegasRow , column = omegasColumn, rowspan = 1, columnspan = 1)
		self.window.omegaTable = scrolledtext.ScrolledText(self.window,  
                                      wrap = tk.WORD,  
                                      width = 8,  
                                      height = 10,  
                                              ) 
		self.window.omegaTable.grid(column = omegasColumn,row=omegasRow+1, pady = 10, padx = 8,rowspan=10)
		self.window.omegaTable.insert('end', '', 'tag-center')
		self.window.omegaTable.delete("1.0", tk.END)
		if hasattr(self.IOResolvent,'Omegas'):
			self.listToOmegas(self.IOResolvent.Omegas)
		else:
			self.listToOmegas(AllSettings['Omegas']['default'])

		self.window.omegaTable.configure(state ='disabled')

		# Single Frequency
		fColumn=2
		fRow=1
		self.window.fL  = tk.Label(self.window, text='Omegas', font=labelFont())
		self.window.fL.grid         (row =fRow , column = fColumn, rowspan = 1, columnspan = 1)
		self.window.fStr    = tk.StringVar()
		self.window.fE         = tk.Entry(self.window, textvariable = self.window.fStr, width = 5)
		self.window.fE.grid         (row =fRow , column = fColumn+1, rowspan = 1, columnspan = 1)
		
		self.window.addSingleFreqB     = tk.Button(self.window,text='Add single omega', width=12, command=lambda: self.addSingleFrequency())
		self.window.addSingleFreqB.grid        (row = fRow+1, column = fColumn, rowspan = 1, columnspan = 2, sticky="W")		

		# Multi Frequency Min
		fMinColumn=2
		fMinRow=3
		self.window.fMinL  = tk.Label(self.window, text='omega_min', font=labelFont())
		self.window.fMinL.grid         (row =fMinRow , column = fMinColumn, rowspan = 1, columnspan = 1)
		self.window.fMinStr    = tk.StringVar()
		self.window.fMinE         = tk.Entry(self.window, textvariable = self.window.fMinStr, width = 5)
		self.window.fMinE.grid         (row =fMinRow , column = fMinColumn+1, rowspan = 1, columnspan = 1)
		# Multi Frequency Step
		fStepColumn=2
		fStepRow=4
		self.window.fStepL  = tk.Label(self.window, text='omega_step', font=labelFont())
		self.window.fStepL.grid         (row =fStepRow , column = fStepColumn, rowspan = 1, columnspan = 1)
		self.window.fStepStr    = tk.StringVar()
		self.window.fStepE         = tk.Entry(self.window, textvariable = self.window.fStepStr, width = 5)
		self.window.fStepE.grid         (row =fStepRow , column = fStepColumn+1, rowspan = 1, columnspan = 1)
		# MultiFrequency Max
		fMaxColumn=2
		fMaxRow=5
		self.window.fMaxL  = tk.Label(self.window, text='omega_max', font=labelFont())
		self.window.fMaxL.grid         (row =fMaxRow , column = fMaxColumn, rowspan = 1, columnspan = 1)
		self.window.fMaxStr    = tk.StringVar()
		self.window.fMaxE         = tk.Entry(self.window, textvariable = self.window.fMaxStr, width = 5)
		self.window.fMaxE.grid         (row =fMaxRow , column = fMaxColumn+1, rowspan = 1, columnspan = 1)
		
		#Add multi Freq button
		self.window.addSingleFreqB     = tk.Button(self.window,text='Add multiple omega', width=12, command=lambda: self.addMultipleFrequencies())
		self.window.addSingleFreqB.grid        (row = fMaxRow+1, column = fColumn, rowspan = 1, columnspan = 2, sticky="W")		

		#Add multi Freq button
		self.window.addSingleFreqB     = tk.Button(self.window,text='Remove all omega', width=12, command=lambda: self.removeAllFrequencies())
		self.window.addSingleFreqB.grid        (row = fMaxRow+2, column = fColumn, rowspan = 1, columnspan = 2, sticky="W")		
		# Forcing Mode (Body or Boundary)
		ForcingModeColumn=4
		ForcingModeRow=1
		self.window.ForcingModeL  = tk.Label(self.window, text='Forcing\nMode', font=labelFont())
		self.window.ForcingModeL.grid         (row =ForcingModeRow , column = ForcingModeColumn, rowspan = 1, columnspan = 1)
		self.window.ForcingModeStr    = tk.StringVar()
		ForcingModeCH = {'Body','Boundary'}
		self.window.ForcingModeM  = tk.OptionMenu(self.window, self.window.ForcingModeStr, *ForcingModeCH,command=self.doNothing)
		self.window.ForcingModeM.grid         (row =ForcingModeRow+1 , column = ForcingModeColumn, rowspan = 1, columnspan = 1)
		if hasattr(self.IOResolvent,'ForcingMode'):
			self.window.ForcingModeStr.set(str(self.IOResolvent.ForcingMode))
		else:
			self.window.ForcingModeStr.set(AllSettings['ForcingMode']['default'])
		
		# Forcing Boundary
		ForcingBoundaryIndicesColumn=5
		ForcingBoundaryIndicesRow=1
		self.window.ForcingBoundaryIndicesL  = tk.Label(self.window, text='Forcing\nBoundary\nIndex', font=labelFont())
		self.window.ForcingBoundaryIndicesL.grid         (row =ForcingBoundaryIndicesRow , column = ForcingBoundaryIndicesColumn, rowspan = 1, columnspan = 3)
		self.window.ForcingBoundaryIndicesStr    = tk.StringVar()
		self.window.ForcingBoundaryIndicesEntry  = tk.Entry(self.window,
								   textvariable = self.window.ForcingBoundaryIndicesStr,
								    width = 5) 
							 
		self.window.ForcingBoundaryIndicesEntry.grid(row =ForcingBoundaryIndicesRow+1,
							     column = ForcingBoundaryIndicesColumn,
							     rowspan = 1, 
							     columnspan = 3)

		self.listToBoundaryIndices(self.IOResolvent.ForcingBoundaryIndices)

		# Define export frame
		NormFrameColumn=4
		NormFrameRow=3
		ExportFrameColumn=1
		ExportFrameRow=1

		self.window.NormFrame = tk.LabelFrame(self.window, text="Resolvent norms",font=labelFont())
		self.window.NormFrame.grid        (row = NormFrameRow, column = NormFrameColumn, rowspan = 10, columnspan = 1, sticky="")
		# Labels of variables for norms
		variableLabelColumn=1
		variableLabelRow=2
		variableLableDict={}
		variableList = mainGUI.param.Case.getExtendedTransportedQuantityList()
		for variable in variableList:
			variableLableDict[variable] = tk.Label(self.window.NormFrame, text=variable, font=labelFont())
			variableLableDict[variable].grid(row = variableLabelRow+1+variableList.index(variable) , column = variableLabelColumn, rowspan = 1, columnspan = 1)
		# Forcing Norm



		forcingNormColumn=variableLabelColumn+1
		forcingNormRow=variableLabelRow
		responseNormRow=forcingNormRow
		responseNormColumn=forcingNormColumn+1
		forcingNormL     = tk.Label(self.window.NormFrame, text='Forcing', font=labelFont())
		forcingNormL.grid         (row = forcingNormRow , column = forcingNormColumn, rowspan = 1, columnspan = 1)
		responseNormL     = tk.Label(self.window.NormFrame, text='Response', font=labelFont())
		responseNormL.grid         (row = responseNormRow , column = responseNormColumn, rowspan = 1, columnspan = 1)

		self.window.forcingNormCheckBoxDict={}
		self.window.forcingNormBoolDict={}
		self.window.responseNormCheckBoxDict={}
		self.window.responseNormBoolDict={}
		for variable in variableList:
			self.window.forcingNormBoolDict[variable]=tk.BooleanVar()
			self.window.forcingNormCheckBoxDict[variable]=tk.Checkbutton(self.window.NormFrame,
				text='',
				variable = self.window.forcingNormBoolDict[variable],
				command=self.doNothing)
			self.window.forcingNormCheckBoxDict[variable].grid(row =forcingNormRow+1+variableList.index(variable) ,
				column = forcingNormColumn,
				rowspan = 1,
				columnspan = 1)
			if variableList.index(variable) in self.IOResolvent.ForcingCoeff:
				self.window.forcingNormBoolDict[variable].set(True)
			else:
				self.window.forcingNormBoolDict[variable].set(False)
			self.window.responseNormBoolDict[variable]=tk.BooleanVar()
			self.window.responseNormCheckBoxDict[variable]=tk.Checkbutton(self.window.NormFrame,
				text='',
				variable = self.window.responseNormBoolDict[variable],
				command=self.doNothing)
			self.window.responseNormCheckBoxDict[variable].grid(row =responseNormRow+1+variableList.index(variable) ,
				column = responseNormColumn,
				rowspan = 1,
				columnspan = 1)
			if variableList.index(variable) in self.IOResolvent.ResponseCoeff:
				self.window.responseNormBoolDict[variable].set(True)
			else:
				self.window.responseNormBoolDict[variable].set(False)


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
		mainGUI.param.IOResolvent.Omegas=self.omegasToList()
		mainGUI.param.IOResolvent.ForcingMode=self.window.ForcingModeStr.get()
		# Arrange list of forcing coefficients
		variableList= self.Case.getExtendedTransportedQuantityList()
		ForcingCoeff=[]
		for key in list(self.window.forcingNormBoolDict.keys()):
			if self.window.forcingNormBoolDict[key].get():
				ForcingCoeff.append(variableList.index(key))
		mainGUI.param.IOResolvent.ForcingCoeff=ForcingCoeff
		# Arrange list of response coefficients
		ResponseCoeff=[]
		for key in list(self.window.responseNormBoolDict.keys()):
			if self.window.responseNormBoolDict[key].get():
				ResponseCoeff.append(variableList.index(key))
		mainGUI.param.IOResolvent.ResponseCoeff=ResponseCoeff
		mainGUI.param.IOResolvent.ForcingBoundaryIndices=self.boundaryIndicesToList()
		
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
		pass

	def doNothing(self,*args):
		pass

	def addSingleFrequency(self):
		'''Adding a single frequency (from the fStr-object) to the frequency list (using the function listToOmegas)'''
		FreqList=list(self.window.fStr.get().split(','))
		for i in range(len(FreqList)):
			FreqList[i]=float(FreqList[i])
		self.listToOmegas(FreqList)

	def addMultipleFrequencies(self):
		''' Writing a the valoes between fMin and fMax to the omegas List in the GUI'''
		from numpy import arange
		fMin=float(self.window.fMinStr.get())
		fStep=float(self.window.fStepStr.get())
		fMax=float(self.window.fMaxStr.get())
		FreqArray=arange(fMin,fMax,fStep)
		# Round all entries in the array to the 14th digit, to correct floating point error
		for i in range(len(FreqArray)):
			FreqArray[i]=round(FreqArray[i],14)
		FreqList=FreqArray.tolist()
		self.listToOmegas(FreqList)

	def listToOmegas(self,listToAdd):
		'''Adding a list of floats to the omegaTable'''
		oldList=self.omegasToList()
		listToAdd=list(set(listToAdd)-set(oldList))
		listToAdd.extend(oldList)
		listToAdd.sort()
		textToInsert=''
		for entry in listToAdd:
			textToInsert += str(entry)+'\n'
		
		self.window.omegaTable.configure(state ='normal')
		self.window.omegaTable.delete("0.0", tk.END)
		self.window.omegaTable.insert(tk.INSERT, textToInsert)
		self.window.omegaTable.configure(state ='disabled')

	def removeAllFrequencies(self):
		'''Clear all frequencies from the table'''	
		textToInsert=''
		self.window.omegaTable.configure(state ='normal')
		self.window.omegaTable.delete("0.0", tk.END)
		self.window.omegaTable.insert(tk.INSERT, textToInsert)
		self.window.omegaTable.configure(state ='disabled')
		text=self.window.omegaTable.get("1.0",tk.END)
	def omegasToList(self):
		''' Function reading all entries in the omega table and returning a list containing all omegas'''
		text=self.window.omegaTable.get("1.0",tk.END)
		exportList=list(text.split('\n'))
		exportList= [i for i in exportList if not i == '']
		for i in range(len(exportList)):
			exportList[i]=float(exportList[i])
		return exportList
	def listToBoundaryIndices(self,inputList):
		'''Function writing a inputList of integers to the Boundary List Entry '''
		tempString = ''
		if inputList:
			for item in inputList:
				tempString+=str(item)+','
			tempString=tempString[:-1]
		self.window.ForcingBoundaryIndicesStr.set(tempString)
	def boundaryIndicesToList(self):
		''' returning the entries of BoundaryIndices as a list of integers'''
		tempString = self.window.ForcingBoundaryIndicesStr.get()
		tempList=tempString.split(',')
		outputList=[]
		for entry in tempList:
			if not entry=='':
				outputList.append(int(entry))
		return outputList
