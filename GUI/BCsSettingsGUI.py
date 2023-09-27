import tkinter as tk
from GUI.GUISettings import labelFont
class BCsSettingsGUI():
	def __init__(self,mainGUI):
		'''Initializing the BCs GUI
		\t Input: 
		\t -mainGUI: The FELiCS main GUI object
		'''
		from dolfinx import (plot,
			#MeshTags,
			)
		#import meshio
		import matplotlib.pyplot as plt
		from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
		import numpy as np
		from matplotlib.patches import Polygon
		from mpl_toolkits.mplot3d import Axes3D
		from mpl_toolkits.mplot3d.art3d import Poly3DCollection
		
		#Local libraries and methods
		from itertools import (
			cycle,
			)
		
		self.BCs=mainGUI.param.BCs
		self.Case=mainGUI.param.Case
		self.workDir=mainGUI.workDir
		
		self.window =tk.Toplevel(mainGUI.window)
		self.window.resizable(True,True)
		self.window.title("Boundary Condition Settings")
		
		#mesh=Mesh.ReadMesh(mainGUI.param.Case.MeshFilePath)
		mesh=self.BCs.getMesh()
		#gmsh=meshio.read(mainGUI.param.Case.MeshFilePath)
		#Load BC file
		BCsFilePathColumn=1
		BCsFilePathRow=1
		self.BCsFilePathB = tk.Button(self.window, text='Load BCs File', command=self.askForBCsFilePath)
		self.BCsFilePathB.grid         (row = BCsFilePathRow , column = BCsFilePathColumn, rowspan = 1, columnspan = 2)
		self.BCsFilePathStr    = tk.StringVar()
		self.BCsFilePathL   = tk.Label(self.window, textvariable=self.BCsFilePathStr)
		self.BCsFilePathL.grid         (row = BCsFilePathRow , column = BCsFilePathColumn+2, rowspan = 1,columnspan=5)
		if hasattr(self.BCs,'BCsFilePath'):
			self.BCsFilePathStr.set(self.BCs.BCsFilePath)
		else:
			self.BCsFilePathStr.set(AllSettings['BCsFilePath']['default'])
		# Create BC file
		BCsFilePathColumn=1
		BCsFilePathRow=2
		self.CreateBCsFileB = tk.Button(self.window, text='Create empty BCs File', command=self.createBCsFile)
		self.CreateBCsFileB.grid         (row = BCsFilePathRow , column = BCsFilePathColumn, rowspan = 1, columnspan = 2)

		if self.BCs.dim == 3:
			boundarykey='triangle'	
		else:
			boundarykey='line'
		# Plot the boundaries
####		borderIDX = gmsh.cells[boundarykey]
####		Lines = gmsh.cell_data[boundarykey]['gmsh:physical']
		# get a list of all kinds of BC indices
		self.BCIDs=self.BCs.getBCIDs()
			
		# Create cycled list of colors
		count = 0
		plotColorsCycle = cycle(['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black'])
		plotColorsCycleElement = []
		# Iterate for every boundary

		meshPlot=plt.figure()
		# Plot the mesh
####		if not self.BCs.dim==3:
####			plot(mesh, linewidth=0.1)
		canvas = FigureCanvasTkAgg(meshPlot, master = self.window)
		canvas.draw()
		if self.BCs.dim==3:
			ax=meshPlot.gca(projection='3d')
		toolbar_frame = tk.Frame(self.window)
		toolbar_frame.grid(row=4,column=1,columnspan=10)
		tBar = NavigationToolbar2Tk( canvas, toolbar_frame )	
		tBar.update()
		canvas.get_tk_widget().grid             (row = 3, column = 1, rowspan = 1, columnspan = 10, pady=(20,20), sticky=     "N")
		for i in self.BCIDs:
####			# Get the lines of the loops boundary
####			LinesLocal=list(borderIDX[Lines==i])
####			# Get their coordinates and plot the lines
####			verts=[]
####			
####			# Append next color element to list (required for re-accessing identical BoundaryLabel colors further below)
			plotColorsCycleElement.append(next(plotColorsCycle))
####			
####			for line in LinesLocal:
####				CoordX=[]
####				CoordY=[]
####				CoordZ=[]
####				for point in line:
####					CoordX.append(mesh.coordinates()[point][0])
####					CoordY.append(mesh.coordinates()[point][1])
####					if self.BCs.dim>2:
####						CoordZ.append(mesh.coordinates()[point][2])
####					
####				if self.BCs.dim==2:
####					plt.plot(CoordX,CoordY,color=plotColorsCycleElement[count])
####				elif self.BCs.dim==3:
####					verts.append(list(zip(CoordX, CoordY, CoordZ)))
####			if self.BCs.dim==3:
####				ax.set_xlabel('x')
####				ax.set_ylabel('y')
####				ax.set_zlabel('z')
####				
####				srf = Poly3DCollection(verts, facecolor=plotColorsCycleElement[count])
####				plt.gca().add_collection3d(srf)
####			count += 1
					
		ExtendedSolutionList=mainGUI.param.Case.getExtendedTransportedQuantityList()
####		# Construct a dictionary of Objects
		MatrixRow=5
####		BoundaryLabels={}
####		#Iterate over all boundaries
####		gmshBCDict=gmsh.field_data
####		# Delete all internal regions in mesh (which are not boundaries)
####		for key in list(gmshBCDict.keys()):
####			if gmshBCDict[key][1]>self.BCs.dim-1:
####				del gmshBCDict[key]
####		BCNames=list(gmshBCDict.keys())
####		# The dictionaries NameToID and IDToName allow to reference quickly between the name of a Boundary and its ID, and back...
####		self.BCNameToID={}
####		self.BCIDToName={}
####		for i,BCName in enumerate(BCNames):
####			BoundaryLabels[BCName]=tk.Label(self.window, text='B'+str(gmshBCDict[BCName][0])+': '+BCName, font=labelFont(), highlightthickness=4,  highlightbackground=plotColorsCycleElement[i])
####			BoundaryLabels[BCName].grid(row =MatrixRow , column = (1+BCNames.index(BCName))*2, rowspan = 1, columnspan = 2)
####			self.BCNameToID[BCName]=gmshBCDict[BCName][0]
####			self.BCIDToName[gmshBCDict[BCName][0]]=BCName 
		BCNames =[str(entry) for entry in list(self.BCs.getBCIDs()) ]
		self.BCNameToID={}
		self.BCIDToName={}
		BoundaryLabels={}
		for i,BCName in enumerate(BCNames):
			BoundaryLabels[BCName] = tk.Label(self.window, text='B: '+BCName, font=labelFont(), highlightthickness=4,  highlightbackground=plotColorsCycleElement[i])
			BoundaryLabels[BCName].grid(row =MatrixRow , column = (1+BCNames.index(BCName))*2, rowspan = 1, columnspan = 2)
			self.BCNameToID[BCName]=self.BCs.getBCIDs()[i]
			self.BCIDToName[BCName]=BCName
		# Make a dictionary for the labels of transported quantities
		UnknownLabels={}
		# Write the labels of transportet quantities to the GUI
		for i,Unknown in enumerate(ExtendedSolutionList):
			UnknownLabels[Unknown]=tk.Label(self.window, text=Unknown, font=labelFont())
			UnknownLabels[Unknown].grid(row =MatrixRow+i+1 , column = 1, rowspan = 1, columnspan =1)
				
		#Define choices for Dropdown Menues	
		coordSysCH = {'Dirichlet', 'Neumann'}
		# Define Dictionary(of Dictionaries) for the Dropdown Menues
		self.BCTypeDropdown={}
		# And define a dictionary of dictionaries for the respective StringVariables...
		self.BCTypeStr = {}
		#Define a Dictionary of dictionary for the values of the BCs and the Strin variable
		self.BCValueEntries = {}
		self.BCValueStr= {}
		#Iterate through Unknowns and Boundaries
		for i,Variable in enumerate(ExtendedSolutionList):
			self.BCTypeDropdown[Variable]={}
			self.BCTypeStr[Variable]={}
			self.BCValueEntries[Variable]={}
			self.BCValueStr[Variable]={}
			for k,BCName in enumerate(BCNames):
				ID=self.BCNameToID[BCName]
				self.BCTypeStr[Variable][ID] = tk.StringVar()
				self.BCTypeDropdown[Variable][ID] = tk.OptionMenu(self.window, self.BCTypeStr[Variable][ID], *coordSysCH,command=self.doNothing)
				self.BCTypeDropdown[Variable][ID].grid(row =MatrixRow+i+1 , column = (1+BCNames.index(BCName))*2, rowspan = 1, columnspan = 1)
				self.BCValueStr[Variable][ID] = tk.StringVar()
				self.BCValueEntries[Variable][ID] = tk.Entry(self.window, textvariable = self.BCValueStr[Variable][ID], width = 5)
				self.BCValueEntries[Variable][ID].grid(row =MatrixRow+i+1 , column = (1+BCNames.index(BCName))*2+1, rowspan = 1, columnspan = 1)
			
		if not self.BCs.BCsFilePath=='':
			self.BCs.importBCsDict(self.Case.getExtendedTransportedQuantityList())
		self.BCsDictToGUI(self.BCs.getBCsDict())
		

		
		# Buttons to close and save and close
		self.CancelB     = tk.Button(self.window,text='Cancel', width=10, command=lambda: self.Cancel())
		self.CancelB.grid        (row = 99, column = 1, rowspan = 1, columnspan = 2, sticky="E")		
		self.SaveNCloseB     = tk.Button(self.window,text='Save BCs to File&Close', width=14, command=lambda: self.SaveNClose(mainGUI))
		self.SaveNCloseB.grid        (row = 99, column = 3, rowspan = 1, columnspan = 2, sticky="E")		

	def ReturnParametersToMain(self,mainGUI):
		'''Passing the settings to the parameters object of the mainGUI
		\t Input: 
		\t -mainGUI: Needed to pass the changes in the parameter file'''
		mainGUI.param.BCs.BCsFilePath=self.BCsFilePathStr.get()
		mainGUI.param.BCs.dim=self.BCs.dim
		mainGUI.param.BCs.setBCsDict(self.BCs.getBCsDict())

	def Cancel(self):
		'''Function closing the current window'''
		self.window.destroy()	

	def SaveNClose(self,mainGUI):
		''' Function when user pushes Save and Close'''
		self.sendGUIBCsToBCsObject()
		self.BCs.exportBCs()
		self.ReturnParametersToMain(mainGUI)
		self.window.destroy()	
		mainGUI.refresh()

	def doNothing(self,*args):
		pass

	def askForBCsFilePath(self):
		'''
		Funciton opening a BCs file in the GUI
		'''
		import os
		tempstring = os.path.relpath(tk.filedialog.askopenfilename(filetypes = (("FELiCS BC files","*.bc"),("all files","*    .*"))),self.workDir)
		self.BCsFilePathStr.set(tempstring)
		self.BCs.BCsFilePath =tempstring
		self.BCs.importBCsDict(self.Case.getExtendedTransportedQuantityList())
		self.BCsDictToGUI(self.BCs.getBCsDict())
	

	def createBCsFile(self):
		'''
		Funciton opening a BCs file in the GUI
		'''
		import os
		tempstring = os.path.relpath(tk.filedialog.asksaveasfilename(filetypes = (("FELiCS BC files","*.bc"),("all files","*    .*"))),self.workDir)
		self.BCsFilePathStr.set(tempstring)
		self.BCs.BCsFilePath =tempstring
		self.sendGUIBCsToBCsObject()
		self.BCs.exportBCs()
		
		#self.BCs.importBCsDict(self.Case.getExtendedTransportedQuantityList())
		#self.BCsDictToGUI(self.BCs.getBCsDict())
	
	def BCsDictToGUI(self,BCsDict):
		'''Passing on the entries in the BCsDict to the GUI
		\tInput:
		\t-BCsDict: Dictionary containing all Boundary conditions to be writen to the GUI'''
		#ExtendedSolutionList=self.Case.getExtendedTransportedQuantityList()
		if not self.BCs.BCsFilePath=='':
			for Variable in list(BCsDict.keys()):
				for BC in BCsDict[Variable]:			
					self.BCTypeStr[Variable][BC['ID']].set(BC['type'])
					self.BCValueStr[Variable][BC['ID']].set(BC['value'])
		else:
			for Variable in list(BCsDict.keys()):
				for BC in BCsDict[Variable]:			
					self.BCTypeStr[Variable][BC['ID']].set('Neumann')
					self.BCValueStr[Variable][BC['ID']].set(0.0)
			

	def sendGUIBCsToBCsObject(self):
		'''Function writing the BCs from the GUI to the BCs object (see setBCsFromGUI function in BCsSettingsClass.py'''
		self.BCs.setBCsFromGUI(self.BCTypeStr,self.BCValueStr)
