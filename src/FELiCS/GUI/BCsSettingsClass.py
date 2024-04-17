from os.path import isfile
from os import system, path

from dolfinx import __version__
from dolfinx.mesh import (
	locate_entities,
	locate_entities_boundary,
	meshtags,
	meshtags_from_entities,
	Mesh,
	)
from dolfinx.io import gmshio
from mpi4py import MPI

from ufl import triangle

import gmsh
from FELiCS.GUI.SettingsClass import Settings
import h5py
import numpy as np
import pdb
from FELiCS.tensorUtils import CoordinateSystem
from ufl import SpatialCoordinate
from FELiCS.functions import printDeprecatedWarning

class FELiCSMesh(Mesh):
	'''
	This class is an extension to the fenics mesh class
	'''
	def __init__(
			self,
			coordinateSystem,
			filename=None,
			gdim=0,
                        m=0,
			inputMesh=None,
			):

		if inputMesh is None:
			gmsh.initialize()
			if __version__.find('0.4') >= 0:
				printDeprecatedWarning("Dolfinx version <0.5.0 is used.")
				from FELiCSGUI.gmsh_helpers import gmsh_model_to_mesh, read_from_msh
				mesh, cell_tags, hi, facet_tags = read_from_msh(filename, cell_data=True, facet_data=True, gdim=gdim)
				self.coordinatesGMSH = extract_gmsh_geometry(gmsh.model)
			else:
				gmsh.open(filename)
				from dolfinx.io import gmshio
				mesh_comm = MPI.COMM_WORLD
				model_rank = 0
				mesh, _, facet_tags = gmshio.model_to_mesh(gmsh.model, mesh_comm, model_rank, gdim=gdim)

			#Mesh.__init__(self, MPI.COMM_WORLD, mesh.topology, mesh.geometry, mesh.ufl_domain())
			try:    #try new version of dolfinx 
				Mesh.__init__(self,  mesh, mesh.ufl_domain())
			except: #use old language 
				printDeprecatedWarning("Mesh module from dolfinx version <0.7.0 is used.")
				Mesh.__init__(self, MPI.COMM_WORLD, mesh.topology, mesh.geometry, mesh.ufl_domain())
				#Mesh.__init__(self, MPI.COMM_WORLD, mesh.topology, mesh.geometry)

			self.dolfinxMesh = mesh
			self.facet_tags = facet_tags
			self.gdim = gdim
			self._ufl_domain = mesh._ufl_domain
			self.calcConnectivity()

			# save the coordinates in gmsh order:
			gmsh.open(filename)


		#
		else:
			try:    #try new version of dolfinx 
				Mesh.__init__(self, inputMesh, inputMesh.ufl_domain())
			except: #use old language 
				printDeprecatedWarning("Mesh module from dolfinx version <0.7.0 is used.")
				Mesh.__init__(self, MPI.COMM_WORLD, inputMesh.topology, inputMesh.geometry, inputMesh.ufl_domain())
			self.gdim = inputMesh.topology.dim
			self.dolfinxMesh = inputMesh
        
		x = SpatialCoordinate(self)
		# Define tensor coordinate system, we always assume the third dimension to be homogenous
		if coordinateSystem =='Cartesian':
		    self.__coordinateSystem = CoordinateSystem(
                                    x, 
                                    coordinateSystem.lower(), 
                                    m = m,
                                    mesh_dims = (1, 1, 0),
                                    )
		elif coordinateSystem =='Cylindrical':
		    self.__coordinateSystem = CoordinateSystem(
                                    x,
                                    "cylindricalfelics", 
                                    m = m,
                                    mesh_dims = (1, 1, 0),
                                    )
		else:
		    printError('Coord. syst not yet implemented in tensor framework.')
		self._coordinates = self.coordinates()


	def saveInFELiCSFormat(self, filename):
		'''
		This function saves the computational mesh in the FELiCS format

		Function arguments:
		- filename: The path where to save the mesh

		Function returns:
		'''
		coordinates = self.coordinates()
		self.calcConnectivity()
		meshCells = self.meshCells
		nDim = coordinates.shape[1]
		coordinateNames = ['x']
		if nDim > 1:
			coordinateNames.append('y')
		if nDim > 2:
			coordinateNames.append('z')
		hf = h5py.File(filename, 'w')
		g1 = hf.create_group('coordinates')
		for i_coordinateName,coordinateName in enumerate(coordinateNames):
			g1.create_dataset(coordinateName,data=coordinates[:,i_coordinateName])
		g2 = hf.create_group('cells')
		g2.create_dataset('triangles',data=np.array(meshCells))
		hf.close()

	def calcConnectivity(self):
		"""
		this method calculates the meshCells array in the fenics representation
		"""
		connectivityCells = self.topology.connectivity(2, 0)


		try:    #try new version of dolfinx 
			self.meshCells = connectivityCells.array.reshape(
			    [self.topology.original_cell_index.shape[0], self.topology.cell_types[0].value])
		except: #use old language. TODO: handle DEPRECATED stuff uniformly
			printDeprecatedWarning("Mesh module from dolfinx version <0.7.0 is used.")
			self.meshCells = connectivityCells.array.reshape(
				[self.topology.original_cell_index.shape[0], self.topology.cell_type.value])

	def cells(self):
		"""
		this methods returns the cell-connectivity information
		"""
		self.calcConnectivity()
		return self.meshCells


	def coordinates(self):
		"""
		This method acts as a getter-method for the vertex-coordinates.
		"""
		return self.geometry.x[:, 0:self.gdim]

	@property
	def coordinateSystem(self):
		return self.__coordinateSystem

class BCsSettingsClass(Settings):
	def __init__(self,MeshFilePath):
		'''Initializng the BCs Class
		\tInput:
		\t\t-MeshFilePath: Path of the file, needed to obtain boundaries
		'''
		super().__init__()
		self._settingsKind = 'BCs'
		SettingsDict=self.getAllSettingsDict()
		self.__BCIDs__=[]
		for key in list(SettingsDict.keys()):
			tempStr='self.'+key+'=\"'+str(SettingsDict[key]['default'])+'\"'
			exec(tempStr)

	def getAllSettingsDict(self):
		'''Function returning all boundary condition settings with default values'''
		SettingsDict={\
			'BCsFilePath':{'datatype':str,'default':''},
			'dim':{'datatype':int,'default':0}
		}
		return SettingsDict

	def exportBCs(self):
		''' Exports the BCDict to file filename'''

		print('Saving BCs...')
		if not self.BCsFilePath == '':
			BCsFile=open(self.BCsFilePath,'w')
			BCsFile.write(str(self.__BCsDict__))

	def importSettings(self,settingFilePath):
		''' Loading BCs parameters from file, without checking the consistency of BCs and mesh. To read the BCs with checking, use importBCsDict() '''
		SettingsDict=self.getAllSettingsDict()
		if not settingFilePath =='':
			file = open(settingFilePath)

			#Read whole file
			lines = file.readlines()
			# Add every line of the file as an attribute to the object
			for line in lines:
				if line.split('=')[0].strip() in list(SettingsDict.keys()):
					exec('self.'+line)
			file.close()

	def initBCsDict(self,VariableList):
		''' Initialize BCsDict '''
		from FELiCS.functions import printWarning
		BCIDList=self.__BCIDs__
		#First define local BCsDict
		BCsDict={}
		for Variable in VariableList:
			BCsDict[Variable]=[]
			for BCID in BCIDList:
				BCsDict[Variable].append({'ID':BCID,'type':'Neumann','value':0.0})
		self.__BCsDict__=BCsDict

	def importBCsDict(self,VariableList):
		''' Import a boundary condition file with checking the consistency of BCs and mesh. To read the BCs without checking use importSettings()'''
		from FELiCS.functions import printWarning
		BCIDList=self.__BCIDs__
		#First define local BCsDict
		BCsDict={}
		filepath=self.BCsFilePath
		if not filepath == '' and isfile(filepath):
			self.BCsFilePath = filepath
			for Variable in VariableList:
				BCsDict[Variable]=[]
				for BCID in BCIDList:
					BCsDict[Variable].append({'ID':BCID,'type':'Neumann','value':0.0})
			# Read BCFile
			BCFile=open(self.BCsFilePath)
			importDict=eval(BCFile.readline())
			# Loop over all variables and IDs and if needed values present in BCFile, copy the contents to the local BCsDict
			for Variable in VariableList:
				if Variable in list(importDict.keys()):
					BCsDict[Variable]=[]
					for BC in importDict[Variable]:
						if BC['ID'] in BCIDList:
							BCsDict[Variable].append(BC)
						else:
							printWarning('Boundary condition of variable '+Variable+' for boundary with ID '+str(BC['ID'])+' not found in file. Choosing homogeneous Neumann instead.')

				else:
					printWarning('Boundary conditions for variable '+Variable+' not found in file. Choosing homogeneous Neumann instead.')
			# Finally, copy local BCsDict to the object
			self.__BCsDict__=BCsDict


	def setBC(self,field,BoundaryID,BCType,BCvalue):
		''' Setting the Boundary condition of a single variable '''
		for BC in self.__BCsDict__[field]:
			if BC['ID']== BoundaryID:
				self.__BCsDict__[field][BoundaryID]['type']=BCType
				self.__BCsDict__[field][BoundaryID]['value']=BCvalue

	def setBCsFromGUI(self,BCTypeStr,BCValueStr):
		'''Function to communicate with the GUI.
		This function writes the BCs set in the GUI to a private object calles __BCsDict__
		\tInput:
		\t\t-BCsTypeStr: A dictionary containing the type of the BC for each variable and boundary
		\t\t-BCsValueStr: A dictionary containing the value of the BCs
		'''
		BCsDict={}
		for Variable in list(BCTypeStr.keys()):
			BCsDict[Variable]=[]
			for ID in list(BCTypeStr[Variable].keys()):
				BCsDict[Variable].append({'ID':ID,'type':BCTypeStr[Variable][ID].get(),'value':float(BCValueStr[Variable][ID].get())})
		self.__BCsDict__=BCsDict

	def getBCsDict(self):
		''' The BCsDict is a private variable of the class.
		This function returns the BCsDict '''
		return self.__BCsDict__
	def setBCsDict(self,BCsDict):
		''' The BCsDict is a private variable of the class.
		This function returns the BCsDict '''
		self.__BCsDict__ = BCsDict

	def readBCInfo(self,MeshFilePath, felicsMesh):
		''' Input: - MeshFilePath
		This function reads both the IDs of the boundary conditions from the mesh and stores them
		in a private list of the class and also the boundary nodes and stores them in __boundaries__'''
		from numpy import unique

		# get a list of all kinds of BC indices
		self.__BCIDs__ = unique(felicsMesh.facet_tags.values)
		self.__boundaries__ = felicsMesh.facet_tags
	def getBCIDs(self):
		''' Returning a list containing all indices of the boundary conditions'''
		return self.__BCIDs__
	def getBoundaries(self):
		''' Returning the boundary nodes'''
		return self.__boundaries__

	def complete(self,BCVariableList):
		''' Checking if all necessary case attributes are present. A Variable List, as well as the index list of boundaries needs to be provided, to check if all necessary boundaries are set. '''
		from FELiCS.functions import printOK
		BCIDList= self.__BCIDs__
		if isfile(self.BCsFilePath):
			EverythingPresent=True
		else:
			printOK('Set boundary conditions file!')
			EverythingPresent=False
		return EverythingPresent

	def readDomainData(
                        self,
                        Meshfile, 
                        gDim, 
                        ExtendedTransportedQuantityList,
                        coordinateSystem,
                        m,
                        ):
		''' Input: Mesfile
		Read all the domain data from the meshfile '''
		self.readMesh(
                            Meshfile, 
                            gDim,
                            coordinateSystem,
                            m
                            )
		self.readBCInfo(Meshfile, self.getMesh())
		self.initBCsDict(ExtendedTransportedQuantityList)
		self.importBCsDict(ExtendedTransportedQuantityList)

	def readMesh(
                    self,
                    MeshFile, 
                    dim, 
                    coordinateSystem,
                    m,
                    ):
		'''
			Reading Meshfile and saving it as private object

			Function Arguments:
			- MeshFile: File of a gmsh-meshfile. File needs to be in .msh format
			- gdim: Geometrical Dimension of the mesh. This argument is needed
			by the gmsh helper-functions, which read in the mesh

			Function returns:

		'''
		if not MeshFile == '' and path.isfile(MeshFile):
			self.__mesh__ = FELiCSMesh(
                                                coordinateSystem,
                                                MeshFile,
                                                dim,
                                                m,
                                                )
			self.dim = self.__mesh__.gdim

	def getMesh(self):
		''' Function is returning the mesh '''
		return self.__mesh__
        
