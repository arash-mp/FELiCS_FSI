import  os
import  h5py
import  numpy               as np
from    scipy               import interpolate
from    FELiCS.Misc.logging import Logger

# Debug
import  matplotlib
import  matplotlib.pyplot   as plt
matplotlib.use("TkAgg")


# Get the logger
logger = Logger.get_logger("felics")

class Reader:
    
    def __init__(
            self, 
            param, 
            FEMSpaces, 
            filePath, 
            GroupName, 
            needInterpolation   = True,
            originalMeshFile    = None,
            isComplex           = False,
            cacheData           = True,
            ):
        """
        Initialize the Reader for loading and interpolating data onto FELiCS mesh.
        
        This class handles reading data from HDF5 files and mapping/interpolating
        them onto the FELiCS mesh. Each Reader instance corresponds to one file.
        
        Parameters
        ----------
        param : object
            Object containing calculation parameters.
        FEMSpaces : object
            Object containing calculation FEM spaces, export FEM spaces and meshes.
        filePath : str
            Path to the HDF5 file to be read.
        GroupName : str
            Name of the group in the file (used for coordinates in main file).
        needInterpolation : bool, optional
            Whether to interpolate data onto the FELiCS mesh, by default True.
        isComplex : bool, optional
            Whether the data is complex-valued, by default False.
        originalMeshFile : str or None, optional
            Path to original mesh file if coordinates not in main file, by default None.
        cacheData : bool, optional
            Whether to load all available data in file to a cache, by default True.
            
        TODOs:
            - Implement interpolation method.
            - Handle mixed function spaces when setting arrays.
            - Bypass the whole group name issue, by opening until we get variable names?

        Raises
        ------
        FileNotFoundError
            If the specified file does not exist.
        """
        self._FEMSpaces         = FEMSpaces
        self._param             = param
        self._ndim              = param.Case.nDim
        self._mapping           = FEMSpaces.mappingObj
        self._mesh              = param.getMesh()
        self._coordnames        = self._mesh.coordinateNames
        self._groupName         = GroupName
        self._needInterpolation = needInterpolation
        self._filePath          = filePath
        self._cacheData         = cacheData
        self._isComplex         = isComplex
        self._rawDataDict       = {}        # Where we store the numpy arrays

        # Validate file existence
        if not os.path.isfile(filePath):
            logger.error(f"File '{filePath}' not found.")
            raise FileNotFoundError(f"File '{filePath}' not found.")
        
        # Cache for loaded data to avoid reloading
        self._availableVars     = self._getListOfAvailableVariables()
        if self._cacheData:
            self._loadArraysFromFile(self._availableVars)
        
        # Get calculation mesh coordinates
        self._calcMeshCoords    = self._FEMSpaces.P2.tabulate_dof_coordinates()[:, 0:self._ndim]
        # NOTE: For P1 we rarely need it, so commented for now
        # self._calcP1MeshCoords  = FEMSpaces.P1.tabulate_dof_coordinates()[:, 0:self._ndim]
        
        # Load import mesh coordinates
        self._load_importMesh_coordinates(needInterpolation, originalMeshFile, filePath)
        
        # Create mapping from import to calculation mesh if we don't interpolate
        if not self._needInterpolation:
            self._import_to_P2calc  = self._mapping._mappingFunc(
                self._importMeshCoords, 
                self._calcMeshCoords
            )
        # NOTE: For P1 we rarely need it, so commented for now
        # self._import_to_P1calc = self._mapping._mappingFunc(
        #     self._importMeshCoords, 
        #     self._calcP1MeshCoords
        # )
        # Mapping check: for debugging purposes
        # assert np.allclose(
        #     self._importMeshCoords[self._import_to_P2calc],
        #     Field.space.sub(icomp).collapse()[0].tabulate_dof_coordinates()[:, 0:self._ndim],
        #     atol=1e-12, rtol=0
        # ), "Coordinate mismatch: import and P2 coords differ beyond tolerance."

    def _getListOfAvailableVariables(self):
        """
        Get list of available variables in the file group.
        
        Returns
        -------
        list
            List of available variable names in the specified group.
        """
        with h5py.File(self._filePath, 'r') as f:
            if self._groupName in f:
                self._availableVars = list(f[self._groupName].keys())
                logger.debug(f"Available variables in file under group '{self._groupName}': {self._availableVars}")
            else:
                logger.warning(f"Group '{self._groupName}' not found in file.")
                self._availableVars = []
        return self._availableVars

    def _loadArraysFromFile(self, listOfVars):
        """
        Load variables in the list into a dictionnary.
        """
        with h5py.File(self._filePath, 'r') as f:
            for varName in listOfVars:
                if varName not in self._rawDataDict:
                    datasetName = f'{self._groupName}/{varName}'
                    if self._isComplex:
                        # For complex variables, load magnitude and angle to make real and imag parts
                        magnitude   = self._load_or_error(f, datasetName + '/magnitude')
                        angle       = self._load_or_error(f, datasetName + '/angle')
                        self._rawDataDict[varName + '_real'] = magnitude * np.cos(angle)
                        self._rawDataDict[varName + '_imag'] = magnitude * np.sin(angle)
                    else:
                        if not self._needInterpolation:
                            # If we don't interpolate, we load a FELiCS file with magnitude only
                            self._rawDataDict[varName] = self._load_or_error(f, datasetName + '/magnitude')
                        else:
                            # Otherwise we load the full variable
                            self._rawDataDict[varName] = self._load_or_error(f, datasetName)

    def _load_importMesh_coordinates(self, needInterpolation, originalMeshFile, filePath):
        """
        Load mesh coordinates based on interpolation requirements and file sources.
        Three scenarios:
            1. Interpolation needed, original mesh file provided: load from original mesh file.
            2. Interpolation needed, no original mesh file: load from main file under group name.
            3. No interpolation needed: load from FELiCS exported mesh file.

        Parameters
        ----------
        needInterpolation : bool
            Whether interpolation is needed.
        originalMeshFile : str or None
            Path to the original mesh file, if available.
        filePath : str
            Path to the main data file.
        """
        if needInterpolation and originalMeshFile is not None:
            # Load coordinates from original mesh file
            mesh_file                               = originalMeshFile
            coord_prefix                            = 'coordinates/' # NOTE: we might want to make this more general
            
            with h5py.File(mesh_file, 'r') as handle:
                self._sizeImportCoord               = handle[coord_prefix + self._coordnames[0]][:].shape
                self._importMeshCoords              = np.zeros((max(self._sizeImportCoord), len(self._coordnames)))
                for i, coord in enumerate(self._coordnames):
                    self._importMeshCoords[:, i]    = self._load_or_error(handle, coord_prefix + coord)
                
        elif needInterpolation and originalMeshFile is None:
            # Assume coordinates are in the main file with same group name
            # First, load them if not already done
            if not all(coord in self._availableVars for coord in self._coordnames):
                listOfVars                          = [coord for coord in self._coordnames if coord not in self._availableVars]
                self._loadArraysFromFile(self, listOfVars)
                
            # Then set as mesh coords
            self._sizeImportCoord                   = self._rawDataDict[self._coordnames[0]].shape
            self._importMeshCoords                  = np.zeros((max(self._sizeImportCoord), len(self._coordnames)))
            for i, coord in enumerate(self._coordnames):
                self._importMeshCoords[:, i]        = self._rawDataDict[coord]
                
        else:
            # Even if we don't interpolate, we still need to load the FELiCS mesh coordinates
            mesh_file                               = os.path.join(self._param.Export.ExportFolder, f"{self._param.Case.AnalysisMode}_mesh.h5")
            coord_prefix                            = 'coordinates/'
            
            with h5py.File(mesh_file, 'r') as handle:
                self._sizeImportCoord               = handle[coord_prefix + self._coordnames[0]][:].shape
                self._importMeshCoords              = np.zeros((max(self._sizeImportCoord), len(self._coordnames)))
                
                # TODO: in FELiCS meshes it's always x,y,z -> need to adapt to coord sys!
                for i, coord in enumerate(['x','y','z'][:self._ndim]):
                    self._importMeshCoords[:, i]    = self._load_or_error(handle, coord_prefix + coord)
        

    def _load_or_error(self, fileHandle, datasetName):
        """
        Helper function to load a dataset from an HDF5 file or raise an error if not found.
        """
        if datasetName in fileHandle:
            return fileHandle[datasetName][:]
        else:
            logger.error(f"Dataset '{datasetName}' not found in the file.")
            raise KeyError(f"Dataset '{datasetName}' not found in the file.")
    
    def _getNameVarsToLoadForField(self, Field):
        """
        This function returns the names of the arrays to be loaded from file
        depending on the configuration.
        It takes into account that vector fields need to load each component
        separately.
        NOTE: Assumes the names in the Field are the same as in the file.
        """
        nameList                = []
        
        # Check what type of field we are loading the data into
        fieldInfo               = Field.describeFunctionSpace()
        
        # If we deal with a mixed space, get the subfields
        if fieldInfo['type'] == "mixed":
            subFields           = Field.getListOfSubFields()
        else:
            subFields           = [Field]

        # Loop over the subfield to check if vectors or scalars are needed
        for subField in subFields:
            subFieldInfo        = subField.describeFunctionSpace()
            
            if subFieldInfo['type'] == "vector":
                # Get component names from the Field
                componentNames  = subField.getComponentsNames()
                
                # Add the components to the nameList
                for compName in componentNames:
                    # Specific case for u_forcing. TODO: make the same pattern for all vector variables!
                    if Field.getName()[0] == 'u':
                        nameList.append(subField.getName()[0] + compName + subField.getName()[1:])
                    else:
                        nameList.append(subField.getName() + compName)

            else:
                # We have a scalar field, so we just add the name
                nameList.append(subField.getName())
                
        return nameList
    
    def _interpolateDataToFELiCSMesh(self, Field, nameList):
        """
        This function interpolates the data from the import mesh to the FELiCS calculation mesh.
        """
        logger.info(f'Interpolating {nameList} on FELiCS mesh.')
        importMesh                  = self._importMeshCoords[:, 0:self._ndim]
        destinationMesh             = self._calcMeshCoords
        interpolatedData            = {}
        
        # Set all fields to interpolate in one array for speed
        dataToInterpolate           = np.zeros((max(self._sizeImportCoord), len(nameList)))
        for i, varName in enumerate(nameList):
            dataToInterpolate[:, i] = self._rawDataDict[varName]
        
        # Debug
        # plt.figure()
        # plt.tricontourf(importMesh[:,0], importMesh[:,1], dataToInterpolate[:, 0], levels=14)
        # plt.colorbar()
        # plt.title(f'{nameList[0]} imported')

        # Use nearest-neighbor interpolation for backup
        backupNearestData           = interpolate.griddata(
            importMesh, 
            dataToInterpolate, 
            destinationMesh, 
            method='nearest'
        )
        
        # Debug
        # plt.figure()
        # plt.tricontourf(destinationMesh[:,0], destinationMesh[:,1], backupNearestData[:, 0], levels=14)
        # plt.colorbar()
        # plt.title(f'{nameList[0]} nearest-neighbor')
        
        # Try linear interpolation and fill in NaNs with nearest-neighbor
        interpolateDataArray        = interpolate.griddata(
            importMesh, 
            dataToInterpolate, 
            destinationMesh, 
            method='linear'
        )
        
        # Use mask to fill in NaNs from linear with nearest-neighbor values
        # Assumes they are NaNs in all variables at the same locations
        mask                                = np.argwhere(np.isnan(interpolateDataArray[:,0]))
        if np.any(mask):
            logger.warning(f"Linear interpolation produced NaNs for {mask.size} points. Filling with nearest-neighbor values.")
            # Fill in the interpolated data
            interpolateDataArray[mask, :]   = backupNearestData[mask, :]
            
            # Debug
            # plt.figure()
            # plt.scatter(destinationMesh[:,0], destinationMesh[:,1], marker='.', c='k', label='All points')
            # plt.scatter(destinationMesh[mask,0], destinationMesh[mask,1], marker='o', c='r', label='NaN points')
            # plt.title('Points with NaNs after linear interpolation')
        
        # Debug
        # plt.figure()
        # plt.tricontourf(destinationMesh[:,0], destinationMesh[:,1], interpolateDataArray[:, 0], levels=14)
        # plt.colorbar()
        # plt.title(f'{nameList[0]} linear interpolated with NaNs filled')
        # plt.show()

        # Put into dictionary
        for i, varName in enumerate(nameList):
            interpolatedData[varName]       = interpolateDataArray[:, i]
        
        return interpolatedData
    
    def _mapDataToFELiCSMesh(self, Field, nameList):
        """
        This function maps the data from the import mesh to the FELiCS calculation mesh.
        """
        mappedData              = {}
        
        # Loop over variables and map the data using the precomputed mapping
        for varName in nameList:
            if varName in self._rawDataDict:
                mappedData[varName] = self._rawDataDict[varName][self._import_to_P2calc]
            else:
                logger.warning(f"Variable '{varName}' not found in cached data.")
                
            # Extra debug:
            # import matplotlib
            # import matplotlib.pyplot as plt
            # matplotlib.use("TkAgg")
            # plt.figure()
            # plt.tricontourf(self._calcMeshCoords[:,0], self._calcMeshCoords[:,1], mappedData[varName], levels=14)
            # plt.colorbar()
            # plt.title(f'{varName} mapped')
            # plt.show()

        return mappedData
    
    def _setArraysToField(self, dictOfArrays, Field):
        """
        This function sets the values of a given field from a dictionary of arrays.
        The keys of the dictionary are the names of the variables.
        """
        
        # Get the Field information
        fieldInfo               = Field.describeFunctionSpace()
        
        # NOTE: We might have to be careful to set the arrays into sub.sub fields if we have mixed type.
        if fieldInfo['type'] == "mixed":
            logger.error("Setting arrays to mixed fields is not implemented yet.")
            raise NotImplementedError("Setting arrays to mixed fields is not implemented yet.")
        
        if fieldInfo['type'] == "vector":
            # Get component names from the Field
            componentNames      = Field.getComponentsNames()
            
            # Loop over components sub-fields
            for icomp, compName in enumerate(componentNames):
                # Specific case for u_forcing. TODO: make the same pattern for all vector variables!
                if Field.getName()[0] == 'u':
                    fullCompName = Field.getName()[0] + compName + Field.getName()[1:]
                else:
                    fullCompName = Field.getName() + compName
                    
                if fullCompName in dictOfArrays:
                    # Get indices of the subfield in the mixed space
                    indicesOfSubField = Field.space.sub(icomp).collapse()[1]
                    Field.function.x.array[indicesOfSubField] = dictOfArrays[fullCompName]
                    
        else:
            # We have a scalar field, so we just add the name
            varName = Field.getName()
            if varName in dictOfArrays:
                Field.function.x.array[:] = dictOfArrays[varName]
        
        return Field
        
    def importInField(
            self,
            Field,
        ):
        """
        This function does the loading of the data from the file into the given Field.
        NOTE: Assumes the names in the Field are the same as in the file!
        """
        
        # Get the list of variables we want in the field
        nameListField       = self._getNameVarsToLoadForField(Field)
        logger.debug(f"Loading {nameListField}...")
        
        # Check which variables are not available in the file and filter the list to only those available
        missingVars         = [name for name in nameListField if name not in self._availableVars]
        nameListField       = [name for name in nameListField if name not in missingVars]
        if missingVars:
            logger.warning(f"Variables {missingVars} not found in file '{self._filePath}'. Set to default values.")
            
        # If we didn't cache the data, we need to load the variables now
        if not self._cacheData and nameListField:
            self._loadArraysFromFile(nameListField)
            
        # Debug that the data matches the import coordinates
        # import matplotlib
        # import matplotlib.pyplot as plt
        # matplotlib.use("TkAgg")
        # plt.figure()
        # plt.tricontourf(self._importMeshCoords[:,0], self._importMeshCoords[:,1], self._rawDataDict['ux'], levels=14)
        # plt.colorbar()
        # plt.title(f'ux imported')
        # plt.show()
            
        # Interpolate or map the data to the FELiCS mesh, if we have some variables to load
        if self._needInterpolation and nameListField:
            dataForField    = self._interpolateDataToFELiCSMesh(Field, nameListField)
        elif nameListField:
            dataForField    = self._mapDataToFELiCSMesh(Field, nameListField)

        # Set the data to the corresponding field entries
        if nameListField:
            Field           = self._setArraysToField(dataForField, Field)

        return Field, missingVars