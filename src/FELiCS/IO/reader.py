import  os
import  h5py
import  numpy               as np
from    scipy               import interpolate
from    FELiCS.Misc.logging import Logger


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
        self._dataCache         = {}        # Where we store the numpy arrays

        # Validate file existence
        if not os.path.isfile(filePath):
            logger.error(f"File '{filePath}' not found.")
            raise FileNotFoundError(f"File '{filePath}' not found.")
        
        # Cache for loaded data to avoid reloading
        self._availableVars     = self._getListOfAvailableVariables()
        if self._cacheData:
            self._cacheVariables(self._availableVars)
        
        # Get calculation mesh coordinates
        self._calcP2MeshCoords  = FEMSpaces.P2.tabulate_dof_coordinates()[:, 0:self._ndim]
        self._calcP1MeshCoords  = FEMSpaces.P1.tabulate_dof_coordinates()[:, 0:self._ndim]
        
        # Load import mesh coordinates
        self._load_importMesh_coordinates(needInterpolation, originalMeshFile, filePath)
        
        # Create mapping from import to calculation mesh
        self._import_to_P2calc  = self._mapping._mappingFunc(
            self._importMeshCoords, 
            self._calcP2MeshCoords
        )
        # NOTE: For P1 we rarely need it, so commented for now
        # self._import_to_P1calc = self._mapping._mappingFunc(
        #     self._importMeshCoords, 
        #     self._calcP1MeshCoords
        # )

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
            else:
                logger.warning(f"Group '{self._groupName}' not found in file.")
                self._availableVars = []
        return self._availableVars

    def _cacheVariables(self, listOfVars):
        """
        Load variables in the list into a dictionnary.
        """
        with h5py.File(self._filePath, 'r') as f:
            for varName in listOfVars:
                if varName not in self._dataCache:
                    datasetName = f'{self._groupName}/{varName}'
                    if self._isComplex:
                        # For complex variables, load magnitude and angle to make real and imag parts
                        magnitude   = self._load_or_error(f, datasetName + '/magnitude')
                        angle       = self._load_or_error(f, datasetName + '/angle')
                        self._dataCache[varName + '_real'] = magnitude * np.cos(angle)
                        self._dataCache[varName + '_imag'] = magnitude * np.sin(angle)
                    else:
                        if not self._needInterpolation:
                            # If we don't interpolate, we load a FELiCS file with magnitude only
                            self._dataCache[varName] = self._load_or_error(f, datasetName + '/magnitude')
                        else:
                            # Otherwise we load the full variable
                            self._dataCache[varName] = self._load_or_error(f, datasetName)

    def _load_importMesh_coordinates(self, needInterpolation, originalMeshFile, filePath):
        """
        Load mesh coordinates based on interpolation requirements and file sources.

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
            mesh_file       = originalMeshFile
            coord_prefix    = 'coordinates/' # NOTE: we might want to make this more general
        elif needInterpolation and originalMeshFile is None:
            # Load coordinates from main file with group name
            mesh_file       = filePath
            coord_prefix    = f'{self._groupName}'
        else:
            # Load coordinates from FELiCS exported mesh
            mesh_file       = os.path.join(self._param.Export.ExportFolder, f"{self._param.Case.AnalysisMode}_mesh.h5")
            coord_prefix    = 'coordinates/'

        with h5py.File(mesh_file, 'r') as handle:
            self._sizeImportCoord   = handle[coord_prefix + self._coordnames[0]][:].shape
            self._importMeshCoords  = np.zeros((max(self._sizeImportCoord), len(self._coordnames)))
            for i, coord in enumerate(self._coordnames):
                self._importMeshCoords[:, i] = self._load_or_error(handle, coord_prefix + coord)

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
                    nameList.append(subField.getName() + compName)
            
            else:
                # We have a scalar field, so we just add the name
                nameList.append(subField.getName())
                
        return nameList
    
    def _interpolateDataToFELiCSMesh(self, Field, nameList):
        """
        This function interpolates the data from the import mesh to the FELiCS calculation mesh.
        """
        interpolatedData        = {}
        
        logger.error("Interpolation not yet implemented.")
        raise NotImplementedError("Interpolation not yet implemented.")
        
        return interpolatedData
    
    def _mapDataToFELiCSMesh(self, Field, nameList):
        """
        This function maps the data from the import mesh to the FELiCS calculation mesh.
        """
        mappedData              = {}
        
        # Loop over variables and map the data using the precomputed mapping
        for varName in nameList:
            if varName in self._dataCache:
                mappedData[varName] = self._dataCache[varName][self._import_to_P2calc]
            else:
                logger.warning(f"Variable '{varName}' not found in cached data.")

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
                
                # Mapping check: for debugging purposes
                assert np.allclose(
                    self._importMeshCoords[self._import_to_P2calc],
                    Field.space.sub(icomp).collapse()[0].tabulate_dof_coordinates()[:, 0:self._ndim],
                    atol=1e-12, rtol=0
                ), "Coordinate mismatch: import and P2 coords differ beyond tolerance."

                fullCompName    = Field.getName() + compName
                if fullCompName in dictOfArrays:
                    # Get indices of the subfield in the mixed space
                    indicesOfSubField = Field.space.sub(icomp).collapse()[1]
                    Field.function.x.array[indicesOfSubField] = dictOfArrays[fullCompName]
                    
                # Extra debug:
                # import matplotlib
                # import matplotlib.pyplot as plt
                # matplotlib.use("TkAgg")
                # comp_coords = Field.space.sub(icomp).collapse()[0].tabulate_dof_coordinates()[:, 0:self._ndim]
                # plt.figure()
                # plt.tricontourf(comp_coords[:,0], comp_coords[:,1], dictOfArrays[fullCompName], levels=14)
                # plt.colorbar()
                # plt.title(f'{fullCompName} mapped')
                # plt.show()
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
            self._cacheVariables(nameListField)
            
        # Debug that the data matches the import coordinates
        # import matplotlib
        # import matplotlib.pyplot as plt
        # matplotlib.use("TkAgg")
        # plt.figure()
        # plt.tricontourf(self._importMeshCoords[:,0], self._importMeshCoords[:,1], self._dataCache['ux'], levels=14)
        # plt.colorbar()
        # plt.title(f'ux imported')
        # plt.show()
            
        # Interpolate or map the data to the FELiCS mesh
        if self._needInterpolation:
            dataForField    = self._interpolateDataToFELiCSMesh(Field, nameListField)
        else:
            dataForField    = self._mapDataToFELiCSMesh(Field, nameListField)

        # Set the data to the corresponding field entries
        Field               = self._setArraysToField(dataForField, Field)

        return Field, missingVars