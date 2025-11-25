# Standard libraries
import  os
import  h5py
from    functools           import cached_property
from    typing              import Optional, List, Dict, Any, Tuple

# Third party libraries
import  numpy               as np
from    scipy.spatial       import Delaunay
from    scipy.interpolate   import LinearNDInterpolator, NearestNDInterpolator

# Local libraries and methods
from    FELiCS.Misc.logging import Logger
from    FELiCS.IO.Mapping   import Mapping

# Debug
# import matplotlib
# import matplotlib.pyplot as plt
# matplotlib.use("TkAgg")

# Get the logger
logger = Logger.get_logger("felics")

class Reader:
    """
    A class for reading and interpolating data from HDF5 files onto FELiCS meshes.

    This class manages the loading of data from HDF5 files, handles interpolation or direct mapping
    to the calculation mesh, and assigns the data to FELiCS fields. It supports caching for performance
    and can handle complex-valued data.

    Attributes
    ----------
    _sourceDir : str
        The source directory for data files.
    _needInterpolation : bool
        Whether to interpolate data onto the FELiCS mesh.
    _felicsMeshFilePath : str or None
        Path to the FELiCS mesh file.
    _isComplex : bool
        Whether the data is complex-valued.
    _cacheData : bool
        Whether to cache loaded data for reuse.
    _sourceKey : tuple or None
        Current source key (filePath, groupName).
    _availableVars : list or None
        List of available variables in the current source.
    _rawDataDict : dict
        Raw data dictionary from HDF5 files.
    _dataForField : dict
        Processed data for fields.
    _calcMeshCoords : numpy.ndarray or None
        Coordinates of the calculation mesh.
    _calcMeshP1Coords : numpy.ndarray or None
        Coordinates of the P1 calculation mesh.
    _mappingImportToCalc : numpy.ndarray or None
        Mapping from import to calculation mesh.
    _fullAxisNames : list or None
        Full axis names.
    _meshAxisNames : list or None
        Mesh axis names.
    _numDoFsCalcMeshP2 : int
        Number of DoFs for P2 calculation mesh.
    _numDoFsCalcMeshP1 : int
        Number of DoFs for P1 calculation mesh.
    _spacesDegrees : list or None
        Degrees of spaces.

    Notes
    -----
    - Caching is per-file and cleared when switching sources.
    - Interpolation uses Delaunay triangulation for efficiency.
    - Retro-compatibility with previous FELiCS files with groupNames is ensured (for now).

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    RuntimeError
        If environment is not bound or source is not set.
    NotImplementedError
        For unsupported field types like mixed function spaces.
    """

    def __init__(
        self,
        sourceDir:          str,
        needInterpolation:  bool = True,
        felicsMeshFilePath: Optional[str] = None,
        isComplex:          bool = False,
        cacheData:          bool = True,
    ) -> None:
        """
        Initialize the Reader instance.

        Parameters
        ----------
        needInterpolation : bool, optional
            Whether to interpolate data onto the FELiCS mesh, by default True.
        felicsMeshFilePath : str or None, optional
            Path to the FELiCS mesh file WHEN it is in a different directory than the source data.
        isComplex : bool, optional
            Whether the data is complex-valued, by default False.
        cacheData : bool, optional
            Whether to cache loaded data for reuse, by default True.
        """
        # Configuration (fixed for this Reader instance)
        self._sourceDir             = sourceDir
        self._needInterpolation     = needInterpolation
        self._felicsMeshFilePath    = felicsMeshFilePath
        self._isComplex             = isComplex
        self._cacheData             = cacheData

        # Per-file caches (reset when source changes)
        self._sourceKey             = None  # (filePath, groupName) TODO: remove groupname from files
        self._availableVars         = None
        self._rawDataDict           = {}
        self._dataForField          = {}
        
        # Per-calculation mesh/FEM space caches
        self._calcMeshCoords        = None
        self._calcMeshP1Coords      = None
        self._mappingImportToCalc   = None
        self._fullAxisNames         = None
        self._meshAxisNames         = None
        self._numDoFsCalcMeshP2     = 0
        self._numDoFsCalcMeshP1     = 0
        self._spacesDegrees         = None

    # --------------------------
    # Cached properties
    # --------------------------
    @cached_property
    def meshAxisNames(self):
        """
        Names of the full spatial dimensions: mesh + spectral.

        Returns
        -------
        list of str
            The full dimension names (e.g., ['x', 'y', 'z']).

        Raises
        ------
        RuntimeError
            If accessed before bind_env is called.
        """
        return self._meshAxisNames
    
    @cached_property
    def fullAxisNames(self):
        """
        Get the names of all spatial dimensions: mesh + spectral.

        Returns
        -------
        list of str
            The full dimension names (e.g., ['x', 'y', 'z']).

        Raises
        ------
        RuntimeError
            If accessed before bind_env is called.
        """
        return self._fullAxisNames

    @cached_property
    def _triangulationImportMesh(self):
        """
        Cached Delaunay triangulation of the import mesh for interpolation.

        Returns
        -------
        scipy.spatial.Delaunay
            The triangulation object for the import mesh.
        """
        nDimMesh        = len(self.meshAxisNames)
        importCoords    = self.importMeshCoords[:, :nDimMesh]
        return Delaunay(importCoords)

    @cached_property
    def calcMeshP2Coords(self):
        """
        Get the coordinates of the calculation mesh.

        Returns
        -------
        numpy.ndarray
            Array of shape (nPoints, nDimMesh) with the coordinates of the calculation mesh.
        """
        return self._calcMeshP2Coords
    
    @cached_property
    def calcMeshP1Coords(self):
        """
        Get the coordinates of the P1 calculation mesh.

        Returns
        -------
        numpy.ndarray
            Array of shape (nPoints, nDimMesh) with the coordinates of the P1 calculation mesh.
        """
        return self._calcMeshP1Coords

    @cached_property
    def mappingImportToCalcP2(self):
        """
        Get the mapping from import mesh to calculation mesh (P2).

        Returns
        -------
        numpy.ndarray or None
            The mapping array if no interpolation is needed, otherwise None.

        Notes
        -----
        """
        if self._needInterpolation:
            return None
        else:
            return self._mappingImportToCalc

    @cached_property
    def importMeshCoords(self):
        """
        Get the coordinates of the import mesh based on the configuration.

        Three scenarios:
        1. Interpolation needed, original mesh file provided: load from original mesh file.
        2. Interpolation needed, no original mesh file: load from main file under group name.
        3. No interpolation needed: load from FELiCS exported mesh file.

        Returns
        -------
        numpy.ndarray
            Array of shape (nPoints, nDimMesh) with the coordinates of the import mesh.
        """

        logger.debug("Associating import mesh to Reader instance.")

        # Number and names of mesh axes
        meshAxisNames                       = self.meshAxisNames
        numMeshAxis                         = len(meshAxisNames)
        logger.debug(f"Getting the {numMeshAxis}D import mesh coordinates ({meshAxisNames}).")

        def _load_coordinates_from_file(file_path: str, prefix: str, coord_names: List[str]) -> np.ndarray:
            """
            Helper to load coordinates from HDF5 file.

            Parameters
            ----------
            file_path : str
                Path to the HDF5 file.
            prefix : str
                Prefix for the dataset paths.
            coord_names : list of str
                Names of the coordinate datasets.

            Returns
            -------
            numpy.ndarray
                Array of coordinates.
            """
            with h5py.File(file_path, "r") as fh:
                first_shape         = fh[prefix + coord_names[0]][:].shape
                coords              = np.zeros((max(first_shape), len(coord_names)))
                for i, name in enumerate(coord_names):
                    coords[:, i]    = fh[prefix + name][:]
            return coords
        
        if self._needInterpolation:
            
            if self._felicsMeshFilePath:
                # Option 1: load from a given FELiCS mesh file
                coordsArray = _load_coordinates_from_file(self._felicsMeshFilePath, "coordinates/", meshAxisNames)
                return coordsArray
            
            else:
                # Option 2: Coords into file containing the data. Load at first file access
                # NOTE: assumes the coordinates are in the same group as the variables
                self._ensure_source_set()
                self._get_list_available_vars()
                missingCoords               = [coord for coord in meshAxisNames if coord not in self._rawDataDict.keys()]
                if missingCoords:
                    self._load_from_h5(missingCoords)
                firstCoordShape             = self._rawDataDict[meshAxisNames[0]].shape
                coordsArray                 = np.zeros((max(firstCoordShape), len(meshAxisNames)))
                for i, coordName in enumerate(meshAxisNames):
                    coordsArray[:, i]       = self._rawDataDict[coordName]
                return coordsArray
        
        else:
            # Option 3: Load the FELiCS mesh
            if self._felicsMeshFilePath is not None:
                # For example when we load the mean flow, the FELiCS mesh is usually in another directory
                meshFile                    = self._felicsMeshFilePath
            else:
                # Default FELiCS mesh file in the source directory (e.g. when loading modes)
                meshFile                    = os.path.join(self._sourceDir, "mesh.h5")
            
            # Check that the file exists
            if not os.path.isfile(meshFile):
                logger.error(f"FELiCS mesh file for import not found: {meshFile}")
                raise FileNotFoundError(meshFile)

            coordsArray = _load_coordinates_from_file(meshFile, "coordinates/", ["x", "y", "z"][:numMeshAxis])
            return coordsArray

    # --------------------------
    # Source management
    # --------------------------
    def _update_data_source(
        self, 
        filePath: str, 
        groupName: str
    ) -> None:
        """
        Switch to a new file/group if different, clearing per-file caches.

        Parameters
        ----------
        filePath : str
            Path to the HDF5 file.
        groupName : str
            Name of the group in the HDF5 file.

        Raises
        ------
        FileNotFoundError
            If the file does not exist.
        """

        # Check that the file exists
        if not os.path.isfile(filePath):
            raise FileNotFoundError(filePath)

        # Check that the source directory fits that of the reader instance
        currentSourceDir         = os.path.dirname(filePath)
        if os.path.abspath(currentSourceDir) != os.path.abspath(self._sourceDir):
            # NOTE: maybe we do something specific in that case?
            logger.warning(
                f"Reader source directory '{currentSourceDir}' does not match "
                f"the Reader instance sourceDir '{self._sourceDir}': expected behavior issues."
            )

        # Set new source if different
        key = (os.path.abspath(filePath), groupName)
        # NOTE: might need completion of the variable list
        if key != self._sourceKey:
            self._sourceKey     = key
            self._availableVars = None
            self._rawDataDict   = {}
            self._dataForField  = {}
            logger.debug(f"Reader: switched source to {key}, cleared per-file caches.")
            
    def _update_calc_mesh(
        self, 
        field: Any
    ) -> None:
        """
        Switch to a new calculation mesh if different from cached one.
        Clearing per-calculation mesh caches.
        
        The calculation mesh is defined as the DoFs in the field.
        Specific treatment needed for mixed spaces with different DOF counts (e.g., P1 vs P2).
        
        This updates the following cached properties:
        - fullAxisNames
        - meshAxisNames
        - calcMeshP2Coords
        - mappingImportToCalcP2
        and the following internal attributes:
        - _numDoFsCalcMeshP2
        - _numDoFsCalcMeshP1

        Parameters
        ----------
        Field : object
            The FELiCS field object.

        Raises
        ------
        NotImplementedError
            If the field type is not supported.
        """
        
        # Get the field type
        fieldType               = field.info["type"]

        # Polynomial degrees of all subfields (not for vector subspaces)
        subFields               = field.getListOfSubFields() if fieldType == "mixed" else [field]
        degrees                 = [subField.info.get("degree", None) for subField in subFields]
        if max(degrees) > 2:
            logger.error("Reader currently only supports P1 and P2 FEM spaces.")
            raise NotImplementedError("Reader currently only supports P1 and P2 FEM spaces.")
        
        # Get the number of DoFs in P2 and P1 subfields if present
        p2SubFieldIndex            = degrees.index(2) if 2 in degrees else None    # Keep only the first occurrence
        p1SubFieldIndex            = degrees.index(1) if 1 in degrees else None    # Keep only the first occurrence
        numDofsForP2               = subFields[p2SubFieldIndex].info['nDofsSpace'] if p2SubFieldIndex is not None else 0
        numDofsForP1               = subFields[p1SubFieldIndex].info['nDofsSpace'] if p1SubFieldIndex is not None else 0

        # If we are dealing with a new space/mesh, clear cached properties
        if numDofsForP2 != self._numDoFsCalcMeshP2 or numDofsForP1 != self._numDoFsCalcMeshP1:
            logger.debug("New calculation mesh / FEM space type, updating cache.")

            # Clear cached properties that depend on calculation meshes
            propertiesToClear       = [
                "mappingImportToCalcP2",
                "fullAxisNames",
                "meshAxisNames",
            ]
            if numDofsForP1 != self._numDoFsCalcMeshP1:
                propertiesToClear.append("calcMeshP1Coords")
            if numDofsForP2 != self._numDoFsCalcMeshP2:
                propertiesToClear.append("calcMeshP2Coords")
            self.clear_cached_properties(cachedProps=propertiesToClear)

            # Update the number of DoFs in calculation mesh and space degrees
            self._numDoFsCalcMeshP2 = numDofsForP2
            self._numDoFsCalcMeshP1 = numDofsForP1
            self._spacesDegrees     = degrees

            # Set the names in cached properties
            felicsMesh              = field.mesh
            self._meshAxisNames     = felicsMesh.meshAxisNames
            self._fullAxisNames     = felicsMesh.axisNames
            
            # If we have sub-spaces, get the indices of the sub-space DoFs
            if fieldType == "mixed":
                # TODO: check that this does not mess up with the DoF indices
                if numDofsForP2 > 0:
                    if field.getListOfSubFields()[p2SubFieldIndex].info['type'] == 'vector':
                        self._calcMeshP2Coords  = field.space.sub(p2SubFieldIndex).sub(0).collapse()[0].tabulate_dof_coordinates()[:, :len(self._meshAxisNames)]
                    else:
                        self._calcMeshP2Coords  = field.space.sub(p2SubFieldIndex).collapse()[0].tabulate_dof_coordinates()[:, :len(self._meshAxisNames)]
                if numDofsForP1 > 0:
                    if field.getListOfSubFields()[p1SubFieldIndex].info['type'] == 'vector':
                        self._calcMeshP1Coords  = field.space.sub(p1SubFieldIndex).sub(0).collapse()[0].tabulate_dof_coordinates()[:, :len(self._meshAxisNames)]
                    else:
                        self._calcMeshP1Coords  = field.space.sub(p1SubFieldIndex).collapse()[0].tabulate_dof_coordinates()[:, :len(self._meshAxisNames)]
            else:
                # Single space
                dofsArray                       = field.space.tabulate_dof_coordinates()[:, :len(self._meshAxisNames)]
                if numDofsForP2 > 0:
                    self._calcMeshP2Coords      = dofsArray
                if numDofsForP1 > 0:
                    self._calcMeshP1Coords      = dofsArray
                    
            # Compute the mapping (only to P2 mesh) if needed
            if not self._needInterpolation:
                numMeshDimensions               = len(self._meshAxisNames)
                self._mappingImportToCalc       = Mapping.calculateMappingFromDofs(
                    self.importMeshCoords[:, :numMeshDimensions],
                    self.calcMeshP2Coords
                )

    def _ensure_source_set(self) -> None:
        """
        Ensure that a data source is set.

        Raises
        ------
        RuntimeError
            If no data source is set.
        """
        if self._sourceKey is None:
            raise RuntimeError("No data source set. Call importInField(..., filePath, groupName) first.")

    def clear_cached_properties(
        self, 
        cachedProps: List[str] = []
    ) -> None:
        """
        Remove all cached @cached_property values from this Reader instance.

        They will be recomputed on the next access.

        Parameters
        ----------
        cachedProps : list of str, optional
            List of property names to clear. If empty, clears default properties.
        """
        
        if not cachedProps:
            # TODO: complete the list if there are some missing
            cachedProps = [
                "fullAxisNames",
                "meshAxisNames",
                "calcMeshP2Coords",
                "importMeshCoords",
                "mappingImportToCalcP2",
            ]
        cleared = []
        for propName in cachedProps:
            if propName in self.__dict__:
                self.__dict__.pop(propName)
                cleared.append(propName)
        if cleared:
            logger.debug(f"Cleared cached properties: {cleared}")
        else:
            logger.debug("No cached properties to clear.")

    # --------------------------
    # HDF5 file handling
    # --------------------------
    def _get_list_available_vars(self) -> None:
        """
        List variable names available in the current group.

        Raises
        ------
        RuntimeError
            If source is not set.
        """
        if self._availableVars is not None:
            return
        filePath, groupName         = self._sourceKey
        
        # Open file and get variable names
        with h5py.File(filePath, "r") as fileHandle:
            if groupName is None:
                self._availableVars = list(fileHandle.keys())
            elif groupName in fileHandle:
                self._availableVars = list(fileHandle[groupName].keys())
            else:
                logger.warning(f"Group '{groupName}' not found in file '{filePath}'.")
                self._availableVars = []
                
        # Remove variables that should not be there (frequency, gains, etc.)
        if "omega" in self._availableVars:
            self._availableVars.remove("omega")
        if "gain" in self._availableVars:
            self._availableVars.remove("gain")
        if "number" in self._availableVars:
            self._availableVars.remove("number")
            
    def _load_from_h5(
        self, 
        varNames: List[str]
    ) -> None:
        """
        Load selected variables from HDF5 file into _rawDataDict.

        Parameters
        ----------
        varNames : list of str
            List of variable names to load.

        Notes
        -----
        NOTE: Maybe we just check if the groups exist instead of relying on _isComplex?
        """
        # NOTE: maybe we just check if the groups exist instead of relying on _isComplex?
        filePath, groupName         = self._sourceKey
        with h5py.File(filePath, "r") as fileHandle:
            for varName in varNames:
                if varName in self._rawDataDict:
                    continue
                
                # Retro-compatibility with old FELiCS files
                if groupName is None:
                    basePath            = varName
                else:
                    basePath            = f"{groupName}/{varName}"
                
                # Check if variable value is in further "magnitude" group (retro-compatibility with old FELiCS files)
                if (basePath + "/magnitude") in fileHandle:
                    self._rawDataDict[varName] = fileHandle[basePath + "/magnitude"][:].squeeze()
                else:
                    self._rawDataDict[varName] = fileHandle[basePath][:].squeeze()

    # --------------------------
    # Interpolation / mapping
    # --------------------------
    def _interpolate_to_calc_mesh(
        self, 
        varNames:   List[str], 
        toP1:       bool = False
    ) -> Dict[str, Any]:
        """
        Interpolate given variables from import mesh to calculation mesh.

        Parameters
        ----------
        varNames : list of str
            List of variable names to interpolate.
        toP1 : bool, optional
            Whether to interpolate to P1 calculation mesh, by default False.
            NOTE: not very elegant

        Returns
        -------
        dict
            Dictionary mapping variable names to interpolated arrays.

        Raises
        ------
        ValueError
            If variable lengths do not match the import mesh.
        """
        numMeshDimensions                    = len(self.meshAxisNames)
        importCoords                = self.importMeshCoords[:, :numMeshDimensions]
        calcCoords                  = self.calcMeshP1Coords if toP1 else self.calcMeshP2Coords
        numImportPoints             = importCoords.shape[0]
        
        # If we have complex data, we have _real and _imag in _rawDataDict
        if self._isComplex:
            varNamesComplex        = []
            for varName in varNames:
                varNamesComplex.append(varName + "_real")
                varNamesComplex.append(varName + "_imag")
        else:
            varNamesComplex        = varNames

        # Sanity check: all variables have same length as import mesh
        for varName in varNamesComplex:
            if self._rawDataDict[varName].shape[0] != numImportPoints:
                raise ValueError(
                    f"Variable '{varName}' has length {len(self._rawDataDict[varName])}, "
                    f"but import mesh has {numImportPoints} points."
                )

        # Stack data columns (N, K)
        values                      = np.column_stack([self._rawDataDict[varName] for varName in varNamesComplex])

        # Use cached triangulation
        linearInterp                = LinearNDInterpolator(self._triangulationImportMesh, values)
        interpolated                = linearInterp(calcCoords)  # shape (M, K)

        # Handle NaNs by nearest-neighbour fill (vectorized)
        if np.isnan(interpolated).any():
            nearestInterp           = NearestNDInterpolator(importCoords, values)
            nearestValues           = nearestInterp(calcCoords)
            nanMask                 = np.isnan(interpolated)
            interpolated[nanMask]   = nearestValues[nanMask]
            
        # Output dict
        output                      = {varName: interpolated[:, i] for i, varName in enumerate(varNamesComplex)}
            
        # If we had complex data, reconstruct complex arrays
        if self._isComplex:
            for varName in varNames:
                output[varName] = output[varName + "_real"] + 1j * output[varName + "_imag"]
                del output[varName + "_real"]
                del output[varName + "_imag"]

        # Return a dict {var_name: array}
        return output

    def _map_to_calc_mesh(
        self, 
        varNames: List[str]
    ) -> Dict[str, Any]:
        """
        Applies precomputed mapping.

        Parameters
        ----------
        varNames : list of str
            List of variable names to map.

        Returns
        -------
        dict
            Dictionary mapping variable names to mapped arrays.
        """
        
        indices                         = self.mappingImportToCalcP2
        return {varName: self._rawDataDict[varName][indices] for varName in varNames}

    # --------------------------
    # Field helpers
    # --------------------------
    def _names_for_field(
        self, 
        field: Any
    ) -> List[str]:
        """
        Determine variable names to load from file for this field.

        Parameters
        ----------
        field : object
            The FELiCS field object.

        Returns
        -------
        list of str
            List of variable names.
        """
        names           = []
        info            = field.info
        
        subFields       = field.getListOfSubFields() if info["type"] == "mixed" else [field]
        
        for subField in subFields:
            subInfo     = subField.info
            
            if subInfo["type"] == "vector":
                for component in subField.getNamesOfSubFields():
                    
                    # Append "_real" and "_imag" if complex
                    if self._isComplex:
                        names.append(component + "_real")
                        names.append(component + "_imag")
                    else:
                        names.append(component)
                    
            elif subInfo["type"] == "scalar":
                if self._isComplex:
                    names.append(subField.name + "_real")
                    names.append(subField.name + "_imag")
                else:
                    names.append(subField.name)
                
            else:
                logger.error(f"Field type '{subInfo['type']}' not supported in Reader yet.")
                raise NotImplementedError("Field type not supported in Reader yet.")
            
        return names

    def _set_arrays_to_field(
        self, 
        arrays:  Dict[str, Any],
        field: Any
    ) -> Any:
        """
        Assign interpolated/mapped arrays to a given field.

        Parameters
        ----------
        arrays : dict
            Dictionary of arrays to assign.
        field : object
            The FELiCS field object.

        Returns
        -------
        object
            The updated field.

        Raises
        ------
        NotImplementedError
            For mixed function spaces.
        """
        # Field info
        info        = field.info
        
        def _assembleComplexArrays(varName: str, arrays: Dict[str, Any]):
            """Helper to assemble complex arrays from real and imaginary parts. 
            Or just return real array."""
            if self._isComplex:
                assembledArray = arrays[varName + "_real"] + 1j * arrays[varName + "_imag"]
            else:
                assembledArray = arrays[varName]
            return assembledArray
        
        # If we load complex data, get list of names without _real/_imag
        if self._isComplex:
            baseNames = [name[:-5] for name in arrays.keys() if name.endswith("_real")]
        else:
            baseNames = list(arrays.keys())
            
        # Set the arrays in FEM depending on field type
        if info["type"] == "mixed":
            # Loop over the subfields
            for iField, subFieldName in enumerate(field.getNamesOfSubFields()):
                # If subfield is a vector, loop over its components
                if info['subspaces'][iField]['type'] == 'vector':
                    subFields                               = field.getListOfSubFields()
                    for jComp, compName in enumerate(subFields[iField].getNamesOfSubFields()):
                        if compName in baseNames:
                            indices                         = field.space.sub(iField).sub(jComp).collapse()[1]
                            field.function.x.array[indices] = _assembleComplexArrays(compName, arrays)
                        else:
                            logger.warning(f"Component '{compName}' not found in loaded arrays for vector subfield '{subFieldName}'. Set to default values.")
                # For a scalar subfield, just set the array
                elif info['subspaces'][iField]['type'] == 'scalar':
                    if subFieldName in baseNames:
                        indices                             = field.space.sub(iField).collapse()[1]
                        field.function.x.array[indices]     = _assembleComplexArrays(subFieldName, arrays)
                    else:
                        logger.warning(f"Subfield '{subFieldName}' not found in loaded arrays for scalar subfield. Set to default values.")
                else:
                    logger.error(f"Subfield type '{info['subspaces'][iField]['type']}' not supported in Reader yet.")
                    raise NotImplementedError("Subfield type not supported in Reader yet.")
                
        elif info["type"] == "vector":
            for i, subFieldName in enumerate(field.getNamesOfSubFields()):
                if subFieldName in baseNames:
                    indices                         = field.space.sub(i).collapse()[1]
                    field.function.x.array[indices] = _assembleComplexArrays(subFieldName, arrays)
                else:
                    logger.warning(f"Component '{subFieldName}' not found in loaded arrays for vector field '{field.name}'. Set to default values.")
        else:
            # Then it's a scalar
            varName                                 = field.name
            if varName in baseNames:
                field.function.x.array[:]           = _assembleComplexArrays(varName, arrays)
            else:
                logger.warning(f"Variable '{varName}' not found in loaded arrays for scalar field. Set to default values.")
        return field

    # --------------------------
    # Main API
    # --------------------------
    def importInField(
        self, 
        field: Any, 
        filePath:   str, 
        groupName:  str
    ) -> Tuple[Any, List[str]]:
        """
        Import data from the specified file/group into the given field.

        Automatically clears per-file caches when file/group changes.

        Parameters
        ----------
        field : object
            The FELiCS field to import data into.
        filePath : str
            Path to the HDF5 file.
        groupName : str
            Name of the group in the HDF5 file.

        Returns
        -------
        tuple
            (updated field, list of missing variables)

        Raises
        ------
        FileNotFoundError
            If the file does not exist.
        NotImplementedError
            If the field type is not supported.
        """
        # Update the calculation mesh and purge cached properties if source changed 
        self._update_calc_mesh(field)
        self._update_data_source(filePath, groupName)
        
        # List of FEM spaces degrees in the field
        subFields               = field.getListOfSubFields() if field.info["type"] == "mixed" else [field]
        degrees                 = [subField.info.get("degree", None) for subField in subFields]
        
        # Check variable availability
        self._get_list_available_vars()
        wantedVars              = self._names_for_field(field)
        missingVars             = [var for var in wantedVars if var not in self._availableVars]
        presentVars             = [var for var in wantedVars if var in self._availableVars]
        if not presentVars:
            logger.warning(f"No variables for field '{field.name}' found in file '{filePath}'. Set to default values.")
            return field, missingVars

        # Load and process data
        if self._cacheData:
            if not self._dataForField:
                logger.debug(f"Cache empty, loading and interpolating/mapping all available data: {self._availableVars}")
                
                # Load all non-coordinate variables for this file
                coordsSet           = set(self.fullAxisNames)
                allVars             = [var for var in self._availableVars if var not in coordsSet]
                self._load_from_h5(allVars)

                # If a mixed space has a P1 variable, we need to interpolate the variables to P1
                # NOTE: This is not optimal if all subfields are P1, but this is a rare case?
                if field.info["type"] == "mixed" and 1 in degrees and not self._needInterpolation:
                    logger.debug("Mixed field with P1 subfield detected.")
                    variablesP1     = [name for name, deg in zip(field.getNamesOfSubFields(), degrees) if deg == 1]
                    # Dict of only P1 variables
                    logger.debug(f"Interpolating P1 variable(s): {variablesP1}")
                    processedP1Data = self._interpolate_to_calc_mesh(variablesP1, toP1=True)
                    
                    # Dict of other variables mapped to calc mesh
                    otherVars       = [var for var in allVars if var not in variablesP1]
                    processedOther  = self._map_to_calc_mesh(otherVars)
                    # Combine both
                    processedData   = {**processedP1Data, **processedOther}

                else:
                    processedData   = (
                        self._interpolate_to_calc_mesh(allVars) if self._needInterpolation else
                        self._map_to_calc_mesh(allVars)
                    )
                    
                # Update cached data
                self._dataForField.update(processedData)

            # Check if there are any new variables to load, missing from cached data
            extraVars           = [var for var in presentVars if var not in self._dataForField]
            if extraVars:
                logger.debug(f"Loading and interpolating/mapping extra variables: {extraVars}")
                self._load_from_h5(extraVars)
                updatedData     = (self._interpolate_to_calc_mesh(extraVars)
                                    if self._needInterpolation else
                                    self._map_to_calc_mesh(extraVars)
                                )
                self._dataForField.update(updatedData)
                
            logger.debug(f"Getting {presentVars} from cached data.")
            arrays              = {var: self._dataForField[var] for var in presentVars}

        else:
            # TODO: test this workflow
            # Nothing in cache yet
            self._load_from_h5(presentVars)
            arrays              = (self._interpolate_to_calc_mesh(presentVars)
                                    if self._needInterpolation else
                                    self._map_to_calc_mesh(presentVars)
                                )

        self._set_arrays_to_field(arrays, field)
        return field, missingVars
