import  os
import  h5py
import  numpy               as np
from    scipy.spatial       import Delaunay
from    scipy.interpolate   import LinearNDInterpolator, NearestNDInterpolator
from    FELiCS.Misc.logging import Logger
from    functools           import cached_property
from    .Mapping            import Mapping

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

    **Initialize the Reader object**

    Parameters
    ----------
    needInterpolation : bool, optional
        Whether to interpolate data onto the FELiCS mesh, by default True.
    originalMeshFile : str or None, optional
        Path to the original mesh file if coordinates are not in the main file, by default None.
    isComplex : bool, optional
        Whether the data is complex-valued, by default False.
    cacheData : bool, optional
        Whether to cache loaded data for reuse, by default True.

    Attributes
    ----------
    _needInterpolation : bool
        Indicates if interpolation is required.
    _originalMeshFile : str or None
        Path to the original mesh file.
    _isComplex : bool
        Indicates if data is complex-valued.
    _cacheData : bool
        Indicates if data caching is enabled.
    _param : object or None
        FELiCS parameter object, set via bind_env.
    _FEMSpaces : object or None
        FELiCS finite element spaces object, set via bind_env.
    _sourceKey : tuple or None
        Current data source as (filePath, groupName).
    _availableVars : list or None
        List of available variables in the current source.
    _rawDataDict : dict
        Dictionary of raw loaded data.
    _dataForField : dict
        Dictionary of processed data for fields.

    Notes
    -----
    - The class relies on external dependencies like param and FEMSpaces, which must be bound before use.
    - Caching is per-file and cleared when switching sources.
    - Interpolation uses Delaunay triangulation for efficiency.

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
        needInterpolation   = True,
        originalMeshFile    = None,
        isComplex           = False,  # NOTE: Should be an attribute of the field, not the reader?
        cacheData           = True,
    ):
        """
        Initialize the Reader instance.

        Parameters
        ----------
        needInterpolation : bool, optional
            Whether to interpolate data onto the FELiCS mesh, by default True.
        originalMeshFile : str or None, optional
            Path to the original mesh file if coordinates are not in the main file, by default None.
        isComplex : bool, optional
            Whether the data is complex-valued, by default False.
        cacheData : bool, optional
            Whether to cache loaded data for reuse, by default True.
        """
        # Configuration (fixed for this Reader instance)
        self._needInterpolation = needInterpolation
        self._originalMeshFile  = originalMeshFile
        self._isComplex         = isComplex
        self._cacheData         = cacheData

        # Environment dependencies (TODO: get rid of these dependencies?)
        self._param             = None
        self._FEMSpaces         = None
        self._mesh              = None

        # Per-file caches (reset when source changes)
        self._sourceKey         = None  # (filePath, groupName)
        self._availableVars     = None
        self._rawDataDict       = {}
        self._dataForField      = {}

    # --------------------------
    # Environment binding
    # --------------------------
    # NOTE: TEMPORARY solution until updates in Field are made
    # TODO: get rid of dependencies to these objects?
    def bind_env(self, param, FEMSpaces):
        """
        Bind the FELiCS environment objects to this Reader.

        Parameters
        ----------
        param : object
            The FELiCS parameter object.
        FEMSpaces : object
            The FELiCS finite element spaces object.
        """
        self._param     = param
        self._FEMSpaces = FEMSpaces
        self._mesh      = param.getMesh()

    # --------------------------
    # Cached mesh and mapping data
    # --------------------------
    @cached_property
    def meshDim(self):
        """
        Get the number of spatial dimensions in the mesh.
        TODO: rename property

        Returns
        -------
        list of str
            The full dimension names (e.g., ['x', 'y', 'z']).

        Raises
        ------
        RuntimeError
            If accessed before bind_env is called.
        """
        if self._mesh is None:
            logger.error("meshDim accessed before bind_env() which sets self._mesh")
            raise RuntimeError("meshDim accessed before bind_env() which sets self._mesh")
        return self._mesh.meshAxisNames
    
    @cached_property
    def fullDim(self):
        """
        Get the full number of spatial dimensions: mesh + spectral.
        TODO: rename property

        Returns
        -------
        list of str
            The full dimension names (e.g., ['x', 'y', 'z']).

        Raises
        ------
        RuntimeError
            If accessed before bind_env is called.
        """
        if self._mesh is None:
            raise RuntimeError("fullDim accessed before bind_env() which sets self._mesh")
        return self._mesh.axisNames

    @cached_property
    def _triangulationImportMesh(self):
        """
        Cached Delaunay triangulation of the import mesh for interpolation.

        Returns
        -------
        scipy.spatial.Delaunay
            The triangulation object for the import mesh.
        """
        nDimMesh = len(self.meshDim)
        importCoords = self.importMeshCoords[:, :nDimMesh]
        return Delaunay(importCoords)

    @cached_property
    def calcMeshCoords(self):
        """
        Get the coordinates of the calculation mesh (P2).
        TODO: get from the Field, but needs to be kept in memory for multiple fields or import files
        NOTE: Concept is 1 calc mesh per Reader instance

        Returns
        -------
        numpy.ndarray
            Array of shape (nPoints, nDimMesh) with the coordinates of the calculation mesh.

        Notes
        -----
        TODO: Might want to add P1 support as well.
        """
        logger.debug("Getting the calculation mesh coordinates.")
        nDimMesh = len(self.meshDim)
        return self._FEMSpaces.P2.tabulate_dof_coordinates()[:, :nDimMesh]

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
        
        # Number and names of mesh axes
        meshAxisNames                   = self.meshDim
        numMeshAxis                     = len(meshAxisNames)
        logger.debug(f"Getting the {numMeshAxis}D import mesh coordinates ({meshAxisNames}).")

        if self._needInterpolation:
            
            if self._originalMeshFile:
                # Option 1: load from original mesh file
                with h5py.File(self._originalMeshFile, "r") as fileHandle:
                    prefix                  = "coordinates/"  # NOTE: we want to get rid of groups in files
                    firstCoordShape         = fileHandle[prefix + meshAxisNames[0]][:].shape
                    coordsArray             = np.zeros((max(firstCoordShape), numMeshAxis))
                    for i, coordName in enumerate(meshAxisNames):
                        coordsArray[:, i]   = fileHandle[prefix + coordName][:]
                return coordsArray
            else:
                # Option 2: load from main file under group name
                # NOTE: we assume the coordinates are in the same group as the variables
                # Build from the first current file that has coords
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
            # Option 3: load from FELiCS exported mesh file
            meshFile = os.path.join(
                self._param.Export.ExportFolder,
                "mesh.h5"
            )
            with h5py.File(meshFile, "r") as fileHandle:
                prefix                  = "coordinates/"
                firstCoordShape         = fileHandle[prefix + "x"][:].shape
                coordsArray             = np.zeros((max(firstCoordShape), numMeshAxis))
                # NOTE: in FELiCS mesh the coordinates are always x,y,z! Needed for ParaView?
                for i, coordName in enumerate(["x", "y", "z"][:numMeshAxis]):
                    coordsArray[:, i]   = fileHandle[prefix + coordName][:]
            return coordsArray

    @cached_property
    def import_to_P2calc(self):
        """
        Get the mapping from import mesh to calculation mesh (P2).

        Returns
        -------
        numpy.ndarray or None
            The mapping array if no interpolation is needed, otherwise None.

        Notes
        -----
        TODO: Might want to add mapping to P1 as well.
        TODO: Use the mapping from Field once it's implemented.
        """
        if self._needInterpolation:
            return None
        else:
            nDimMesh = len(self.meshDim)
        return Mapping.calculateMappingFromDofs(
            self.importMeshCoords[:, :nDimMesh],
            self.calcMeshCoords
        )

    # --------------------------
    # Source management
    # --------------------------
    def _update_data_source(self, filePath, groupName):
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
        if not os.path.isfile(filePath):
            raise FileNotFoundError(filePath)

        # Set new source if different
        key = (os.path.abspath(filePath), groupName)
        # NOTE: might need completion of the variable list
        if key != self._sourceKey:
            self._sourceKey     = key
            self._availableVars = None
            self._rawDataDict   = {}
            self._dataForField  = {}
            logger.debug(f"Reader: switched source to {key}, cleared per-file caches.")

    def _ensure_source_set(self):
        """
        Ensure that a data source is set.

        Raises
        ------
        RuntimeError
            If no data source is set.
        """
        if self._sourceKey is None:
            raise RuntimeError("No data source set. Call importInField(..., filePath, groupName) first.")

    def clear_cached_properties(self):
        """
        Remove all cached @cached_property values from this Reader instance.

        They will be recomputed on the next access.
        
        TODO: Complete the method with all cached properties.
        """
        cachedProps = [
            "fullDim",
            "meshDim",
            "calcMeshCoords",
            "importMeshCoords",
            "import_to_P2calc",
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
    def _get_list_available_vars(self):
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
        with h5py.File(filePath, "r") as fileHandle:
            if groupName in fileHandle:
                self._availableVars = list(fileHandle[groupName].keys())
            else:
                logger.warning(f"Group '{groupName}' not found in file '{filePath}'.")
                self._availableVars = []

    def _load_from_h5(self, varNames):
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
                basePath            = f"{groupName}/{varName}"
                if self._isComplex:  # NOTE: not tested yet
                    # TODO: make _real and _imag in file
                    magnitude       = fileHandle[basePath + "/magnitude"][:].squeeze()
                    angle           = fileHandle[basePath + "/angle"][:].squeeze()
                    self._rawDataDict[varName + "_real"] = magnitude * np.cos(angle)
                    self._rawDataDict[varName + "_imag"] = magnitude * np.sin(angle)
                else:
                    if not self._needInterpolation and (basePath + "/magnitude") in fileHandle:
                        # If we don't interpolate, we load a FELiCS file with magnitude only
                        self._rawDataDict[varName] = fileHandle[basePath + "/magnitude"][:].squeeze()
                    else:
                        # Otherwise we load the full variable
                        self._rawDataDict[varName] = fileHandle[basePath][:].squeeze()

    # --------------------------
    # Interpolation / mapping
    # --------------------------
    def _interpolate_to_calc_mesh(self, varNames):
        """
        Interpolate given variables from import mesh to calculation mesh.

        Parameters
        ----------
        varNames : list of str
            List of variable names to interpolate.

        Returns
        -------
        dict
            Dictionary mapping variable names to interpolated arrays.

        Raises
        ------
        ValueError
            If variable lengths do not match the import mesh.
        """
        nMeshDim                    = len(self.meshDim)
        importCoords                = self.importMeshCoords[:, :nMeshDim]
        calcCoords                  = self.calcMeshCoords
        numImportPoints             = importCoords.shape[0]

        # Sanity check: all variables have same length as import mesh
        for varName in varNames:
            if self._rawDataDict[varName].shape[0] != numImportPoints:
                raise ValueError(
                    f"Variable '{varName}' has length {len(self._rawDataDict[varName])}, "
                    f"but import mesh has {numImportPoints} points."
                )

        # Stack data columns (N, K)
        values                      = np.column_stack([self._rawDataDict[varName] for varName in varNames])

        # Use cached triangulation
        linearInterp                = LinearNDInterpolator(self._triangulationImportMesh, values)
        interpolated                = linearInterp(calcCoords)  # shape (M, K)

        # Handle NaNs by nearest-neighbour fill (vectorized)
        if np.isnan(interpolated).any():
            nearestInterp           = NearestNDInterpolator(importCoords, values)
            nearestValues           = nearestInterp(calcCoords)
            nanMask                 = np.isnan(interpolated)
            interpolated[nanMask]   = nearestValues[nanMask]

        # Return a dict {var_name: array}
        return {varName: interpolated[:, i] for i, varName in enumerate(varNames)}

    def _map_to_calc_mesh(self, varNames):
        """
        Map given variables using precomputed index mapping.

        Parameters
        ----------
        varNames : list of str
            List of variable names to map.

        Returns
        -------
        dict
            Dictionary mapping variable names to mapped arrays.
        """
        
        indices = self.import_to_P2calc
        
        # If we have complex data, we have _real and _imag in _rawDataDict
        # In which case we reconstruct the complex array here
        if self._isComplex:
            output = {}
            for varName in varNames:
                realPart                = self._rawDataDict.get(varName + "_real", None)
                imagPart                = self._rawDataDict.get(varName + "_imag", None)
                if realPart is not None and imagPart is not None:
                    complexArray        = realPart + 1j * imagPart
                    output[varName]     = complexArray[indices]
                else:
                    logger.warning(f"Variable '{varName}' not found as complex in raw data. Skipping.")
            return output
        
        else:
            return {varName: self._rawDataDict[varName][indices] for varName in varNames}

    # --------------------------
    # Field helpers
    # --------------------------
    def _names_for_field(self, field):
        """
        Determine variable names to load from file for this Field.

        Parameters
        ----------
        Field : object
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
                    names.append(component)
                    
            elif subInfo["type"] == "scalar":
                names.append(subField.getName())
                
            else:
                logger.error(f"Field type '{subInfo['type']}' not supported in Reader yet.")
                raise NotImplementedError("Field type not supported in Reader yet.")
            
        return names

    def _set_arrays_to_field(self, arrays, Field):
        """
        Assign interpolated/mapped arrays to a given Field.

        Parameters
        ----------
        arrays : dict
            Dictionary of arrays to assign.
        Field : object
            The FELiCS field object.

        Returns
        -------
        object
            The updated Field.

        Raises
        ------
        NotImplementedError
            For mixed function spaces.
        """
        # Field info
        info        = Field.describeFunctionSpace()
            
        # Set the arrays in FEM depending on field type
        if info["type"] == "mixed":
            # Loop over the subfields
            for iField, subFieldName in enumerate(Field.getNamesOfSubFields()):
                # If subfield is a vector, loop over its components
                if info['subspaces'][iField]['type'] == 'vector':
                    subFields                               = Field.getListOfSubFields()
                    for jComp, compName in enumerate(subFields[iField].getNamesOfSubFields()):
                        if compName in arrays:
                            indices                         = Field.space.sub(iField).sub(jComp).collapse()[1]
                            Field.function.x.array[indices] = arrays[compName]
                        else:
                            logger.warning(f"Component '{compName}' not found in loaded arrays for vector subfield '{subFieldName}'. Leaving unchanged.")
                # For a scalar subfield, just set the array
                elif info['subspaces'][iField]['type'] == 'scalar':
                    if subFieldName in arrays:
                        indices                             = Field.space.sub(iField).collapse()[1]
                        Field.function.x.array[indices]     = arrays[subFieldName]
                    else:
                        logger.warning(f"Subfield '{subFieldName}' not found in loaded arrays for scalar subfield. Leaving unchanged.")
                else:
                    logger.error(f"Subfield type '{info['subspaces'][iField]['type']}' not supported in Reader yet.")
                    raise NotImplementedError("Subfield type not supported in Reader yet.")
                
        elif info["type"] == "vector":
            for i, subFieldName in enumerate(Field.getNamesOfSubFields()):
                if subFieldName in arrays:
                    indices                         = Field.space.sub(i).collapse()[1]
                    Field.function.x.array[indices] = arrays[subFieldName]
                else:
                    logger.warning(f"Component '{subFieldName}' not found in loaded arrays for vector field '{Field.getName()}'. Leaving unchanged.")
        else:
            # Then it's a scalar
            varName                                 = Field.getName()
            if varName in arrays:
                Field.function.x.array[:]           = arrays[varName]
        return Field

    # --------------------------
    # Main API
    # --------------------------
    def importInField(self, Field, filePath, groupName):
        """
        Import data from the specified file/group into the given Field.

        Automatically clears per-file caches when file/group changes.

        Parameters
        ----------
        Field : object
            The FELiCS field to import data into.
        filePath : str
            Path to the HDF5 file.
        groupName : str
            Name of the group in the HDF5 file.

        Returns
        -------
        tuple
            (updated_Field, list_of_missing_variables)

        Raises
        ------
        RuntimeError
            If environment is not bound.
        """
        if self._param is None or self._FEMSpaces is None:
            raise RuntimeError("Call bind_env(param, FEMSpaces) before importInField().")

        self._update_data_source(filePath, groupName)
        self._get_list_available_vars()

        wantedVars  = self._names_for_field(Field)
        missingVars = [var for var in wantedVars if var not in self._availableVars]
        presentVars = [var for var in wantedVars if var in self._availableVars]
        if not presentVars:
            logger.warning(f"No variables for Field '{Field.getName()}' found in file '{filePath}'. Set to default values.")
            return Field, missingVars

        # Load and process data
        if self._cacheData:
            if not self._dataForField:
                logger.debug("Cache empty, loading and interpolating/mapping all available data.")
                # Load and interpolate/map all non-coordinate variables for this file
                coordsSet       = set(self.fullDim)
                allVars         = [var for var in self._availableVars if var not in coordsSet]
                self._load_from_h5(allVars)
                processedData   = (self._interpolate_to_calc_mesh(allVars)
                                    if self._needInterpolation else
                                    self._map_to_calc_mesh(allVars)
                                )
                # Update cached data
                self._dataForField.update(processedData)
            else:
                logger.debug("Using cached data.")

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
            arrays              = {var: self._dataForField[var] for var in presentVars}

        else:
            # TODO: test this workflow
            # Nothing in cache yet
            self._load_from_h5(presentVars)
            arrays              = (self._interpolate_to_calc_mesh(presentVars)
                                    if self._needInterpolation else
                                    self._map_to_calc_mesh(presentVars)
                                )

        self._set_arrays_to_field(arrays, Field)
        return Field, missingVars
