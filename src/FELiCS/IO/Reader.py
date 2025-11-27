# Standard libraries
import  os
import  h5py
from    functools           import cached_property
from    typing              import Optional, List, Dict, Any, Tuple
from    pathlib             import Path

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
    Reader for importing and interpolating data from HDF5 files onto FELiCS meshes.

    This class manages loading of data from HDF5 files, interpolation or direct mapping
    to a calculation mesh, and assignment of the resulting arrays to FELiCS fields.
    It supports optional caching for performance and can handle complex-valued data.

    **Initialize the Reader object**

    Parameters
    ----------
    sourceDir : str, optional
        Base directory where the data files are located. This is used as the
        default directory when resolving relative HDF5 file paths.
    needInterpolation : bool, optional
        If True, data are interpolated from the import mesh onto the calculation
        mesh using a Delaunay-based interpolator. If False, a precomputed
        mapping between mesh DoFs is used instead. Default is True.
    felicsMeshFilePath : str or None, optional
        Path to the FELiCS mesh file when the mesh is stored in a different
        directory than the source data. If None, a default ``mesh.h5`` in
        ``sourceDir`` is assumed when interpolation is not needed.
    isComplex : bool, optional
        If True, the Reader expects complex-valued data stored as separate
        ``_real`` and ``_imag`` datasets and reconstructs complex arrays.
        Default is False.
    cacheData : bool, optional
        If True, processed data (after interpolation or mapping) are cached
        per file/group to speed up repeated imports into different fields.
        Default is True.

    Attributes
    ----------
    _sourceDir : str
        Source directory for data files.
    _needInterpolation : bool
        Flag indicating whether interpolation onto the FELiCS mesh is required.
    _felicsMeshFilePath : str or None
        Path to the FELiCS mesh file, if explicitly provided.
    _isComplex : bool
        Flag indicating whether data are treated as complex-valued.
    _cacheData : bool
        Flag indicating whether processed data are cached per file/group.
    _sourceKey : tuple of (str, str) or None
        Current data source identifier as ``(filePath, groupName)``.
    _availableVars : list of str or None
        List of variable names available in the current HDF5 group.
    _rawDataDict : dict
        Raw data loaded from HDF5 datasets, keyed by variable name.
    _dataForField : dict
        Processed (interpolated or mapped) data ready to be assigned to fields.
    _calcMeshCoords : numpy.ndarray or None
        Coordinates of the current calculation mesh DoFs (if available).
    _calcMeshP1Coords : numpy.ndarray or None
        Coordinates of the P1 calculation mesh DoFs (if available).
    _mappingImportToCalc : numpy.ndarray or None
        Precomputed mapping from import mesh indices to calculation mesh DoFs,
        used when interpolation is disabled.
    _fullAxisNames : list of str or None
        Names of all axes (mesh + spectral) for the current FELiCS mesh.
    _meshAxisNames : list of str or None
        Names of the spatial mesh axes for the current FELiCS mesh.
    _numDoFsCalcMeshP2 : int
        Number of DoFs in the P2 calculation mesh for the current field.
    _numDoFsCalcMeshP1 : int
        Number of DoFs in the P1 calculation mesh for the current field.
    _spacesDegrees : list of int or None
        Polynomial degrees of the FEM spaces of all subfields for the current field.

    Notes
    -----
    - Caching is maintained per file/group and is cleared whenever the data
      source changes.
    - Interpolation uses a Delaunay triangulation of the import mesh and
      performs a nearest-neighbour fill where linear interpolation produces
      NaN values.
    - The implementation is retro-compatible with older FELiCS HDF5 layouts,
      including optional grouping and magnitude datasets.
    - Currently only P1 and P2 Lagrange FEM spaces are supported.

    Raises
    ------
    FileNotFoundError
        If a specified mesh or data file does not exist.
    RuntimeError
        If methods that require a bound source are called before any source
        is set.
    NotImplementedError
        If unsupported field or subfield types (e.g., mixed spaces beyond
        the supported structure) are encountered.
    """
    def __init__(
        self,
        sourceDir:          str  = "",
        needInterpolation:  bool = True,
        felicsMeshFilePath: Optional[str] = None,
        isComplex:          bool = False,
        cacheData:          bool = True,
    ) -> None:
        """
        Initialize the Reader instance.

        Parameters
        ----------
        sourceDir : str, optional
            Base directory where the HDF5 data files are located.
            Relative file paths passed to :meth:`importInField` are
            resolved with respect to this directory. Default is ``""``.
        needInterpolation : bool, optional
            If True, variables are interpolated from the import mesh to
            the calculation mesh using a Delaunay-based interpolator.
            If False, a mapping between mesh DoFs is used instead.
            Default is True.
        felicsMeshFilePath : str or None, optional
            Path to a FELiCS mesh file to be used as the import mesh when
            it is not located in the same directory as the HDF5 data.
            If None, a default ``mesh.h5`` is assumed in ``sourceDir``
            when interpolation is disabled. Default is None.
        isComplex : bool, optional
            If True, the Reader expects complex-valued variables stored as
            ``<name>_real`` and ``<name>_imag`` datasets and reconstructs
            complex arrays from them. Default is False.
        cacheData : bool, optional
            If True, processed data (after interpolation or mapping) are
            cached per file/group and reused across subsequent calls to
            :meth:`importInField`. Default is True.
        """        # Configuration (fixed for this Reader instance)
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
        Names of the spatial mesh dimensions.

        Returns
        -------
        list of str
            Names of the mesh axes, e.g. ``['x', 'y', 'z']`` for a 3D mesh.

        Raises
        ------
        RuntimeError
            If accessed before any field is bound and the mesh environment
            has not been initialized.
        """
        return self._meshAxisNames
    
    @cached_property
    def fullAxisNames(self):
        """
        Names of all spatial dimensions (mesh + spectral).

        This includes both the mesh coordinates and any additional
        spectral or parametric axes associated with the FELiCS mesh.

        Returns
        -------
        list of str
            The full list of axis names, e.g. ``['x', 'y', 'z']`` or
            ``['x', 'y', 'z', 'k']``.

        Raises
        ------
        RuntimeError
            If accessed before any field is bound and the mesh environment
            has not been initialized.
        """
        return self._fullAxisNames

    @cached_property
    def _triangulationImportMesh(self):
        """
        Cached Delaunay triangulation of the import mesh.

        The triangulation is built in the space of mesh coordinates
        (excluding spectral axes) and reused for interpolation of multiple
        variables.

        Returns
        -------
        scipy.spatial.Delaunay
            Delaunay triangulation of the import mesh coordinates.
        """
        nDimMesh        = len(self.meshAxisNames)
        importCoords    = self.importMeshCoords[:, :nDimMesh]
        return Delaunay(importCoords)

    @cached_property
    def calcMeshP2Coords(self):
        """
        Coordinates of the P2 calculation mesh.

        These coordinates correspond to the DoFs of the highest-degree
        (P2) scalar or vector subspace present in the field.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(n_points, n_dim_mesh)`` containing the
            coordinates of the P2 calculation mesh DoFs.
        """
        return self._calcMeshP2Coords
    
    @cached_property
    def calcMeshP1Coords(self):
        """
        Coordinates of the P1 calculation mesh.

        These coordinates correspond to the DoFs of the P1 scalar or
        vector subspace present in the field.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(n_points, n_dim_mesh)`` containing the
            coordinates of the P1 calculation mesh DoFs.
        """
        return self._calcMeshP1Coords

    @cached_property
    def mappingImportToCalcP2(self):
        """
        Mapping from import mesh indices to P2 calculation mesh DoFs.

        When interpolation is disabled, this mapping is used to transfer
        variables defined on the import mesh directly to the P2 calculation
        mesh. When interpolation is enabled, no mapping is used.

        Returns
        -------
        numpy.ndarray or None
            One-dimensional array of indices mapping each P2 DoF to a
            point in the import mesh, or None if interpolation is enabled.
        """
        if self._needInterpolation:
            return None
        else:
            return self._mappingImportToCalc

    @cached_property
    def importMeshCoords(self):
        """
        Coordinates of the import mesh.

        Depending on the configuration, the import mesh is obtained in one
        of the following ways:

        1. If interpolation is needed and ``felicsMeshFilePath`` is set,
           coordinates are loaded from the external FELiCS mesh file.
        2. If interpolation is needed and no external mesh file is given,
           coordinates are loaded from datasets co-located with the
           variables in the main HDF5 file.
        3. If interpolation is not needed, coordinates are loaded from a
           FELiCS mesh file (either the explicit path or a default
           ``mesh.h5`` in ``sourceDir``).

        Returns
        -------
        numpy.ndarray
            Array of shape ``(n_points, n_dim_mesh)`` containing the
            coordinates of the import mesh.

        Raises
        ------
        FileNotFoundError
            If the expected FELiCS mesh file does not exist.
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
        Switch to a new HDF5 file/group and clear per-file caches if needed.

        If the requested ``(filePath, groupName)`` differs from the current
        source, the internal caches of available variables, raw data, and
        processed data are reset.

        Parameters
        ----------
        filePath : str
            Absolute or relative path to the HDF5 file containing the data.
        groupName : str
            Name of the HDF5 group in which the variables are stored.
            Use ``None`` to access datasets at the root level.

        Raises
        ------
        FileNotFoundError
            If the specified file does not exist.
        """

        # Check that the file exists
        if not os.path.isfile(filePath):
            raise FileNotFoundError(filePath)

        # Check that the source directory fits that of the reader instance
        currentSourceDir         = os.path.dirname(filePath)
        if os.path.abspath(currentSourceDir) != os.path.abspath(self._sourceDir):
            # NOTE: maybe we do something specific in that case?
            #logger.warning(
            #    f"Reader source directory '{currentSourceDir}' does not match "
            #    f"the Reader instance sourceDir '{self._sourceDir}': expected behavior issues."
            #)
            pass

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
        Ensure that a data source has been configured.

        This checks whether ``_sourceKey`` has been set by a previous call
        to :meth:`importInField` or a related method.

        Raises
        ------
        RuntimeError
            If no data source is currently set.
        """
        if self._sourceKey is None:
            raise RuntimeError("No data source set. Call importInField(..., filePath, groupName) first.")

    def clear_cached_properties(
        self, 
        cachedProps: List[str] = []
    ) -> None:
        """
        Clear cached properties from this Reader instance.

        Any attributes corresponding to ``@cached_property`` results are
        removed from ``__dict__`` so that they are recomputed upon the next
        access.

        Parameters
        ----------
        cachedProps : list of str, optional
            Names of cached properties to clear explicitly. If an empty
            list is given, a default set of mesh-related properties
            (``fullAxisNames``, ``meshAxisNames``, ``calcMeshP2Coords``,
            ``importMeshCoords``, ``mappingImportToCalcP2``) is cleared.
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
        Discover available variable names in the current HDF5 group.

        The method populates ``_availableVars`` based on the current
        ``_sourceKey``. Some known non-field datasets (e.g. 'omega',
        'gain', 'number') are removed from the list.

        Raises
        ------
        RuntimeError
            If no data source has been set prior to calling this method.
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
        Load selected variables from the current HDF5 source.

        For each requested variable name, the corresponding dataset is read
        from the HDF5 file and stored in ``_rawDataDict``. Some legacy
        layouts store the data inside an additional ``magnitude`` subgroup,
        which is handled transparently.

        Parameters
        ----------
        varNames : list of str
            Names of the variables to load from the HDF5 file. Variables
            already present in ``_rawDataDict`` are skipped.

        Notes
        -----
        The interpretation of complex data (``_real`` and ``_imag`` suffixes)
        is handled at a later stage; here only the raw real-valued datasets
        are loaded.
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
        Interpolate variables from the import mesh to the calculation mesh.

        The interpolation is performed using a linear interpolator based on
        the Delaunay triangulation of the import mesh. Where linear
        interpolation yields NaN values, a nearest-neighbour interpolator is
        used to fill in missing values. Complex-valued variables are
        reconstructed from their ``_real`` and ``_imag`` components.

        Parameters
        ----------
        varNames : list of str
            Names of the variables to interpolate. For complex data, the
            corresponding ``_real`` and ``_imag`` datasets are used.
        toP1 : bool, optional
            If True, variables are interpolated to the P1 calculation mesh
            coordinates; otherwise to the P2 calculation mesh. Default is False.

        Returns
        -------
        dict
            Dictionary mapping each requested variable name to its
            interpolated array on the chosen calculation mesh.

        Raises
        ------
        ValueError
            If any variable length does not match the number of points in
            the import mesh.
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
        Map variables from the import mesh to the P2 calculation mesh.

        This method uses a precomputed index mapping between the import
        mesh and the P2 calculation mesh (see
        :attr:`mappingImportToCalcP2`) to transfer variables without
        interpolation.

        Parameters
        ----------
        varNames : list of str
            Names of the variables to map using the precomputed indices.

        Returns
        -------
        dict
            Dictionary mapping each variable name to its mapped array on
            the P2 calculation mesh.
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
        Determine the HDF5 variable names corresponding to a given field.

        For mixed and vector fields, this method expands the field into its
        scalar components and returns the corresponding dataset names. For
        complex-valued data, ``_real`` and ``_imag`` suffixes are appended
        as needed.

        Parameters
        ----------
        field : object
            FELiCS field whose components should be mapped to variable
            names in the HDF5 file.

        Returns
        -------
        list of str
            List of variable names that should be loaded for this field.

        Raises
        ------
        NotImplementedError
            If the field or subfield type is not supported by the Reader.
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
        Assign processed arrays to a FELiCS field.

        The arrays are written into the underlying FEM function associated
        with the field, taking into account mixed spaces, vector components
        and scalar subfields. For complex-valued data, real and imaginary
        parts are combined before assignment.

        Parameters
        ----------
        arrays : dict
            Dictionary mapping variable names to numpy arrays on the
            calculation mesh. For complex data, arrays may be provided as
            ``<name>_real`` and ``<name>_imag`` entries.
        field : object
            FELiCS field to be populated with the provided data.

        Returns
        -------
        object
            The same field instance, with its underlying FEM function
            updated in-place.

        Raises
        ------
        NotImplementedError
            If the field or subfield type is not supported by the Reader.
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
        Import data from an HDF5 file/group into a FELiCS field.

        This is the main public entry point of the Reader. It performs the
        following steps:

        1. Normalizes the file path (adding ``.h5`` extension if needed).
        2. Updates the cached calculation mesh for the given field.
        3. Switches the internal data source and clears per-file caches
           if the file/group has changed.
        4. Determines which variables are required for the field and checks
           their availability in the file.
        5. Loads variables from HDF5 and either interpolates them or maps
           them to the calculation mesh, optionally using per-file caching.
        6. Assigns the processed arrays to the field.

        Parameters
        ----------
        field : object
            FELiCS field object to which the imported data will be written.
        filePath : str
            Path to the HDF5 file (with or without ``.h5``/``.fel`` extension).
        groupName : str
            Name of the group within the HDF5 file containing the variables.
            Use ``None`` for root-level datasets.

        Returns
        -------
        tuple
            A tuple ``(field, missingVars)`` where:
            - ``field`` is the updated FELiCS field object,
            - ``missingVars`` is a list of variable names that were expected
              for the field but not found in the file.

        Raises
        ------
        FileNotFoundError
            If the specified file does not exist.
        NotImplementedError
            If the field type or its subspaces are not supported.
        """

        # Add ".h5" to the file path it not specified
        # TODO Sophie: ask Simon
        fullFilePath = Path(filePath)
        if fullFilePath.suffix != ".h5" and fullFilePath.suffix != ".fel":
            fullFilePath = str(fullFilePath.with_name(fullFilePath.stem + ".h5"))

        # Update the calculation mesh and purge cached properties if source changed 
        self._update_calc_mesh(field)
        self._update_data_source(fullFilePath, groupName)
        
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
