#  ___________________________________   _______________________________________________________
# /-----------------------------------\ /-------------------------------------------------------\
# |   (         (                (     |  This source code is part of FELiCS                     |
# |   )\ )     ) )        (      )\ )  |  (F)inite (E)lement (Li)nearized (C)ombustion (S)olver  |  
# |  (()/(  (  (()/( (    )\   (()/(   |                                                         |  
# |  /(_)) )\  /(_)))\  (((_)  /(_))   |  Licensed under the GNU GPLv3                           |
# |  (_)_)((_) (_)) ((_) )\___ (_))    |                                                         |
# |  | __|| __|| |   (_)((/ __|/ __|   |  (C) 2018-2025: The FELiCS Developers (www.felics.eu)   |
# |  | _| | _| | |__ | | | (__ \__ \   |  Visit          www.felics.eu                           |
# |  |_|  |___||____||_|  \___||___/   |  Contact        info@felics.eu                          |
# \___________________________________/ \_______________________________________________________/
#
# Standard libraries
from    functools   import cached_property
import  h5py
import  os
from    pathlib     import Path
import  time
from    typing      import Optional, List, Dict, Any, Tuple


# Third party libraries
import  numpy               as np
from    scipy.spatial       import Delaunay
from    scipy.interpolate   import (
    LinearNDInterpolator,
    NearestNDInterpolator,
)

# Local libraries and methods
from    FELiCS.Misc.logging import Logger, log_and_raise
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
        sourceDir:          str  = "",
        needInterpolation:  bool = True,
        felicsMeshFilePath: Optional[str] = None,
        cacheData:          bool = True
    ) -> None:
        """
        Initialize the Reader instance.

        Parameters
        ----------
        needInterpolation : bool, optional
            Whether to interpolate data onto the FELiCS mesh, by default True.
        felicsMeshFilePath : str or None, optional
            Path to the FELiCS mesh file WHEN it is in a different directory than the source data.
        cacheData : bool, optional
            Whether to cache loaded data for reuse, by default True.
        """
        # Configuration (fixed for this Reader instance)
        self._sourceDir             = sourceDir
        self._needInterpolation     = needInterpolation
        self._felicsMeshFilePath    = felicsMeshFilePath
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
    def mesh_axis_names(
        self
    ):
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
    def full_axis_names(
        self
    ):
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
    def _triangulation_import_mesh(
        self
    ):
        """
        Cached Delaunay triangulation of the import mesh for interpolation.

        Returns
        -------
        scipy.spatial.Delaunay
            The triangulation object for the import mesh.
        """
        nDimMesh        = len(self.mesh_axis_names)
        importCoords    = self.import_mesh_coords[:, :nDimMesh]
        return Delaunay(importCoords)

    @cached_property
    def calc_mesh_p2_coords(
        self
    ):
        """
        Get the coordinates of the calculation mesh.

        Returns
        -------
        numpy.ndarray
            Array of shape (nPoints, nDimMesh) with the coordinates of the calculation mesh.
        """
        return self._calcMeshP2Coords
    
    @cached_property
    def calc_mesh_p1_coords(
        self
    ):
        """
        Get the coordinates of the P1 calculation mesh.

        Returns
        -------
        numpy.ndarray
            Array of shape (nPoints, nDimMesh) with the coordinates of the P1 calculation mesh.
        """
        return self._calcMeshP1Coords

    @cached_property
    def mapping_import_to_calc_p2(
        self
    ):
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
    def import_mesh_coords(
        self
    ):
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
        mesh_axis_names     = self.mesh_axis_names
        numMeshAxis         = len(mesh_axis_names)
        logger.debug(f"Getting the {numMeshAxis}D import mesh coordinates ({mesh_axis_names}).")

        def _load_coordinates_from_file(
            file_path: str,
            prefix: str,
            coord_names: List[str]
        ) -> np.ndarray:
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
            with h5py.File(
                file_path,
                "r",
            ) as fh:
                first_shape         = fh[prefix + coord_names[0]][:].shape
                coords              = np.zeros((max(first_shape), len(coord_names)))
                for i, name in enumerate(coord_names):
                    coords[:, i]    = fh[prefix + name][:]
            return coords
        
        if self._needInterpolation:
            if self._felicsMeshFilePath:

                # Option 1: load from a given FELiCS mesh file
                logger.debug(f"Executing Option 1: load from a given FELiCS mesh file. {self._felicsMeshFilePath}")
                coordsArray = _load_coordinates_from_file(
                    self._felicsMeshFilePath,
                    "coordinates/",
                    mesh_axis_names,
                )
                return coordsArray
            
            else:

                # Option 2: Coords into file containing the data. Load at first file access
                # NOTE: assumes the coordinates are in the same group as the variables
                logger.debug("Executing Option 2: loading from main file under group name.")
                self._ensure_source_set()
                self._get_list_available_vars()
                missingCoords               = [coord for coord in mesh_axis_names if coord not in self._rawDataDict.keys()]
                if missingCoords:
                    self._load_from_h5(missingCoords)
                firstCoordShape             = self._rawDataDict[mesh_axis_names[0]].shape
                coordsArray                 = np.zeros((max(firstCoordShape), len(mesh_axis_names)))
                for i, coordName in enumerate(mesh_axis_names):
                    coordsArray[:, i]       = self._rawDataDict[coordName]
                return coordsArray
        
        else:
            logger.debug("Executing Option 3: No interpolation needed: load from FELiCS exported mesh file.")

            # Option 3: Load the FELiCS mesh
            if self._felicsMeshFilePath is not None:

                # For example when we load the mean flow, the FELiCS mesh is usually in another directory
                logger.debug(f"Using {self._felicsMeshFilePath}")
                meshFile                    = self._felicsMeshFilePath
            else:

                # Default FELiCS mesh file in the source directory (e.g. when loading modes)
                logger.debug(f"self._felicsMeshFilePath is None, therefore using {os.path.join(self._sourceDir, 'mesh.h5')}")
                meshFile                    = os.path.join(
                    self._sourceDir,
                    "mesh.h5",
                )
            
            # Check that the file exists
            if not os.path.isfile(meshFile):
                logger.error(f"FELiCS mesh file for import not found: {meshFile}")
                raise FileNotFoundError(meshFile)

            coordsArray = _load_coordinates_from_file(
                meshFile,
                "coordinates/",
                ["x", "y", "z"][:numMeshAxis],
            )
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
        subFields               = field.get_list_of_sub_fields() if fieldType == "mixed" else [field]
        degrees                 = [subField.info.get(
            "degree",
            None,
        ) for subField in subFields]
        if max(degrees) > 2:
            log_and_raise(logger, "Reader currently only supports P1 and P2 FEM spaces.", NotImplementedError)

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
            self._meshAxisNames     = felicsMesh.mesh_axis_names
            self._fullAxisNames     = felicsMesh.axis_names
            
            # If we have sub-spaces, get the indices of the sub-space DoFs
            if fieldType == "mixed":

                # TODO: check that this does not mess up with the DoF indices
                if numDofsForP2 > 0:
                    if field.get_list_of_sub_fields()[p2SubFieldIndex].info['type'] == 'vector':
                        self._calcMeshP2Coords  = field.space.sub(p2SubFieldIndex).sub(0).collapse()[0].tabulate_dof_coordinates()[:, :len(self._meshAxisNames)]
                    else:
                        self._calcMeshP2Coords  = field.space.sub(p2SubFieldIndex).collapse()[0].tabulate_dof_coordinates()[:, :len(self._meshAxisNames)]
                if numDofsForP1 > 0:
                    if field.get_list_of_sub_fields()[p1SubFieldIndex].info['type'] == 'vector':
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
                self._mappingImportToCalc       = Mapping.calculate_mapping_from_dofs(
                    self.import_mesh_coords[:, :numMeshDimensions],
                    self.calc_mesh_p2_coords,
                )

    def _ensure_source_set(
        self
    ) -> None:
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
                "mappingImportToCalcP2"
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
    def _get_list_available_vars(
        self
    ) -> None:
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
        with h5py.File(
            filePath,
            "r",
        ) as fileHandle:
            if groupName is None:
                self._availableVars = list(fileHandle.keys())
            elif groupName in fileHandle:
                self._availableVars = list(fileHandle[groupName].keys())
            else:
                logger.warning(f"Group '{groupName}' not found in file '{filePath}'.")
                logger.warning("Switching to root group for variable listing.")
                logger.debug(f"fileHandle contains the following keys: {list(fileHandle.keys())}")
                self._availableVars = list(fileHandle.keys())
                logger.debug(f"now overwriting sourceKey: changing groupName to None")
                self._sourceKey = (filePath, None)
                if self._availableVars is []:
                    logger.error(f"No variables found in file '{filePath}'.")
                    raise RuntimeError(f"No variables found in file '{filePath}'.")
                else:
                    logger.debug(f"Found variables: {self._availableVars}")
                
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
        """
        filePath, groupName         = self._sourceKey
        with h5py.File(
            filePath,
            "r",
        ) as fileHandle:
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
        numMeshDimensions           = len(self.mesh_axis_names)
        importCoords                = self.import_mesh_coords[:, :numMeshDimensions]
        calcCoords                  = self.calc_mesh_p1_coords if toP1 else self.calc_mesh_p2_coords
        numImportPoints             = importCoords.shape[0]
        
        varNamesComplex             = varNames

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
        logger.debug("Now performing linear interpolation using chached Delaunay triangulation.")

        # if len(self.meshAxisNames) == 3:
        #     logger.warning("PRELIMINARY FIX: FOR 3D interpolation nearest-neighbour interpolation is used.")
        #     nearestInterp           = NearestNDInterpolator(self._triangulationImportMesh, values)
        #     interpolated           = nearestInterp(calcCoords)
        # else:
        start_inter = time.time()
        linearInterp                = LinearNDInterpolator(
            self._triangulation_import_mesh,
            values,
        )
        interpolated                = linearInterp(calcCoords)  # shape (M, K)
        time_inter = time.time() - start_inter
        logger.info(f"Linear interpolation took {time_inter:.1f} seconds. (Mesh size: {calcCoords.shape})")

        # Handle NaNs by nearest-neighbour fill (vectorized)
        if np.isnan(interpolated).any():
            logger.debug("NaNs detected in interpolation result; applying nearest-neighbour fill.")
            nearestInterp           = NearestNDInterpolator(
                importCoords,
                values,
            )
            nearestValues           = nearestInterp(calcCoords)
            nanMask                 = np.isnan(interpolated)
            interpolated[nanMask]   = nearestValues[nanMask]
            
        # Output dict
        output                      = {varName: interpolated[:, i] for i, varName in enumerate(varNamesComplex)}
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
        
        indices                         = self.mapping_import_to_calc_p2
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
        
        subFields       = field.get_list_of_sub_fields() if info["type"] == "mixed" else [field]
        
        for subField in subFields:
            subInfo     = subField.info
            
            if subInfo["type"] == "vector":
                for component in subField.get_names_of_sub_fields():
                    names.append(component)
                    
            elif subInfo["type"] == "scalar":
                names.append(subField.name)
                
            else:
                log_and_raise(logger, f"Field type '{subInfo['type']}' not supported in Reader yet.", NotImplementedError)
        return names

    def _set_arrays_to_field(
        self,
        arrays:  Dict[str, Any],
        field: Any,
        typeVars: List[Tuple[str, str]]
    ) -> Any:
        """
        Assign interpolated/mapped arrays to a given field.

        Parameters
        ----------
        arrays : dict
            Dictionary of arrays to assign.
        field : object
            The FELiCS field object.
        typeVars : list of tuple
            List of (variable name, type) tuples indicating if variable is real or complex.

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
        
        def _assemble_complex_arrays(
            varName: str,
            arrays: Dict[str, Any],
            typeVars,
        ) -> Any:
            """Helper to assemble complex arrays from real and imaginary parts. 
            Or just return real array."""
            
            # Check in typeVars if variable is complex
            currentVarType      = typeVars[[v[0] for v in typeVars].index(varName)][1]
            
            if currentVarType == "complex":
                assembledArray  = arrays[varName + "_real"] + 1j * arrays[varName + "_imag"]
            else:
                assembledArray  = arrays[varName]
            return assembledArray
        
        # List of names available in arrays (without _real/_imag)
        baseNames = list({var[:-5] if var.endswith(("_real", "_imag")) else var for var in list(arrays.keys())})
            
        # Set the arrays in FEM depending on field type
        if info["type"] == "mixed":

            # Loop over the subfields
            for iField, subFieldName in enumerate(field.get_names_of_sub_fields()):

                # If subfield is a vector, loop over its components
                if info['subspaces'][iField]['type'] == 'vector':
                    subFields                               = field.get_list_of_sub_fields()
                    for jComp, compName in enumerate(subFields[iField].get_names_of_sub_fields()):
                        if compName in baseNames:
                            indices                         = field.space.sub(iField).sub(jComp).collapse()[1]
                            field.function.x.array[indices] = _assemble_complex_arrays(
                                compName,
                                arrays,
                                typeVars,
                            )
                        else:
                            logger.warning(f"Component '{compName}' not found in loaded arrays for vector subfield '{subFieldName}'. Set to default values.")

                # For a scalar subfield, just set the array
                elif info['subspaces'][iField]['type'] == 'scalar':
                    if subFieldName in baseNames:
                        indices                             = field.space.sub(iField).collapse()[1]
                        field.function.x.array[indices]     = _assemble_complex_arrays(
                            subFieldName,
                            arrays,
                            typeVars,
                        )
                    else:
                        logger.warning(f"Subfield '{subFieldName}' not found in loaded arrays for scalar subfield. Set to default values.")
                else:
                    log_and_raise(logger, f"Subfield type '{info['subspaces'][iField]['type']}' not supported in Reader yet.", NotImplementedError)
        elif info["type"] == "vector":
            for i, subFieldName in enumerate(field.get_names_of_sub_fields()):
                if subFieldName in baseNames:
                    indices                         = field.space.sub(i).collapse()[1]
                    field.function.x.array[indices] = _assemble_complex_arrays(
                        subFieldName,
                        arrays,
                        typeVars,
                    )
                else:
                    logger.warning(f"Component '{subFieldName}' not found in loaded arrays for vector field '{field.name}'. Set to default values.")
        else:

            # Then it's a scalar
            varName                                 = field.name
            if varName in baseNames:
                field.function.x.array[:]           = _assemble_complex_arrays(
                    varName,
                    arrays,
                    typeVars,
                )
            else:
                logger.warning(f"Variable '{varName}' not found in loaded arrays for scalar field. Set to default values.")
        return field
    
    
    def _check_variable_availability_and_type(
        self,
        field: Any
    ) -> Tuple[List[str], List[str]]:
        """
        Check which variables needed for the field are available in the data source.
        And check which ones are real or complex.

        Parameters
        ----------
        field : object
            The FELiCS field object.

        Returns
        -------
        tuple
            (list of present variable names, list of missing variable names)
        """
        
        # Get all the variables available in the file
        self._get_list_available_vars()
        
        # Determine variable types (real/complex)
        typeVars                = []
        for varName in self._availableVars:
            if varName.endswith("_real") or (varName[:-5] + "_imag") in self._availableVars:
                typeVars.append((varName[:-5], "complex"))
            else:
                typeVars.append((varName, "real"))
                
        # Remove duplicated entries for complex variables
        def unique_list(seq):
            seen = set()
            seen_add = seen.add
            return [x for x in seq if not (x in seen or seen_add(x))]
        typeVars                = unique_list(typeVars)
        
        # Get wanted variable names for this field
        wantedVars              = self._names_for_field(field)
        
        # Determine which variables are present or missing (not considering _real/_imag suffixes)
        missingVars             = [var for var in wantedVars if var not in [x[0] for x in typeVars]]
        presentVars             = [var for var in wantedVars if var in [x[0] for x in typeVars]]
        
        if not presentVars:
            logger.warning(f"No variables for field '{field.name}' found in file. Set to default values.")
            return presentVars, typeVars, missingVars
        
        return presentVars, typeVars, missingVars

    # --------------------------
    # Main API
    # --------------------------
    def import_in_field(
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

        # Add ".h5" to the file path it not specified
        fullFilePath            = Path(filePath)
        if fullFilePath.suffix != ".h5" and fullFilePath.suffix != ".fel":
            fullFilePath        = str(fullFilePath.with_name(fullFilePath.stem + ".h5"))

        # Update the calculation mesh and purge cached properties if source changed 
        self._update_calc_mesh(field)
        self._update_data_source(
            fullFilePath,
            groupName,
        )
        
        # List of FEM spaces degrees in the field
        subFields               = field.get_list_of_sub_fields() if field.info["type"] == "mixed" else [field]
        degrees                 = [
            subField.info.get(
                "degree",
                None,
            ) for subField in subFields
        ]
        
        # Check variable availability
        presentVars, typeVars, missingVars = self._check_variable_availability_and_type(
            field
        )

        # Load and process data
        if self._cacheData:
            
            # Exclude the coordinates from caching variables
            coordsSet           = set(self.full_axis_names)
            allVars             = [var for var in self._availableVars if var not in coordsSet]
            
            if not self._dataForField:
                logger.debug(f"Cache empty, loading and interpolating/mapping all available data: {self._availableVars}")
                
                # Load all non-coordinate variables for this file
                self._load_from_h5(allVars)

                # If a mixed space has a P1 variable, we need to interpolate the variables to P1
                # NOTE: This is not optimal if all subfields are P1, but this is a rare case?
                if field.info["type"] == "mixed" and 1 in degrees and not self._needInterpolation:
                    logger.debug("Mixed field with P1 subfield detected.")

                    # Dict of only P1 variables
                    variablesP1     = [name for name, deg in zip(
                    field.get_names_of_sub_fields(),
                    degrees,
                    ) if deg == 1]

                    # Add "_real" and "_imag" suffixes for complex variables
                    variablesP1     = [var for var in allVars if (var in variablesP1) or (var.endswith("_real") and var[:-5] in variablesP1) or (var.endswith("_imag") and var[:-5] in variablesP1)]
                    logger.debug(f"Interpolating P1 variable(s): {variablesP1}")
                    processedP1Data = self._interpolate_to_calc_mesh(
                        variablesP1,
                        toP1=True,
                    )
                    
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
                
            logger.debug(f"Getting {presentVars} from cached data.")
            arrays                  = {var: self._dataForField[var] for var in allVars}

        else:
            # TODO: test this workflow
            # Nothing in cache yet
            self._load_from_h5(presentVars)
            arrays                  = (self._interpolate_to_calc_mesh(presentVars)
                if self._needInterpolation else
                self._map_to_calc_mesh(presentVars)
            )

        self._set_arrays_to_field(
            arrays,
            field,
            typeVars,
        )
        return field, missingVars
