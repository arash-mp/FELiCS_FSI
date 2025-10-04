import  os
import  h5py
import  numpy               as np
# from    scipy               import interpolate
from    scipy.spatial       import Delaunay
from    scipy.interpolate   import LinearNDInterpolator, NearestNDInterpolator
from    FELiCS.Misc.logging import Logger
from    functools           import cached_property

# Debug
# import  matplotlib
# import  matplotlib.pyplot   as plt
# matplotlib.use("TkAgg")

# Get the logger
logger = Logger.get_logger("felics")

class Reader:
    
    def __init__(
            self,
            needInterpolation   = True,
            originalMeshFile    = None,
            isComplex           = False, # NOTE: Should be an attribute of the field, not the reader?
            cacheData           = True,
            ):
        """
        Initialize the Reader for loading and interpolating data onto FELiCS mesh.
        
        This class handles reading data from HDF5 files and mapping/interpolating
        them onto the FELiCS mesh. Each Reader instance corresponds to one file.
        
        Parameters
        ----------
        needInterpolation : bool, optional
            Whether to interpolate data onto the FELiCS mesh, by default True.
        isComplex : bool, optional
            Whether the data is complex-valued, by default False.
        originalMeshFile : str or None, optional
            Path to original mesh file if coordinates not in main file, by default None.
        cacheData : bool, optional
            Whether to load all available data in file to a cache, by default True.
            
        TODOs:
            - Handle mixed function spaces when setting arrays.

        Raises
        ------
        FileNotFoundError
            If the specified file does not exist.
        """
        
        # config (fixed for this Reader instance)
        self._needInterpolation = needInterpolation
        self._originalMeshFile  = originalMeshFile
        self._isComplex         = isComplex
        self._cacheData         = cacheData
        
        # environment TODO: get rid of these dependences?
        self._param             = None
        self._FEMSpaces         = None
        
        # per-file caches (reset when source changes)
        self._source_key        = None          # (filePath, groupName)
        self._availableVars     = None
        self._rawDataDict       = {}
        self._dataForField      = {}
        
    # --------------------------
    # Environment binding
    # --------------------------
    # NOTE: TEMPORARY solution until updates in Field are made
    # TODO: get rid of dependences to these objects?
    def bind_env(self, param, FEMSpaces):
        self._param     = param
        self._FEMSpaces = FEMSpaces

    # --------------------------
    # Cached mesh and mapping data
    # --------------------------
    @cached_property
    def ndim(self):
        if self._param is None:
            raise RuntimeError("ndim accessed before bind_env()")
        return self._param.Case.nDim

    @cached_property
    def coordnames(self):
        if self._param is None:
            raise RuntimeError("coordnames accessed before bind_env()")
        return self._param.getMesh().coordinateNames
    
    @cached_property
    def _triangulationImportMesh(self):
        """Cached Delaunay triangulation of the import mesh for interpolation."""
        import_coords = self.importMeshCoords[:, :self.ndim]
        return Delaunay(import_coords)
    
    @cached_property
    def calcMeshCoords(self):
        """
        Get the coordinates of the calculation mesh (P2).
        # TODO: might want to add P1 as well
        Returns
        -------
        numpy.ndarray
            Array of shape (nPoints, nDim) with the coordinates of the calculation mesh.
        """
        logger.debug("Getting the calculation mesh coordinates.")
        return self._FEMSpaces.P2.tabulate_dof_coordinates()[:, : self.ndim]
        
    @cached_property
    def importMeshCoords(self):
        """
        Get the coordinates of the import mesh based on the configuration.
        Three scenarios:
            1. Interpolation needed, original mesh file provided: load from original mesh file.
            2. Interpolation needed, no original mesh file: load from main file under group name.
            3. No interpolation needed: still load from FELiCS exported mesh file.
        Returns
        -------
        numpy.ndarray
            Array of shape (nPoints, nDim) with the coordinates of the import mesh.
        """
        logger.debug("Getting the import mesh coordinates.")
        if self._needInterpolation:
            if self._originalMeshFile:
                # Option 1: load from original mesh file
                with h5py.File(self._originalMeshFile, "r") as h:
                    pref            = "coordinates/"   # NOTE: we want to get rid of groups in files
                    n               = len(self.coordnames)
                    size0           = h[pref + self.coordnames[0]][:].shape
                    arr             = np.zeros((max(size0), n))
                    for i, c in enumerate(self.coordnames):
                        arr[:, i]   = h[pref + c][:]
                return arr
            else:
                # Option 2: load from main file under group name
                # NOTE: we assume the coordinates are in the same group as the variables
                # Build from the first current file that has coords
                self._ensure_source_set()
                self._get_list_available_vars()
                missing             = [c for c in self.coordnames if c not in self._rawDataDict.keys()]
                if missing: 
                    self._load_from_h5(missing)
                size0               = self._rawDataDict[self.coordnames[0]].shape
                arr                 = np.zeros((max(size0), len(self.coordnames)))
                for i, c in enumerate(self.coordnames):
                    arr[:, i]       = self._rawDataDict[c]
                return arr
        else:
            # Option 3: load from FELiCS exported mesh file
            mesh_file               = os.path.join(
                self._param.Export.ExportFolder, 
                f"{self._param.Case.AnalysisMode}_mesh.h5" # TODO: call all FELiCS meshes the same way
            )
            with h5py.File(mesh_file, "r") as h:
                pref                = "coordinates/"
                size0               = h[pref + "x"][:].shape
                arr                 = np.zeros((max(size0), self.ndim))
                # NOTE: in FELiCS mesh the coordinates are always x,y,z! Change?
                for i, c in enumerate(["x","y","z"][:self.ndim]):
                    arr[:, i]       = h[pref + c][:]
            return arr
            
    @cached_property
    def import_to_P2calc(self):
        """
        Get the mapping from import mesh to calculation mesh (P2).
        TODO: might want to add mapping to P1 as well
        TODO: use the mapping from Field once it's implemented
        Returns
        -------
        numpy.ndarray
            Array of shape (nPoints, nDim) with the coordinates of the import mesh.
        """
        if self._needInterpolation: 
            return None
        return self._FEMSpaces.mappingObj._mappingFunc(
            self.importMeshCoords[:, :self.ndim],
            self.calcMeshCoords
        )
            
    # --------------------------
    # Source management
    # --------------------------
    def _update_data_source(self, filePath, groupName):
        """Switch to new file/group if different, clearing per-file caches."""
        
        if not os.path.isfile(filePath):
            raise FileNotFoundError(filePath)
        
        # Set new source if different
        key                     = (os.path.abspath(filePath), groupName)
        # NOTE: might need completion of the variable list
        if key != self._source_key:
            self._source_key    = key
            self._availableVars = None
            self._rawDataDict   = {}
            self._dataForField  = {}
            logger.debug(f"Reader: switched source to {key}, cleared per-file caches.")
        
    def _ensure_source_set(self):
        if self._source_key is None:
            raise RuntimeError("No data source set. Call importInField(..., filePath, groupName) first.")
        
    def clear_cached_properties(self):
        """
        Remove all cached @cached_property values from this Reader instance.
        They will be recomputed on the next access.
        """
        cached_props = [
            "ndim",
            "coordnames",
            "calcMeshCoords",
            "importMeshCoords",
            "import_to_P2calc",
        ]
        cleared = []
        for name in cached_props:
            if name in self.__dict__:
                self.__dict__.pop(name)
                cleared.append(name)
        if cleared:
            logger.debug(f"Cleared cached properties: {cleared}")
        else:
            logger.debug("No cached properties to clear.")

    # --------------------------
    # HDF5 file handling
    # --------------------------
    def _get_list_available_vars(self):
        """List variable names available in the current group."""
        if self._availableVars is not None:
            return
        filePath, groupName         = self._source_key
        with h5py.File(filePath, "r") as f:
            if groupName in f:
                self._availableVars = list(f[groupName].keys())
            else:
                logger.warning(f"Group '{groupName}' not found in file '{filePath}'.")
                self._availableVars = []
                
    def _load_from_h5(self, var_names):
        """Load selected variables from HDF5 file into _rawDataDict."""
        # NOTE: maybe we just check if the groups exist instead on relying on _isComplex?
        filePath, groupName = self._source_key
        with h5py.File(filePath, "r") as f:
            for v in var_names:
                if v in self._rawDataDict:
                    continue
                base = f"{groupName}/{v}"
                if self._isComplex:  # NOTE: not tested yet
                    mag = f[base + "/magnitude"][:].squeeze()  # TODO: make _real and _imag in file
                    ang = f[base + "/angle"][:].squeeze()
                    self._rawDataDict[v + "_real"] = mag * np.cos(ang)
                    self._rawDataDict[v + "_imag"] = mag * np.sin(ang)
                else:
                    if not self._needInterpolation and (base + "/magnitude") in f:
                        # If we don't interpolate, we load a FELiCS file with magnitude only
                        self._rawDataDict[v] = f[base + "/magnitude"][:].squeeze()
                    else:
                        # Otherwise we load the full variable
                        self._rawDataDict[v] = f[base][:].squeeze()

    # --------------------------
    # Interpolation / mapping
    # --------------------------
    def _interpolate_to_calc_mesh(self, var_names):
        """Interpolate given variables from import mesh to calculation mesh."""
        import_coords   = self.importMeshCoords[:, :self.ndim]
        calc_coords     = self.calcMeshCoords
        n_import        = import_coords.shape[0]
        
        # Sanity check: all variables have same length as import mesh
        for v in var_names:
            if self._rawDataDict[v].shape[0] != n_import:
                raise ValueError(
                    f"Variable '{v}' has length {len(self._rawDataDict[v])}, "
                    f"but import mesh has {n_import} points."
                )
    
        # Stack data columns (N, K)
        values          = np.column_stack([self._rawDataDict[v] for v in var_names])

        # Use cached triangulation
        linear_interp   = LinearNDInterpolator(self._triangulationImportMesh, values)
        interpolated    = linear_interp(calc_coords)   # shape (M, K)
        
        # Handle NaNs by nearest-neighbour fill (vectorized)
        if np.isnan(interpolated).any():
            nearest_interp          = NearestNDInterpolator(import_coords, values)
            nearest_values          = nearest_interp(calc_coords)
            nan_mask                = np.isnan(interpolated)
            interpolated[nan_mask]  = nearest_values[nan_mask]

        # Return a dict {var_name: array}
        return {v: interpolated[:, i] for i, v in enumerate(var_names)}

    def _map_to_calc_mesh(self, var_names):
        """Map given variables using precomputed index mapping."""
        idx     = self.import_to_P2calc
        return {v: self._rawDataDict[v][idx] for v in var_names}
    
    # --------------------------
    # Field helpers
    # --------------------------
    def _names_for_field(self, Field):
        """Determine variable names to load from file for this Field."""
        names   = []
        info    = Field.describeFunctionSpace()
        subs    = Field.getListOfSubFields() if info["type"] == "mixed" else [Field]
        for sub in subs:
            si  = sub.describeFunctionSpace()
            if si["type"] == "vector":
                for c in sub.getComponentsNames():
                    # NOTE: this is an annoying work around for u_forcing
                    # TODO: make the same pattern for all vector variables!
                    if Field.getName()[0] == "u":
                        names.append(sub.getName()[0] + c + sub.getName()[1:])
                    else:
                        names.append(sub.getName() + c)
            else:
                # Then it's a scalar
                names.append(sub.getName())
        return names

    def _set_arrays_to_field(self, arrays, Field):
        """Assign interpolated/mapped arrays to a given Field."""
        info = Field.describeFunctionSpace()
        if info["type"] == "mixed":
            # TODO: implement :P
            raise NotImplementedError("Setting arrays to mixed fields is not implemented yet.")
        if info["type"] == "vector":
            for i, c in enumerate(Field.getComponentsNames()):
                # NOTE: this is an annoying work around for u_forcing
                # TODO: make the same pattern for all vector variables!
                full = (Field.getName()[0] + c + Field.getName()[1:]) if Field.getName()[0] == "u" \
                    else (Field.getName() + c)
                if full in arrays:
                    idx = Field.space.sub(i).collapse()[1]
                    Field.function.x.array[idx] = arrays[full]
        else:
            # Then it's a scalar
            v = Field.getName()
            if v in arrays:
                Field.function.x.array[:] = arrays[v]
        return Field
   
   
    # --------------------------
    # Main API
    # --------------------------
    def importInField(self, Field, filePath, groupName):
        """
        Import data from the specified file/group into the given Field.
        Automatically clears per-file caches when file/group changes.
        """
        if self._param is None or self._FEMSpaces is None:
            raise RuntimeError("Call bind_env(param, FEMSpaces) before importInField().")

        self._update_data_source(filePath, groupName)
        self._get_list_available_vars()

        wanted  = self._names_for_field(Field)
        missing = [n for n in wanted if n not in self._availableVars]
        present = [n for n in wanted if n in self._availableVars]
        if not present:
            logger.warning(f"No variables for Field '{Field.getName()}' \
                found in file '{filePath}'. Set to default values.")
            return Field, missing

        # Load and process data
        if self._cacheData:
            if not self._dataForField:
                logger.debug("Cache empty, loading and interpolating/mapping all available data.")
                # Load and interpolate/map all non-coordinate variables for this file
                coords      = set(self.coordnames)
                all_vars    = [v for v in self._availableVars if v not in coords]
                self._load_from_h5(all_vars)
                ready       = (self._interpolate_to_calc_mesh(all_vars)
                                if self._needInterpolation else
                                self._map_to_calc_mesh(all_vars)
                            )
                # Update cached data
                self._dataForField.update(ready)

            extra           = [v for v in present if v not in self._dataForField]
            if extra:
                self._load_from_h5(extra)
                upd         = (self._interpolate_to_calc_mesh(extra)
                                if self._needInterpolation else
                                self._map_to_calc_mesh(extra)
                            )
                self._dataForField.update(upd)
            logger.debug("Reader: using cached data for field assignment.")
            arrays          = {v: self._dataForField[v] for v in present}
            
        else:
            # Nothing in cache yet
            self._load_from_h5(present)
            arrays          = (self._interpolate_to_calc_mesh(present)
                                if self._needInterpolation else
                                self._map_to_calc_mesh(present)
                            )

        self._set_arrays_to_field(arrays, Field)
        return Field, missing