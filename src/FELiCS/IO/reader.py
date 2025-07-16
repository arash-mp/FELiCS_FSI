import  os
import  h5py
import  numpy               as np
from    scipy               import interpolate
from    FELiCS.Misc.logging import Logger
from    FELiCS.IO.Mapping   import Mapping
from    FELiCS.Fields.Field import Field


# Get the logger
logger = Logger.get_logger("felics")

class Reader:
    def __init__(self, param, FEMSpaces):
        """

        Function arguments:
        - FEMSpaces: Object containing the Cacluation FEM-Spaces, the export
            FEM-Spaces and the corresponding meshes.
        - param: Object, containing the parameters of the calculation.

        Function returns:
        """
        self._FEMSpaces     = FEMSpaces
        self._param         = param

    def import_meanflow_from_file(self, felics_mesh, nameListMean):
        """
        This function does the loading of the mean flow data from the hdf5 file
        corresponding to fields already interpolated on the FELiCS mesh.
        TODO: cleanup
        """
        meanflow        = h5py.File(self._param.FlowInput.MeanFlowFilePath, 'r')
        exportMesh      = h5py.File(f"{self._param.Case.AnalysisMode}_mesh.h5", 'r')
        coordNameList   = ['x', 'y', 'z']

        # arrange the export mesh coordinates in an array
        coordArray = np.zeros((exportMesh['coordinates/x'][:].shape[0], self._param.Case.nDim))
        coordArray[:, 0] = exportMesh['coordinates/x'][:]
        for i in range(self._param.Case.nDim - 1):
            coordArray[:, i + 1] = exportMesh[f'coordinates/{coordNameList[i + 1]}'][:]

        # map the coordinates of the P2-mesh to the export mesh coordinates
        coordinatesOfP2Mesh = self._FEMSpaces.P2.tabulate_dof_coordinates()[:,0:self._param.Case.nDim]
        mapping = Mapping(self._FEMSpaces)
        indexMappingArray = mapping._mappingFunc(coordArray, coordinatesOfP2Mesh)

        # get the velocity components
        velocityComponents = self._param.getInternalVelocityComponents()
        for index, comp in enumerate(velocityComponents):
            velocityComponents[index] = f'u{comp}'

        # loop over nameListMean and create the corresponding fields
        fieldDict = {}
        notInFileList = []
        nameListMean = nameListMean
        for name in nameListMean:
            if name[0] == 'u' and not (name == 'ut' or name == 'ut_forcing'):
                fieldDict[name] = Field(self._FEMSpaces.FunctionSpaceVectorVelocity, felics_mesh, name = name)
            else:
                fieldDict[name] = Field(self._FEMSpaces.P2, felics_mesh, name=name)

            if name == 'u':
                for index, component in enumerate(velocityComponents):
                    indicesOfSubField = fieldDict[name].space.sub(index).collapse()[1]
                    fieldDict[name].function.x.array[indicesOfSubField] = \
                        meanflow[f'meanflow/{component}/magnitude'][:][indexMappingArray]
            else:
                if name in list(meanflow[f'meanflow'].keys()):
                    fieldDict[name].function.x.array[:] = \
                        meanflow[f'meanflow/{component}/magnitude'][:][indexMappingArray]
                else:
                    notInFileList.append(name)
        
        # set the constant value 0 for the ut field if it exists
        if 'ut' in list(fieldDict.keys()):
            fieldDict['ut'].setConstant(0.)

        # set the resulting fieldDict
        return fieldDict, notInFileList

    def import_interpolate_meanflow_from_file(self, felics_mesh, nameListMean, ScalarFunctionSpace, VectorFunctionSpace):
        """
        This method reads mean flow field data from a specified file, interpolates it onto the 
        FELiCS mesh, and stores the results in the object's field dictionary.
    
        Notes
        -----
        Vector fields like velocity ('u') and momentum ('rhou') are handled specially,
        with components being assembled into vector function spaces.
        TODO: Handle vector fields in a more automatic way based on "nameListMean".
        Raises
        ------
        ValueError
            If the coordinate system is not recognized.
        """
        
        param               = self._param
        filePath            = param.FlowInput.MeanFlowFilePath
        velocityComponents  = param.BoundaryCondition.VelocityComponents
        coordSys            = param.Case.CoordinateSystem
        
        # Get the list of mean flow fields to be read
        nameListMeanOriginal = nameListMean.copy()
        
        # TODO: Include in the nameListMean if variables are vectors and deal with components inside reader!
        # Check if we ask for "u" and add corresponding veloctiy components
        # (we always load 3D velocity components, even if the case is 2D)
        if 'u' in nameListMean:
            nameListMean.remove('u')
            for component in reversed(velocityComponents):
                nameListMean.insert(0, 'u' + component)
        if 'rhou' in nameListMean:
            nameListMean.remove('rhou')
            for component in reversed(velocityComponents):
                nameListMean.insert(0, 'rhou' + component)
        if 'u_forcing_i' in nameListMean:
            nameListMean.remove('u_forcing_i')
            for component in reversed(velocityComponents):
                nameListMean.insert(0, 'u' + component +'_forcing_i')
        if 'u_forcing_r' in nameListMean:
            nameListMean.remove('u_forcing_r')
            for component in reversed(velocityComponents):
                nameListMean.insert(0, 'u' + component +'_forcing_r')
        
        # Define coordinate names
        if coordSys == 'Cartesian':
            list_of_coords = ['MeanFlow/x', 'MeanFlow/y']
            if param.Case.nDim == 3:
                list_of_coords.append('MeanFlow/z')
        elif coordSys == 'Cylindrical':
            list_of_coords = ['MeanFlow/x', 'MeanFlow/r']
            if param.Case.nDim == 3:
                list_of_coords.append('MeanFlow/t')
        else:
            logger.error(f"Names of coordinate system '{coordSys}' not recognized.")
            raise ValueError(f"Names of coordinate system '{coordSys}' not recognized.")
            
        # Add "MeanFlow/" prefix to the names
        nameListMean = [f'MeanFlow/{name}' for name in nameListMean]
        
        self._interpFlowDict, notInFileList = self.read_and_interpolate_on_felics_mesh(
            filePath,
            nameListMean,
            list_of_coords
        )
        list_of_available_vars = list(self._interpFlowDict.keys())
        
        fieldDict = {}
        # Define fields in the fieldDict
        for name in nameListMeanOriginal:
            # TODO: deal with vectors in an automatic way! using the namelist
            if name == "u" or name == "rhou":
                # Create a vector field for the velocity
                fieldDict[name] = Field(
                    VectorFunctionSpace, 
                    felics_mesh, 
                    name=name
                )
                # Loop over components to set function coeff values
                for icomp, comp in enumerate(velocityComponents):
                    if name+comp in list_of_available_vars:
                        dofIDX = VectorFunctionSpace.sub(icomp).collapse()[1]
                        # Set the function values
                        fieldDict[name].function.sub(icomp).x.array[dofIDX] = \
                            np.array(self._interpFlowDict[name+comp])
            elif name == "u_forcing_r" or name == "u_forcing_i":
                # Create a vector field for the velocity
                fieldDict[name] = Field(
                    VectorFunctionSpace, 
                    felics_mesh, 
                    name=name
                )
                # Loop over components to set function coeff values
                for icomp, comp in enumerate(velocityComponents):
                    if name[0]+comp+name[1:] in list_of_available_vars:
                        dofIDX = VectorFunctionSpace.sub(icomp).collapse()[1]
                        fieldDict[name].function.sub(icomp).x.array[dofIDX] = \
                            np.array(self._interpFlowDict[name[0]+comp+name[1:]])
       
            else:
                if name in list_of_available_vars:
                    # Create a scalar field for the other quantities
                    fieldDict[name] = Field(ScalarFunctionSpace, felics_mesh, name=name)
                    fieldDict[name].function.x.array[:] = np.array(self._interpFlowDict[name])
        
        del self._interpFlowDict
        return fieldDict, notInFileList


    def read_and_interpolate_on_felics_mesh(
        self, 
        filename, 
        list_of_variables, 
        list_of_coords, 
        destinationSpace=None
    ):
        """
        Interpolate a field on the FELiCS mesh.

        This method is not yet implemented.
        """
        
        list_of_missing_vars                = []
        
        # Default space is P2 scalar space
        if destinationSpace is None:
            destinationSpace                = self._FEMSpaces.P2
        
        # Open file and check for variables
        with h5py.File(filename, 'r') as h5file:
            # Check which variables are present
            for var in list_of_variables:
                if var not in h5file:
                    list_of_missing_vars.append(var)
                    logger.warning(f"Variable '{var}' not found in file '{filename}', set to default (see fieldProperties.py).")
            list_of_available_vars          = [var for var in list_of_variables if var not in list_of_missing_vars]
                    
            # Get the number of input data points
            coord_shape     = h5file[list_of_coords[0]][:].shape
            N_input_data    = max(coord_shape) if isinstance(coord_shape, tuple) else coord_shape
                    
            # Read the input mesh
            input_mesh                      = np.zeros((N_input_data, len(list_of_coords)))
            for i, coord in enumerate(list_of_coords):
                input_mesh[:, i]            = h5file[coord][:]
                
            # Read the input data
            input_data                      = np.zeros((N_input_data, len(list_of_available_vars)))
            for i, var in enumerate(list_of_available_vars):
                input_data[:, i]            = h5file[var][:]
        
        # Interpolation on FELiCS mesh
        interpolated_data = self.interpolate_on_felics_mesh(
            input_mesh,
            input_data,
            destinationSpace
        )
        
        # List of variables without prefix
        list_of_available_vars              = [var.split('/')[-1] for var in list_of_available_vars]
        list_of_missing_vars                = [var.split('/')[-1] for var in list_of_missing_vars]
        
        # Populating dictionary with the interpolated fields
        dict_interpolated_fields            = {}
        for i, var in enumerate(list_of_available_vars):
            dict_interpolated_fields[var]   = interpolated_data[:, i]
        
        return dict_interpolated_fields, list_of_missing_vars
        
        
    def interpolate_on_felics_mesh(self, input_mesh, input_data, destinationSpace):
        """
        Interpolate data from input mesh to destination space.
        NOTE: At this stage we only interpolate on "single" spaces, 
        i.e. not on vector spaces or mixed spaces.
        """
        
        # Get the Dofs corresponding to destination space
        nDim                    = self._param.Case.nDim    
        destinationMesh         = destinationSpace.tabulate_dof_coordinates()[:, np.arange(nDim)]
        
        # First attempt to interpolate with nearest neighbor
        interp_data_nearest     = interpolate.griddata(
            input_mesh, 
            input_data, 
            destinationMesh, 
            method='nearest'
        )
        
        # Attempt linear interpolation, fallback to nearest neighbor if it fails
        try:
            interp_data             = interpolate.griddata(
                input_mesh, 
                input_data, 
                destinationMesh, 
                method='linear'
            )
            
            # Replace NaNs with nearest neighbor values
            nan_mask                = np.isnan(interp_data)
            interp_data[nan_mask]   = interp_data_nearest[nan_mask]
            
        except Exception as e:
            logger.warning(f"Linear interpolation failed: {str(e)}. Using nearest neighbor interpolation instead.")
            interp_data             = interp_data_nearest

        return interp_data
        