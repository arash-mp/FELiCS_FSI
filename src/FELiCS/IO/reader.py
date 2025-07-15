import  os
import  h5py
import  numpy               as np
from    scipy               import interpolate
from    FELiCS.Misc.logging import Logger


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
        