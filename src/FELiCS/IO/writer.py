import  os
import  h5py
import  numpy               as np
from    FELiCS.Fields.Field import Field
from    FELiCS.Misc.logging import Logger

# Get the logger
logger = Logger.get_logger("felics")

class Writer:
    
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
        self._exportMesh    = FEMSpaces.exportMesh

        # self._exportZeroScalarField = Field(self._FEMSpaces.P1Export, self._exportMesh, name="exportZeroScalar")
        # self._exportZeroVectorField = Field(
        #         self._FEMSpaces.FunctionSpaceVectorVelocityExport,
        #         self._exportMesh, 
        #         name="exportZeroVector"
        #         )
        
        