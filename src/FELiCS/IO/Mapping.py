import  numpy as np
from 	FELiCS.Misc.logging         import Logger

# Get the logger
logger = Logger.get_logger("felics")

class Mapping:
    """
    Class for computing degree-of-freedom (DoF) mappings between different FEM spaces.

    This class constructs index mappings between high-order and low-order function
    spaces used in finite element calculations and result exports. It is essential
    for transferring solution data between calculation and export spaces.

    **Initialize the Mapping object**

    Parameters
    ----------
    FEMSpaces : object
        An object that provides access to multiple FEM function spaces with
        methods for tabulating degree-of-freedom coordinates.

    Attributes
    ----------
    P2CalcDofCoordinates : ndarray
        DoF coordinates for the P2 calculation space.
    P1ExportDofCoordinates : ndarray
        DoF coordinates for the P1 export space.
    VectorP2CalcDofCoordinates : ndarray
        DoF coordinates for the vector-valued P2 calculation space.
    VectorP1ExportCalcDofCoordinates : ndarray
        DoF coordinates for the vector-valued P1 export space.
    VMixedVectorDofCoords : ndarray
        DoF coordinates for the vector component of the mixed function space.
    P2CalcToP1ExportIndecies : ndarray
        Mapping indices from P2 calculation space to P1 export space.
    VectorCalcToP1ExportIndecies : ndarray
        Mapping indices from vector-valued calculation space to P1 export space.
    """

    def __init__(self, FEMSpaces):
        """
        Initialize the Mapping class with given FEM function spaces.

        Parameters
        ----------
        FEMSpaces : object
            An object with attributes representing various function spaces,
            each supporting the `tabulate_dof_coordinates()` method.
        """
        self.P2CalcDofCoordinates               = FEMSpaces.P2.tabulate_dof_coordinates()
        self.P1ExportDofCoordinates             = FEMSpaces.P1Export.tabulate_dof_coordinates()
        
        self.VectorP2CalcDofCoordinates         = FEMSpaces.FunctionSpaceVectorVelocity.tabulate_dof_coordinates()
        self.VectorP1ExportCalcDofCoordinates   = FEMSpaces.FunctionSpaceVectorVelocityExport.tabulate_dof_coordinates()

        # NOTE: Names are confusing here, since we only use the sub(0) to get the vector part of the mixed space
        logger.debug('Calculating mappings from Calculation Space to Export Space.')
        self.VMixedVectorDofCoords              = FEMSpaces.VMixed.sub(0).collapse()[0].tabulate_dof_coordinates()
        self.P2CalcToP1ExportIndecies           = self._mappingFunc(self.P2CalcDofCoordinates, self.P1ExportDofCoordinates)
        self.VectorCalcToP1ExportIndecies       = self._mappingFunc(self.VMixedVectorDofCoords, self.P1ExportDofCoordinates)
        # logger.debug('Calculating mapping from Export Space to Calculation Space.')
        # self.P1exportToP2CalcIndecies           = self._mappingFunc(self.P1ExportDofCoordinates, self.P2CalcDofCoordinates)
        # self.P1exportToVectorCalcIndecies       = self._mappingFunc(self.P1ExportDofCoordinates, self.VMixedVectorDofCoords)

        # print('Calculating mapping Vector from Vector-Calc-P2 Space to Vector-P1-Export Space...')
        # self.VectorCalcToP1VectorExportIndecies = self._mappingFunc(self.VMixedVectorDofCoords, self.VectorP1ExportCalcDofCoordinates)

        # print('Calculating mapping Vector from Vector-Export-P1 Space to P1-Export Space...')
        # self.VectorExportToP1ExportIndices = self._mappingFunc(self.P1ExportDofCoordinates, self.VectorP1ExportCalcDofCoordinates)


    def _mappingFunc(self, inputCoords, outputCoords):
        """
        This function provides a mapping between two sets of coordinates.
        It assumes that the two sets contain the same points, possibly in 
        different orders.
    
        Parameters
        ----------
        inputCoords : ndarray
            DoF coordinates of mesh 1.
        outputCoords : ndarray
            DoF coordinates of mesh 2.

        Returns
        -------
        mapping : ndarray
            Index mapping array from calculation mesh to export mesh.

        Notes
        -----
        - This is a fairly slow implementation and could be improved.
        - Coordinates are rounded to 11 decimals for sorting to ensure numerical stability.
        - Sorting is performed lexicographically by x, y, then z.
        - Output is an index array aligning outputMesh to inputMesh order.
        """

        # copy coordinates to avoid changing the original arrays
        inputMesh               = np.copy(inputCoords)
        outputMesh              = np.copy(outputCoords)

        # append indices as last column
        inputMesh               = np.append(inputMesh, np.arange(len(inputMesh)).reshape(len(inputMesh),1), axis=1).round(11)
        outputMesh              = np.append(outputMesh, np.arange(len(outputMesh)).reshape(len(outputMesh),1), axis=1).round(11)

        # sort by x,y,z
        inputMeshSorted         = inputMesh[np.lexsort((inputMesh[:,2], inputMesh[:,1], inputMesh[:,0]))].astype(int)
        outputMeshSorted        = outputMesh[np.lexsort((outputMesh[:, 2], outputMesh[:, 1], outputMesh[:, 0]))].astype(int)

        # find indices of sorted inputMesh in outputMesh
        index_array             = np.vstack((outputMeshSorted[:, -1], inputMeshSorted[:, -1])).T
        mapping                 = index_array[index_array[:, 0].argsort()][:, 1]

        return mapping