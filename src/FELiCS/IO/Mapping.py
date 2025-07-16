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

        logger.debug('Calculating mapping from P2 Calculation Space to P1 Export Space.')
        self.VMixedVectorDofCoords              = FEMSpaces.VMixed.sub(0).collapse()[0].tabulate_dof_coordinates()
        self.P2CalcToP1ExportIndecies           = self._mappingFunc(self.P2CalcDofCoordinates, self.P1ExportDofCoordinates)
        self.VectorCalcToP1ExportIndecies       = self._mappingFunc(self.VMixedVectorDofCoords, self.P1ExportDofCoordinates)

        # print('Calculating mapping Vector from Vector-Calc-P2 Space to Vector-P1-Export Space...')
        # self.VectorCalcToP1VectorExportIndecies = self._mappingFunc(self.VMixedVectorDofCoords, self.VectorP1ExportCalcDofCoordinates)

        # print('Calculating mapping Vector from Vector-Export-P1 Space to P1-Export Space...')
        # self.VectorExportToP1ExportIndices = self._mappingFunc(self.P1ExportDofCoordinates, self.VectorP1ExportCalcDofCoordinates)


    def _mappingFunc(self, exportMeshDOFCoordinates, calcMeshDOFCoordinates):
        """
        Generate a mapping from the calculation mesh to the export mesh.

        Parameters
        ----------
        exportMeshDOFCoordinates : ndarray
            DoF coordinates for the export mesh.
        calcMeshDOFCoordinates : ndarray
            DoF coordinates for the calculation mesh.

        Returns
        -------
        mapping : ndarray
            Index mapping array from calculation mesh to export mesh.

        Notes
        -----
        - Coordinates are rounded to 11 decimals for sorting to ensure numerical stability.
        - Sorting is performed lexicographically by x, y, then z.
        - Output is an index array aligning calcMesh to exportMesh order.
        """
    #    old docstring, for review purposes: (please delete after reviewing the docstrings)
    #    """
    #    This function provides a mapping from the calculation mesh to the export
    #    mesh. It adds a column with indices to both coordinate matrices each,
    #    sorts the matrices by x, y and z coordinate and then takes the last
    #    column of indices to build a new matrix. Sorting this new matrix by the
    #    column with indices of the calculation mesh will reveal the wanted
    #    mapping to the export mesh in the other column.
    #    Note that the coordinates are rounded to 11 decimals for the sorting
    #    process. This does not affect the actual coordinates and is only
    #    relevant within this function.
    #    """

        # copy coordinates to avoid changing the original arrays
        exportMesh          = np.copy(exportMeshDOFCoordinates)
        calcMesh            = np.copy(calcMeshDOFCoordinates)

        # append indices as last column
        exportMesh          = np.append(exportMesh,np.arange(len(exportMesh)).reshape(len(exportMesh),1),axis=1).round(11)
        calcMesh            = np.append(calcMesh,np.arange(len(calcMesh)).reshape(len(calcMesh),1),axis=1).round(11)

        # sort by x,y,z
        exportMeshSorted    = exportMesh[np.lexsort((exportMesh[:,2],exportMesh[:,1],exportMesh[:,0]))].astype(int)
        calcMeshSorted      = calcMesh[np.lexsort((calcMesh[:, 2], calcMesh[:, 1], calcMesh[:, 0]))].astype(int)

        # find indices of sorted exportMesh in calcMesh
        index_array         = np.vstack((calcMeshSorted[:, -1], exportMeshSorted[:, -1])).T
        mapping             = index_array[index_array[:, 0].argsort()][:, 1]

        return mapping

    # @property
    # def MappingVector(self):
    #     return self._indexVector
