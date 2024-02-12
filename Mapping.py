
import numpy as np
import pdb

from functions import (
						printDebug,
						printError,
					)

class Mapping:
    """
    This class contains mappings

    Attributes:
    -
    """
    def __init__(self, FEMSpaces):

        self.P2CalcDofCoordinates = FEMSpaces.P2.tabulate_dof_coordinates()
        self.P1ExportDofCoordinates = FEMSpaces.P1Export.tabulate_dof_coordinates()
        #pdb.set_trace()
        self.VectorP2CalcDofCoordinates = FEMSpaces.FunctionSpaceVectorVelocity.tabulate_dof_coordinates()
        self.VectorP1ExportCalcDofCoordinates = FEMSpaces.FunctionSpaceVectorVelocityExport.tabulate_dof_coordinates()

        self.VMixedVectorDofCoords = FEMSpaces.VMixed.sub(0).collapse()[0].tabulate_dof_coordinates()
        printDebug(True, '-- Calculating mapping Vector from P2 Calculation Space to P1 Export Space...')
        self.P2CalcToP1ExportIndecies = self._mappingFunc(self.P2CalcDofCoordinates, self.P1ExportDofCoordinates)
        printDebug(True, '-- Calculating mapping Vector from VMixed Vector-P2-Sub Space to P1 Export Space...')
        self.VectorCalcToP1ExportIndecies = self._mappingFunc(self.VMixedVectorDofCoords, self.P1ExportDofCoordinates)

        # print('Calculating mapping Vector from Vector-Calc-P2 Space to Vector-P1-Export Space...')
        # self.VectorCalcToP1VectorExportIndecies = self._mappingFunc(self.VMixedVectorDofCoords, self.VectorP1ExportCalcDofCoordinates)

        # print('Calculating mapping Vector from Vector-Export-P1 Space to P1-Export Space...')
        # self.VectorExportToP1ExportIndices = self._mappingFunc(self.P1ExportDofCoordinates, self.VectorP1ExportCalcDofCoordinates)


    def _mappingFunc(self, exportMeshDOFCoordinates, calcMeshDOFCoordinates):
        """
        This function provides a mapping from the calculation mesh to the export
        mesh. It adds a column with indices to both coordinate matrices each,
        sorts the matrices by x, y and z coordinate and then takes the last
        column of indices to build a new matrix. Sorting this new matrix by the
        column with indices of the calculation mesh will reveal the wanted
        mapping to the export mesh in the other column.
        Note that the coordinates are rounded to 11 decimals for the sorting
        process. This does not affect the actual coordinates and is only
        relevant within this function.
        """

        # copy coordinates to avoid changing the original arrays
        exportMesh = np.copy(exportMeshDOFCoordinates)
        calcMesh = np.copy(calcMeshDOFCoordinates)

        # append indices as last column
        exportMesh = np.append(exportMesh,np.arange(len(exportMesh)).reshape(len(exportMesh),1),axis=1).round(11)
        calcMesh = np.append(calcMesh,np.arange(len(calcMesh)).reshape(len(calcMesh),1),axis=1).round(11)

        # sort by x,y,z
        exportMeshSorted = exportMesh[np.lexsort((exportMesh[:,2],exportMesh[:,1],exportMesh[:,0]))].astype(np.int)
        calcMeshSorted = calcMesh[np.lexsort((calcMesh[:, 2], calcMesh[:, 1], calcMesh[:, 0]))].astype(np.int)

        # find indices of sorted exportMesh in calcMesh
        index_array = np.vstack((calcMeshSorted[:, -1], exportMeshSorted[:, -1])).T
        mapping = index_array[index_array[:, 0].argsort()][:, 1]

        return mapping

    # @property
    # def MappingVector(self):
    #     return self._indexVector
