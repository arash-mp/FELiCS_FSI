
import numpy as np
import pdb

import time
import sys

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
        print('Calculating mapping Vector from P2 Calculation Space to P1 Export Space...')
        self.P2CalcToP1ExportIndecies = self._mappingFunc(self.P2CalcDofCoordinates, self.P1ExportDofCoordinates)
        print('Calculating mapping Vector from VMixed Vector-P2-Sub Space to P1 Export Space...')
        self.VectorCalcToP1ExportIndecies = self._mappingFunc(self.VMixedVectorDofCoords, self.P1ExportDofCoordinates)

        # print('Calculating mapping Vector from Vector-Calc-P2 Space to Vector-P1-Export Space...')
        # self.VectorCalcToP1VectorExportIndecies = self._mappingFunc(self.VMixedVectorDofCoords, self.VectorP1ExportCalcDofCoordinates)

        # print('Calculating mapping Vector from Vector-Export-P1 Space to P1-Export Space...')
        # self.VectorExportToP1ExportIndices = self._mappingFunc(self.P1ExportDofCoordinates, self.VectorP1ExportCalcDofCoordinates)


    def _mappingFunc(self, exportMeshDOFCoordinates, calcMeshDOFCoordinates):
        s_t = time.time()
        # mapping from CalcMesh to exportMesh needs to be done:
        exportMesh = np.copy(exportMeshDOFCoordinates)
        exportMeshNorm = np.linalg.norm(exportMesh,axis=1)
        calcMesh = np.copy(calcMeshDOFCoordinates)
        calcMeshNorm = np.linalg.norm(calcMesh,axis=1)
        exportMesh = np.append(exportMesh,np.arange(len(exportMesh)).reshape(len(exportMesh),1),axis=1)
        calcMesh = np.append(calcMesh,np.arange(len(calcMesh)).reshape(len(calcMesh),1),axis=1)
        #exportMesh = exportMesh[np.lexsort((exportMesh[:,2],exportMesh[:,1],exportMesh[:,0]))].astype(np.int)
        calcMesh = calcMesh[np.lexsort((calcMesh[:, 2], calcMesh[:, 1], calcMesh[:, 0]))].astype(np.int)
        mappingtest = exportMesh[:,-1][calcMesh[:,-1]]

        #print(time.time() - s_t, file=sys.stderr)

        mapping = np.zeros(exportMesh.shape[0], dtype=int)
        del_indices = []
        for index, coordinate in enumerate(calcMeshDOFCoordinates):
            ind= np.isclose(coordinate, exportMesh[:,:-1])
            cur_index = np.where(ind.all(axis=1))[0]
            mapping[index] = exportMesh[cur_index,-1]
            del_indices.append(cur_index)
            if len(del_indices) > 100:
                exportMesh = np.delete(exportMesh, del_indices,axis=0)
                del_indices = []
        #print(mapping, file=sys.stderr)
        #print(mappingtest, file=sys.stderr)
        #print(time.time() - s_t, file=sys.stderr)

        return mapping

    # @property
    # def MappingVector(self):
    #     return self._indexVector
