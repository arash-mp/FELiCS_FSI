
import numpy as np
import pdb

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
    	# mapping from CalcMesh to exportMesh needs to be done:
    	mapping = np.zeros(exportMeshDOFCoordinates.shape[0], dtype=int)
    	for index, coordinate in enumerate(calcMeshDOFCoordinates):
    		mapping[index] = np.where(np.isclose(coordinate, exportMeshDOFCoordinates).all(axis=1) == True)[0]

    	return mapping

    # @property
    # def MappingVector(self):
    #     return self._indexVector
