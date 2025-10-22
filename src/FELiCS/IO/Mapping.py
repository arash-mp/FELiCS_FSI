import  numpy as np
from 	FELiCS.Misc.logging         import Logger

# Get the logger
logger = Logger.get_logger("felics")

class Mapping:
    #TODO: Docstrings

    @staticmethod
    def calculateMappingFromSpaces(
            inputSpace,
            outputSpace):

        # get numpy arrays of dof coordinates
        inputDofs  = inputSpace.tabulate_dof_coordinates()
        outputDofs = outputSpace.tabulate_dof_coordinates()

        inputMesh  = np.copy(inputDofs)
        outputMesh = np.copy(outputDofs)

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



    @staticmethod
    def calculateMappingFromDofs(
            inputDofs,
            outputDofs):

        # get numpy arrays of dof coordinates
        inputMesh  = np.copy(inputDofs)
        outputMesh = np.copy(outputDofs)

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


