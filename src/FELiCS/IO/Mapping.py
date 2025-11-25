import  numpy as np
from 	FELiCS.Misc.logging         import Logger

# Get the logger
logger = Logger.get_logger("felics")

class Mapping:
    """
    A class for calculating mappings between finite element spaces or DOF coordinates.

    This class provides static methods to compute index mappings between input and output
    spaces or DOF arrays, useful for data transfer in FEM simulations.
    """

    @staticmethod
    def calculateMappingFromSpaces(
            inputSpace,
            outputSpace):
        """
        Calculate the mapping indices from input space to output space.

        Parameters
        ----------
        inputSpace : object
            The input finite element space.
        outputSpace : object
            The output finite element space.

        Returns
        -------
        numpy.ndarray
            Array of indices mapping input DOFs to output DOFs.
        """

        # Get numpy arrays of dof coordinates
        inputDofs           = inputSpace.tabulate_dof_coordinates()
        outputDofs          = outputSpace.tabulate_dof_coordinates()

        return Mapping.calculateMappingFromDofs(inputDofs, outputDofs)

    @staticmethod
    def calculateMappingFromDofs(
            inputDofs,
            outputDofs):
        """
        Calculate the mapping indices from input DOF coordinates to output DOF coordinates.

        Parameters
        ----------
        inputDofs : numpy.ndarray
            Array of input DOF coordinates.
        outputDofs : numpy.ndarray
            Array of output DOF coordinates.

        Returns
        -------
        numpy.ndarray
            Array of indices mapping input DOFs to output DOFs.
        """

        # Get numpy arrays of dof coordinates
        inputMesh           = np.copy(inputDofs)
        outputMesh          = np.copy(outputDofs)

        # Append indices as last column
        inputMesh           = np.append(inputMesh, np.arange(len(inputMesh))[:, None], axis=1).round(11)
        outputMesh          = np.append(outputMesh, np.arange(len(outputMesh))[:, None], axis=1).round(11)

        # Sort by x,y,z
        inputMeshSorted     = inputMesh[np.lexsort((inputMesh[:,2], inputMesh[:,1], inputMesh[:,0]))].astype(int)
        outputMeshSorted    = outputMesh[np.lexsort((outputMesh[:, 2], outputMesh[:, 1], outputMesh[:, 0]))].astype(int)

        # Find indices of sorted inputMesh in outputMesh
        index_array         = np.vstack((outputMeshSorted[:, -1], inputMeshSorted[:, -1])).T
        mapping             = index_array[index_array[:, 0].argsort()][:, 1]

        return mapping


