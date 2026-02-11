#  ___________________________________   _______________________________________________________
# /-----------------------------------\ /-------------------------------------------------------\
# |   (         (                (     |  This source code is part of FELiCS                     |
# |   )\ )     ) )        (      )\ )  |  (F)inite (E)lement (Li)nearized (C)ombustion (S)olver  |  
# |  (()/(  (  (()/( (    )\   (()/(   |                                                         |  
# |  /(_)) )\  /(_)))\  (((_)  /(_))   |  Licensed under the GNU GPLv3                           |
# |  (_)_)((_) (_)) ((_) )\___ (_))    |                                                         |
# |  | __|| __|| |   (_)((/ __|/ __|   |  (C) 2018-2025: The FELiCS Developers (www.felics.eu)   |
# |  | _| | _| | |__ | | | (__ \__ \   |  Visit          www.felics.eu                           |
# |  |_|  |___||____||_|  \___||___/   |  Contact        info@felics.eu                          |
# \___________________________________/ \_______________________________________________________/
#
import  numpy as np
from 	FELiCS.Misc.logging         import Logger

# Get the logger
logger = Logger.get_logger("felics")

class Mapping:
    """
    Utility class for constructing index mappings between FEM spaces or DOF coordinates.

    This class contains static methods that compute permutation arrays used to transfer
    values between meshes or finite element spaces with identical layouts but differing
    internal DOF orderings.

    **Initialize the Mapping object**

    This class is not intended to be instantiated; all functionality is provided through
    static methods.

    Notes
    -----
    - The mapping assumes that input and output meshes have identical geometric DOF locations, possibly in different order.
    - Coordinates are rounded to 11 decimal places before comparison to avoid small floating-point inconsistencies.
    """

    @staticmethod
    def calculate_mapping_from_spaces(
            inputSpace,
    outputSpace,
    ):
        """
        Compute a DOF index mapping between two finite element spaces.

        This method extracts the coordinates of DOFs from both spaces and delegates the
        mapping computation to :meth:`calculateMappingFromDofs`.

        Parameters
        ----------
        inputSpace : object
            Finite element space providing the input DOF coordinates. Must implement
            ``tabulate_dof_coordinates()``.
        outputSpace : object
            Finite element space providing the output DOF coordinates. Must implement
            ``tabulate_dof_coordinates()``.

        Returns
        -------
        numpy.ndarray
            One-dimensional array of indices such that ``output[index] = input[mapping[index]]``.

        Raises
        ------
        ValueError
            If the DOF coordinate arrays from the two spaces are incompatible in dimension.
        """

        # Get numpy arrays of dof coordinates
        inputDofs           = inputSpace.tabulate_dof_coordinates()
        outputDofs          = outputSpace.tabulate_dof_coordinates()

        return Mapping.calculate_mapping_from_dofs(
        inputDofs,
        outputDofs,
        )

    @staticmethod
    def calculate_mapping_from_dofs(
            inputDofs,
    outputDofs,
    ):
        """
        Compute a DOF index mapping from coordinate arrays.

        The algorithm constructs a bijection between two sets of DOF coordinates by:
        1. Copying the coordinate arrays,
        2. Appending the original DOF indices,
        3. Sorting by spatial coordinates,
        4. Aligning sorted entries,
        5. Extracting the permutation that maps input DOFs to output DOFs.

        Parameters
        ----------
        inputDofs : numpy.ndarray
            Array of shape ``(n_dofs, dim)`` containing input DOF coordinates.
        outputDofs : numpy.ndarray
            Array of shape ``(n_dofs, dim)`` containing output DOF coordinates.

        Returns
        -------
        numpy.ndarray
            One-dimensional integer array representing the mapping from input DOF indices
            to output DOF indices.

        Raises
        ------
        ValueError
            If the coordinate arrays have different lengths or incompatible shapes.

        Notes
        -----
        - Coordinates are rounded to 11 decimal places before matching.
        - The mapping requires exact geometric correspondence of DOFs between input and
          output meshes.
        """

        # Get numpy arrays of dof coordinates
        inputMesh           = np.copy(inputDofs)
        outputMesh          = np.copy(outputDofs)

        # Append indices as last column
        inputMesh           = np.append(
        inputMesh,
        np.arange(len(inputMesh))[:, None],
        axis=1,
        ).round(11)
        outputMesh          = np.append(
        outputMesh,
        np.arange(len(outputMesh))[:, None],
        axis=1,
        ).round(11)

        # Sort by x,y,z
        inputMeshSorted     = inputMesh[np.lexsort((inputMesh[:,2], inputMesh[:,1], inputMesh[:,0]))].astype(int)
        outputMeshSorted    = outputMesh[np.lexsort((outputMesh[:, 2], outputMesh[:, 1], outputMesh[:, 0]))].astype(int)

        # Find indices of sorted inputMesh in outputMesh
        index_array         = np.vstack((outputMeshSorted[:, -1], inputMeshSorted[:, -1])).T
        mapping             = index_array[index_array[:, 0].argsort()][:, 1]

        return mapping


