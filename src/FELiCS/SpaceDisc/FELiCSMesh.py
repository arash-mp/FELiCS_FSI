import  gmsh
import  h5py
import  numpy                   as np
from    FELiCS.Misc.tensorUtils import CoordinateSystem
from    FELiCS.Misc.functions   import printDeprecatedWarning, printError
from    mpi4py                  import MPI
from    ufl                     import SpatialCoordinate
from    dolfinx                 import __version__
from    dolfinx.mesh            import Mesh 
from    FELiCS.Misc.logging     import Logger


comm = MPI.COMM_WORLD
# Get the logger
logger = Logger.get_logger("felics")

class FELiCSMesh:
    """
    Wrapper class for the DOLFINx mesh, with support for custom coordinate systems and FELiCS-specific utilities.

    This class allows loading and saving of mesh data, extraction of cell and coordinate information, 
    and computation of connectivity for use in FEM simulations. It supports both Cartesian and Cylindrical 
    coordinate systems and includes functionality to convert and export mesh data in FELiCS format.

    **Initialize the FELiCSMesh object**

    Parameters
    ----------
    coordinateSystem : str
        Name of the coordinate system ('Cartesian' or 'Cylindrical').
    meshFileName : str, optional
        Path to the mesh file to load.
    gdim : int, optional
        Geometric dimension of the mesh. If None, it is inferred from the Gmsh model.
    m : int, optional
        Parameter used in the coordinate system configuration.
    inputMesh : dolfinx.mesh.Mesh, optional
        Existing DOLFINx mesh object to wrap instead of reading from a file.

    Attributes
    ----------
    dolfinxMesh : dolfinx.mesh.Mesh
        The wrapped DOLFINx mesh object.
    facet_tags : dolfinx.mesh.MeshTags
        Boundary facet tags parsed from the mesh.
    gdim : int
        Geometric dimension of the mesh.
    coordinateSystemName : str
        Name of the selected coordinate system.
    _coordinates : numpy.ndarray
        Cached array of mesh vertex coordinates.
    """

    def __init__(self, coordinateSystemName, meshFileName=None, gdim = None, m=0, inputMesh=None):
        """
        Initializes the FELiCSMesh object, loading a mesh from file or using an existing mesh.

        Parameters
        ----------
        coordinateSystem : str
            Coordinate system type ('Cartesian' or 'Cylindrical').
        meshFileName : str, optional
            Path to the mesh file.
        gdim : int, optional
            Geometric dimension of the mesh.
        m : int, optional
            A parameter used when constructing the coordinate system.
        inputMesh : dolfinx.mesh.Mesh, optional
            An existing DOLFINx mesh object.
        """

        if inputMesh is None:
            # Initialize gmsh and suppress its output
            gmsh.initialize()
            gmsh.option.setNumber("General.Terminal", 0)  # Disable console log
            gmsh.option.setNumber("General.Verbosity", 0) # Disable all logging
            
            logger.info(f"Opening mesh file: {meshFileName}")
            gmsh.open(meshFileName)
            
            # Get mesh statistics
            nodes           = gmsh.model.mesh.getNodes()
            elements        = gmsh.model.mesh.getElements()
            num_nodes       = len(nodes[0])
            num_elements    = sum(len(elements[1][i]) for i in range(len(elements[1])))
            logger.info(f"Mesh contains {num_nodes} nodes and {num_elements} elements")

            from dolfinx.io import gmshio
            mesh_comm = MPI.COMM_WORLD
            model_rank = 0
            if gdim == None:
                gdim = gmsh.model.getDimension()
            mesh, _, facet_tags = gmshio.model_to_mesh(gmsh.model, mesh_comm, model_rank, gdim=gdim)
                
            self.dolfinxMesh = Mesh(mesh, mesh.ufl_domain())._cpp_object
            self._cpp_object = self.dolfinxMesh._cpp_object

            self.facet_tags = facet_tags
            self.gdim = gdim
            self._ufl_domain = mesh._ufl_domain
            self.calcConnectivity()
            # save the coordinates in gmsh order:
            gmsh.open(meshFileName)
        else:
            #self.dolfinxMesh = Mesh(inputMesh, inputMesh.ufl_domain())
            self.dolfinxMesh = inputMesh 
            self.gdim        = inputMesh.topology.dim
            self._cpp_object = self.dolfinxMesh._cpp_object
        x = SpatialCoordinate(self.dolfinxMesh)
        # Define tensor coordinate system, we always assume the third dimension to be homogenous
        self.coordinateSystemName = coordinateSystemName
        if coordinateSystemName =='Cartesian':
            self.__coordinateSystem = CoordinateSystem(
                                    x, 
                                    coordinateSystemName.lower(), 
                                    m = m,
                                    gdim = self.gdim,
                                    )
        elif coordinateSystemName =='Cylindrical':
            self.__coordinateSystem = CoordinateSystem(
                                    x,
                                    "cylindricalfelics", 
                                    m = m,
                                    gdim = self.gdim,
                                    )
        else:
            printError('Coord. syst not yet implemented in tensor framework.')
        self._coordinates = self.coordinates()

    def saveInFELiCSFormat(self, filename):
        """
        Saves the computational mesh in the FELiCS HDF5-based format.

        Parameters
        ----------
        filename : str
            Path to the output HDF5 file where the mesh will be saved.

        Notes
        -----
        - Only vertex coordinates and triangular cell connectivity are saved.
        - Currently, degrees of freedom (DoFs) for boundary conditions are not included.
        """
        
        # TODO: save the DoFs corresponding to the different BCs
        
        coordinates = self.coordinates()
        self.calcConnectivity()
        meshCells = self.meshCells
        nDim = coordinates.shape[1]
        coordinateNames = ['x']
        if nDim > 1:
            coordinateNames.append('y')
        if nDim > 2:
            coordinateNames.append('z')
        hf = h5py.File(filename, 'w')
        g1 = hf.create_group('coordinates')
        for i_coordinateName,coordinateName in enumerate(coordinateNames):
            g1.create_dataset(coordinateName,data=coordinates[:,i_coordinateName])
        g2 = hf.create_group('cells')
        g2.create_dataset('triangles',data=np.array(meshCells))
        hf.close()

    def calcConnectivity(self):
        """
        Calculates and updates the internal mesh cell connectivity array.

        Notes
        -----
        The meshCells attribute is reshaped from the mesh's topology connectivity.
        This is essential for writing mesh data or querying cell connectivity.
        """

        connectivityCells = self.dolfinxMesh.topology.connectivity(2, 0)
        topology          = self.dolfinxMesh.topology
        self.meshCells    = connectivityCells.array.reshape(
                            [topology.original_cell_index.shape[0], topology.cell_type.value])

    def cells(self):
        """
        Returns the cell connectivity array of the mesh.

        Returns
        -------
        numpy.ndarray
            Array of mesh cells with vertex indices.
        """

        self.calcConnectivity()
        return self.meshCells

    def coordinates(self):
        """
        Retrieves the vertex coordinates of the mesh.

        Returns
        -------
        numpy.ndarray
            Array of vertex coordinates, truncated to the mesh's geometric dimension.
        """

        return self.dolfinxMesh.geometry.x[:, 0:self.gdim]
    
    def getBCInfo(self):
        """
        Retrieves boundary condition tags and corresponding boundary facets.

        Returns
        -------
        numpy.ndarray
            Unique boundary condition IDs.
        dolfinx.mesh.MeshTags
            Mesh tags object representing the boundary facets.
        """

        from numpy import unique
        return unique(self.facet_tags.values), self.facet_tags 

    @property
    def coordinateSystem(self):
        """
        Coordinate system object associated with the mesh.

        Returns
        -------
        CoordinateSystem
            Tensor-based coordinate system initialized for the mesh.
        """

        return self.__coordinateSystem
    
    @property
    def coordinateNames(self):
        """
        Names of the coordinates in the mesh's coordinate system.
        Currently implemented systems are:
            - Cartesian:    ['x', 'y', 'z']
            - Cylindrical:  ['x', 'r', 'theta']

        Returns
        -------
        list of str
            List of coordinate names (e.g., ['x', 'y', 'z']).
        """
        if self.coordinateSystemName =='Cartesian':
            coordinateNames = ['x']
            if self.gdim > 1:
                coordinateNames.append('y')
            if self.gdim > 2:
                coordinateNames.append('z')
        
        elif self.coordinateSystemName =='Cylindrical':
            coordinateNames = ['x']
            if self.gdim > 1:
                coordinateNames.append('r')
            if self.gdim > 2:
                coordinateNames.append('theta')
                
        else:
            logger.error('Coord. syst not yet implemented in tensor framework.')
            raise NotImplementedError('Coord. syst not yet implemented in tensor framework.')

        return coordinateNames
