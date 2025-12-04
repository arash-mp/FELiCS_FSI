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
import  gmsh
import  h5py
import  numpy                   as np
from    FELiCS.Misc.tensorUtils import CoordinateSystem
from    mpi4py                  import MPI
from    ufl                     import SpatialCoordinate
from    dolfinx.mesh            import Mesh, refine 
from    FELiCS.Misc.logging     import Logger


comm = MPI.COMM_WORLD
# Get the logger
logger = Logger.get_logger("felics")

class FELiCSMesh:
    """
    Wrapper class for the DOLFINx mesh, with support for custom coordinate systems
    and FELiCS-specific utilities.

    This class allows loading and saving of mesh data, extraction of cell and
    coordinate information, and computation of connectivity for use in FEM
    simulations. It supports both Cartesian and Cylindrical coordinate systems
    and provides functionality to convert and export mesh data in a FELiCS
    HDF5-based format.

    **Initialize the FELiCSMesh object**

    Parameters
    ----------
    coordinateSystemName : str
        Name of the coordinate system (``"Cartesian"`` or ``"Cylindrical"``).
    meshFileName : str, optional
        Path to the Gmsh mesh file to load. Ignored if ``inputMesh`` is given.
    gdim : int, optional
        Geometric dimension of the mesh. If ``None`` and a mesh file is read,
        it is inferred from the Gmsh model.
    m : int, optional
        Spectral wave number used in the associated :class:`CoordinateSystem`.
        If nonzero, the logical dimension ``dim`` is set to ``gdim + 1``.
    inputMesh : dolfinx.mesh.Mesh, optional
        Existing DOLFINx mesh object to wrap instead of reading from a file.
        In this case, ``meshFileName`` and Gmsh are not used.

    Attributes
    ----------
    dolfinxMesh : dolfinx.mesh.Mesh
        The wrapped DOLFINx mesh object.
    facet_tags : dolfinx.mesh.MeshTags or None
        Boundary facet tags parsed from the Gmsh model when a mesh file is
        read. Not set when the mesh is provided via ``inputMesh``.
    gdim : int
        Geometric dimension of the mesh.
    dim : int
        Logical dimension of the system, including a possible spectral
        dimension (``dim == gdim`` if ``m == 0``).
    coordinateSystemName : str
        Name of the selected coordinate system.
    coordinateSystem : CoordinateSystem
        Tensor-based coordinate system initialized for the mesh.
    meshCells : numpy.ndarray
        Array containing the cell connectivity in FELiCS ordering; populated
        by :meth:`calcConnectivity`.
    """
    def __init__(self, coordinateSystemName, meshFileName=None, gdim = None, m=0, inputMesh=None):
        """
        Initialize the FELiCSMesh object, loading a mesh from file or using an
        existing DOLFINx mesh.

        Parameters
        ----------
        coordinateSystemName : str
            Name of the coordinate system (``"Cartesian"`` or ``"Cylindrical"``).
        meshFileName : str, optional
            Path to the Gmsh mesh file. Used only if ``inputMesh`` is ``None``.
        gdim : int, optional
            Geometric dimension of the mesh. If ``None`` when reading from a
            file, the dimension is obtained from the Gmsh model.
        m : int, optional
            Spectral wave number passed to the associated
            :class:`CoordinateSystem`. If nonzero, a spectral dimension is
            assumed and ``self.dim`` is set to ``gdim + 1``.
        inputMesh : dolfinx.mesh.Mesh, optional
            Existing DOLFINx mesh object. If provided, Gmsh is not used and
            the mesh is wrapped directly.
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
            logger.error('Coord. syst not yet implemented in tensor framework.')
            raise NotImplementedError('Coord. syst not yet implemented in tensor framework.')

        # update dimension if m!=0, i.e. if there is a spectral dimension
        if m!=0:
            self.dim = gdim + 1
        else:
            self.dim = self.gdim

        self._coordinates = self.coordinates()

    @property
    def exportMesh(self):
        """
        Lazy-loaded refined export mesh for exporting simulation results.

        On first access, a uniformly refined P1 mesh is created from this
        instance and wrapped in an :class:`ExportMesh`. Subsequent accesses
        return the cached export mesh.

        Returns
        -------
        ExportMesh
            An :class:`ExportMesh` instance associated with this base mesh.
        """
        if not hasattr(self, '_exportMesh') or self._exportMesh is None:
            self._exportMesh = ExportMesh(self)
        return self._exportMesh

    def setTrueDimension(self, dim):
        """
        Set the true dimension of the system, updating the coordinate system.

        This must be called if a spectral dimension is included, in which case
        the true dimension is higher than the geometrical dimension of the
        mesh.

        Parameters
        ----------
        dim : int
            True dimension of the system, including a possible spectral
            dimension.
        """

        self.dim = dim
        self.__coordinateSystem.setTrueDimension(dim)

    def saveInFELiCSFormat(self, filename):
        """
        Save the refined (export) mesh in the FELiCS HDF5-based format.

        The method uses the lazily constructed :attr:`exportMesh` to obtain a
        refined P1 mesh, then writes its vertex coordinates and triangular cell
        connectivity to an HDF5 file.

        Parameters
        ----------
        filename : str
            Path to the output HDF5 file where the mesh will be saved.

        """
        
        # TODO: save the DoFs corresponding to the different BCs

        # Compute additional export mesh properties
        coordinates     = self.exportMesh.coordinates()
        self.exportMesh.calcConnectivity()
        meshCells       = self.exportMesh.meshCells
        nDim            = coordinates.shape[1]
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
        Compute and update the internal mesh cell connectivity array.

        The connectivity is extracted from the DOLFINx topology and stored in
        :attr:`meshCells` as an ``(n_cells, n_vertices_per_cell)`` array.

        """

        connectivityCells = self.dolfinxMesh.topology.connectivity(self.gdim, 0)
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
        Retrieve boundary condition tags and corresponding boundary facets.

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
    def axisNames(self):
        """
        Names of the axes in the coordinate system.
        These include both spectral and mesh dimensions.
        Currently implemented systems are:
        
        - Cartesian:    ['x', 'y', 'z']
        - Cylindrical:  ['x', 'r', 't']

        Returns
        -------
        list of str
            List of coordinate names (e.g., ['x', 'y', 'z']).
        """
        if self.coordinateSystemName =='Cartesian':
            coordinateNames = ['x']
            if self.dim > 1:
                coordinateNames.append('y')
            if self.dim > 2:
                coordinateNames.append('z')
        
        elif self.coordinateSystemName =='Cylindrical':
            coordinateNames = ['x']
            if self.dim > 1:
                coordinateNames.append('r')
            if self.dim > 2:
                coordinateNames.append('t')
                
        else:
            logger.error('Coord. syst not yet implemented in tensor framework.')
            raise NotImplementedError('Coord. syst not yet implemented in tensor framework.')

        return coordinateNames
    
    @property
    def meshAxisNames(self):
        """
        Names of the axes in the mesh's coordinate system.
        Currently implemented systems are:

        - Cartesian:    ['x', 'y', 'z']
        - Cylindrical:  ['x', 'r', 't']

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
                coordinateNames.append('t')
                
        else:
            logger.error('Coord. syst not yet implemented in tensor framework.')
            raise NotImplementedError('Coord. syst not yet implemented in tensor framework.')

        return coordinateNames


class ExportMesh(FELiCSMesh):
    """
    Refined export mesh built from a base :class:`FELiCSMesh`.

    This subclass creates a uniformly refined P1 mesh from an existing
    :class:`FELiCSMesh` instance. It is primarily used for exporting
    simulation results on a finer mesh than the one used for computation.
    """
    
    def __init__(self, base_mesh: FELiCSMesh):
        """
        Construct an ExportMesh by refining a base FELiCSMesh.

        The underlying DOLFINx mesh of ``base_mesh`` is uniformly refined using
        :func:`dolfinx.mesh.refine`, and the resulting mesh is wrapped as a new
        :class:`FELiCSMesh` instance with the same coordinate system and
        geometric dimension.

        Parameters
        ----------
        base_mesh : FELiCSMesh
            The base mesh from which the refined export mesh is created.
        """
        logger.debug("Initializing refined P1 export mesh.")
        
        # Refine the mesh
        base_mesh.dolfinxMesh.topology.create_entities(1)
        refine_tuple        = refine(base_mesh.dolfinxMesh)
        exportMesh_dolfinx  = refine_tuple[0]

        # Create the new ExportMesh instance
        super().__init__(
            coordinateSystemName    = base_mesh.coordinateSystemName,
            inputMesh               = exportMesh_dolfinx,
            gdim                    = base_mesh.gdim,
        )
