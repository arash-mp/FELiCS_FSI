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

# Get the logger
logger = Logger.get_logger("felics")

class FELiCSMesh(Mesh):
    '''
    This class is an extension to the fenics mesh class
    '''
    def __init__(self, coordinateSystem, filename=None, gdim=0, m=0, inputMesh=None):
        if inputMesh is None:
            # Initialize gmsh and suppress its output
            gmsh.initialize()
            gmsh.option.setNumber("General.Terminal", 0)  # Disable console log
            gmsh.option.setNumber("General.Verbosity", 0) # Disable all logging
            
            if __version__.find('0.4') >= 0:
                logger.info("DEPRECATED: Dolfinx version <0.5.0 is used.")
                from FELiCSGUI.gmsh_helpers import read_from_msh, extract_gmsh_geometry #, gmsh_model_to_mesh
                mesh, cell_tags, hi, facet_tags = read_from_msh(filename, cell_data=True, facet_data=True, gdim=gdim)
                self.coordinatesGMSH = extract_gmsh_geometry(gmsh.model)
            else:
                logger.info(f"Opening mesh file: {filename}")
                gmsh.open(filename)
                
                # Get mesh statistics
                nodes           = gmsh.model.mesh.getNodes()
                elements        = gmsh.model.mesh.getElements()
                num_nodes       = len(nodes[0])
                num_elements    = sum(len(elements[1][i]) for i in range(len(elements[1])))
                logger.info(f"Mesh contains {num_nodes} nodes and {num_elements} elements")

                from dolfinx.io import gmshio
                mesh_comm = MPI.COMM_WORLD
                model_rank = 0
                mesh, _, facet_tags = gmshio.model_to_mesh(gmsh.model, mesh_comm, model_rank, gdim=gdim)
                
            try:    #try new version of dolfinx 
                super().__init__(mesh, mesh.ufl_domain())
                newMesh = Mesh(mesh, mesh.ufl_domain())
            except: #use old language 
                logger.info("DEPRECATED: Mesh module from dolfinx version <0.7.0 is used.")
                super().__init__(MPI.COMM_WORLD, mesh.topology, mesh.geometry, mesh.ufl_domain())
            #   #Mesh.__init__(self, MPI.COMM_WORLD, mesh.topology, mesh.geometry)
            
            try: 
                self.dolfinxMesh  = mesh
                self._ccp_object  = mesh._cpp_object
            except:
                logger.info("DEPRECATED: Mesh module from dolfinx version <0.7.0 is used.")
                self.dolfinxMesh  = self
                self._cpp_object  = mesh
                
            self.facet_tags = facet_tags
            self.gdim = gdim
            self._ufl_domain = mesh._ufl_domain
            self.calcConnectivity()
            # save the coordinates in gmsh order:
            gmsh.open(filename)
        else:
            try:    #try new version of dolfinx 
                super().__init__(inputMesh, inputMesh.ufl_domain())
            except: #use old language 
                logger.info("DEPRECATED: Mesh module from dolfinx version <0.7.0 is used.")
                super().__init__(MPI.COMM_WORLD, inputMesh.topology, inputMesh.geometry, inputMesh.ufl_domain())
            self.gdim = inputMesh.topology.dim
            try: 
                self.dolfinxMesh = inputMesh
                self._ccp_object  = inputMesh._cpp_object
            except:
                self.dolfinxMesh = self 
                logger.info("DEPRECATED: Mesh module from dolfinx version <0.7.0 is used.")
                self._cpp_object  = inputMesh
        x = SpatialCoordinate(self)
        # Define tensor coordinate system, we always assume the third dimension to be homogenous
        if coordinateSystem =='Cartesian':
            self.__coordinateSystem = CoordinateSystem(
                                    x, 
                                    coordinateSystem.lower(), 
                                    m = m,
                                    mesh_dims = (1, 1, 0),
                                    )
        elif coordinateSystem =='Cylindrical':
            self.__coordinateSystem = CoordinateSystem(
                                    x,
                                    "cylindricalfelics", 
                                    m = m,
                                    mesh_dims = (1, 1, 0),
                                    )
        else:
            printError('Coord. syst not yet implemented in tensor framework.')
        self._coordinates = self.coordinates()

    def saveInFELiCSFormat(self, filename):
        '''
        This function saves the computational mesh in the FELiCS format

        Function arguments:
        - filename: The path where to save the mesh

        Function returns:
        '''
        
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
        this method calculates the meshCells array in the fenics representation
        """
        connectivityCells = self.topology.connectivity(2, 0)
        try:    #try new version of dolfinx 
            self.meshCells = connectivityCells.array.reshape(
                [self.topology.original_cell_index.shape[0], self.topology.cell_types[0].value])
        except: #use old language. TODO: handle DEPRECATED stuff uniformly
            logger.info("DEPRECATED: Mesh module from dolfinx version <0.7.0 is used.")
            self.meshCells = connectivityCells.array.reshape(
                [self.topology.original_cell_index.shape[0], self.topology.cell_type.value])

    def cells(self):
        """
        this methods returns the cell-connectivity information
        """
        self.calcConnectivity()
        return self.meshCells

    def coordinates(self):
        """
        This method acts as a getter-method for the vertex-coordinates.
        """
        return self.geometry.x[:, 0:self.gdim]

    @property
    def coordinateSystem(self):
        return self.__coordinateSystem
