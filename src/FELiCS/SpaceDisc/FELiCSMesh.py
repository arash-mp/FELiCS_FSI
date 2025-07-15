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
    '''
    This class is a wrapper to the fenics mesh class
    '''
    def __init__(self, coordinateSystemName, meshFileName=None, gdim = None, m=0, inputMesh=None):
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
                                    mesh_dims = (1, 1, 0),
                                    )
        elif coordinateSystemName =='Cylindrical':
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
        connectivityCells = self.dolfinxMesh.topology.connectivity(2, 0)
        topology          = self.dolfinxMesh.topology
        self.meshCells    = connectivityCells.array.reshape(
                            [topology.original_cell_index.shape[0], topology.cell_type.value])

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
        return self.dolfinxMesh.geometry.x[:, 0:self.gdim]
    
    def getBCInfo(self):
        """ 
        This function provides both the IDs of the boundary conditions 
        and also the boundary nodes.
        """
        from numpy import unique
        return unique(self.facet_tags.values), self.facet_tags 

    @property
    def coordinateSystem(self):
        return self.__coordinateSystem
