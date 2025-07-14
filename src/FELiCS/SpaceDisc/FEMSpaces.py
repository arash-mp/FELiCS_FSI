from    dolfinx.fem                 import functionspace
from    dolfinx.mesh                import refine
from    basix.ufl                   import element, mixed_element
#from    ufl                         import finiteelement, MixedElement, triangle, VectorElement, tetrahedron
from    ufl                         import triangle, tetrahedron
from    FELiCS.SpaceDisc.FELiCSMesh import FELiCSMesh
from    FELiCS.IO.Mapping           import Mapping
from 	FELiCS.Misc.logging         import Logger

# Get the logger
logger = Logger.get_logger("felics")


def getElementShape(meshDim):
    if meshDim == 2:
        return "triangle"
    elif meshDim == 3:
        return "tetrahedron"
    else:
        raise Exception("nDim is neither 2 or 3!")

def getElementType():
    return "CG"


def getFELiCSSpace(mesh, order = 2, dim = 1):

    elementShape = getElementShape(mesh.gdim)
    elementType  = getElementType()

    if dim == 1:
        return functionspace(
                    mesh.dolfinxMesh,
                    element(elementType, elementShape, order)
                    )
    else:
        return  functionspace(
                     mesh.dolfinxMesh,
                     (elementType, order,
                     (dim,))
                     )


class FEMSpaces():
    """
    Class containing all Finite Element (FE) spaces used in the simulation.

    This class organizes and initializes various FE spaces such as scalar, vector,
    mixed, and export function spaces based on the mesh and problem configuration.

    **Initialize the FEMSpaces object**

    Parameters
    ----------
    param : ParameterClass
        Configuration object containing simulation parameters including Case, 
        BoundaryCondition, Numerics, and Export settings.
    mesh : FELiCSMesh
        Computational mesh object used to define the function spaces.
        
    Attributes
    ----------
    element_shape : ufl.Cell
        Shape of the element (triangle or tetrahedron).
    elementTypeStr : str
        Type of finite element ('CG' for continuous Galerkin).
    _nVelocityComponents : int
        Number of velocity components defined by the case.
    exportMesh : FELiCSMesh
        Refined mesh used for exporting solution fields.
    FunctionSpaceVectorVelocity : dolfinx.fem.FunctionSpace
        Velocity vector space on the main mesh.
    FunctionSpaceVectorVelocityExport : dolfinx.fem.FunctionSpace
        Velocity vector space on the export mesh.
    FunctionSpaceVectorVelocityP1 : dolfinx.fem.FunctionSpace
        First-order vector space for velocity.
    VMixed : dolfinx.fem.FunctionSpace
        Mixed finite element space for all variables.
    VMixedExport : dolfinx.fem.FunctionSpace
        Export version of the mixed space.
    P1, P2 : dolfinx.fem.FunctionSpace
        First- and second-order scalar function spaces on the main mesh.
    P1Export, P2Export : dolfinx.fem.FunctionSpace
        First- and second-order scalar function spaces on the export mesh.
    mappingObj : Mapping
        Mapping object that links the function spaces to coordinate mappings.

    """
    def __init__(self, param, mesh, degree=2):
        logger.info('Defining FEM-spaces.')

        self.element_shape  = getElementShape(mesh.gdim)
        self.elementTypeStr = getElementType()


        ## create export mesh
        # Refine the mesh 
        logger.debug('Defining refined P1 export mesh.')
        refine_tuple        = refine(mesh.dolfinxMesh)
        exportMesh_dolfinx  = refine_tuple[0]
        # create new FELiCSMesh
        exportMesh          = FELiCSMesh(param.Case.CoordinateSystem,inputMesh=exportMesh_dolfinx)
        exportMesh.gdim     = param.Case.nDim
        # save mesh
        meshfileName        = f'{param.Case.AnalysisMode}_mesh.h5'
        exportMesh.saveInFELiCSFormat(f'{param.Export.ExportFolder}/{meshfileName}')
        self.exportMesh     = exportMesh
        

        ## create vector spaces
        self._nVelocityComponents   = param.BoundaryCondition.nVelocityComponents # dimension of velocity vector
        # Get the order of polynomials for velocity components
        if 'u' in param.getTransportedQuantityList():
            velocityOrder   = param.Numerics.PolynomialOrder['u']
        else:
            velocityOrder   = 2
        # Define FEM spaces for the velocity vector
        self.FunctionSpaceVectorVelocity       = getFELiCSSpace(mesh, order = velocityOrder, dimension = self._nVelocityComponents)
        self.FunctionSpaceVectorVelocityExport = getFELiCSSpace(exportMesh, order = 1, dimension = self._nVelocityComponents)
        self.FunctionSpaceVectorVelocityP1     = getFELiCSSpace(mesh, order = 1, dimension = self._nVelocityComponents)


        ## create scalar spaces
        # Get function spaces for first order and second order elements.
        self.P1 = getFELiCSSpace(mesh, order = 1)
        self.P2 = getFELiCSSpace(mesh, order = 2)
        # Get function spaces for first order and second order elements on the export mesh.
        self.P1Export = getFELiCSSpace(exportMesh, order = 1)
        self.P2Export = getFELiCSSpace(exportMesh, order = 2)


        ### create VMixed Space
        # Prepare list to contain all the finite element spaces
        self.MixedList = []
        # create elements for mixed list 
        for name in param.getTransportedQuantityList():
            logger.debug("Adding FEM space of order "+
                         f"{param.Numerics.PolynomialOrder[name]} for "+
                         f"{name}")
           
            if name == 'u':
                FE = element(
                    self.elementTypeStr,
                    self.element_shape,
                    param.Numerics.PolynomialOrder['u'],
                    shape = (self._nVelocityComponents,)
                    )   
            else:
                FE = element(
                    self.elementTypeStr,
                    self.element_shape,
                    param.Numerics.PolynomialOrder[name],
                    )
            self.MixedList.append(FE)        
        logger.debug('Mixed finite element list: ' + str(self.MixedList))
        # Create a element of the mixed function space
        MixedFE = mixed_element(self.MixedList)
        # Create a function space containing of mixed elements on the FEM mesh
        self.VMixed = functionspace(mesh.dolfinxMesh,MixedFE)
        # Create a function space containing of mixed elements on the export mesh
        self.VMixedExport = functionspace(exportMesh.dolfinxMesh, MixedFE)

        # create mapping
        self.mappingObj = Mapping(self)

    def addCustomScalarSpaceToMixedSpace(self, mesh, order):
        """
        Add a scalar finite element of a specified polynomial order to the mixed function space.

        This method appends a new scalar finite element to the internal list of
        elements that make up the mixed function space and rebuilds both the primary
        and export mixed spaces.

        Notes
        -----
        This method is most likely deprecated, it is not used in the codebase.

        Parameters
        ----------
        mesh : FELiCSMesh
            The mesh on which the new mixed function space is defined.
        order : int
            Polynomial order of the scalar finite element to be added.
        """
        # Create a element of the mixed function space
        MixedFE = MixedElement(self.MixedList)
        # Create a function space containing of mixed elements on the given mesh
        self.VMixed = FunctionSpace(mesh.dolfinxMesh,MixedFE)
#       # Create a function space containing of mixed elements on the export mesh
        self.VMixedExport = FunctionSpace(self.exportMesh.dolfinxMesh,MixedFE)

        #self.mappingObj = Mapping(self)

    def _projectField2allFEMSpaces(self, field, nfluctvar, nDim):
        """
        Project a given FEM field to all components of the mixed function space.

        This method is primarily used in fluctuation analyses, where a forcing
        or limiter field must be projected to each component of the mixed space,
        including both scalar and vector subspaces.

        Parameters
        ----------
        field : dolfinx.fem.Function
            The finite element function to be projected.
        nfluctvar : int
            Number of fluctuation variables in the system.
        nDim : int
            Number of spatial dimensions, usually equal to the number of velocity components.

        Returns
        -------
        fieldVMixed : dolfinx.fem.Function
            A function in the mixed space with the input field projected into all subspaces.

        Notes
        -----
        Originally used for projecting forcing/response fields in resolvent analyses.
        Supports fields of arbitrary polynomial order. This method is most likely deprecated, it is not used in the codebase.
        """
        
        # NOTE: (Simon) I think that this function is not required anymore

        from dolfinx.fem import Function

        fieldVMixed = Function(self.VMixed)
        for i in range(self.VMixed.num_sub_spaces):
            
            if self.VMixed.sub(i).num_sub_spaces > 0:   # Vector field, we go through sub-space

                for j in range(self.VMixed.sub(i).num_sub_spaces):
                    space_ii, map_ii = self.VMixed.sub(i).sub(j).collapse()
                    targetFunction = Function(space_ii)
                    targetFunction.interpolate(field)
                    fieldVMixed.x.array[map_ii] = targetFunction.x.array
                    
            else:
                space_i, map_i = self.VMixed.sub(i).collapse() 
                targetFunction = Function(space_i)
                targetFunction.interpolate(field)
                fieldVMixed.x.array[map_i] = targetFunction.x.array

        return fieldVMixed


