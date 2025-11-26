from    dolfinx.fem                 import functionspace
from    basix.ufl                   import element, mixed_element
#from    ufl                         import finiteelement, MixedElement, triangle, VectorElement, tetrahedron
from    ufl                         import triangle, tetrahedron
from 	FELiCS.Misc.logging         import Logger

# Get the logger
logger = Logger.get_logger("felics")


def getElementShape(meshDim):
    """
    Get the element shape (triangle or tetrahedron) based on mesh dimension.

    Parameters
    ----------
    meshDim : int
        Geometric dimension of the mesh (2 or 3).

    Returns
    -------
    str
        'triangle' for 2D, 'tetrahedron' for 3D.

    Raises
    ------
    Exception
        If meshDim is not 2 or 3.
    """
    if meshDim == 2:
        return "triangle"
    elif meshDim == 3:
        return "tetrahedron"
    else:
        raise Exception("nDim is neither 2 or 3!")

def getElementType():
    """
    Get the default element type string for finite elements.

    Returns
    -------
    str
        The element type string ('CG' for continuous Galerkin).
    """
    return "CG"


def createFunctionSpace(mesh, degree = 2, dim = 1):
    """
    Create a finite element function space for the given mesh, order, and dimension. 
    ElementType is "CG".

    Parameters
    ----------
    mesh : FELiCSMesh
        The mesh object.
    degree : int, optional
        Polynomial degree of the element (default is 2).
    dim : int, optional
        Number of components (default is 1 for scalar fields).

    Returns
    -------
    dolfinx.fem.FunctionSpace
        The created function space.
    """

    elementShape = getElementShape(mesh.gdim)
    elementType  = getElementType()

    if dim == 1:
        return functionspace(
                    mesh.dolfinxMesh,
                    element(elementType, elementShape, degree)
                    )
    else:
        return  functionspace(
                     mesh.dolfinxMesh,
                     (elementType, degree,
                     (dim,))
                     )


class FEMSpaces():
    """
    Class containing all Finite Element (FE) spaces used in the simulation.

    Organizes and initializes various FE spaces such as scalar, vector, mixed, and export function spaces
    based on the mesh and problem configuration.

    **Initialize the FEMSpaces object**

    Parameters
    ----------
    param : ParameterClass
        Configuration object containing simulation parameters including Case, BoundaryCondition, Numerics, and Export settings.
    mesh : FELiCSMesh
        Computational mesh object used to define the function spaces.
    degree : int, optional
        Polynomial degree for the function spaces (default is 2).

    Attributes
    ----------
    element_shape : ufl.Cell
        Shape of the element (triangle or tetrahedron).
    elementTypeStr : str
        Type of finite element ('CG' for continuous Galerkin).
    _nVelocityComponents : int
        Number of velocity components defined by the case.
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

        # get names for state vector (to put in the mixed space)
        self.stateVectorNames = param.Case.StateVectorVariables

        ## Create vector spaces
        self._nVelocityComponents   = param.BoundaryCondition.nVelocityComponents # dim of velocity vector
        # Get the order of polynomials for velocity components
        if 'u' in param.getTransportedQuantityList():
            velocityOrder   = param.Numerics.PolynomialOrder['u']
        else:
            velocityOrder   = 2
        # Define FEM spaces for the velocity vector
        self.FunctionSpaceVectorVelocity       = createFunctionSpace(mesh, degree = velocityOrder, dim = self._nVelocityComponents)
        self.FunctionSpaceVectorVelocityExport = createFunctionSpace(mesh.exportMesh, degree = 1, dim = self._nVelocityComponents)
        self.FunctionSpaceVectorVelocityP1     = createFunctionSpace(mesh, degree = 1, dim = self._nVelocityComponents)

        ## create scalar spaces
        # Get function spaces for first order and second order elements.
        self.P1 = createFunctionSpace(mesh, degree = 1)
        self.P2 = createFunctionSpace(mesh, degree = 2)
        # Get function spaces for first order and second order elements on the export mesh.
        self.P1Export = createFunctionSpace(mesh.exportMesh, degree = 1)
        self.P2Export = createFunctionSpace(mesh.exportMesh, degree = 2)


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
        self.VMixedExport = functionspace(mesh.exportMesh.dolfinxMesh, MixedFE)

        # Set state vector names to the VMixed space
        self.VMixed.stateVectorNames = self.stateVectorNames
        self.VMixedExport.stateVectorNames = self.stateVectorNames


    def addCustomScalarSpaceToMixedSpace(self, mesh, order):
        """
        Add a scalar finite element of a specified polynomial order to the mixed function space.

        This method appends a new scalar finite element to the internal list of elements that make up the mixed
        function space and rebuilds both the primary and export mixed spaces.

        Parameters
        ----------
        mesh : FELiCSMesh
            The mesh on which the new mixed function space is defined.
        order : int
            Polynomial order of the scalar finite element to be added.

        Notes
        -----
        This method is most likely deprecated and is not used in the codebase.
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
            Number of spatial dims, usually equal to the number of velocity components.

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


