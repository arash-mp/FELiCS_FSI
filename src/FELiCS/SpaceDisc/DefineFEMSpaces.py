from    dolfinx.fem                 import FunctionSpace, VectorFunctionSpace
from    dolfinx.mesh                import refine
from    ufl                         import FiniteElement, MixedElement, triangle, VectorElement, tetrahedron
from    FELiCS.SpaceDisc.FELiCSMesh import FELiCSMesh
from    FELiCS.IO.Mapping           import Mapping
from 	FELiCS.Misc.logging         import Logger

# Get the logger
logger = Logger.get_logger("felics")

class FEMSpacesClass():
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

        if param.Case.nDim==2:
            self.element_shape = triangle
        elif param.Case.nDim==3:
            self.element_shape = tetrahedron
        else:
            raise Exception("nDim is neither 2 or 3!")
        
        self.elementTypeStr         = 'CG'
        self._nVelocityComponents   = param.BoundaryCondition.nVelocityComponents

        # Refine the mesh to get the exportMesh:
        logger.debug('Defining refined P1 export mesh.')
        exportMesh_dolfinx  = refine(mesh.dolfinxMesh)
        exportMesh          = FELiCSMesh(param.Case.CoordinateSystem,inputMesh=exportMesh_dolfinx)
        exportMesh.gdim     = param.Case.nDim
        meshfileName        = f'{param.Case.AnalysisMode}_mesh.h5'
        exportMesh.saveInFELiCSFormat(f'{param.Export.ExportFolder}/{meshfileName}')
        self.exportMesh     = exportMesh
        
        # Get the order of polynomials for velocity components
        if 'u' in param.getTransportedQuantityList():
            velocityOrder   = param.Numerics.PolynomialOrder['u']
        else:
            velocityOrder   = 2
        
        # Define FEM spaces for the velocity vector
        self.FunctionSpaceVectorVelocity = VectorFunctionSpace(
            mesh.dolfinxMesh,
            (self.elementTypeStr, velocityOrder),
            dim = self._nVelocityComponents,
            )
        self.FunctionSpaceVectorVelocityExport = VectorFunctionSpace(
            exportMesh.dolfinxMesh,
            (self.elementTypeStr, 1),
            dim = self._nVelocityComponents
            )

        # Prepare list to contain all the finite element spaces
        self.MixedList = []

        # Then a scalar for all the remaining quantities
        self.FunctionSpaceList          = []
        self.FunctionSpaceListExport    = []
        for name in param.getTransportedQuantityList():
            logger.debug("Adding FEM space of order "+
                         f"{param.Numerics.PolynomialOrder[name]} for "+
                         f"{name}")
            
            if name == 'u':
                FE = VectorElement(
                    self.elementTypeStr,
                    self.element_shape,
                    param.Numerics.PolynomialOrder['u'],
                    dim = self._nVelocityComponents,
                    )   
            else:
                FE = FiniteElement(
                    self.elementTypeStr,
                    self.element_shape,
                    param.Numerics.PolynomialOrder[name],
                    )
            self.MixedList.append(FE)        
        logger.debug('Mixed finite element list: ' + str(self.MixedList))

        # Create a element of the mixed function space
        MixedFE = MixedElement(self.MixedList)

        # Create a function space containing of mixed elements on the FEM mesh
        self.VMixed = FunctionSpace(mesh.dolfinxMesh,MixedFE)

        self.FunctionSpaceVectorVelocityP1 = VectorFunctionSpace(
            mesh.dolfinxMesh,
            (self.elementTypeStr, 1),
            dim=self._nVelocityComponents,
            )
        self.FunctionSpaceListExport = []

#       # Get function spaces for first order and second order elements.
        self.P1 = FunctionSpace(
            mesh.dolfinxMesh,
            FiniteElement(self.elementTypeStr, self.element_shape, 1)
            )
        self.P2 = FunctionSpace(
            mesh.dolfinxMesh,
            FiniteElement(self.elementTypeStr, self.element_shape, 2)
            )

#       # Create a function space containing of mixed elements on the export mesh
        self.VMixedExport = FunctionSpace(exportMesh.dolfinxMesh, MixedFE)

        # Get function spaces for first order and second order elements on the export mesh.
        self.P1Export = FunctionSpace(
            exportMesh.dolfinxMesh,
            FiniteElement(self.elementTypeStr, self.element_shape, 1)
            )
        self.P2Export = FunctionSpace(
            exportMesh.dolfinxMesh,
            FiniteElement(self.elementTypeStr, self.element_shape, 2)
            )

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


