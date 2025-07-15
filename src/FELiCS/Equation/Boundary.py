import  json
import  numpy as np
from    enum                         import Enum
from    petsc4py.PETSc               import ScalarType
from 	FELiCS.Misc.logging          import Logger
# Get the logger
logger = Logger.get_logger("felics")



class BoundaryHandler():
    """
    Boundary condition manager for FEM simulations.

    Handles the initialization, validation, and construction of boundary condition
    objects for each boundary defined in the mesh and specified in a JSON file.

    **Initialize the BoundaryHandler object**

    Parameters
    ----------
    variables : list of tuples
        Each tuple contains the name of a variable and its components (e.g., [("u", ["x", "y"])]).
    mesh : dolfinx.Mesh
        The mesh object with tagged facets defining boundaries.
    BCsFilePath : str
        Path to the JSON file containing boundary condition specifications.

    Attributes
    ----------
    facet_tags : dolfinx.mesh.meshtags
        Facet tags from the mesh that identify boundary regions.
    IDs : numpy.ndarray
        Unique boundary IDs found in the mesh.
    variables : list of tuples
        Stored list of simulation variables.
    boundaryList : list
        List of boundary condition objects created from the boundary JSON file.

    Notes
    -----
    Raises a `ValueError` if:
    - The boundary ID in the file does not exist in the mesh.
    - The boundary condition name is invalid.
    """
        
    def __init__(self, variables, mesh, BCsFilePath):

        ###### initialize list of boundaries ########

        #1. get info from felics mesh and store variables
        self.facet_tags             = mesh.facet_tags
        self.IDs                    = np.unique(self.facet_tags.values)
        self.variables              = variables

        #2. read bc file 
        # TODO Sophie: 
        # 1. throw error if file is not there
        # 2. throw error if file is empty or not in the correct format
        logger.debug(f"Reading boundary conditions from '{BCsFilePath}'")
        file    = open(BCsFilePath) 
        BCsInfo = json.load(file)
        

        #3. create boundary object for each boundary
        self.boundaryList = []
        for ID in BCsInfo:
            id_int = int(ID)
            if id_int not in self.IDs: # check if the boundary ID exists in mesh file
                logger.error(f"The given boundary ID '{ID}' in your boundary file does not exist. The boundary IDs given from your mesh file are: '{self.IDs}'. ")
                raise ValueError("One of the given boundary IDs does not exist. Please read the FELiCS error message for details.")
            info = BCsInfo[ID]
            name = info["name"].lower() # make the info in the bc file non-case sensitive
            if   name  == "custom":
                bc = Custom(id_int, BCsInfo[ID], self)
            elif name  == "zerodirichlet":
                bc = ZeroDirichlet(id_int, BCsInfo[ID], self )
            elif name  == "wall":
                bc = Wall(id_int, BCsInfo[ID], self)
            elif name  == "symmetry":
                bc = Symmetry(id_int, BCsInfo[ID], self)
            elif name  == "none":
                bc = BoundaryCondition(id_int, BCsInfo[ID], self) # this is the default boundary condition, nothing is done for any variable
            else:
                logger.error(f"The boundary condition with name '{name}' does not exist. Please choose from the following list: [custom, zeroDirichlet, wall, symmetry, none]. The names are not case sensitive. ")
                raise ValueError("One of the set boudary conditions does not exist. Please read the FELiCS error message for details.")
            self.boundaryList.append(bc)


    def getListOfBoundaries(self):
        """
        Return the list of all boundary condition objects.

        Returns
        -------
        list
        List of boundary condition objects associated with each mesh boundary.
        """

        return self.boundaryList


    def getListOfDirichletBCsForDolfinx(self, functionSpace):
        """
        Construct the list of Dirichlet boundary conditions for Dolfinx.

        Loops through each variable and boundary, and extracts Dirichlet boundary conditions
        for scalar and vector components as required.

        Parameters
        ----------
        functionSpace : dolfinx.fem.FunctionSpace
            The function space from which subspaces will be extracted for boundary condition application.

        Returns
        -------
        list
            List of DirichletBC objects from dolfinx.fem.dirichletbc.

        Notes
        -----
        This function supports mixed function spaces and variables with multiple components.
        Logs a debug message each time a boundary condition is added.
        """

        from dolfinx.fem      import dirichletbc, locate_dofs_topological
        BCs = []
        #TODO: give a good description of what is done here:
        for boundary in self.boundaryList:
            for var in self.variables:
                index = self.variables.index(var)
                if len(var[1])==0 and boundary.types[index][0] == BoundaryType.DIRICHLET:
                    value = boundary.values[index][0]
                    space = functionSpace.sub(index)
                    dofs  = locate_dofs_topological(space, 1, self.facet_tags.indices[self.facet_tags.values==boundary.ID])
                    BCs.append(dirichletbc(value, dofs, space))
                    logger.debug("Adding Dirichlet BC for "+var[0]+ " in equation "+str(index)+" with value "+str(value)+" on boundary with index "+str(boundary.ID))
                else: 
                    for index2 in range(len(var[1])):
                        if boundary.types[index][index2] == BoundaryType.DIRICHLET:
                            value = boundary.values[index][index2]
                            space = functionSpace.sub(index).sub(index2)
                            dofs  = locate_dofs_topological(space, 1, self.facet_tags.indices[self.facet_tags.values==boundary.ID])
                            BCs.append(dirichletbc(value, dofs, space))
                            logger.debug("Adding Dirichlet BC for "+var[0]+var[1][index2] + " in equation "+str(index)+" with value "+str(value)+" on boundary with index "+str(boundary.ID))

        return BCs
    

    def getListOfNonlinearBoundaries(self):
        """
        Placeholder for returning nonlinear boundary condition objects.

        Notes
        -----
        To be implemented for base flow computations.
        """
        #TODO Sophie: fill out later for base flow computations
        pass


    def getListOfNonlinearDirichletBCsForDolfinx(self):
        """
        Placeholder for returning nonlinear DirichletBCs for Dolfinx.

        Notes
        -----
        To be implemented for base flow computations.
        """
        #TODO Sophie: fill out later for base flow computations
        pass


class BoundaryType(Enum):
    """
    Enum for different types of boundary conditions.

    Attributes
    ----------
    NONE : int
        No boundary condition.
    DIRICHLET : int
        Dirichlet boundary condition.
    NEUMANN : int
        Neumann boundary condition (currently equivalent to NONE).

    Notes
    -----
    The NEUMANN type currently has no distinct behavior; this may change in future updates.
    """

    NONE      = 0
    DIRICHLET = 1
    NEUMANN   = 2 # equal to "none" at the moment; this should be changed in the future, and also communicated really well
    #MIXED     = 3


class BoundaryCondition():
    """
    Generic boundary condition class.

    Serves both as the base class for all boundary conditions and
    as the 'none' type condition when no boundary constraints are imposed.

    **Initialize the BoundaryCondition object**

    Parameters
    ----------
    boundaryID : int
        Identifier for the boundary.
    boundaryInfo : dict
        Dictionary with boundary specifications from the JSON file.
    boundaryHandler : BoundaryHandler
        Reference to the parent BoundaryHandler.

    Attributes
    ----------
    ID : int
        Boundary ID.
    info : dict
        Boundary configuration details.
    bH : BoundaryHandler
        Reference to the boundary handler.
    types : list of lists
        Boundary condition types for each variable/component.
    values : list of lists
        Boundary values for each variable/component.
    """


    # This class serves two functions:
    # 1. it is the parent class of all boundary conditions (all variables are initialized, all types are "none", all values are "0")
    # 2. it is the boundary condition "None" for all variables
    def __init__(self, boundaryID, boundaryInfo, boundaryHandler):
        self.ID   = boundaryID
        self.info = boundaryInfo
        self.bH   = boundaryHandler

        # initialize types and values lists with "NONE" and "0"
        self.types         = []
        self.values        = []
        for var in self.bH.variables:
            name, components = var
            comp_types  = []
            comp_values = []
            if len(components) > 0:
                for comp in components:
                    comp_types.append(BoundaryType.NONE)
                    comp_values.append(0)
            else:
                comp_types.append(BoundaryType.NONE)
                comp_values.append(0)
            self.types.append(comp_types)
            self.values.append(comp_values)



class Custom(BoundaryCondition):
    """
    Custom boundary condition with user-defined specifications.

    Interprets the "specifics" field in the boundary JSON file to assign boundary
    conditions to individual variable components.

    **Initialize the Custom object**

    Parameters
    ----------
    boundaryID : int
        Identifier for the boundary.
    boundaryInfo : dict
        Dictionary with boundary specifications including "specifics".
    boundaryHandler : BoundaryHandler
        Reference to the parent BoundaryHandler.

    Notes
    -----
    Assumes valid keys and values exist under "specifics". Raises may occur if
    structure or content is invalid (to be implemented).
    """

    def __init__(self, boundaryID, boundaryInfo, boundaryHandler):
        super().__init__(boundaryID, boundaryInfo, boundaryHandler)
        self.name = "custom"
        # This boundary condition takes what is given under "specifics". 

        # TODO Sophie: write warning if no specifics are there, and say that everything has been set to "None" (which basically means no boundary conditions) 
        specs = self.info["specifics"]

        for spec in specs:
            # 1. read specs: get variable name, type and value
            # TODO Sophie: catch "KeyError" if spec type or spec variable does not exist or if the value is not a number; give out easy to understand error message
            var    = spec["variable"]
            bcType = BoundaryType[spec["type"].upper()]
            value  = spec["value"]
            if len(spec["variable"])>1:
                comp = spec["variable"][1]
            else:
                comp = ""

            # 2. get index of specific variable
            for v in self.bH.variables:
                if  v[0] == var or v[0] == var[:-1]: #var can be e.g. ux, uy or rhoux, rhouy; var can also be e.g. rho or p; thus both has to be checked 
                    index1 = self.bH.variables.index(v)
                    if len(v[1])>1:
                        index2 = v[1].index(comp)
                    else:
                        index2 = 0

            # 3. set boundary condition
            self.types[index1][index2]  = bcType
            self.values[index1][index2] = ScalarType(np.real(value) + 1j*np.imag(value))



class ZeroDirichlet(BoundaryCondition):
    """
    Homogeneous Dirichlet condition for all variables.

    Sets all variable components to have Dirichlet type with zero value.

    **Initialize the ZeroDirichlet object**

    Parameters
    ----------
    boundaryID : int
        Identifier for the boundary.
    boundaryInfo : dict
        Dictionary with boundary specifications.
    boundaryHandler : BoundaryHandler
        Reference to the parent BoundaryHandler.
    """

    def __init__(self, boundaryID, boundaryInfo, boundaryHandler):
        super().__init__(boundaryID, boundaryInfo, boundaryHandler)
        self.name = "zeroDirichlet"
        # This boundary condition sets a zero dirichlet condition for every variable.

        for i in range(len(self.types)):
            for j in range(len(self.types[i])):
                self.types[i][j]  = BoundaryType.DIRICHLET
                self.values[i][j] = ScalarType(0.+0.j) 



class Wall(BoundaryCondition):
    """
    Wall boundary condition enforcing zero velocity.

    Applies Dirichlet(0) for all components of velocity-type variables
    (e.g., "u" or "rhou"). All other variables are left unconstrained.

    **Initialize the Wall object**

    Parameters
    ----------
    boundaryID : int
        Identifier for the boundary.
    boundaryInfo : dict
        Dictionary with boundary specifications.
    boundaryHandler : BoundaryHandler
        Reference to the parent BoundaryHandler.

    Notes
    -----
    Only handles velocity conditions. Future extensions may include temperature or other constraints.
    """

    def __init__(self, boundaryID, boundaryInfo, boundaryHandler):
        super().__init__(boundaryID, boundaryInfo, boundaryHandler)
        self.name = "wall"
        # This boundary condition, at the moment, sets only the velocity components to zero, 
        # all other variables have no boundary condition ("none").
        # TODO Sophie: add "attribute":  e.g. "adiabatic", "isothermal"
      
        for var in self.bH.variables:
            if var[0] == "u" or var[0] == "rhou":
                index_u = self.bH.variables.index(var)
                for comp in var[1]:
                    index_comp = var[1].index(comp)
                    self.types[index_u][index_comp]  = BoundaryType.DIRICHLET
                    self.values[index_u][index_comp] = ScalarType(0.+0.j) 
                break


class Symmetry(BoundaryCondition):
    """
    Symmetry boundary condition for specified variables.

    Uses the "specifics" field from the JSON file and validates:
    - Only 'dirichlet' or 'neumann' types
    - Only zero-valued conditions
    - All variables are defined

    **Initialize the Symmetry object**

    Parameters
    ----------
    boundaryID : int
        Identifier for the boundary.
    boundaryInfo : dict
        Dictionary with boundary specifications including "specifics".
    boundaryHandler : BoundaryHandler
        Reference to the parent BoundaryHandler.

    Notes
    -----
    Performs stricter validation than the 'Custom' class to ensure
    boundary conditions conform to symmetry constraints.
    """

    def __init__(self, boundaryID, boundaryInfo, boundaryHandler):
        super().__init__(boundaryID, boundaryInfo, boundaryHandler)
        self.name = "symmetry"
        # This boundary condition takes what is given under "specifics". 
        # There can only be the types "dirichlet" or "neumann", and only the value 0.
        # Specifics have to be given for every variable.
        # Difference to boundary condition "Custom": it is checked if all conditions 
        # for a symmetry BC are met. In "Custom", there can appear more general BC combinations.
 
        # TODO Sophie: write error message if no "specifics" are there and stop FELiCS
        # TODO Sophie: also write error message if not all variables are specified
        # TODO Sophie: also write error message if not all types are "dirichlet" or "neumann" or if not all values are "0".
        specs = self.info["specifics"]

        for spec in specs:
            # 1. read specs: get variable name, type and value
            # TODO Sophie: catch "KeyError" if spec type does not exist and give out easy to understand error message
            var    = spec["variable"]
            bcType = BoundaryType[spec["type"].upper()]
            value  = spec["value"]
            if len(spec["variable"])>1:
                comp = spec["variable"][1]
            else:
                comp = ""

            # 2. get index of specific variable
            for v in self.bH.variables:
                if  v[0] == var or v[0] == var[:-1]: #var can be e.g.  ux, uy or rhoux, rhouy; var can also be e.g. rho or p; thus both has to be checked 
                    index1 = self.bH.variables.index(v)
                    if len(v[1])>1:
                        index2 = v[1].index(comp)
                    else:
                        index2 = 0

            # 3. set boundary condition
            self.types[index1][index2]  = bcType
            self.values[index1][index2] = ScalarType(np.real(value) + 1j*np.imag(value))



