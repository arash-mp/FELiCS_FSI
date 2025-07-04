import  json
import  numpy as np
from    enum                         import Enum
from    petsc4py.PETSc               import ScalarType
from 	FELiCS.Misc.logging          import Logger
# Get the logger
logger = Logger.get_logger("felics")



class BoundaryHandler():
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
        # TODO Sophie: check if ID is in self.Ids and throw an error message if not (also naming all the ids that are there)
        self.boundaryList = []
        for ID in BCsInfo:
            id_int = int(ID)
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
                # TODO Sophie: write error message if the name is not recognized and stop FELiCS (give list of possible boundary condition names)
                pass
            self.boundaryList.append(bc)


    def getListOfBoundaries(self):
        return self.boundaryList


    def getListOfDirichletBCsForDolfinx(self,functionSpace):
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
        #TODO Sophie: fill out later for base flow computations
        pass


    def getListOfNonlinearDirichletBCsForDolfinx(self):
        #TODO Sophie: fill out later for base flow computations
        pass


class BoundaryType(Enum):
    NONE      = 0
    DIRICHLET = 1
    NEUMANN   = 2 # equal to "none" at the moment; this should be changed in the future, and also communicated really well
    #MIXED     = 3


class BoundaryCondition():
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
    def __init__(self, boundaryID, boundaryInfo, boundaryHandler):
        super().__init__(boundaryID, boundaryInfo, boundaryHandler)
        self.name = "zeroDirichlet"
        # This boundary condition sets a zero dirichlet condition for every variable.

        for i in range(len(self.types)):
            for j in range(len(self.types[i])):
                self.types[i][j]  = BoundaryType.DIRICHLET
                self.values[i][j] = ScalarType(0.+0.j) 



class Wall(BoundaryCondition):
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



