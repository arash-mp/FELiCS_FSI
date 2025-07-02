import  json
import  numpy as np
from    enum                         import Enum
from 	FELiCS.Misc.logging          import Logger

# Get the logger
logger = Logger.get_logger("felics")



class BoundaryHandler():
    def __init__(self, variables, mesh, BCsFilePath):

        ###### initialize list of boundaries ########

        #1. get info from felics mesh and store variables
        self.facet_tags             = mesh.facet_tags
        self.Ids                    = np.unique(self.facet_tags.values)
        self.variables              = variables

        #2. read bc file 
        logger.debug(f"Reading boundary conditions from '{BCsFilePath}'")
        file    = open(BCsFilePath) 
        BCsInfo = json.load(file)
        
        #3. create boundary object for each boundary
        self.boundaryList = []
        for ID in BCsInfo:
            info = BCsInfo[ID]
            if   info["name"]   == "custom":
                bc = Custom(ID, BCsInfo[ID], self)
            elif info["name"]   == "zeroDirichlet":
                bc = ZeroDirichlet(ID, BCsInfo[ID], self )
            elif info["name"]   == "wall":
                bc = Wall(ID, BCsInfo[ID], self)
            elif info["name"]   == "symmetry":
                bc = Symmetry(ID, BCsInfo[ID], self)
            elif info["name"]   != "none":
                # TODO Sophie: write error message if the name is not recognized and stop FELiCS (give list of possible boundary condition names)
                pass
            self.boundaryList.append(bc)

    def getListOfBoundaries(self):
        return self.boundaryList

    def getListOfDirichletBCsForDolfinx(self,functionSpace):
        from dolfinx.fem      import dirichletbc, locate_dofs_topological
        from petsc4py.PETSc   import ScalarType
        BCs = []
        #TODO: give a good description of what is done here:
        for boundary in self.boundaryList:
            for var in self.variables:
                index = self.variables.index(var)
                if len(var[1])==1:
                    value = boundary.values[index][0]
                    BCs.append(dirichletbc(ScalarType(value), 
                                     locate_dofs_topological(functionSpace.sub(index), 1, self.facet_tags.indices[self.facet_tags.values==boundary.ID]), 
                                     functionSpace.sub(index)))
                    logger.debug("Adding Dirichlet BC for "+var[0]+ " in equation "+str(index)+" with value "+str(value)+" on boundary with index "+str(boundary.ID))
                else: 
                    for index2 in range(len(var[1])):
                        value = boundary.values[index][index2]
                        BCs.append(dirichletbc(ScalarType(value), 
                                         locate_dofs_topological(functionSpace.sub(index).sub(index2), 1, self.facet_tags.indices[self.facet_tags.values==boundary.ID]), 
                                         functionSpace.sub(index).sub(index2)))
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
            var    = spec["variable"][0]
            # TODO Sophie: catch "KeyError" if spec type does not exist and give out easy to understand error message
            bcType = BoundaryType[spec["type"].upper()]
            value  = spec["value"]
            if len(spec["variable"])>1:
                comp = spec["variable"][1]
            else:
                comp = ""
            for v in self.bH.variables:
                if v[0] == var:
                    index1 = self.bH.variables.index(v)
                    if len(v[1])>1:
                        index2 = v[1].index(comp)
                    else:
                        index2 = 0
            self.types[index1][index2]  = bcType
            self.values[index1][index2] = value



class ZeroDirichlet(BoundaryCondition):
    def __init__(self, boundaryID, boundaryInfo, boundaryHandler):
        super().__init__(boundaryID, boundaryInfo, boundaryHandler)
        self.name = "zeroDirichlet"
        # This boundary condition sets a zero dirichlet condition for every variable.

        for i in range(len(self.types)):
            for j in range(len(self.types[i])):
                self.types[i][j]  = BoundaryType.DIRICHLET
                self.values[i][j] = 0.


class Wall(BoundaryCondition):
    def __init__(self, boundaryID, boundaryInfo, boundaryHandler):
        super().__init__(boundaryID, boundaryInfo, boundaryHandler)
        self.name = "wall"
        # This boundary condition, at the moment, sets only the velocity components to zero. 
        # TODO: add "attribute":  e.g. "adiabatic", "isothermal"
      
        for var in self.bH.variables:
            if var[0] == "u" or var[0] == "rhou":
                index_u = self.bH.variables.index(var)
                for comp in var[1]:
                    index_comp = var[1].index(comp)
                    self.types[index_u][index_comp]  = BoundaryType.DIRICHLET
                    self.values[index_u][index_comp] = 0. 
                break


class Symmetry(BoundaryCondition):
    def __init__(self, boundaryID, boundaryInfo, boundaryHandler):
        super().__init__(boundaryID, boundaryInfo, boundaryHandler)
        self.name = "symmetry"
        # This boundary condition takes what is given under "specifics". 
        # There can only be the types "dirichlet" or "neumann", and only the value 0.
        # Specifics have to be given for every variable.
 
        # TODO Sophie: write error message if no "specifics" are there and stop FELiCS
        # TODO Sophie: also write error message if not all variables are specified
        # TODO Sophie: also write error message if not all types are "dirichlet" or "neumann" or if not all values are "0".
        specs = self.info["specifics"]

        for spec in specs:
            var    = spec["variable"][0]
            # TODO Sophie: catch "KeyError" if spec type does not exist and give out easy to understand error message
            bcType = BoundaryType[spec["type"].upper()]
            value  = spec["value"]
            if len(spec["variable"])>1:
                comp = spec["variable"][1]
            else:
                comp = ""
            for v in self.bH.variables:
                if v[0] == var:
                    index1 = self.bH.variables.index(v)
                    if len(v[1])>1:
                        index2 = v[1].index(comp)
                    else:
                        index2 = 0
            self.types[index1][index2]  = bcType
            self.values[index1][index2] = value




