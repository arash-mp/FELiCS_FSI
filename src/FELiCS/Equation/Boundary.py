import  json
from 	FELiCS.Misc.logging          import Logger

# Get the logger
logger = Logger.get_logger("felics")



class BoundaryHandler():
    def __init__(self, mesh, BCsFilePath):
        ###### initialize list of boundaries ########
        #1. get info from felics mesh
        self.Ids, self.boundaries = mesh.getBCInfo()
        #2. read bc file 
        logger.debug(f"Reading boundary conditions from '{BCsFilePath}'")
        file    = open(BCsFilePath) 
        BCsInfo = json.load(file)
        #3. create boundary object for each boundary
        boundaryList = []
        for ID in BCsInfo:
            info = BCsInfo[ID]
            if   info["name"]   == "custom":
                bc = Custom(ID)
            elif info["name"]   == "zeroDirichlet":
                bc = ZeroDirichlet(ID)
            elif info["name"]   == "wall":
                bc = Wall(ID)
            elif info["name"]   == "symmetry":
                bc = Symmetry(ID)
                 
            boundaryList.append(bc)


    def getListOfBoundaries():
        pass

    def getListOfDirichletBCsForDolfinx():
        pass


    def getListOfNonlinearBoundaries():
        #TODO Sophie: fill out later for base flow computations
        pass


    def getListOfNonlinearDirichletBCsForDolfinx():
        #TODO Sophie: fill out later for base flow computations
        pass



class BoundaryCondition():
    def __init__(self, boundaryID):
        self.boundaryID = boundaryID
        self.listOfVariables = []
        self.name =""

        self.ds =0.
        self.n = 0.
        self.types = []
        self.values = []

    def appendDirichletBCsForDolfinx(self, listOfDirichletBCs):
        pass

    def isNeumann(self, variable):
        return False



class Custom(BoundaryCondition):
    def __init__(self, boundaryID):
        self.name = "custom"


class ZeroDirichlet(BoundaryCondition):
    def __init__(self, boundaryID):
        self.name = "zeroDirichlet"


class Wall(BoundaryCondition):
    def __init__(self, boundaryID):
        self.name = "wall"


class Symmetry(BoundaryCondition):
    def __init__(self, boundaryID):
        self.name = "symmetry"



