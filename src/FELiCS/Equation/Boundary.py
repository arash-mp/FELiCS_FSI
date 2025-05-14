class Boundary():

    def __init__(self, boundaryID):
        self.boundaryID = boundaryID
        self.listOfVariables = []
        self.name =""

        self.ds =0.
        self.n = 0.
        self.types = []
        self.values = []


    def appendDirichletBC(self):
        pass


    def isNeumann(self, variable):
        return False


    def setBCs_matrix(self, petscMatrix):
        pass

    def setBCs_vector(self, petscVector):
        pass




