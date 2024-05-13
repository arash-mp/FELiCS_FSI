from dolfinx import Function


class Field:


    def __init__(self, FEMSpace, mesh):
        self.space = FEMSpace
        self.mesh  = mesh

        self.function = Function(FEMSpace)



    def getListOfSingleFunctions(self):#
        ##ToDo: return single functions if e.g. mixed Function
        pass


