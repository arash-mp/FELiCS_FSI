from .Field import Field

class MeanField(Field):

    def __init__(FEMSpace, mesh):
        super().__init(FEMSpace, mesh)


    def importFromFile(self, fileName):
        ##ToDo
        pass
