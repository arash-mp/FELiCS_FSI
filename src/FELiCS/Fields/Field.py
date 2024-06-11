from dolfinx.fem import Function


class Field:


    def __init__(self, FEMSpace, mesh):
        self.space = FEMSpace
        self.mesh  = mesh

        self.function = Function(FEMSpace)



    def getListOfSingleFields(self):
        listOfFields = []

        numberOfSubSpaces = self.space.num_sub_spaces

        if numberOfSubSpaces == 0: 
            listOfFields.append(self)
            return listOfFields

        for i in range(numberOfSubSpaces):
            space, mapping            = self.space.sub(i).collapse()
            field                     = Field(space, self.mesh)
            field.function.x.array[:] = self.function.x.array[mapping]
            listOfFields.append(field)

        return listOfFields


    def getCoefficientArray(self):
        import numpy as np
        array = np.empty(len(self.function.x.array[:]),dtype=complex)
        array[:] = self.function.x.array[:]
        return array

    def setCoefficientArray(self, array): 
        self.function.x.array[:] = array[:]

    def getPetscVector(self):
        from petsc4py import PETSc
        return PETSc.Vec().createWithArray(self.getCoefficientArray())

    def conjugate(self):
        import numpy as np
        self.setCoefficientArray(np.conj(self.getCoefficientArray()))

