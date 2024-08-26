from dolfinx.fem import Function, petsc


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



    def setListOfSingleFields(self, listOfFields):
        numberOfSubSpaces = self.space.num_sub_spaces

        if numberOfSubSpaces == 0:
            try:
                self.setCoefficientArray(listOfFields[0].getCoefficientArray())
            except:
                print("ERROR")
                #TODO: Throw error!
            return 

        for i in range(numberOfSubSpaces):
            try:
                space, mapping                 = self.space.sub(i).collapse()
                self.function.x.array[mapping] = listOfFields[i].getCoefficientArray()
            except:
                print("ERROR")
                #TODO: Throw error!

        return 





    def getCoefficientArray(self):
        import numpy as np
        array = np.empty(len(self.function.x.array[:]),dtype=complex)
        array[:] = self.function.x.array[:]
        return array

    def setCoefficientArray(self, array): 
        self.function.x.array[:] = array[:]

    def setConstantValue(self, value): 
        self.function.x.array[:] = value

    def getPetscVector(self):
        from petsc4py import PETSc
        return PETSc.Vec().createWithArray(self.getCoefficientArray())

    def conjugate(self):
        import numpy as np
        self.setCoefficientArray(np.conj(self.getCoefficientArray()))

    def setBoundaryConditions(self, bcs):
        petscArray = self.getPetscVector()
        petsc.set_bc(petscArray,bcs)
        self.setCoefficientArray(petscArray.getArray())

