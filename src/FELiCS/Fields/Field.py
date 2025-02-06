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


    def evaluateUflExpression(self, ufl_expression, bcs=[], restartSolver=False):
        # evaluates an ufl expression by 
        from FELiCS.Solvers.LinearSolver import LinearSolver
        from FELiCS.Equation.WeakForm import WeakForm
        import ufl 
        import dolfinx
        ## create petsc solver and save it as attribute to the corresponding space - to use the LU-decomposition later 
        if not hasattr(self.space, 'FEMWeightSolver') and not restartSolver:
            test_FEM   = ufl.TestFunctions(self.space)
            trial_FEM  = ufl.TrialFunctions(self.space)
            matrix_ufl = WeakForm()
            i=0
            for test in test_FEM:
                try:
                    j=0
                    for subTest in test:
                        matrix_ufl.add(ufl.conj(subTest)*trial_FEM[i][j]*ufl.dx)
                        j+=1
                except:
                    matrix_ufl.add(ufl.conj(test)*trial_FEM[i]*ufl.dx)
                i+=1
            try:
                matrix_ufl.setCorrectMeshObject(self.mesh)
            except:
                pass
            matrix = petsc.assemble_matrix(dolfinx.fem.form(matrix_ufl.lhs), bcs=bcs)
            matrix.assemble()
            self.space.FEMWeightSolver = LinearSolver.createEquationSystemSolver(matrix)

        ## assemble rhs and solve equation system
        expr_ufl = WeakForm(ufl_expression)
        try:
            expr_ufl.setCorrectMeshObject(self.mesh)
        except:
            pass
        petscVec = petsc.assemble_vector(dolfinx.fem.form(-expr_ufl.rhs))
        petscVec.assemble()
        petsc.set_bc(petscVec, bcs)
        self.setCoefficientArray(LinearSolver.solveEquationSystemWithPredefinedSolver(self.space.FEMWeightSolver, petscVec))
        self.setBoundaryConditions(bcs)


    def evaluateUflTensorExpression(self, ufl_expression, bcs=[], restartSolver=False):
        # evaluates an ufl expression by 
        from FELiCS.Solvers.LinearSolver import LinearSolver
        from FELiCS.Equation.WeakForm import WeakForm
        import ufl 
        import dolfinx

        from   FELiCS.Misc.tensorUtils import (
            Tensor,
            iDot,
            iConj,
        )
        ## create petsc solver and save it as attribute to the corresponding space - to use the LU-decomposition later 
        if not hasattr(self.space, 'FEMWeightSolver') and not restartSolver:
            test_FEM   = ufl.TestFunctions(self.space)
            trial_FEM  = ufl.TrialFunctions(self.space)
            matrix_ufl = WeakForm()
            coordinateSystem = self.mesh.coordinateSystem
            J_hat = coordinateSystem.J_hat
            i=0
            for test in test_FEM:
                iTest = Tensor(test, coordinateSystem, containsTestFunction=True)
                iFluc = Tensor(trial_FEM[i], coordinateSystem, containsFluctuation=True)
                if iTest.order  == 1:
                    matrix_ufl.add( ( iDot(iFluc, iConj(iTest)) ).ufl_tens*J_hat*ufl.dx)
                elif iTest.order == 0:
                    matrix_ufl.add( ( iFluc * iConj(iTest)).ufl_tens*J_hat*ufl.dx)
                else:
                    ##LOGGING TODO (Sophie): throw error 
                    print("ERROR, in 'Field.evaluateUflExpression'")
                i+=1
            try:
                matrix_ufl.setCorrectMeshObject(self.mesh)
            except:
                pass
            matrix = petsc.assemble_matrix(dolfinx.fem.form(matrix_ufl.lhs), bcs=bcs)
            matrix.assemble()
            self.space.FEMWeightSolver = LinearSolver.createEquationSystemSolver(matrix)

        ## assemble rhs and solve equation system
        expr_ufl = WeakForm(ufl_expression)
        try:
            expr_ufl.setCorrectMeshObject(self.mesh)
        except:
            pass
        petscVec = petsc.assemble_vector(dolfinx.fem.form(-expr_ufl.rhs))
        petscVec.assemble()
        petsc.set_bc(petscVec, bcs)
        self.setCoefficientArray(LinearSolver.solveEquationSystemWithPredefinedSolver(self.space.FEMWeightSolver, petscVec))
        self.setBoundaryConditions(bcs)


    def smoothUflTensorExpression(self, ufl_expression, smoothFactor, bcs=[], restartSolver=False):
        # evaluates an ufl expression by 
        from FELiCS.Solvers.LinearSolver import LinearSolver
        from FELiCS.Equation.WeakForm import WeakForm
        import ufl 
        import dolfinx

        from   FELiCS.Misc.tensorUtils import (
            Tensor,
            iDot,
            iConj,
        )
        ## create petsc solver and save it as attribute to the corresponding space - to use the LU-decomposition later 
        if not hasattr(self.space, 'FEMSmoothSolver') and not restartSolver:
            test_FEM   = ufl.TestFunctions(self.space)
            trial_FEM  = ufl.TrialFunctions(self.space)
            matrix_ufl = WeakForm()
            coordinateSystem = self.mesh.coordinateSystem
            J_hat = coordinateSystem.J_hat
            i=0
            for test in test_FEM:
                iTest = Tensor(test, coordinateSystem, containsTestFunction=True)
                iFluc = Tensor(trial_FEM[i], coordinateSystem, containsFluctuation=True)
                if iTest.order  == 1:
                    matrix_ufl.add( ( iDot(iFluc, iConj(iTest)) ).ufl_tens*J_hat*ufl.dx)
                    matrix_ufl.add((smoothFactor* iInner(iGrad(iFluc),iGrad(iConj(iTest)))).ufl_tens*J_hat*dx)
                elif iTest.order == 0:
                    matrix_ufl.add( ( iFluc * iConj(iTest)).ufl_tens*J_hat*ufl.dx)
                    matrix_ufl.add((smoothFactor* iDot(iGrad(iFluc),iGrad(iConj(iTest)) )).ufl_tens*J_hat*dx)
                else:
                    ##LOGGING TODO (Sophie): throw error 
                    print("ERROR, in 'Field.smoothTensorUflExpression'")
                i+=1
            try:
                matrix_ufl.setCorrectMeshObject(self.mesh)
            except:
                pass
            matrix = petsc.assemble_matrix(dolfinx.fem.form(matrix_ufl.lhs), bcs=bcs)
            matrix.assemble()
            self.space.FEMSmoothSolver = LinearSolver.createEquationSystemSolver(matrix)

        ## assemble rhs and solve equation system
        expr_ufl = WeakForm(ufl_expression)
        try:
            expr_ufl.setCorrectMeshObject(self.mesh)
        except:
            pass
        petscVec = petsc.assemble_vector(dolfinx.fem.form(-expr_ufl.rhs))
        petscVec.assemble()
        petsc.set_bc(petscVec, bcs)
        self.setCoefficientArray(LinearSolver.solveEquationSystemWithPredefinedSolver(self.space.FEMSmoothSolver, petscVec))
        self.setBoundaryConditions(bcs)





    ### dunder methods for overloading arithmetic operators ###
    def __add__(self, other):
        ## overrides '+'
        ## returns newly created Field with a coefficient array, which is the sum of two given coefficientarrays
        if isinstance(other, Field):
            result = Field(self.space, self.mesh)
            result.setCoefficientArray(self.getCoefficientArray() + other.getCoefficientArray())
            return result 
        return NotImplemented



