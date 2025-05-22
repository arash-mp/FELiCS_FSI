from dolfinx.fem import Function, petsc


class Field:
    """
    Class representing a finite element field.

    This class provides methods for handling finite element fields, including 
    coefficient manipulation, boundary conditions application, and expression 
    evaluation. It supports complex-valued fields and operates within a 
    tensorial framework.

    The class interacts with `dolfinx.fem.Function` for finite element operations 
    and includes utilities for working with `PETSc` vectors and UFL expressions.

    **Initialize the Field object**

    Parameters
    ----------
    FEMSpace : dolfinx.fem.FunctionSpace
        The finite element function space.
    mesh : dolfinx.mesh.Mesh
        The computational mesh associated with the function space.
    """

    def __init__(self, FEMSpace, mesh):
        self.space = FEMSpace
        self.mesh  = mesh

        self.function = Function(FEMSpace)



    def getListOfSingleFields(self):
        """
        Get a list of single-component fields.

        Returns
        -------
        list
            A list containing individual scalar fields if the function space has 
            multiple subspaces. Otherwise, returns a list containing only this field.

        Notes
        -----
        - If the space has multiple subspaces, each subspace is extracted as an 
          individual field.
        """
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
        """
        Set the field coefficients from a list of single-component fields.

        Parameters
        ----------
        listOfFields : list
            A list of `Field` objects representing individual subspaces.

        Notes
        -----
        - If the space has multiple subspaces, their coefficients are mapped back.
        - Throws an error if the input list does not match the expected size.
        """
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
        """
        Get the coefficient array of the field.

        Returns
        -------
        numpy.ndarray
            A complex-valued array representing the field coefficients.
        """
        import numpy as np
        array = np.empty(len(self.function.x.array[:]),dtype=complex)
        array[:] = self.function.x.array[:]
        return array

    def setCoefficientArray(self, array):
        """
        Set the coefficient array of the field.

        Parameters
        ----------
        array : numpy.ndarray
            A complex-valued array of coefficients.
        """
        self.function.x.array[:] = array[:]

    def setConstantValue(self, value): 
        """
        Set all coefficients to a constant value.

        Parameters
        ----------
        value : complex or float
            The constant value to assign to all coefficients.
        """
        self.function.x.array[:] = value

    def getPetscVector(self):
        """
        Convert the field to a PETSc vector.

        Returns
        -------
        petsc4py.PETSc.Vec
            A PETSc vector created from the coefficient array.
        """
        from petsc4py import PETSc
        return PETSc.Vec().createWithArray(self.getCoefficientArray())

    def conjugate(self):
        """
        Compute the complex conjugate of the field.
        """
        import numpy as np
        self.setCoefficientArray(np.conj(self.getCoefficientArray()))

    def setBoundaryConditions(self, bcs):
        """
        Apply boundary conditions to the field.

        Parameters
        ----------
        bcs : list
            A list of Dirichlet boundary conditions to be applied.
        """
        petscArray = self.getPetscVector()
        petsc.set_bc(petscArray,bcs)
        self.setCoefficientArray(petscArray.getArray())


    def evaluateUflExpression(self, ufl_expression, bcs=[], restartSolver=False):
        """
        Evaluate a UFL expression and update the field accordingly.

        Parameters
        ----------
        ufl_expression : ufl.Form
            The UFL expression to evaluate.
        bcs : list, optional
            A list of boundary conditions to apply.
        restartSolver : bool, optional
            Whether to restart the solver instead of reusing an existing one.
        """
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
        expr_ufl = WeakForm()
        expr_ufl.add(ufl_expression)
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
        """
        Evaluate a tensor-based UFL expression and update the field accordingly.

        Parameters
        ----------
        ufl_expression : ufl.Form
            The UFL expression to be evaluated.
        bcs : list, optional
            A list of boundary conditions to be applied.
        restartSolver : bool, optional
            Whether to restart the solver instead of reusing an existing one.

        Notes
        -----
        - This method assembles and solves a tensor-based weak form.
        - Uses a predefined solver if available to improve performance.
        - The weak form includes integration over the computational domain.
        """
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
        expr_ufl = WeakForm()
        expr_ufl.add(ufl_expression)
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
        """
        Smooth a tensor-based UFL expression using a diffusion-like approach.

        Parameters
        ----------
        ufl_expression : ufl.Form
            The UFL expression to be smoothed.
        smoothFactor : float
            A smoothing factor controlling the influence of the gradient term.
        bcs : list, optional
            A list of boundary conditions to apply.
        restartSolver : bool, optional
            Whether to restart the solver instead of reusing an existing one.

        Notes
        -----
        - This method applies a smoothing operation by adding a gradient-based 
          regularization term to the weak form.
        - It is particularly useful for regularizing noisy numerical solutions.
        """
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
        """
        Overload the `+` operator for adding two Field objects.

        Parameters
        ----------
        other : Field
            Another Field object.

        Returns
        -------
        Field
            A new Field object with the summed coefficient arrays.

        Raises
        ------
        NotImplementedError
            If `other` is not a Field object.
        """
        ## overrides '+'
        ## returns newly created Field with a coefficient array, which is the sum of two given coefficientarrays
        if isinstance(other, Field):
            result = Field(self.space, self.mesh)
            result.setCoefficientArray(self.getCoefficientArray() + other.getCoefficientArray())
            return result 
        return NotImplemented



