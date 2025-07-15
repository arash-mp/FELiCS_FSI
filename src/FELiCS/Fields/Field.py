from dolfinx.fem             import Function, petsc

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
    mesh:      FELiCS.SpaceDisc.FELiCSMesh
        FEliCS mesh associated with the function space.
    name:      optional, list of tuples; e.g. [("name",[])] for a scalar space, or [("name",["x","y","z"]] for a vector space, a list of both types for a mixed space
    m:         optional, integer
        Wave number. If this is set, the Field is assumed to have one spectral spatial dimension, regardless the value of m.
    """

    def __init__(self, FEMSpace, mesh, name=[], m=None):
        self.space    = FEMSpace
        self.mesh     = mesh

        # if the name is in the wrong format, re-format
        # TODO Sophie: reformat if there is only a list of variables for a mixed space?
        # TODO Sophie: throw warning / error if the number of names does not coincide with the number of subspaces
        if isinstance (name, list) and len(name)>0 and isinstance(name[0],tuple):
            self.name  = name
        elif isinstance(name, str):
            self.name  = [(name,[])]
        else:
            self.name = name

        # handle spectral dimension and wave number
        if m != None:
            self.hasSpectralDimension = True
            self.m = m
        else:
            self.hasSpectralDimension = False
            self.m = 0

        # initialize function
        self.function = Function(FEMSpace)

    def getName(self):
        if len(self.name) == 0:
            return ""
        elif len(self.name) == 1 and len(self.name[0][1])==0:
            return self.name[0][0]
        else:
            return self.name

    def getComponentsNames(self):
        ## This is  a workaround for now, to use for the retreat.
        ## TODO Sophie: make this independent of the coordinate system, and also usable for mixed function spaces.
        numberOfSubSpaces = self.space.num_sub_spaces
        if numberOfSubSpaces == 0:
            return []
        elif self.mesh.coordinateSystemName == "Cartesian" and numberOfSubSpaces == 2:
            return ["x","y"]
        elif self.mesh.coordinateSystemName == "Cartesian" and numberOfSubSpaces == 3:
            return ["x","y","z"]
        elif self.mesh.coordinateSystemName == "Cylindrical" and numberOfSubSpaces == 2:
            return ["x","r"]
        elif self.mesh.coordinateSystemName == "Cylindrical" and numberOfSubSpaces == 3:
            return ["x","r","t"]

    def getTensor(self):
        from FELiCS.Misc.tensorUtils import Tensor
        # TODO Sophie: handle Tensors of mixed functions (later)
        return Tensor(self.function, self.mesh.coordinateSystem, m = self.m, hasSpectralDimension = self.hasSpectralDimension)

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
            # if len(self.name) == numberOfSubSpaces:  
            #     field.name = [self.name[i]]
            # elif len(self.name)==1 and len(self.name[0]) == numberOfSubSpaces:
                # field.name = [(self.name[0][0] + self.name[0][1][i], [])]
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
                self.name = listOfFields[0].name
            except:
                print("ERROR")
                #TODO: Throw error!
            return 

        self.name  = [None]*numberOfSubSpaces
        for i in range(numberOfSubSpaces):
            try:
                space, mapping                 = self.space.sub(i).collapse()
                self.function.x.array[mapping] = listOfFields[i].getCoefficientArray()
                # TODO Sophie: how to recognize velocity names? ux, uy, uz => [(u, [x,y,z]) ?
                if len(listOfFields[i].name) > 0: #and len(listOfFields[i].name[0][1] == 0:
                    self.name[i]               = listOfFields[i].name[0] 
                      
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
        import numpy as np
        if np.isscalar(array):
            self.function.x.array[:] = array
        else:
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

        
    def getGradientField(self):
        # TODO Sophie: throw error if Field is not scalar    
        from ufl import TestFunction, dx
        from FELiCS.Misc.tensorUtils import iGrad, iConj, iDot, Tensor
        from FELiCS.SpaceDisc.FEMSpaces import getFELiCSSpace
        # Create  a Field for the gradient
        # The space must be a vector vor a scalar field
        # TODO Sophie: handle order (get it from function?)
        order = 2
        dim = self.mesh.gdim 
        if self.hasSpectralDimension:
            dim += 1
        gradientSpace = getFELiCSSpace(self.mesh, order=order, dim=dim)
        gradientField = Field(gradientSpace, self.mesh)

        coordSystem = self.mesh.coordinateSystem
        J_hat = coordSystem.J_hat

        v = TestFunction(gradientSpace)
        # Sophie: m and hasSpectralDimension may not be needed for test function
        v_tens = Tensor(v, CoordSys=coordSystem, m = self.m, hasSpectralDimension=self.hasSpectralDimension)
        expression = iDot(iGrad(self.getTensor()), iConj(v_tens)).ufl_tens * J_hat* dx

        gradientField.evaluateUflTensorExpression(expression)    

        # TODO Sophie: check how to destroy all petsc objects, also matrix
        gradientSpace.FEMWeightSolver.destroy()

        return gradientField

    def exportH5(self, fileName, mesh = None):
        # for now this is a dummy method that we use for the scripting part of the retreat.
        # fileName: WITHOUT SUFFIX, but WITH PATH
        # mesh: FELiCSMesh
        # What it should be: export the field in h5 format, using the fileName (which should contain the whole path)
        # There should also be an optional possibility to give the mesh, for scripting (or use it inside FELiCS as such?)
        from dolfinx.io import XDMFFile
        from mpi4py import MPI
        import numpy as np

        if mesh != None:
            with XDMFFile(MPI.COMM_WORLD, fileName+".xdmf", "w") as xdmf:
                xdmf.write_mesh(mesh.dolfinxMesh)
                xdmf.write_function(self._function)

        # this is only a dummy for the scripting
        np.save(fileName+".npy", self.getCoefficientArray())


    def importH5(self, fileName, meshFileName = None):
        # for now this is a dummy method that we use for the scripting part of the retreat.
        # fileName: WITHOUT SUFFIX, but WITH PATH
        # mesh: FELiCSMesh
        # What it should be: import the field in h5 format, using the fileName (which should contain the whole path)
        # There should also be an optional possibility to give the mesh, for scripting (or use it inside FELiCS as such?)
        from dolfinx.io import XDMFFile
        from mpi4py import MPI
        import numpy as np

        #if meshFileName != None:
        #    with XDMFFile(MPI.COMM_WORLD, fileName+".xdmf", "r") as xdmf:
        #        mesh          = xdmf.read_mesh(meshFileName)

        # this is only a dummy for the scripting
        self.setCoefficientArray(np.load(fileName+".npy"))


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
        from FELiCS.Equation.UflDecorator import UflDecorator
        import ufl
        import dolfinx
        ## create petsc solver and save it as attribute to the corresponding space - to use the LU-decomposition later 
        if not hasattr(self.space, 'FEMWeightSolver') and not restartSolver:
            test_FEM   = ufl.TestFunctions(self.space)
            trial_FEM  = ufl.TrialFunctions(self.space)
            matrix_ufl = UflDecorator()
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
        
            # matrix = petsc.assemble_matrix(dolfinx.fem.form(matrix_ufl.lhs), bcs=bcs)
            # matrix.assemble()
            matrix = matrix_ufl.getAssembledMatrix(self.mesh, bcs)
            self.space.FEMWeightSolver = LinearSolver.createEquationSystemSolver(matrix)

        ## assemble rhs and solve equation system
        expr_ufl = UflDecorator(ufl_expression)
        petscVec = expr_ufl.getAssembledVector(self.mesh, bcs)
        self.setCoefficientArray(LinearSolver.solveEquationSystemWithPredefinedSolver(self.space.FEMWeightSolver, petscVec))


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
        from FELiCS.Equation.UflDecorator import UflDecorator
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
            matrix_ufl = UflDecorator()
            coordinateSystem = self.mesh.coordinateSystem
            J_hat = coordinateSystem.J_hat
            i=0
            for test in test_FEM:
                iTest = Tensor(test, coordinateSystem, hasSpectralDirection=True)
                iFluc = Tensor(trial_FEM[i], coordinateSystem, hasSpectralDirection=True)
                if iTest.order  == 1:
                    matrix_ufl.add( ( iDot(iFluc, iConj(iTest)) ).ufl_tens*J_hat*ufl.dx)
                elif iTest.order == 0:
                    matrix_ufl.add( ( iFluc * iConj(iTest)).ufl_tens*J_hat*ufl.dx)
                else:
                    ##LOGGING TODO (Sophie): throw error 
                    print("ERROR, in 'Field.evaluateUflExpression'")
                i+=1
            matrix = matrix_ufl.getAssembledMatrix(self.mesh, bcs) 
            self.space.FEMWeightSolver = LinearSolver.createEquationSystemSolver(matrix)

        ## assemble rhs and solve equation system
        expr_ufl = UflDecorator(ufl_expression)
        petscVec = expr_ufl.getAssembledVector(self.mesh, bcs)
        self.setCoefficientArray(LinearSolver.solveEquationSystemWithPredefinedSolver(self.space.FEMWeightSolver, petscVec))


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
        from FELiCS.Equation.UflDecorator import UflDecorator
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
            matrix_ufl = UflDecorator()
            coordinateSystem = self.mesh.coordinateSystem
            J_hat = coordinateSystem.J_hat
            i=0
            for test in test_FEM:
                iTest = Tensor(test, coordinateSystem, hasSpectralDirection=True)
                iFluc = Tensor(trial_FEM[i], coordinateSystem, hasSpectralDirection=True)
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
        expr_ufl = UflDecorator(ufl_expression)
        petscVec = expr_ufl.getAssembledVector(self.mesh, bsc)
        self.setCoefficientArray(LinearSolver.solveEquationSystemWithPredefinedSolver(self.space.FEMWeightSolver, petscVec))



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
        # TODO Sophie: raise error / not implemented if fields are not defined on the same space
        if isinstance(other, Field):
            result = Field(self.space, self.mesh)
            result.setCoefficientArray(self.getCoefficientArray() + other.getCoefficientArray())
            return result 
        return NotImplemented

    def __mul__(self, other):
        import numpy as np
        ## overrides '*'
        ## returns newly created Field with a coefficient array, which is the product of two given coefficientarrays, or the product of its coefficientarray with a scalar value
        # TODO Sophie: raise error / not implemented if fields are not defined on the same space
        if isinstance(other, Field):
            result = Field(self.space, self.mesh)
            result.setCoefficientArray(self.getCoefficientArray()*other.getCoefficientArray())
            return result
        elif np.isscalar(other):
            result = Field(self.space, self.mesh)
            result.setCoefficientArray(self.getCoefficientArray()*other)
            return result
        return NotImplemented

    def __truediv__(self, other):
        import numpy as np
        ## overrides '*'
        ## returns newly created Field with a coefficient array, which is the division of two given coefficientarrays, or the division of its coefficientarray with a scalar value
        # TODO Sophie: raise error / not implemented if fields are not defined on the same space
        if isinstance(other, Field):
            result = Field(self.space, self.mesh)
            result.setCoefficientArray(self.getCoefficientArray()/other.getCoefficientArray())
            return result
        elif np.isscalar(other):
            result = Field(self.space, self.mesh)
            result.setCoefficientArray(self.getCoefficientArray()/other)
            return result
        return NotImplemented

    def describeFunctionSpace(self):
        """
        Describe the structure of the function space and its subspaces.

        Returns
        -------
        dict
            A dictionary containing:
            - 'type': str - 'scalar', 'vector', or 'mixed'
            - 'num_subspaces': int - Number of subspaces
            - 'subspaces': list - List of subspace descriptions for mixed spaces
            - 'value_size': int - Total number of components
            - 'description': str - Human-readable description

        Examples
        --------
        For a scalar field:
        {'type': 'scalar', 'num_subspaces': 0, 'value_size': 1, 'description': 'Single scalar space'}

        For a vector field:
        {'type': 'vector', 'num_subspaces': 3, 'value_size': 3, 'description': 'Single vector space with 3 components'}

        For a mixed field:
        {'type': 'mixed', 'num_subspaces': 2, 'subspaces': [...], 'value_size': 4, 'description': 'Mixed space with 2 subspaces'}
        """
        num_subspaces = self.space.num_sub_spaces
        value_size = self.space.value_size
        
        result = {
            'num_subspaces': num_subspaces,
            'value_size': value_size,
            'subspaces': []
        }
        
        if num_subspaces == 0:
            # Single space (scalar or vector)
            result['type']              = 'scalar'
            result['description']       = 'Single scalar space'
            
        else:
            if num_subspaces == value_size:
                # Single vector space
                result['type']          = 'vector'
                result['description']   = f'Single vector space with {value_size} components'
            
            else:
                result['type']          = 'mixed'
                subspace_descriptions   = []
                
                for i in range(num_subspaces):
                    subspace = self.space.sub(i)
                    sub_value_size = subspace.value_size
                    sub_num_subspaces = subspace.num_sub_spaces
                    
                    if sub_value_size == 1:
                        subspace_type = 'scalar'
                        subspace_desc = f'Subspace {i}: scalar'
                    else:
                        subspace_type = 'vector'
                        subspace_desc = f'Subspace {i}: vector ({sub_value_size} components)'
                    
                    subspace_info = {
                        'index': i,
                        'type': subspace_type,
                        'value_size': sub_value_size,
                        'num_sub_subspaces': sub_num_subspaces,
                        'description': subspace_desc
                    }
                    
                    subspace_descriptions.append(subspace_info)
            
                result['subspaces']     = subspace_descriptions
                result['description']   = f'Mixed space with {num_subspaces} subspaces: ' + \
                                  ', '.join([sub['description'].split(': ')[1] for sub in subspace_descriptions])
        
        return result



