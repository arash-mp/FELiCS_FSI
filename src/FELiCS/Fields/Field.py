from dolfinx.fem             import Function, petsc
from FELiCS.IO.reader        import Reader
import basix

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

    def __init__(self, FEMSpace, mesh, name=[], isStateVector=False, m=None):
        self.space    = FEMSpace
        self.mesh     = mesh

        self.info     = self.describeFunctionSpace()
        # if the name is in the wrong format, re-format
        # TODO Sophie: reformat if there is only a list of variables for a mixed space?
        # TODO Sophie: throw warning / error if the number of names does not coincide with the number of subspaces
        if isinstance (name, list) and len(name)>0 and isinstance(name[0],tuple):
            self.name  = name
        elif isinstance(name, str):
            self.name  = [(name,[])]
        else:
            self.name = name

        self.isStateVector = isStateVector

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
        if self.name == None or len(self.name) == 0:
            return ""
        elif len(self.name) == 1: 
            return self.name[0][0]
        else:
            return self.name

    def getComponentsNames(self):
        ## This is  a workaround for now, to use for the retreat.
        ## TODO Sophie: make this independent of the coordinate system, and also usable for mixed function spaces.
        # TODO: make this independent of the variable name
        # space_info = self.describeFunctionSpace
        if isinstance (self.name, list) and len(self.name)>0 and isinstance(self.name[0],tuple) and \
            self.getName() in ["u", "rhou", "u_forcing_r", "u_forcing_i"]:
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
        if len(self.name)>1:
            componentsNames = []
            for n in self.name:
                componentsNames.append(n[0])
            return componentsNames
        else:
            return []

    def getTensor(self):
        from FELiCS.Misc.tensorUtils import Tensor
        # TODO Sophie: handle Tensors of mixed functions (later)
        return Tensor(self.function, self.mesh.coordinateSystem, m = self.m, mayHaveSpectralDimension = self.hasSpectralDimension)

    def isReal(self):
        import numpy as np
        return np.linalg.norm(np.imag(self.getCoefficientArray()))==0

    def getListOfSubFields(self):
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
        listOfFields            = []
        
        # Get the info about the field we are examining
        infoSpaceField          = self.describeFunctionSpace()
        numberOfSubSpaces       = self.space.num_sub_spaces

        if numberOfSubSpaces == 0: 
            listOfFields.append(self)
            return listOfFields
        
        # If we have a vector space, we need to create a new list of names
        if infoSpaceField['type'] == 'vector':
            name                = [None]*numberOfSubSpaces
        else:
            name                = self.name
        
        # Loop over sub-fields
        for i in range(numberOfSubSpaces):
            # For a vector space we give the name of the components
            if infoSpaceField['type'] == 'vector':
                comp    = self.getComponentsNames()
                name[i] = [(self.name[0][0]+comp[i],[])]
    
            # transfer content
            space, mapping            = self.space.sub(i).collapse()
            field                     = Field(space, self.mesh, name=name[i])
            field.function.x.array[:] = self.function.x.array[mapping]
            if len(self.name) == numberOfSubSpaces:  
                field.name = [self.name[i]]
            elif len(self.name)==1 and self.name[0] == None:
                field.name = []
            elif len(self.name)==1 and len(self.name[0]) == numberOfSubSpaces:
                field.name = [(self.name[0][0] + self.name[0][1][i], [])]
            listOfFields.append(field)

        return listOfFields

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
        num_subspaces                   = self.space.num_sub_spaces
        value_size                      = self.space.value_size
        nDofsMesh                       = len(self.mesh._coordinates)

        result = {
            'num_subspaces':            num_subspaces,
            'value_size':               value_size,
            'nDofsMesh':                nDofsMesh,
            'subspaces':                []
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
                    subspace            = self.space.sub(i)
                    sub_value_size      = subspace.value_size
                    sub_num_subspaces   = subspace.num_sub_spaces

                    if sub_value_size == 1:
                        subspace_type   = 'scalar'
                        subspace_desc   = f'Subspace {i}: scalar'
                    else:
                        subspace_type   = 'vector'
                        subspace_desc   = f'Subspace {i}: vector ({sub_value_size} components)'
                    
                    subspace_info = {
                        'index':        i,
                        'type':         subspace_type,
                        'value_size':   sub_value_size,
                        'num_sub_subspaces': sub_num_subspaces,
                        'description':  subspace_desc
                    }
                    
                    subspace_descriptions.append(subspace_info)
            
                result['subspaces']     = subspace_descriptions
                result['description']   = f'Mixed space with {num_subspaces} subspaces: ' + \
                                  ', '.join([sub['description'].split(': ')[1] for sub in subspace_descriptions])
        
        return result

    def importData(
            self,
            reader,
            importFilePath,
            groupName = None
        ):
        
        # Just call the reader function
        self, notInFile = reader.importInField(
            self,
            importFilePath,
            groupName
        )
        
        return self, notInFile






    def setListOfSubFields(self, listOfFields):
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
        # Get info about the field we are examining
        numberOfSubSpaces       = self.space.num_sub_spaces
        infoSpaceField          = self.describeFunctionSpace()
        
        # Prepare names for field
        # If we are assembling a vector we need only one name
        if infoSpaceField['type'] == 'vector':
            self.name = [listOfFields[0].name[0][0][:-1], self.getComponentsNames()]
        # Maybe we are silly and "assemble" a single scalar...    
        elif infoSpaceField['type'] == 'scalar':
            self.name = listOfFields[0].name
        # For a mixed space we keep the name as it was given
        elif infoSpaceField['type'] == 'mixed':
            self.name  = [None]*numberOfSubSpaces
        
        if numberOfSubSpaces == 0:
            # TODO: check that the functionSpaces are the same!
            self.setCoefficientArray(listOfFields[0].getCoefficientArray())
            self.name = listOfFields[0].name
            
        else:
            for i in range(numberOfSubSpaces):
                # TODO: check that the functionSpaces are the same!
                space, mapping                  = self.space.sub(i).collapse()
                self.function.x.array[mapping]  = listOfFields[i].getCoefficientArray()
                
                # Define the names of sub-fields if we deal with a Mixed field
                if infoSpaceField['type'] == 'mixed':
                    self.name[i]                = listOfFields[i].name

    def getSize(self):
        """
        Get the length of the coefficient array from the underlying function.

        Returns
        -------
        integer
            length of coefficient array 
        """
        return len(self.function.x.array[:])


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


    def getRealCoefficientArray(self):
        """
        Get the real part of the coefficient array of the field.

        Returns
        -------
        numpy.ndarray
            A float-valued array representing the field coefficients.
        """
        import numpy as np
        array = np.empty(len(self.function.x.array[:]),dtype=float)
        array[:] = np.real(self.function.x.array[:])
        return array

    def getImagCoefficientArray(self):
        """
        Get the imaginary part of the coefficient array of the field.

        Returns
        -------
        numpy.ndarray
            A float-valued array representing the field coefficients.
        """
        import numpy as np
        array = np.empty(len(self.function.x.array[:]),dtype=float)
        array[:] = np.imag(self.function.x.array[:])
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
        order = self.space.element.basix_element.degree
        dim = self.mesh.gdim 
        if self.hasSpectralDimension:
            dim += 1
        gradientSpace = getFELiCSSpace(self.mesh, order=order, dim=dim)
        gradientField = Field(gradientSpace, self.mesh)

        coordSystem = self.mesh.coordinateSystem
        J_hat = coordSystem.J_hat

        v = TestFunction(gradientSpace)
        # Sophie: m and hasSpectralDimension may not be needed for test function
        v_tens = Tensor(v, CoordSys=coordSystem, m = self.m, mayHaveSpectralDimension=self.hasSpectralDimension)
        expression = iDot(iGrad(self.getTensor()), iConj(v_tens)).ufl_tens * J_hat* dx

        gradientField.evaluateUflTensorExpression(expression)    

        # TODO Sophie: check how to destroy all petsc objects, also matrix
        gradientSpace.FEMWeightSolver.destroy()

        return gradientField

    def calculateL2Norm(self):
        import ufl
        import numpy as np
        from FELiCS.Misc.tensorUtils import iConj, iDot
        from FELiCS.Equation.UflDecorator import UflDecorator
        norm_squared = 0.
        list1 = self.getListOfSubFields()
        J_hat = self.mesh.coordinateSystem.J_hat
        norm_ufl = UflDecorator()
        for field in list1:
            fieldTens = field.getTensor()
            norm_ufl += (iDot(fieldTens, iConj(fieldTens))).ufl_tens*J_hat*ufl.dx        
        norm_squared = norm_ufl.getAssembledScalar(self.mesh)
        return np.sqrt(norm_squared)


    def getVorticityField(self):
        # TODO Sophie: throw error if Field is not vector
        # TODO: Add vorticity 3D field

        dim = self.space.num_sub_spaces
        if dim < 2:
            # Error message
            pass
        componentGradient = []
        velocityComponents = self.getListOfSubFields()

        for field in velocityComponents:
            componentGradient.append(field.getGradientField())
        
        
        # 2D field -> scalar vorticity field
        if dim == 2:
            dvdx = componentGradient[1].getListOfSubFields()[0]
            dudy = componentGradient[0].getListOfSubFields()[1]
            vorticityField = dvdx - dudy
            
            return vorticityField
        # 3D field -> vector vorticity field
        # TODO: Still needs to be tested
        if dim == 3:
            dwdy = componentGradient[2].getListOfSubFields()[1]
            dvdz = componentGradient[1].getListOfSubFields()[2]

            dudz = componentGradient[0].getListOfSubFields()[2]
            dwdx = componentGradient[2].getListOfSubFields()[0]

            dvdx = componentGradient[1].getListOfSubFields()[0]
            dudy = componentGradient[0].getListOfSubFields()[1]

            vorticity_x = dwdy - dvdz
            vorticity_y = dudz - dwdx
            vorticity_z = dvdx - dudy

            vorticityField = Field(self.space, self.mesh)
            vorticityField.setListOfSubFields([vorticity_x, vorticity_y, vorticity_z])


    def exportToH5(self, writer, fileName=None):
        # Sophie: This will be the final method
        if fileName == None:
            fileName = self.getName()
        writer.exportFieldToH5(self, fileName)


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

        # # this is only a dummy for the scripting
        if fileName == "function_values_2d":
            data = np.load(fileName+".npy").reshape(2,-1).T
            field1, field2 = self.getListOfSubFields()
            field1.setCoefficientArray(data[:,0])  
            field2.setCoefficientArray(data[:,1])
            self.setListOfSubFields([field1, field2])
        else:
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
        if not hasattr(self.space, 'FEMWeightSolver') or restartSolver:
            test_FEM   = ufl.TestFunctions(self.space)
            trial_FEM  = ufl.TrialFunctions(self.space)
            matrix_ufl = UflDecorator()
            coordinateSystem = self.mesh.coordinateSystem
            J_hat = coordinateSystem.J_hat
            i=0
            for test in test_FEM:
                iTest = Tensor(test, coordinateSystem, mayHaveSpectralDimension=True)
                iFluc = Tensor(trial_FEM[i], coordinateSystem, mayHaveSpectralDimension=True)
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
            iGrad,
        )
        ## create petsc solver and save it as attribute to the corresponding space - to use the LU-decomposition later 
        if not hasattr(self.space, 'FEMSmoothSolver') or restartSolver:
            test_FEM   = ufl.TestFunctions(self.space)
            trial_FEM  = ufl.TrialFunctions(self.space)
            matrix_ufl = UflDecorator()
            coordinateSystem = self.mesh.coordinateSystem
            J_hat = coordinateSystem.J_hat
            i=0
            for test in test_FEM:
                iTest = Tensor(test, coordinateSystem, mayHaveSpectralDimension=True)
                iFluc = Tensor(trial_FEM[i], coordinateSystem, mayHaveSpectralDimension=True)
                if iTest.order  == 1:
                    matrix_ufl.add( ( iDot(iFluc, iConj(iTest)) ).ufl_tens*J_hat*ufl.dx)
                    matrix_ufl.add((smoothFactor* iInner(iGrad(iFluc),iGrad(iTest))).ufl_tens*J_hat*ufl.dx)
                elif iTest.order == 0:
                    matrix_ufl.add( ( iFluc * iConj(iTest)).ufl_tens*J_hat*ufl.dx)
                    matrix_ufl.add((smoothFactor* iDot(iGrad(iFluc),iGrad(iConj(iTest)) )).ufl_tens*J_hat*ufl.dx)
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
        petscVec = expr_ufl.getAssembledVector(self.mesh, bcs)
        self.setCoefficientArray(LinearSolver.solveEquationSystemWithPredefinedSolver(self.space.FEMSmoothSolver, petscVec))


    def smooth(self, smoothFactor, bcs=[], restartSolver=False):
        """
        Smooth the field using a diffusion-like approach.

        Parameters
        ----------
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
            iGrad,
        )
        ## create petsc solver and save it as attribute to the corresponding space - to use the LU-decomposition later 
        if not hasattr(self.space, 'FEMSmoothSolver') or restartSolver:
            test_FEM   = ufl.TestFunctions(self.space)
            trial_FEM  = ufl.TrialFunctions(self.space)
            matrix_ufl = UflDecorator()
            coordinateSystem = self.mesh.coordinateSystem
            J_hat = coordinateSystem.J_hat
            i=0
            for test in test_FEM:
                iTest = Tensor(test, coordinateSystem, mayHaveSpectralDimension=True)
                iFluc = Tensor(trial_FEM[i], coordinateSystem, mayHaveSpectralDimension=True)
                if iTest.order  == 1:
                    matrix_ufl.add( ( iDot(iFluc, iConj(iTest)) ).ufl_tens*J_hat*ufl.dx)
                    matrix_ufl.add((smoothFactor* iInner(iGrad(iFluc),iGrad(iTest))).ufl_tens*J_hat*ufl.dx)
                elif iTest.order == 0:
                    matrix_ufl.add( ( iFluc * iConj(iTest)).ufl_tens*J_hat*ufl.dx)
                    matrix_ufl.add((smoothFactor* iDot(iGrad(iFluc),iGrad(iConj(iTest)) )).ufl_tens*J_hat*ufl.dx)
                else:
                    ##LOGGING TODO (Sophie): throw error 
                    print("ERROR, in 'Field.smoothTensorUflExpression'")
                i+=1
            matrix = matrix_ufl.getAssembledMatrix(self.mesh, bcs=bcs)
            self.space.FEMSmoothSolver = LinearSolver.createEquationSystemSolver(matrix)

        ## assemble rhs and solve equation system
        expr_ufl     = UflDecorator()
        listOfFields = self.getListOfSingleFields()
        test_FEM     = ufl.TestFunctions(self.space)
        coordinateSystem = self.mesh.coordinateSystem
        J_hat = coordinateSystem.J_hat
        for i in range(len(listOfFields)):
            iTest = Tensor(test_FEM[i], coordinateSystem, mayHaveSpectralDimension=True)
            field     = listOfFields[i]
            expr_ufl += (iDot(field.getTensor(), iConj(iTest))).ufl_tens*J_hat*ufl.dx
        petscVec = expr_ufl.getAssembledVector(self.mesh, bcs)
        self.setCoefficientArray(LinearSolver.solveEquationSystemWithPredefinedSolver(self.space.FEMSmoothSolver, petscVec))



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
        import numpy as np
        ## overrides '+'
        ## returns newly created Field with a coefficient array, which is the sum of two given coefficientarrays
        # TODO Sophie: raise error / not implemented if fields are not defined on the same space
        if isinstance(other, Field):
            result = Field(self.space, self.mesh)
            result.setCoefficientArray(self.getCoefficientArray() + other.getCoefficientArray())
            return result 
        elif np.isscalar(other):
            result = Field(self.space, self.mesh)
            result.setCoefficientArray(self.getCoefficientArray() + other)
            return result
        return NotImplemented
    

    def __sub__(self, other):
        """
        Overload the `-` operator for adding two Field objects.

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
            result.setCoefficientArray(self.getCoefficientArray() - other.getCoefficientArray())
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
        ## overrides '/'
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



