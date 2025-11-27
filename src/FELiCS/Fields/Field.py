from    dolfinx.fem             import Function, petsc
from    FELiCS.Misc.logging     import Logger

# Get the logger
logger = Logger.get_logger("felics")

class Field:
    """
    Class representing a finite element field.

    This class provides methods for handling finite element fields, including 
    coefficient manipulation, boundary condition application, and expression 
    evaluation. It supports complex-valued fields and can optionally represent
    a state vector in a mixed formulation. The class also supports an optional
    spectral dimension via a wave number.

    The class interacts with ``dolfinx.fem.Function`` for finite element operations 
    and includes utilities for working with PETSc vectors and UFL expressions.

    **Initialize the Field object**

    Parameters
    ----------
    FEMSpace : dolfinx.fem.FunctionSpace
        The finite element function space.
    mesh : FELiCS.SpaceDisc.FELiCSMesh
        FELiCS mesh associated with the function space.
    name : str, optional
        Name of the field. If not provided, a default name is chosen based on
        the field type (scalar, vector, or mixed).
    isStateVector : bool, optional
        If True, the field is interpreted as a state vector in a mixed space,
        and subfield names may be taken from ``space.stateVectorNames`` where
        available. Default is False.
    m : int, optional
        Wave number associated with a spectral spatial dimension. If this is
        set (even to zero), the field is assumed to have one spectral spatial
        dimension, and ``hasSpectralDimension`` is set to True.
    """
    def __init__(self, FEMSpace, mesh, name=None, isStateVector=False, m=None):
        self.space    = FEMSpace
        self.mesh     = mesh

        # get Info
        self.info     = self.describeFunctionSpace()

        if name is None or not isinstance(name, str):
            if self.info['type']=='scalar':
                self._name = 'scalarField'
            elif self.info['type']=='vector':
                self._name = 'vectorField'
            elif self.info['type']=='mixed':
                self._name = 'mixedField'
            else: 
                self._name = 'field'
                logger.warning("Field type not recognized. Using default name 'field'.")

            if not isinstance(name, str) and name is not None:
                logger.warning("Field initialized a 'name' not being a string. Using default names based on field type.")
                
        else:
            self._name = name

        self._namesOfSubFields = []
        
        self.isStateVector = isStateVector

        # handle spectral dimension and wave number
        if m is not None:
            self.hasSpectralDimension = True
            self.m = m
        else:
            self.hasSpectralDimension = False
            self.m = 0

        # initialize function
        self.function = Function(FEMSpace)

    @property 
    def name(self):
        return self._name
  
    @name.setter
    def name(self, name):
        self._name = name

    def getNamesOfSubFields(self):
        """
        Return the list of subfield names for vector or mixed spaces.

        For scalar fields, an empty list is returned and a warning is logged.
        For vector fields, names are generated from the base field name and
        the mesh axis names (e.g. ``u_x``, ``u_y``). For mixed fields, names
        are taken from previously set values, from ``stateVectorNames`` (when
        ``isStateVector`` is True), or reasonable defaults such as
        ``scalar1``, ``vector1``, etc.

        Returns
        -------
        list of str
            List of subfield names. Empty for pure scalar fields.
        """
        # Default value
        subFieldNames           = []
        
        if self.info['type'] == 'scalar':
            logger.warning("getNamesOfSubFields() called for single scalar field. Returning empty list.")
            return subFieldNames
        
        elif self.info['type'] == 'vector':
            axisNames           = self.mesh.axisNames
            numSubSpaces        = self.info['num_subspaces']
            
            # Check that the number of axis names is sufficient
            if len(axisNames) < numSubSpaces:
                logger.error(f"Not enough axis names {axisNames} in the coordinate system for the vector field with {numSubSpaces} components.")
                raise ValueError("Not enough axis names in the coordinate system for the vector field.")

            # If the vector was not given before, we set a default
            if self._name is None or not isinstance(self._name, str):
                self._name       = 'vectorField'

            for i in range(numSubSpaces):
                subFieldNames.append(self._name+axisNames[i])
        
        elif self.info['type'] == 'mixed':
            # If it was not set before and is a state vector, use state vector names
            # NOTE: assumes stateVectorNames is a list of tuples 
            if self.isStateVector and not self._namesOfSubFields:
                for name in self.space.stateVectorNames:
                    subFieldNames.append(name[0])
                return subFieldNames
            
            # If it was set before, return the stored names
            elif self._namesOfSubFields:
                return self._namesOfSubFields
            
            # Otherwise, set default names
            else:
                counter_scalars     = 1
                counter_vectors     = 1
                for i in range(self.info['num_subspaces']):
                    if self.info['subspaces'][i]['type'] == 'scalar':
                        subFieldNames.append(f'scalar{counter_scalars}')
                        counter_scalars += 1
                    elif self.info['subspaces'][i]['type'] == 'vector':
                        subFieldNames.append(f'vector{counter_vectors}')
                        counter_vectors += 1
                    else:
                        logger.error("Subspace type neither scalar nor vector.")
                        raise ValueError("Subspace type neither scalar nor vector.") 

        return subFieldNames

    def setNamesOfSubFields(self, nameList):
        """
        Assign explicit names to the subfields of a mixed or vector field.

        Parameters
        ----------
        nameList : list of str
            List of subfield names. The length must match the number of
            subspaces in the underlying function space; otherwise, a warning
            is logged and the names are not changed.
        """
        if len(nameList) != self.info['num_subspaces']:
            logger.warning("setNamesOfSubFields() called with a list of names that does not match the number of subspaces. No names were set.")
            return
        self._namesOfSubFields = nameList 


    def getTensor(self):
        """
        Wrap the field as a tensor-aware object.

        Returns
        -------
        FELiCS.Misc.tensorUtils.Tensor
            Tensor wrapper around the underlying finite element function,
            constructed with the mesh coordinate system and the field's
            spectral settings (``m`` and ``hasSpectralDimension``).

        Notes
        -----
        - Mixed functions are not yet fully supported and are treated as a
          single tensor object.
        """
        from FELiCS.Misc.tensorUtils import Tensor
        # TODO Sophie: handle Tensors of mixed functions (later)
        return Tensor(self.function, self.mesh.coordinateSystem, m = self.m, mayHaveSpectralDimension = self.hasSpectralDimension)

    def isReal(self):
        """
        Check whether the field is (numerically) real-valued.

        Returns
        -------
        bool
            True if the imaginary part of the coefficient array has zero
            norm, False otherwise.
        """
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
        
        # Get the number of subspaces
        numberOfSubSpaces       = self.info['num_subspaces']

        # If there are no subspaces, return a list with one entry (this field)
        if numberOfSubSpaces == 0: 
            listOfFields.append(self)
            return listOfFields
      
        # Else get the names of the subfields
        namesOfSubFields        = self.getNamesOfSubFields()
        
        # Loop over sub-fields and append to list
        for i in range(numberOfSubSpaces):
            # transfer content
            space, mapping            = self.space.sub(i).collapse()
            field                     = Field(space, self.mesh, name=namesOfSubFields[i])
            field.setCoefficientArray(self.getCoefficientArray()[mapping])
            listOfFields.append(field)

        return listOfFields



    def setListOfSubFields(self, listOfFields, name=[]):
        """
        Set the field coefficients from a list of single-component fields.

        Parameters
        ----------
        listOfFields : list
            A list of `Field` objects representing individual subspaces.
        name: string, optional
            The name of the new field that contains all the given fields as subfields.

        Notes
        -----
        - This only works if the "listOfFields" contains fields with the correct spaces, 
          which are subspaces of the space which which this field has been initialized.
        - If the space has multiple subspaces, their coefficients are mapped back.
        - Throws an error if the input list does not match the expected size.
        """

        # Get info about the field we are examining
        numberOfSubSpaces       = self.info['num_subspaces']
        
        if numberOfSubSpaces == 0:
            # TODO: check that the functionSpaces are the same!
            self.setCoefficientArray(listOfFields[0].getCoefficientArray())
            if not name or not isinstance(name, str):
                self._name = listOfFields[0].name
            else:
                self._name = name
            
        else:
            namesOfSubFields = []
            for i in range(numberOfSubSpaces):
                # TODO: check that the functionSpaces are the same!
                space, mapping                  = self.space.sub(i).collapse()
                self.function.x.array[mapping]  = listOfFields[i].getCoefficientArray()
                namesOfSubFields.append(listOfFields[i].name)
                self.setNamesOfSubFields(namesOfSubFields)
                if isinstance(name,str):
                    self._name = name
                


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
        # TODO: remove after the writer has been updated. The input mesh in field HAS TO BE A FELiCSMesh, and CANNOT be a dolfinx mesh
        try:
            nDofsMesh                       = len(self.mesh._coordinates)
        except:
            nDofsMesh                       = self.space.dofmap.index_map.size_global

        result = {
            'num_subspaces':            num_subspaces,
            'value_size':               value_size,
            'nDofsMesh':                nDofsMesh,
            'subspaces':                []
        }

        # Single space (scalar)
        # This handles a single scalar or a VMixed with a single scalar inside
        if num_subspaces == 0 or (num_subspaces == 1 and value_size == 1):
            result['type']              = 'scalar'
            result['description']       = 'Single scalar space'
            result['degree']            = self.space.ufl_element().degree
            result['nDofsSpace']        = self.space.dofmap.index_map.size_global
            
        else:
            # Single vector space
            if num_subspaces == value_size:
                result['type']          = 'vector'
                result['description']   = f'Single vector space with {value_size} components'
                result['degree']        = self.space.ufl_element().degree
                result['nDofsSpace']    = self.space.sub(0).dofmap.index_map.size_global

            # Mixed space
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
                        subspace_nDofs  = len(subspace.collapse()[1])
                    else:
                        subspace_type   = 'vector'
                        subspace_desc   = f'Subspace {i}: vector ({sub_value_size} components)'
                        subspace_nDofs  = len(subspace.sub(0).collapse()[1])
                    
                    subspace_info = {
                        'index':        i,
                        'type':         subspace_type,
                        'value_size':   sub_value_size,
                        'num_sub_subspaces': sub_num_subspaces,
                        'description':  subspace_desc,
                        'degree':       subspace.ufl_element().degree,
                        'nDofsSpace':   subspace_nDofs
                    }
                    
                    subspace_descriptions.append(subspace_info)
            
                result['subspaces']     = subspace_descriptions
                result['description']   = f'Mixed space with {num_subspaces} subspaces: ' + \
                                  ', '.join([sub['description'].split(': ')[1] for sub in subspace_descriptions])
        
        return result

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
        """
        Compute the gradient of a scalar field as a new Field.

        The gradient is obtained in a vector-valued function space with
        dimension equal to the geometric dimension of the mesh (plus an
        additional spectral dimension if present). The result is computed
        via a tensor-based weak form and projection.

        Returns
        -------
        Field
            A new Field instance representing the gradient of the original
            field.

        """
        # TODO Sophie: throw error if Field is not scalar    
        from ufl import TestFunction, dx
        from FELiCS.Misc.tensorUtils import iGrad, iConj, iDot, Tensor
        from FELiCS.SpaceDisc.FEMSpaces import createFunctionSpace
        # Create  a Field for the gradient
        # The space must be a vector vor a scalar field
        # TODO Sophie: handle order (get it from function?)
        degree = self.space.element.basix_element.degree
        dim = self.mesh.gdim 
        if self.hasSpectralDimension:
            dim += 1
        gradientSpace = createFunctionSpace(self.mesh, degree=degree, dim=dim)
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
        """
        Compute the L2 norm of the field (and its subfields, if mixed).

        For mixed or vector-valued fields, the norm is computed by summing
        the contributions of each scalar subfield in the mixed decomposition.

        Returns
        -------
        float
            The L2 norm of the field.
        """
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
        """
        Compute the vorticity field associated with a velocity field.

        For a 2D velocity field, this returns a scalar vorticity field
        (ω_z = ∂v/∂x − ∂u/∂y). For a 3D velocity field, a vector-valued
        vorticity field is constructed component-wise using the curl of
        the velocity.

        Returns
        -------
        Field or None
            Vorticity field as a scalar (2D) or vector (3D) Field. Returns
            None if the number of components is less than 2.

        Notes
        -----
        - 3D behaviour is marked as experimental in the implementation and
          may not be fully tested.
        - No explicit error is raised if the field is not a velocity-type
          vector; it is the caller's responsibility to ensure consistency.
        """
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
        """
        Export the field to an HDF5 file using a FELiCS writer.

        Parameters
        ----------
        writer : object
            Writer object providing an ``exportFieldToH5(field, fileName)``
            method.
        fileName : str, optional
            Base name for the exported dataset or file. If None, the field's
            ``name`` attribute is used.
        """
        # Sophie: This will be the final method
        if fileName == None:
            fileName = self.name
        writer.exportFieldToH5(self, fileName)


    def importData(
            self,
            reader,
            importFilePath,
            groupName = None
        ):
        """
        Import data into the field using a FELiCS reader.

        This is a convenience wrapper around the reader's
        :meth:`importInField` method.

        Parameters
        ----------
        reader : object
            Reader object providing an ``importInField(field, filePath, groupName)``
            method.
        importFilePath : str
            Path to the file containing the data to be imported.
        groupName : str or None, optional
            Optional group name inside the file from which to read.

        Returns
        -------
        Field
            The updated field (self).
        list of str
            List of variable names that were not found in the file.
        """
        # Just call the reader function
        self, notInFile = reader.importInField(
            self,
            importFilePath,
            groupName
        )
        
        return self, notInFile


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
            iInner,
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

        This method constructs a smoothing operator that combines an
        identity term with a gradient-based regularization term, and then
        applies it to the current field by solving the corresponding weak
        form.

        Parameters
        ----------
        smoothFactor : float
            Smoothing factor controlling the influence of the gradient term.
        bcs : list, optional
            List of boundary conditions to apply.
        restartSolver : bool, optional
            If True, rebuild the smoothing solver even if one already exists.

        Notes
        -----
        - Particularly useful for regularizing noisy numerical solutions.
        - The same solver is reused between calls unless ``restartSolver`` is True.
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
            iInner,
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

    def plot(self):
        """
        Plotting function for debugging purposes. This function can be used, to check if a
        field looks as expected and rule out e.g. import problems.

        Notes
        -----
        - This method provides a simple visualization of the field.
        """
        import matplotlib.pyplot as plt
        from matplotlib.tri import Triangulation
        import numpy as np

        if self.space.num_sub_spaces > 1:
            raise NotImplementedError("Plotting is only implemented for scalar fields. Use getListOfSingleFields() to get subfields. These can then be plotted individually with the same method.")
        
        FieldsList = self.getListOfSubFields()

        # Create figure outside the loop
        fig, axes = plt.subplots(1, 1, figsize=(6, 6))

        # Create the plot
        phi = self.getCoefficientArray()
        dof_coordinates = self.space.tabulate_dof_coordinates()

        x = dof_coordinates[:, 0]
        y = dof_coordinates[:, 1]
        triang_scalar = Triangulation(x, y)

        # Create contour plot of phi
        contour = axes.tricontourf(triang_scalar, phi, levels=20, cmap='RdBu_r', alpha=0.7)
        # contour_lines = axes.tricontour(triang_scalar, phi, levels=10, colors='black', alpha=0.5, linewidths=0.5)

        # Add colorbar for phi
        cbar = plt.colorbar(contour, ax=axes, label=r'$\phi$ '+ self.name)

        # Set labels and title
        axes.set_xlabel('x')
        axes.set_ylabel('y')
        if self.name != "":
            axes.set_title(self.name)
        else:
            axes.set_title('Scalar Field ')
        axes.set_aspect('equal')
        axes.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()
        

    ### dunder methods for overloading arithmetic operators ###
    def __add__(self, other):
        """
        Overload the ``+`` operator for adding fields or scalars.

        Parameters
        ----------
        other : Field or scalar
            Another Field object defined on the same space, or a scalar
            value to be added to all coefficients.

        Returns
        -------
        Field
            A new Field object with the summed coefficient arrays, or the
            original field shifted by a scalar.

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
        """
        Overload the ``*`` operator for pointwise multiplication.

        Parameters
        ----------
        other : Field or scalar
            Another Field defined on the same space, or a scalar value.

        Returns
        -------
        Field
            A new Field whose coefficient array is the pointwise product
            of the two fields, or the field scaled by the scalar.
        """
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
        """
        Overload the ``/`` operator for pointwise division.

        Parameters
        ----------
        other : Field or scalar
            Another Field defined on the same space, or a scalar value.

        Returns
        -------
        Field
            A new Field whose coefficient array is the pointwise division
            ``self / other``.

        Notes
        -----
        - Division by a Field is performed coefficient-wise; it is the
          caller's responsibility to avoid division by zero.

        """
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



