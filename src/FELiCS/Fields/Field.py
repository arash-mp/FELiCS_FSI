#  ___________________________________   _______________________________________________________
# /-----------------------------------\ /-------------------------------------------------------\
# |   (         (                (     |  This source code is part of FELiCS                     |
# |   )\ )     ) )        (      )\ )  |  (F)inite (E)lement (Li)nearized (C)ombustion (S)olver  |  
# |  (()/(  (  (()/( (    )\   (()/(   |                                                         |  
# |  /(_)) )\  /(_)))\  (((_)  /(_))   |  Licensed under the GNU GPLv3                           |
# |  (_)_)((_) (_)) ((_) )\___ (_))    |                                                         |
# |  | __|| __|| |   (_)((/ __|/ __|   |  (C) 2018-2025: The FELiCS Developers (www.felics.eu)   |
# |  | _| | _| | |__ | | | (__ \__ \   |  Visit          www.felics.eu                           |
# |  |_|  |___||____||_|  \___||___/   |  Contact        info@felics.eu                          |
# \___________________________________/ \_______________________________________________________/
#

import numpy as np
from petsc4py import PETSc
# TODO Sophie: throw error if Field is not scalar    
from ufl import TestFunction, dx
from FELiCS.Misc.tensorUtils import i_grad, i_conj, i_dot,  i_inner, Tensor
from FELiCS.SpaceDisc.FEMSpaces import create_function_space
import ufl
from FELiCS.Equation.UflDecorator import UflDecorator
import os
from FELiCS.Solvers.LinearSolver import LinearSolver
import dolfinx
# evaluates an ufl expression by 

import  matplotlib.pyplot           as plt
from    matplotlib.tri              import Triangulation
from    dolfinx                     import fem
from    basix.ufl                   import element
from    dolfinx.fem             import Function, petsc
from    FELiCS.Misc.logging     import Logger
from    mpl_toolkits.axes_grid1 import make_axes_locatable 

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
    def __init__(self,
    FEMSpace,
    mesh,
    name=None,
    isStateVector=False,
    m=None,
    ):
        self.space    = FEMSpace
        self.mesh     = mesh

        # get Info
        self.info     = self.describe_function_space()

        if name is None or not isinstance(
        name,
        str,
        ):
            if self.info['type']=='scalar':
                self._name = 'scalarField'
            elif self.info['type']=='vector':
                self._name = 'vectorField'
            elif self.info['type']=='mixed':
                self._name = 'mixedField'
            else: 
                self._name = 'field'
                logger.warning("Field type not recognized. Using default name 'field'.")

            if not isinstance(
            name,
            str,
            ) and name is not None:
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
    def name(self,
    name,
    ):
        self._name = name

    def get_names_of_sub_fields(self):
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
            axis_names           = self.mesh.axis_names
            numSubSpaces        = self.info['num_subspaces']
            
            # Check that the number of axis names is sufficient
            if len(axis_names) < numSubSpaces:
                logger.error(f"Not enough axis names {axis_names} in the coordinate system for the vector field with {numSubSpaces} components.")
                raise ValueError("Not enough axis names in the coordinate system for the vector field.")

            # If the vector was not given before, we set a default
            if self._name is None or not isinstance(
            self._name,
            str,
            ):
                self._name       = 'vectorField'

            for i in range(numSubSpaces):
                subFieldNames.append(self._name+axis_names[i])
        
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

    def set_names_of_sub_fields(self,
    nameList,
    ):
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


    def get_tensor(self):
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
        # TODO Sophie: handle Tensors of mixed functions (later)
        return Tensor(
        self.function,
        self.mesh.coordinate_system,
        m = self.m,
        mayHaveSpectralDimension = self.hasSpectralDimension,
        )

    def is_real(self):
        """
        Check whether the field is (numerically) real-valued.

        Returns
        -------
        bool
            True if the imaginary part of the coefficient array has zero
            norm, False otherwise.
        """
        return np.linalg.norm(np.imag(self.get_coefficient_array()))==0

    def get_list_of_sub_fields(self):
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
        namesOfSubFields        = self.get_names_of_sub_fields()
        
        # Loop over sub-fields and append to list
        for i in range(numberOfSubSpaces):
            # transfer content
            space, mapping            = self.space.sub(i).collapse()
            # NOTE: the sub-fields inherit the wave number from the field.
            field                     = Field(
            space,
            self.mesh,
            name=namesOfSubFields[i],
            m=self.m,
            )
            field.set_coefficient_array(self.get_coefficient_array()[mapping])
            listOfFields.append(field)

        return listOfFields



    def set_list_of_sub_fields(self,
    listOfFields,
    name=[],
    ):
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
            self.set_coefficient_array(listOfFields[0].get_coefficient_array())
            if not name or not isinstance(
            name,
            str,
            ):
                self._name = listOfFields[0].name
            else:
                self._name = name
            
        else:
            namesOfSubFields = []
            for i in range(numberOfSubSpaces):
                # TODO: check that the functionSpaces are the same!
                space, mapping                  = self.space.sub(i).collapse()
                self.function.x.array[mapping]  = listOfFields[i].get_coefficient_array()
                namesOfSubFields.append(listOfFields[i].name)
                self.set_names_of_sub_fields(namesOfSubFields)
                if isinstance(
                name,
                str,
                ):
                    self._name = name
                


    def describe_function_space(self):
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


    def get_coefficient_array(self):
        """
        Get the coefficient array of the field.

        Returns
        -------
        numpy.ndarray
            A complex-valued array representing the field coefficients.
        """
        array = np.empty(
        len(self.function.x.array[:]),
        dtype=complex,
        )
        array[:] = self.function.x.array[:]
        return array


    def get_real_coefficient_array(self):
        """
        Get the real part of the coefficient array of the field.

        Returns
        -------
        numpy.ndarray
            A float-valued array representing the field coefficients.
        """
        array = np.empty(
        len(self.function.x.array[:]),
        dtype=float,
        )
        array[:] = np.real(self.function.x.array[:])
        return array

    def get_imag_coefficient_array(self):
        """
        Get the imaginary part of the coefficient array of the field.

        Returns
        -------
        numpy.ndarray
            A float-valued array representing the field coefficients.
        """
        array = np.empty(
        len(self.function.x.array[:]),
        dtype=float,
        )
        array[:] = np.imag(self.function.x.array[:])
        return array



    def set_coefficient_array(self,
    array,
    ):
        """
        Set the coefficient array of the field.

        Parameters
        ----------
        array : numpy.ndarray
            A complex-valued array of coefficients.
        """
        if np.isscalar(array):
            self.function.x.array[:] = array
        else:
            self.function.x.array[:] = array[:]

    def set_constant_value(self,
    value,
    ): 
        """
        Set all coefficients to a constant value.

        Parameters
        ----------
        value : complex or float
            The constant value to assign to all coefficients.
        """
        self.function.x.array[:] = value

    def get_petsc_vector(self):
        """
        Convert the field to a PETSc vector.

        Returns
        -------
        petsc4py.PETSc.Vec
            A PETSc vector created from the coefficient array.
        """
        return PETSc.Vec().createWithArray(self.get_coefficient_array())

    def conjugate(self):
        """
        Compute the complex conjugate of the field.
        """
        self.set_coefficient_array(np.conj(self.get_coefficient_array()))

    def set_boundary_conditions(self,
    bcs,
    ):
        """
        Apply boundary conditions to the field.

        Parameters
        ----------
        bcs : list
            A list of Dirichlet boundary conditions to be applied.
        """
        petscArray = self.get_petsc_vector()
        petsc.set_bc(
        petscArray,
        bcs,
        )
        self.set_coefficient_array(petscArray.getArray())

        
    def get_gradient_field(self):
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
        # Create  a Field for the gradient
        # The space must be a vector vor a scalar field
        # TODO Sophie: handle order (get it from function?)
        degree = self.space.element.basix_element.degree
        dim = self.mesh.gdim 
        if self.hasSpectralDimension:
            dim += 1
        gradientSpace = create_function_space(
        self.mesh,
        degree=degree,
        dim=dim,
        )
        gradientField = Field(
        gradientSpace,
        self.mesh,
        )

        coordSystem = self.mesh.coordinate_system
        J_hat = coordSystem.J_hat

        v = TestFunction(gradientSpace)
        # Sophie: m and hasSpectralDimension may not be needed for test function
        v_tens = Tensor(
        v,
        CoordSys=coordSystem,
        m = self.m,
        mayHaveSpectralDimension=self.hasSpectralDimension,
        )
        expression = i_dot(
        i_grad(self.get_tensor()),
        i_conj(v_tens),
        ).ufl_tens * J_hat* dx

        gradientField.evaluate_ufl_tensor_expression(expression)    

        # TODO Sophie: check how to destroy all petsc objects, also matrix
        gradientSpace.FEMWeightSolver.destroy()

        return gradientField

    def calculate_l2_norm(self):
        """
        Compute the L2 norm of the field (and its subfields, if mixed).

        For mixed or vector-valued fields, the norm is computed by summing
        the contributions of each scalar subfield in the mixed decomposition.

        Returns
        -------
        float
            The L2 norm of the field.
        """
        norm_squared = 0.
        list1 = self.get_list_of_sub_fields()
        J_hat = self.mesh.coordinate_system.J_hat
        norm_ufl = UflDecorator()
        for field in list1:
            fieldTens = field.get_tensor()
            norm_ufl += (i_dot(
            fieldTens,
            i_conj(fieldTens),
            )).ufl_tens*J_hat*ufl.dx        
        norm_squared = norm_ufl.get_assembled_scalar(self.mesh)
        return np.sqrt(norm_squared)


    def get_vorticity_field(self):
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
        velocityComponents = self.get_list_of_sub_fields()

        for field in velocityComponents:
            componentGradient.append(field.get_gradient_field())
        
        
        # 2D field -> scalar vorticity field
        if dim == 2:
            dvdx = componentGradient[1].get_list_of_sub_fields()[0]
            dudy = componentGradient[0].get_list_of_sub_fields()[1]
            vorticityField = dvdx - dudy
            vorticityField.name = "vorticity"
            
            return vorticityField
        # 3D field -> vector vorticity field
        # TODO: Still needs to be tested
        if dim == 3:
            dwdy = componentGradient[2].get_list_of_sub_fields()[1]
            dvdz = componentGradient[1].get_list_of_sub_fields()[2]

            dudz = componentGradient[0].get_list_of_sub_fields()[2]
            dwdx = componentGradient[2].get_list_of_sub_fields()[0]

            dvdx = componentGradient[1].get_list_of_sub_fields()[0]
            dudy = componentGradient[0].get_list_of_sub_fields()[1]

            vorticity_x = dwdy - dvdz
            vorticity_y = dudz - dwdx
            vorticity_z = dvdx - dudy

            vorticityField = Field(
            self.space,
            self.mesh,
            name="vorticity",
            )
            vorticityField.set_list_of_sub_fields([vorticity_x, vorticity_y, vorticity_z])
            
            return vorticityField


    def export_to_h5(self,
    writer,
    fileName=None,
    ):
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
        writer.export_field_to_h5(
        self,
        fileName,
        )


    def import_data(
            self,
    reader,
    importFilePath = None,
    groupName = None,
    
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
        # determine filepath if not given
        # TODO Sophie: these two lines should be in the reader
        if importFilePath==None:
            importFilePath = os.path.join(
            reader._sourceDir,
            self.name,
            )
        # Just call the reader function
        self, notInFile = reader.import_in_field(
        self,
        importFilePath,
        groupName,
        )
        
        return self, notInFile


    def evaluate_ufl_expression(self,
    ufl_expression,
    bcs=[],
    restartSolver=False,
    ):
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
        ## create petsc solver and save it as attribute to the corresponding space - to use the LU-decomposition later 
        if not hasattr(
        self.space,
        'FEMWeightSolver',
        ) and not restartSolver:
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
            matrix = matrix_ufl.get_assembled_matrix(
            self.mesh,
            bcs,
            )
            self.space.FEMWeightSolver = LinearSolver.create_equation_system_solver(matrix)

        ## assemble rhs and solve equation system
        expr_ufl = UflDecorator(ufl_expression)
        petscVec = expr_ufl.get_assembled_vector(
        self.mesh,
        bcs,
        )
        self.set_coefficient_array(LinearSolver.solve_equation_system_with_predefined_solver(
        self.space.FEMWeightSolver,
        petscVec,
        ))


    def evaluate_ufl_tensor_expression(self,
    ufl_expression,
    bcs=[],
    restartSolver=False,
    ):
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
        ## create petsc solver and save it as attribute to the corresponding space - to use the LU-decomposition later 
        if not hasattr(
        self.space,
        'FEMWeightSolver',
        ) or restartSolver:
            test_FEM   = ufl.TestFunctions(self.space)
            trial_FEM  = ufl.TrialFunctions(self.space)
            matrix_ufl = UflDecorator()
            coordinate_system = self.mesh.coordinate_system
            J_hat = coordinate_system.J_hat
            i=0
            for test in test_FEM:
                iTest = Tensor(
                test,
                coordinate_system,
                mayHaveSpectralDimension=True,
                )
                iFluc = Tensor(
                trial_FEM[i],
                coordinate_system,
                mayHaveSpectralDimension=True,
                )
                if iTest.order  == 1:
                    matrix_ufl.add( ( i_dot(
                    iFluc,
                    i_conj(iTest),
                    ) ).ufl_tens*J_hat*ufl.dx)
                elif iTest.order == 0:
                    matrix_ufl.add( ( iFluc * i_conj(iTest)).ufl_tens*J_hat*ufl.dx)
                else:
                    ##LOGGING TODO (Sophie): throw error 
                    print("ERROR, in 'Field.evaluateUflExpression'")
                i+=1
            matrix = matrix_ufl.get_assembled_matrix(
            self.mesh,
            bcs,
            ) 
            self.space.FEMWeightSolver = LinearSolver.create_equation_system_solver(matrix)

        ## assemble rhs and solve equation system
        expr_ufl = UflDecorator(ufl_expression)
        petscVec = expr_ufl.get_assembled_vector(
        self.mesh,
        bcs,
        )
        self.set_coefficient_array(LinearSolver.solve_equation_system_with_predefined_solver(
        self.space.FEMWeightSolver,
        petscVec,
        ))


    def smooth_ufl_tensor_expression(self,
    ufl_expression,
    smoothFactor,
    bcs=[],
    restartSolver=False,
    ):
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
        ## create petsc solver and save it as attribute to the corresponding space - to use the LU-decomposition later 
        if not hasattr(
        self.space,
        'FEMSmoothSolver',
        ) or restartSolver:
            test_FEM   = ufl.TestFunctions(self.space)
            trial_FEM  = ufl.TrialFunctions(self.space)
            matrix_ufl = UflDecorator()
            coordinate_system = self.mesh.coordinate_system
            J_hat = coordinate_system.J_hat
            i=0
            for test in test_FEM:
                iTest = Tensor(
                test,
                coordinate_system,
                mayHaveSpectralDimension=True,
                )
                iFluc = Tensor(
                trial_FEM[i],
                coordinate_system,
                mayHaveSpectralDimension=True,
                )
                if iTest.order  == 1:
                    matrix_ufl.add( ( i_dot(
                    iFluc,
                    i_conj(iTest),
                    ) ).ufl_tens*J_hat*ufl.dx)
                    matrix_ufl.add((smoothFactor* i_inner(
                    i_grad(iFluc),
                    i_grad(iTest),
                    )).ufl_tens*J_hat*ufl.dx)
                elif iTest.order == 0:
                    matrix_ufl.add( ( iFluc * i_conj(iTest)).ufl_tens*J_hat*ufl.dx)
                    matrix_ufl.add((smoothFactor* i_dot(
                    i_grad(iFluc),
                    i_grad(i_conj(iTest)),
                    )).ufl_tens*J_hat*ufl.dx)
                else:
                    ##LOGGING TODO (Sophie): throw error 
                    print("ERROR, in 'Field.smoothTensorUflExpression'")
                i+=1
            try:
                matrix_ufl.setCorrectMeshObject(self.mesh)
            except:
                pass
            matrix = petsc.assemble_matrix(
            dolfinx.fem.form(matrix_ufl.lhs),
            bcs=bcs,
            )
            matrix.assemble()
            self.space.FEMSmoothSolver = LinearSolver.create_equation_system_solver(matrix)

        ## assemble rhs and solve equation system
        expr_ufl = UflDecorator(ufl_expression)
        petscVec = expr_ufl.get_assembled_vector(
        self.mesh,
        bcs,
        )
        self.set_coefficient_array(LinearSolver.solve_equation_system_with_predefined_solver(
        self.space.FEMSmoothSolver,
        petscVec,
        ))


    def smooth(self,
    smoothFactor,
    bcs=[],
    restartSolver=False,
    ):
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
        ## create petsc solver and save it as attribute to the corresponding space - to use the LU-decomposition later 
        if not hasattr(
        self.space,
        'FEMSmoothSolver',
        ) or restartSolver:
            test_FEM   = ufl.TestFunctions(self.space)
            trial_FEM  = ufl.TrialFunctions(self.space)
            matrix_ufl = UflDecorator()
            coordinate_system = self.mesh.coordinate_system
            J_hat = coordinate_system.J_hat
            i=0
            for test in test_FEM:
                iTest = Tensor(
                test,
                coordinate_system,
                mayHaveSpectralDimension=True,
                )
                iFluc = Tensor(
                trial_FEM[i],
                coordinate_system,
                mayHaveSpectralDimension=True,
                )
                if iTest.order  == 1:
                    matrix_ufl.add( ( i_dot(
                    iFluc,
                    i_conj(iTest),
                    ) ).ufl_tens*J_hat*ufl.dx)
                    matrix_ufl.add((smoothFactor* i_inner(
                    i_grad(iFluc),
                    i_grad(iTest),
                    )).ufl_tens*J_hat*ufl.dx)
                elif iTest.order == 0:
                    matrix_ufl.add( ( iFluc * i_conj(iTest)).ufl_tens*J_hat*ufl.dx)
                    matrix_ufl.add((smoothFactor* i_dot(
                    i_grad(iFluc),
                    i_grad(i_conj(iTest)),
                    )).ufl_tens*J_hat*ufl.dx)
                else:
                    ##LOGGING TODO (Sophie): throw error 
                    print("ERROR, in 'Field.smoothTensorUflExpression'")
                i+=1
            matrix = matrix_ufl.get_assembled_matrix(
            self.mesh,
            bcs=bcs,
            )
            self.space.FEMSmoothSolver = LinearSolver.create_equation_system_solver(matrix)

        ## assemble rhs and solve equation system
        expr_ufl     = UflDecorator()
        listOfFields = self.get_list_of_sub_fields()
        test_FEM     = ufl.TestFunctions(self.space)
        coordinate_system = self.mesh.coordinate_system
        J_hat = coordinate_system.J_hat
        for i in range(len(listOfFields)):
            iTest = Tensor(
            test_FEM[i],
            coordinate_system,
            mayHaveSpectralDimension=True,
            )
            field     = listOfFields[i]
            expr_ufl += (i_dot(
            field.get_tensor(),
            i_conj(iTest),
            )).ufl_tens*J_hat*ufl.dx
        petscVec = expr_ufl.get_assembled_vector(
        self.mesh,
        bcs,
        )
        self.set_coefficient_array(LinearSolver.solve_equation_system_with_predefined_solver(
        self.space.FEMSmoothSolver,
        petscVec,
        ))

    def plot(self,
    xlim=None,
    ylim=None,
    plotType="real",
    clim=None,
    ):
        """
        Plotting function for debugging purposes. This function can be used, to check if a
        field looks as expected and rule out e.g. import problems.

        Notes
        -----
        - This method provides a simple visualization of the field.
        """

        if self.space.num_sub_spaces > 1:
            raise NotImplementedError("Plotting is only implemented for scalar fields. Use getListOfSingleFields() to get subfields. These can then be plotted individually with the same method.")

        # ---- METHOD 1: Directly use dof coordinates and tricontourf ----
        # NOTE: With this method we loose mesh connectivity information.
        # phi_vertex      = self.getCoefficientArray()
        # dof_coordinates = self.space.tabulate_dof_coordinates()
        # x               = dof_coordinates[:, 0]
        # y               = dof_coordinates[:, 1]
        # triang          = Triangulation(x, y)
        # ---- End of METHOD 1 ----

        # ---- METHOD 2: Interpolate to P1 space and use tricontourf ----
        mesh            = self.mesh.dolfinxMesh      # dolfinx.mesh.Mesh
        u_h             = self.function              # fem.Function in Vh
        tdim            = mesh.topology.dim

        # Ensure cell->vertex connectivity
        mesh.topology.create_connectivity(
        tdim,
        0,
        )

        # --- 1) Interpolate to P1 space on same mesh ---
        V1              = fem.functionspace(
        mesh,
        element(
        "CG",
        "triangle",
        1,
        ),
        )
        u1              = fem.Function(V1)
        u1.interpolate(u_h)   # works if u_h is scalar-valued; see note below for vectors
        
        # Get the values corresponding to plot type
        if plotType == "imag":
            phi_vertex      = np.imag(u1.x.array)
            cmap            = "seismic"
        elif plotType == "magnitude":
            phi_vertex      = np.abs(u1.x.array)
            cmap            = "magma"
        else:
            phi_vertex      = np.real(u1.x.array)
            cmap            = "seismic"
        # Now u1.x.array has one value per vertex, in the same ordering as geometry.x

        # --- 2) Build triangulation from the mesh ---
        cells_to_vertices   = mesh.topology.connectivity(
        tdim,
        0,
        ).array
        triangles           = cells_to_vertices.reshape(
        -1,
        3,
        )
        coords              = mesh.geometry.x
        x                   = coords[:, 0]
        y                   = coords[:, 1]
        triang              = Triangulation(
        x,
        y,
        triangles=triangles,
        )
        # ---- End of METHOD 2 ----
        
        # ---- Plotting ----
        fig, axes       = plt.subplots()
        if clim is None:
            if plotType == "magnitude":
                clim = (0, np.max(phi_vertex))
            else:
                clim = (-0.5*np.max(np.abs(phi_vertex)), 0.5*np.max(np.abs(phi_vertex)))
        contour = axes.tripcolor(
        triang,
        phi_vertex,
        shading='gouraud',
        cmap=cmap,
        vmin=clim[0],
        vmax=clim[1],
        )

        # Set labels and title
        axes.set_xlabel('x')
        axes.set_ylabel('y')
        if self.name != "":
            title = self.name
        else:
            title = "scalar_field"
        if plotType == "imag":
            title += "_imag"
        elif plotType == "magnitude":
            title += "_magnitude"
        else:
            title += "_real"

        axes.set_title(title)
        axes.set_aspect('equal')
        axes.grid(
        True,
        alpha=0.3,
        )
        
        # Add colorbar for phi
        divider         = make_axes_locatable(axes)
        colorbar_axes   = divider.append_axes(
        "right",
        size="2%",
        pad=0.5,
        ) 
        cbar            = plt.colorbar(
        contour,
        label=title,
        cax=colorbar_axes,
        )
        cbar.formatter.set_powerlimits((0, 0))
        cbar.update_ticks()
        if xlim is not None:
            axes.set_xlim(xlim)
        if ylim is not None:
            axes.set_ylim(ylim)

        plt.tight_layout()
        plt.show()
        

    ### dunder methods for overloading arithmetic operators ###
    def __add__(self,
    other,
    ):
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
        ## overrides '+'
        ## returns newly created Field with a coefficient array, which is the sum of two given coefficientarrays
        # TODO Sophie: raise error / not implemented if fields are not defined on the same space
        if isinstance(
        other,
        Field,
        ):
            result = Field(
            self.space,
            self.mesh,
            )
            result.set_coefficient_array(self.get_coefficient_array() + other.get_coefficient_array())
            return result 
        elif np.isscalar(other):
            result = Field(
            self.space,
            self.mesh,
            )
            result.set_coefficient_array(self.get_coefficient_array() + other)
            return result
        return NotImplemented

    def __iadd__(self,
    other,
    ):
        """
        Overload the ``+=`` operator for adding fields or scalars.

        Parameters
        ----------
        other : Field or scalar
            Another Field object defined on the same space, or a scalar
            value to be added to all coefficients.

        Returns
        -------
        self
        """
        ## overrides '+'
        ## returns newly created Field with a coefficient array, which is the sum of two given coefficientarrays
        # TODO Sophie: raise error / not implemented if fields are not defined on the same space
        if isinstance(
        other,
        Field,
        ):
            self.set_coefficient_array(self.get_coefficient_array() + other.get_coefficient_array())
            return self
        elif np.isscalar(other):
            self.set_coefficient_array(self.get_coefficient_array() + other)
            return self
        return NotImplemented
 

    def __sub__(self,
    other,
    ):
        """
        Overload the `-` operator for adding two Field objects.

        Parameters
        ----------
        other : Field or scalar
            Another Field object defined on the same space, or a scalar
            value to be added to all coefficients.

        Returns
        -------
        Field
            A new Field object with the substracted coefficient arrays.

        Raises
        ------
        NotImplementedError
            If `other` is not a Field object.
        """
        ## overrides '-'
        ## returns newly created Field with a coefficient array, which is the sum of two given coefficientarrays
        # TODO Sophie: raise error / not implemented if fields are not defined on the same space
        if isinstance(
        other,
        Field,
        ):
            result = Field(
            self.space,
            self.mesh,
            )
            result.set_coefficient_array(self.get_coefficient_array() - other.get_coefficient_array())
            return result 
        elif np.isscalar(other):
            result = Field(
            self.space,
            self.mesh,
            )
            result.set_coefficient_array(self.get_coefficient_array() - other)
            return result 
        return NotImplemented

    def __isub__(self,
    other,
    ):
        """
        Overload the `-=` operator for adding two Field objects.

        Parameters
        ----------
        other : Field
            Another Field object.

        Returns
        -------
        self

        Raises
        ------
        NotImplementedError
            If `other` is not a Field object.
        """
        ## overrides '-='
        ## returns newly created Field with a coefficient array, which is the sum of two given coefficientarrays
        # TODO Sophie: raise error / not implemented if fields are not defined on the same space
        if isinstance(
        other,
        Field,
        ):
            self.set_coefficient_array(self.get_coefficient_array() - other.get_coefficient_array())
            return self
        elif np.isscalar(other):
            self.set_coefficient_array(self.get_coefficient_array() - other)
            return self
        return NotImplemented


    def __mul__(self,
    other,
    ):
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
        ## overrides '*'
        ## returns newly created Field with a coefficient array, which is the product of two given coefficientarrays, or the product of its coefficientarray with a scalar value
        # TODO Sophie: raise error / not implemented if fields are not defined on the same space
        if isinstance(
        other,
        Field,
        ):
            result = Field(
            self.space,
            self.mesh,
            )
            result.set_coefficient_array(self.get_coefficient_array()*other.get_coefficient_array())
            return result
        elif np.isscalar(other):
            result = Field(
            self.space,
            self.mesh,
            )
            result.set_coefficient_array(self.get_coefficient_array()*other)
            return result
        return NotImplemented


    def __imul__(self,
    other,
    ):
        """
        Overload the ``*=`` operator for pointwise multiplication.

        Parameters
        ----------
        other : Field or scalar
            Another Field defined on the same space, or a scalar value.

        Returns
        -------
        self
        """
        ## overrides '*='
        ## returns newly created Field with a coefficient array, which is the product of two given coefficientarrays, or the product of its coefficientarray with a scalar value
        # TODO Sophie: raise error / not implemented if fields are not defined on the same space
        if isinstance(
        other,
        Field,
        ):
            self.set_coefficient_array(self.get_coefficient_array()*other.get_coefficient_array())
            return self
        elif np.isscalar(other):
            self.set_coefficient_array(self.get_coefficient_array()*other)
            return self
        return NotImplemented


    def __truediv__(self,
    other,
    ):
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
        ## overrides '/'
        ## returns newly created Field with a coefficient array, which is the division of two given coefficientarrays, or the division of its coefficientarray with a scalar value
        # TODO Sophie: raise error / not implemented if fields are not defined on the same space
        if isinstance(
        other,
        Field,
        ):
            result = Field(
            self.space,
            self.mesh,
            )
            result.set_coefficient_array(self.get_coefficient_array()/other.get_coefficient_array())
            return result
        elif np.isscalar(other):
            result = Field(
            self.space,
            self.mesh,
            )
            result.set_coefficient_array(self.get_coefficient_array()/other)
            return result
        return NotImplemented


    def __itruediv__(self,
    other,
    ):
        """
        Overload the ``/=`` operator for pointwise division.

        Parameters
        ----------
        other : Field or scalar
            Another Field defined on the same space, or a scalar value.

        Returns
        -------
        self

        Notes
        -----
        - Division by a Field is performed coefficient-wise; it is the
          caller's responsibility to avoid division by zero.

        """
        ## overrides '/='
        ## returns newly created Field with a coefficient array, which is the division of two given coefficientarrays, or the division of its coefficientarray with a scalar value
        # TODO Sophie: raise error / not implemented if fields are not defined on the same space
        if isinstance(
        other,
        Field,
        ):
            self.set_coefficient_array(self.get_coefficient_array()/other.get_coefficient_array())
            return self
        elif np.isscalar(other):
            self.set_coefficient_array(self.get_coefficient_array()/other)
            return self
        return NotImplemented




