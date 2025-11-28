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
from ufl import (
                lhs,
                rhs,
                )
from dolfinx.fem import (
                form,
                assemble_scalar,
                                )
from dolfinx.fem.petsc import (
                assemble_matrix,
                assemble_vector,
                set_bc
                )

from 	FELiCS.Misc.logging          import Logger

# Get the logger
logger = Logger.get_logger("felics")


class UflDecorator():
    """
    Decorator for handling UFL (Unified Form Language) expressions.

    This class wraps and extends functionality for UFL expressions. It provides methods 
    for symbolic manipulation, analysis, and assembling into finite element matrices, 
    vectors, or scalars. 

    **Initialize the UflDecorator object**

    Parameters
    ----------
    inputUfl : ufl.Form or None, optional
        Initial UFL expression to wrap. If None, the decorator starts without an expression.

    Attributes
    ----------
    lhs : ufl.Form
        The left-hand side form of the expression.
    rhs : ufl.Form
        The right-hand side form of the expression.
    form : ufl.Form
        The complete form from the expression.
    expression : ufl.Form
        The underlying UFL expression.

    Notes
    -----
    Operator overloading is supported (+, +=, -, -=).
    """

    def __init__(self, inputUfl = None):
        """
        Initialize the UflDecorator.

        Parameters
        ----------
        inputUfl : ufl.Form or None, optional
            The initial UFL expression to store. If not provided, the object is empty.
        """

        if inputUfl != None:
            self._expression = inputUfl

    def add(self, input_ufl):
        """
        Add a UFL expression to the current expression.

        Parameters
        ----------
        input_ufl : ufl.Form
            The UFL expression to be added.

        Notes
        -----
        If no expression has been set yet, this initializes the expression.
        """
        if hasattr(self,'_expression'):
            self._expression += input_ufl
        else:
            self._expression = input_ufl

    def subtract(self,input_ufl):
        """
        Subtract a UFL expression from the current expression.

        Parameters
        ----------
        input_ufl : ufl.Form
            The UFL expression to be subtracted.

        Notes
        -----
        If no expression has been set yet, this initializes the expression to the negative of the input.
        """
        if hasattr(self,'_expression'):
            self._expression -= input_ufl
        else:
            self._expression = -input_ufl

    @property
    def expression(self):
        """
        Return the current UFL expression.

        Returns
        -------
        ufl.Form
            The stored UFL expression.
        """
        return self._expression

    @property
    def lhs(self):
        """
        Return the left-hand side of the UFL expression.

        Returns
        -------
        ufl.Form
            The left-hand side form of the expression.
        """
        return lhs(self._expression)

    @property
    def rhs(self):
        """
        Return the right-hand side of the UFL expression.

        Returns
        -------
        ufl.Form
            The right-hand side form of the expression.
        """
        return rhs(self._expression)

    @property
    def form(self):
        """
        Return the general form representation of the UFL expression.

        Returns
        -------
        ufl.Form
            The UFL form of the expression.
        """
        return form(self._expression)


    def getAssembledMatrix(self, mesh, bcs = []):
        """
        Assemble the left-hand side form into a matrix.

        Parameters
        ----------
        mesh : FELiCS.SpaceDisc.FELiCSMesh
            The mesh on which to assemble.
        bcs : list of dolfinx.fem.dirichletbc.DirichletBC, optional
            List of boundary conditions to apply.

        Returns
        -------
        petsc4py.PETSc.Mat
            The assembled PETSc matrix.
        """
        self._setCorrectMeshObject(mesh)
        matrix = assemble_matrix(form(self.lhs), bcs=bcs)
        matrix.assemble()
        return matrix

    def getAssembledVector(self, mesh, bcs = []):
        """
        Assemble the right-hand side form into a vector.

        Parameters
        ----------
        mesh : FELiCS.SpaceDisc.FELiCSMesh
            The mesh on which to assemble.
        bcs : list of dolfinx.fem.dirichletbc.DirichletBC, optional
            List of boundary conditions to apply.

        Returns
        -------
        petsc4py.PETSc.Vec
            The assembled PETSc vector.
        """
        self._setCorrectMeshObject(mesh)
        vector = assemble_vector(form(-self.rhs))
        vector.assemble()
        set_bc(vector, bcs)
        return vector

    def getAssembledScalar(self, mesh):
        """
        Assemble the expression into a scalar.

        Parameters
        ----------
        mesh : FELiCS.SpaceDisc.FELiCSMesh
            The mesh on which to assemble.

        Returns
        -------
        float
            The scalar result of the assembly.

        """
        self._setCorrectMeshObject(mesh)
        return assemble_scalar(form(self._expression))

    def lhsIsZero(self):
        """
        Check whether the left-hand side of the expression is effectively zero.

        Returns
        -------
        bool
            True if lhs is zero or has insufficient arguments; otherwise, False.
        """
        if self.isZero():
            return True
        else:
            temp=lhs(self._expression)
            if len(temp.arguments())<2:
                return True
            else:
                return False

    def rhsIsZero(self):
        """
        Check whether the right-hand side of the expression is effectively zero.

        Returns
        -------
        bool
            True if rhs is zero or has insufficient arguments; otherwise, False.
        """
        if self.isZero():
            return True
        else:
            temp=rhs(self._expression)
            if len(temp.arguments())<1:
                return True
            else:
                return False

    def isZero(self):
        """
        Check whether the expression has been initialized.

        Returns
        -------
        bool
            True if the expression is not set; False otherwise.
        """
        if hasattr(self,'__weakForm__'):
            return False
        else:
            return True
        
    def printExpression(self):
        """
        Print a summary of the current UFL expression.

        Prints the expression, its arguments, and number of arguments.
        """
        # TODO: change this format? uses standard print.
        if hasattr(self,'_expression'):
            print('###### ufl expression:')
            print(self._expression)
            print('###### arguments:')
            print(self._expression.arguments())
            print('###### Number of arguments:')
            print(len(self._expression.arguments()))
        else:
            print('Ufl expression is zero.')


    def _setCorrectMeshObject(self, mesh):
        """
        Internal workaround to ensure mesh compatibility in UFL expressions.

        

        Parameters
        ----------
        mesh : dolfinx.mesh.Mesh
            The mesh object to be patched into the UFL expression.

        Notes
        -----
        Mesh compatibility workaround
        Sophie: This is a weird work-around, because somehow the wrong mesh object is given to the UFL-form 
        when using a newer version of dolfinx (version >= 0.6.*).
        I will try and understand why that is (probably has something to do with the class FelicsMesh?), 
        but for now this works fine. 
        """
        sd = self._expression.subdomain_data()
        domain, = list(sd.keys())  # Assuming single domain
        domain._ufl_cargo = mesh._cpp_object
 

    def __add__(self, other):
        """
        Add another UFL expression or UflDecorator to this one.

        Parameters
        ----------
        other : UflDecorator or ufl.Form
            The expression to add.

        Returns
        -------
        UflDecorator
            A new UflDecorator with the combined expression.

        Notes
        -----
        Does not currently validate that `other` is a valid UFL expression.
        """
        if hasattr(self, "_expression"):
            result = UflDecorator(self._expression)
        else: 
            result = UflDecorator()
        if isinstance(other, UflDecorator):
            result.add(other._expression)
        else:
            result.add(other)
        return result

    def __iadd__(self, other):
        """
        In-place addition of another UFL expression or UflDecorator.

        Parameters
        ----------
        other : UflDecorator or ufl.Form
            The expression to add.

        Returns
        -------
        UflDecorator
            The modified instance.

        Notes
        -----
        Does not currently validate that `other` is a valid UFL expression.
        """
        if isinstance(other, UflDecorator):
            self.add(other._expression)
        else:
            self.add(other)
        return self

    def __sub__(self, other):
        """
        Subtract another UFL expression or UflDecorator from this one.

        Parameters
        ----------
        other : UflDecorator or ufl.Form
            The expression to subtract.

        Returns
        -------
        UflDecorator
            A new UflDecorator with the updated expression.

        Notes
        -----
        Does not currently validate that `other` is a valid UFL expression.
        """
        if hasattr(self, "_expression"):
            result = UflDecorator(self._expression)
        else: 
            result = UflDecorator()
        if isinstance(other, UflDecorator):
            result.subtract(other._expression)
        else:
            result.subtract(other)
        return result

    def __isub__(self, other):
        """
        In-place subtraction of another UFL expression or UflDecorator.

        Parameters
        ----------
        other : UflDecorator or ufl.Form
            The expression to subtract.

        Returns
        -------
        UflDecorator
            The modified instance.

        Notes
        -----
        Does not currently validate that `other` is a valid UFL expression.
        """
        if isinstance(other, UflDecorator):
            self.subtract(other._expression)
        else:
            self.subtract(other)
        return self


