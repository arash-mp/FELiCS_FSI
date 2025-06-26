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

    def __init__(self, inputUfl = None):
        if inputUfl != None:
            self._expression = inputUfl

    def add(self, input_ufl):
        if hasattr(self,'_expression'):
            self._expression += input_ufl
        else:
            self._expression = input_ufl

    def subtract(self,input_ufl):
        if hasattr(self,'_expression'):
            self._expression -= input_ufl
        else:
            self._expression = -input_ufl

    @property
    def expression(self):
        return self._expression

    @property
    def lhs(self):
        return lhs(self._expression)

    @property
    def rhs(self):
        return rhs(self._expression)

    @property
    def form(self):
        return form(self._expression)


    def __add__(self, other):
        # Note (Sophie): This method does not check if "other" is actually an ufl expression. 
        # This could be added to avoid mistakes.
        result = UlfDecorator(self._expression)
        if isinstance(other, UflDecorator):
            result.add(other._expression)
        else:
            result.add(other)
        return result

    def __iadd__(self, other):
        # Note (Sophie): This method does not check if "other" is actually an ufl expression. 
        # This could be added to avoid mistakes.
        if isinstance(other, UflDecorator):
            self.add(other._expression)
        else:
            self.add(other)
        return self

    def __sub__(self, other):
        # Note (Sophie): This method does not check if "other" is actually an ufl expression. 
        # This could be added to avoid mistakes.
        result = UlfDecorator(self._expression)
        if isinstance(other, UflDecorator):
            result.substract(other._expression)
        else:
            result.substract(other)
        return result

    def __isub__(self, other):
        # Note (Sophie): This method does not check if "other" is actually an ufl expression. 
        # This could be added to avoid mistakes.
        if isinstance(other, UflDecorator):
            self.substract(other._expression)
        else:
            self.substract(other)
        return self

    def lhsIsZero(self):
        if self.isZero():
            return True
        else:
            temp=lhs(self._expression)
            if len(temp.arguments())<2:
                return True
            else:
                return False

    def rhsIsZero(self):
        if self.isZero():
            return True
        else:
            temp=rhs(self._expression)
            if len(temp.arguments())<1:
                return True
            else:
                return False

    def isZero(self):
        if hasattr(self,'__weakForm__'):
            return False
        else:
            return True
    def analyseExpression(self):
        if hasattr(self,'_expression'):
            print('###### ufl expression:')
            print(self._expression)
            print('###### arguments:')
            print(self._expression.arguments())
            print('###### Number of arguments:')
            print(len(self._expression.arguments()))
        else:
            print('Ufl expression is zero.')


    def getAssembledMatrix(self, mesh, bcs = []):
        self._setCorrectMeshObject(mesh)
        matrix = assemble_matrix(form(self.lhs), bcs=bcs)
        matrix.assemble()
        return matrix

    def getAssembledVector(self, mesh, bcs = []):
        self._setCorrectMeshObject(mesh)
        vector = assemble_vector(form(-self.rhs))
        vector.assemble()
        set_bc(vector, bcs)
        return vector

    def getAssembledScalar(self, mesh):
        self._setCorrectMeshObject(mesh)
        return assemble_scalar(form(ufl_expression))

    def _setCorrectMeshObject(self, mesh):
        # Mesh compatibility workaround
        # Sophie: This is a weird work-around, because somehow the wrong mesh object is given to the UFL-form 
        # when using a newer version of dolfinx (version >= 0.6.*).
        # I will try and understand why that is (probably has something to do with the class FelicsMesh?), 
        # but for now this works fine. 
        try:
            sd = self._expression.subdomain_data()
            domain, = list(sd.keys())  # Assuming single domain
            domain._ufl_cargo = mesh._cpp_object._cpp_object
        except:
            logger.info("DEPRECATED: Mesh module from dolfinx version <0.7.0 is used.")
 

