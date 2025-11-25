import numpy as np
from mpi4py              import MPI
from petsc4py.PETSc      import ScalarType
from ufl                 import (
    dx,TestFunction,TrialFunction,conj,SpatialCoordinate, Dx,dot, inner, grad
)
from dolfinx.fem         import (
    Function,
    functionspace,
    dirichletbc,
    Constant,
    form,
    locate_dofs_topological,
    petsc
)

class UnitTestHelper():
    def __init__(self):
        pass
    
    @staticmethod
    def checkVectorExpressionAlignment(tensor_expr, valid_expr, tol=1.e-14):
        res1 = petsc.assemble_vector(form(tensor_expr))
        res2 = petsc.assemble_vector(form(valid_expr))
        res1.assemble()
        res2.assemble()
        
        print("   - Tensor expression vector norm: ", np.linalg.norm(res1.getArray()))
        print("   - Valid expression vector norm : ", np.linalg.norm(res2.getArray()))
        print("   - Difference norm              : ", np.linalg.norm(res1.getArray()-res2.getArray()))
        assert np.linalg.norm(res1.getArray())> 0
        assert np.linalg.norm(res2.getArray())> 0
        assert np.linalg.norm(res1.getArray()-res2.getArray()) < tol
        pass
    
    @staticmethod
    def getVecGrad(func_vector: Function, coordinateSystemName, m, r):
        t11 = Dx(func_vector[0],0)
        t12 = Dx(func_vector[0],1)
        t21 = Dx(func_vector[1],0)
        t22 = Dx(func_vector[1],1)
        t31 = Dx(func_vector[2],0)
        t32 = Dx(func_vector[2],1)
        if coordinateSystemName == "cartesian": 
            t13 = 1j*m*func_vector[0]
            t23 = 1j*m*func_vector[1]
            t33 = 1j*m*func_vector[2]
        elif coordinateSystemName == "cylindricalfelics": 
            t13 = 1j*m*func_vector[0]/r
            t23 = 1j*m*func_vector[1]/r - func_vector[2]/r
            t33 = 1j*m*func_vector[2]/r + func_vector[1]/r
        return [[t11, t12, t13], [t21, t22, t23], [t31, t32, t33]]
    
    @staticmethod
    def getScalarGrad(func_scalar: Function, coordinateSystemName, m, r):
        s1 = Dx(func_scalar,0)
        s2 = Dx(func_scalar,1)
        if coordinateSystemName == "cartesian": 
            s3 = 1j*m*func_scalar
        elif coordinateSystemName == "cylindricalfelics": 
            s3 = 1j*m*func_scalar/r
        return [s1, s2, s3]
        
