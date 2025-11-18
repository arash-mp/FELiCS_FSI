from FELiCS.Misc.tensorUtils import*

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
from dolfinx             import mesh
from FELiCS.Fields.Field import Field
from FELiCS.Fields.Mode  import Mode

class RandomFieldsHandle():
    # TODO: write docstrings
    def __init__(self,dim_vector, ):
        # 1. Define mesh 
        self.mesh = mesh.create_rectangle(comm=MPI.COMM_WORLD,
                            points=((0.0, 0.0), (1.0, 1.0)), n=(10, 10),
                            cell_type=mesh.CellType.triangle,
                            ghost_mode=mesh.GhostMode.none)
        # 2. create function spaces & test/trial functions
        self.dim_vector = dim_vector
        self.space_scalar = functionspace(self.mesh, ("CG", 2))
        self.space_vector = functionspace(self.mesh, ("CG", 2,(dim_vector,)))

        self.test_scalar  = TestFunction(self.space_scalar)
        self.test_vector  = TestFunction(self.space_vector)

        self.trial_scalar = TrialFunction(self.space_scalar)
        self.trial_vector = TrialFunction(self.space_vector)
        
        # 3. create random conplex number
        self.random_float = np.random.random() + 1j*np.random.random()
        
    def createDolfinxFunction(self, dim):
        if dim == 1:
            func_sclar = Function(self.space_scalar)
            length = len(func_sclar.x.array[:])
            func_sclar.x.array[:]  = np.random.rand(length) + 1j*np.random.rand(length)
            return func_sclar
        elif dim == self.dim_vector:
            func_vector = Function(self.space_vector)
            length = len(func_vector.x.array[:])
            func_vector.x.array[:] = np.random.rand(length)+ 1j*np.random.rand(length)
            return func_vector
        else:
            raise ValueError("Dimension is not consistent with the class definition.")
    pass
    
    def createFELiCSField(self):
        # TODO: Add implementation
        pass
    
    def smoothing(self, func, smoothFactor):
        if func.function_space == self.space_scalar:
            # Smooth scalar functions
            rhs = func * conj(self.test_scalar)*dx
            lhs = conj(self.test_scalar)*self.trial_scalar * dx + \
                smoothFactor * inner(grad(self.trial_scalar),grad(self.test_scalar))*dx
            problem = petsc.LinearProblem(lhs, rhs, bcs=[],
                                            petsc_options={"ksp_type": "preonly",
                                                            "pc_type": "lu"})
            temp_func = problem.solve()
            return temp_func
        elif func.function_space == self.space_vector:
            # Smooth vector functions
            rhs = inner(func,  self.test_vector)*dx
            lhs = inner(self.trial_vector, self.test_vector)*dx + smoothFactor*inner(grad(self.trial_vector),grad(self.test_vector))*dx
            problem = petsc.LinearProblem(lhs, rhs, bcs=[], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
            temp_func = problem.solve()
            return temp_func
        pass
    pass
    