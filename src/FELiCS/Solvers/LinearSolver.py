#Standard libraries
import time
import os
import sys
import multiprocessing

# Third party libraries
import numpy as np

from dolfinx.fem import (
        Function,
        dirichletbc,
        form,
)

from dolfinx.fem.petsc import (
    assemble_vector,
)

from ufl import (
    dx,
    TestFunctions,
    SpatialCoordinate,
)

from functools import partial

#Local libraries and methods
from FELiCS.Misc.functions import (
    printError,
    printWarning,
    printDebug,
    )




class LinearSolver:


    @staticmethod
    def solveGeneralEigenproblem(
                A, 
                B, 
                sigma, 
                nev, 
                tol=1.e-12, 
                max_it=200, 
                adjoint=False, 
                isForEigenProblem=True, 
                ):

        """
        Solves the generalized eigenvalue problem (GEVP) 
        using the SLEPc and PETSc libraries.

        Function arguments:
        - A, B      Matrices defining the GEVP (A-wB)q = 0
        - sigma     Eigenvalue guess
        - nev       Number of eigenvalues to compute
        - tol           (optional) precision of GEVP
        - max_it    (optional) maximum number of iterations 
        - adjoint   (bool, optional) if set to "True", compute the adjoint GEVP alongside the direct one
        - isForEigenProblem, isIncompressible       (optional) flags that can be set to determine the kind of matrices  

        Function returns:
        - eigVals           numpy array; computed eigenvalues
        - eigVecs           numpy array, 2D; computed eigenvectors
        - eigVecs_adjoint   numpy array, 2D, only returned if "adjoint=True"; computed left-hand-side eigenvectors
        - error             numpy array; for each eigenproblem solution, relative residuum of the eigenproblem
        """

        from slepc4py import SLEPc

        # finish assembling matrices 
        A.assemble()
        B.assemble()

        if adjoint:
            A.hermitianTranspose()
            B.hermitianTranspose()
            guess = np.conj(sigma)
        else:
            guess = sigma

        # create eigenproblem solver 
        eps = SLEPc.EPS().create()
        eps.setOperators(A,B)
        eps.setProblemType(SLEPc.EPS.ProblemType.GNHEP)     # general non-Hermitian eigenproblem with semi-definite B
        
        eps.setTolerances(tol=tol,max_it=max_it)
        
        #eps.setType(SLEPc.EPS.Type.KRYLOVSCHUR) # is standard, does not need to be set
        eps.getST().setType(SLEPc.ST.Type.SINVERT)
        eps.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)
        eps.setTarget(guess)
        eps.setDimensions(nev=nev)
        
        eps.getST().getKSP().getPC().setType('lu')
        eps.getST().getKSP().getPC().setFactorSolverType('mumps')
        
        eps.setFromOptions()
        eps.setUp()
        
        eps.solve()
        
        dim = A.getSize()[0]
        eigVals,  eigVecs  = np.empty(nev,complex), np.empty([nev,dim],complex)
        vec_real, vec_imag = A.getVecs()
        error = np.empty(nev)
        for i in range(nev):
            try:
                eigVals[i]   = eps.getEigenpair(i,vec_real,vec_imag)
                eigVecs[i,:] = vec_real.getArray() + 1j * vec_imag.getArray()
                error[i] = eps.computeError(i, SLEPc.EPS.ErrorType.RELATIVE)
                #printDebug(True,f"SLEPc error relative: {error[i]}")
            
            except:
                printWarning("Could not access eigenpair nb ", nev+1, "!")
        
        eps.getST().getKSP().getPC().destroy()
        eps.getST().getKSP().destroy()
        eps.getST().destroy()
        eps.destroy()
        return eigVals, eigVecs, error
    

    @staticmethod
    def solveSVDOfResolvent(
                resolventOperator, 
                nev,
                tol    = 1.e-16,
                max_it = 200):

        from slepc4py import SLEPc
        from petsc4py import PETSc
        
        # R is matrix free, needs extra class
        R = PETSc.Mat().createPython(resolventOperator.getSize())
        R.setPythonContext(resolventOperator)
        R.setUp()
        
        # create eigenproblem solver 
        eps = SLEPc.EPS().create()
        eps.setOperators(R)
        eps.setDimensions(nev=nev)
        eps.setTolerances(tol=tol,max_it=max_it)
        eps.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_MAGNITUDE)
        eps.getST().getKSP().getPC().setType('none')
        eps.setFromOptions()
        eps.setUp()
        
        eps.solve()
        
        dim = resolventOperator.getSize()[0]
        eigVals,  eigVecs  = np.empty(nev,complex), np.empty([dim,nev],complex)
        vec_real, vec_imag = resolventOperator.getVecs()
        for i in range(nev):
            try:
                eigVals[i] = eps.getEigenpair(i,vec_real,vec_imag)
                eigVecs[:,i] = vec_real.getArray() + 1j * vec_imag.getArray()
                printDebug(True,f"SLEPc error relative: {eps.computeError(i, SLEPc.EPS.ErrorType.RELATIVE)}")
                #printDebug(True,f"SLEPc error absolute: {eps.computeError(i, SLEPc.EPS.ErrorType.ABSOLUTE)}")
            except:
                printWarning("Could not access eigenpair nb ", nev+1, "!")
        
        eps.destroy()
        R.destroy()
        return eigVals, eigVecs



    @staticmethod
    def solveEquationSystem(
        A,
        b):

        """__solveEquationSystem  
        Solves a linear equation system Ax=b, using the PETSc libraries.
        
        Parameters
        ----------
        A : PETSc matrix
            matrix of the linear syste
        b : PETSc vector
            rhs of the linear system
        
        Returns
        -------
        x:  numpy array
            solution of the linear equation system
        """
        from petsc4py import PETSc
        solution,dummy = A.createVecs()
        
        solver = PETSc.KSP().create()
        solver.setOperators(A)
        solver.setType(PETSc.KSP.Type.PREONLY)
        solver.getPC().setType(PETSc.PC.Type.LU)
        solver.getPC().setFactorSolverType('mumps')
        
        solver.solve(b, solution)
        
        x = solution.getArray()
        
        solver.destroy()
        dummy.destroy()
        
        return x 


    @staticmethod
    def solveTransposeEquationSystem(
        A,
        b):

        """
        Solves a linear equation system A^T x=b, using the PETSc libraries.
        
        Parameters
        ----------
        A : PETSc matrix
            matrix of the linear syste
        b : PETSc vector
            rhs of the linear system
        
        Returns
        -------
        x:  numpy array
            solution of the linear equation system
        """
        from petsc4py import PETSc
        solution,dummy = A.createVecs()
        
        solver = PETSc.KSP().create()
        solver.setOperators(A)
        solver.setType(PETSc.KSP.Type.PREONLY)
        solver.getPC().setType(PETSc.PC.Type.LU)
        solver.getPC().setFactorSolverType('mumps')
        
        solver.solveTranspose(b, solution)
        
        x = solution.getArray()
        
        solver.destroy()
        dummy.destroy()
        
        return x 


    @staticmethod
    def createEquationSystemSolver(
        A):

        """__solveEquationSystem  
        Creates a KSP petsc solver to solve a linear equation system. This is useful if several linear equation systems with the same matrix are solved, 
        since it stores the preconditioner and the calculation time is significantly reduced.
        To solve the equation system use the method "solveEquationSystemWithPredefinedSolver".
        
        Parameters
        ----------
        A : PETSc matrix
            matrix of the linear system
        
        Returns
        -------
        solver:  PETSc KSP solver
                 solver for the given matrix 
        """
        from petsc4py import PETSc
        
        solver = PETSc.KSP().create()
        solver.setOperators(A)
        solver.setType(PETSc.KSP.Type.PREONLY)
        solver.getPC().setType(PETSc.PC.Type.LU)
        solver.getPC().setFactorSolverType('mumps')
        
        return solver 


    @staticmethod
    def solveEquationSystemWithPredefinedSolver(
        solver,
        b):

        """__solveEquationSystem  
        Solves a linear equation system Ax=b, using the PETSc libraries.
        
        Parameters
        ----------
        solver : PETSc KSP solver
                 created with the method "createEquationSystemSolver"
        b : PETSc vector
            rhs of the linear system
        
        Returns
        -------
        x:  numpy array
            solution of the linear equation system
        """
        from petsc4py import PETSc
        solution = b.copy()
        
        solver.solve(b, solution)
        
        x = solution.getArray()
        
        return x 






class ResolventOperator(object):

        """
        This class serves as a "matrix-free" representation of the Resolvent operator multiplicated with its Hermitian transposed, to conduct the singular value decomposition of the system. 
        It contains a method called "mult", which is called by the eigenvalue solver, and returns a matrix vector product of the represented matrix. 

        Private attributes:
        Protected attributes:
        Public attributes:
        """

        def __init__(self,
                OP,
                Qf,
                Pu,
                Cr,
                B_forcing,
                B_response):


                from petsc4py import PETSc

                self._size = Qf.getSize()

                self._Z1, self._Z2  = OP.getVecs()
                self._Y1, self._Y2  = Qf.getVecs()

                self._Pu = Pu
                self._Cr = Cr
                self._Bf = B_forcing
                self._Br = B_response

                # create KSP1: This is a solver for the System OP*x=y.
                self._ksp1 = PETSc.KSP().create()
                self._ksp1.setOperators(OP)
                self._ksp1.setType(PETSc.KSP.Type.PREONLY)
                self._ksp1.getPC().setType(PETSc.PC.Type.LU)
                self._ksp1.getPC().setFactorSolverType('mumps')
                self._ksp1.setUp()
                # create KSP2: This is a solver for the System conj(OP)*x=y. It will later be used to solve the transposed system, thus effectively solving OP^H *x=y, which is the Hermitian transpose of the system.
                # TODO Sophie: unfortuantely there is no "solveHermitianTranspose" in the petsc4py (yet?). Thus we have to do an additional LU decomposistion.... Change as soon as this is included in the petsc4py! 
                OP_H = OP.copy()   #create a new matrix, s.t. the original one will not be overwritten 
                OP_H.conjugate()
                OP_H.assemble()
                self._ksp2 = PETSc.KSP().create()
                self._ksp2.setOperators(OP_H)
                self._ksp2.setType(PETSc.KSP.Type.PREONLY)
                self._ksp2.getPC().setType(PETSc.PC.Type.LU)
                self._ksp2.getPC().setFactorSolverType('mumps')
                self._ksp2.setUp()
                # create KSP3: This is a solver for the System Qf*x=y (will later be used to solve the transposed system).
                #Qf.conjugate() #=> is this needed?
                self._ksp3 = PETSc.KSP().create()
                self._ksp3.setOperators(Qf)
                self._ksp3.setType(PETSc.KSP.Type.PREONLY)
                self._ksp3.getPC().setType(PETSc.PC.Type.LU)
                self._ksp3.getPC().setFactorSolverType('mumps')
                self._ksp3.setUp()

        def getSize(self):
                return self._size

        def getVecs(self):
                return self._Y1, self._Y2

        def mult(self, mat, X, Y):
                # returns Y=mat*X 
                # mat = (Qf^T)^-1 * Pu^T * Bf^T * (OP^H)^-1 * Cr^T * Br * Cr * OP^-1 * Bf * Pu

                self._Pu.mult            (X,        self._Z1)  #Z1 = Pu*X
                self._Bf.mult            (self._Z1, self._Z2)  #Z2 = Bf*Z1
                self._ksp1.solve         (self._Z2, self._Z1)  #Z1 = OP^-1 * Z2
    
                self._Cr.mult            (self._Z1, self._Z2)  #Z2 = Cr*Z1
                self._Br.mult            (self._Z2, self._Z1)  #Z1 = Br*Z2
                self._Cr.multTranspose   (self._Z1, self._Z2)  #Z2 = Cr^T * Z1
                self._ksp2.solveTranspose(self._Z2, self._Z1)  #Z1 = (OP^H)^-1 * Z2

                self._Bf.multTranspose   (self._Z1, self._Z2)  #Z2 = Bf^T * Z1
                self._Pu.multTranspose   (self._Z2, self._Y1)  #Y1 = Pu^T * Z2
                self._ksp3.solveTranspose(self._Y1, Y)         #Y  = (Qf^T)^-1 * Y1

                return Y

        def getKSP(self):
                return self._ksp1

        def destroySelf(self):
                """Clean-up the disc space to avoid memory leaks. Should be called if several resolvent SVDs are done one after the other.""" 
                self._ksp1.getPC().destroy()
                self._ksp2.getPC().destroy()
                self._ksp3.getPC().destroy()
                self._ksp1.destroy()
                self._ksp2.destroy()
                self._ksp3.destroy()





