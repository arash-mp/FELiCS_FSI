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
#Standard libraries
# import time
# import os
# import sys
# import multiprocessing

# Third party libraries
import numpy as np

# from dolfinx.fem import (
#         Function,
#         dirichletbc,
#         form,
# )

# from dolfinx.fem.petsc import (
#     assemble_vector,
# )

# from ufl import (
#     dx,
#     TestFunctions,
#     SpatialCoordinate,
# )

# from functools import partial

#Local libraries and methods
from FELiCS.Misc.logging import Logger

# Get the logger
logger = Logger.get_logger("felics")

class LinearSolver:
    """
    Linear algebra utilities using PETSc and SLEPc.

    Provides static methods for solving generalized eigenvalue problems (GEVP),
    singular value decompositions (SVD), and linear systems efficiently using
    PETSc and SLEPc.

    **Initialize the LinearSolver object**

    Parameters
    ----------
    None

    Notes
    -----
    Consists of only static methods that don't need an instance of this class. Call the methods via "LinearSolver.method()".

    Examples
    --------
    Solving a linear system:

    >>> from petsc4py import PETSc
    >>> A = PETSc.Mat().create()
    >>> b = PETSc.Vec().create()
    >>> x = LinearSolver.solveEquationSystem(A, b)

    Solving a generalized eigenvalue problem:

    >>> eigVals, eigVecs, error = LinearSolver.solveGeneralEigenproblem(
    ...     A, B, sigma=0.0, nev=5)
    """

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
        Solve the generalized eigenvalue problem (GEVP) using SLEPc and PETSc.

        Parameters
        ----------
        A : PETSc.Mat
            Matrix defining the eigenvalue problem.
        B : PETSc.Mat
            Second matrix in the generalized eigenvalue problem.
        sigma : float
            Initial guess for the eigenvalue.
        nev : int
            Number of eigenvalues to compute.
        tol : float, optional
            Tolerance for convergence, by default 1.e-12.
        max_it : int, optional
            Maximum number of iterations, by default 200.
        adjoint : bool, optional
            Compute adjoint eigenvalues and eigenvectors if True, by default False.
        isForEigenProblem : bool, optional
            Flag indicating whether this is for an eigenproblem, by default True.

        Returns
        -------
        eigVals : numpy.ndarray
            Computed eigenvalues.
        eigVecs : numpy.ndarray
            Computed eigenvectors.
        error : numpy.ndarray
            Relative residuals for the eigenvalue problem.
        """

        from slepc4py import SLEPc

        # finish assembling matrices 
        A.assemble()
        B.assemble()

        if adjoint:
            # A.hermitianTranspose()
            # B.hermitianTranspose()
            A_adj = A.copy()
            A_adj.hermitianTranspose()
            B_adj = B.copy()
            B_adj.hermitianTranspose()
            guess = np.conj(sigma)
        else:
            guess = sigma

        # create eigenproblem solver 
        eps = SLEPc.EPS().create()
        if adjoint:
            eps.setOperators(A_adj,B_adj)
        else:
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
                logger.warning("Could not access eigenpair nb ", nev+1, "!")
        
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
        """
        Solve the singular value decomposition (SVD) of a resolvent operator.

        Parameters
        ----------
        resolventOperator : ResolventOperator
            Resolvent operator object for which the SVD is computed.
        nev : int
            Number of singular values to compute.
        tol : float, optional
            Tolerance for convergence, by default 1.e-16.
        max_it : int, optional
            Maximum number of iterations, by default 200.

        Returns
        -------
        eigVals : numpy.ndarray
            Singular values of the operator.
        eigVecs : numpy.ndarray
            Singular vectors corresponding to the singular values.
        """

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
                logger.debug(f"SLEPc error relative: {eps.computeError(i, SLEPc.EPS.ErrorType.RELATIVE)}")
            except:
                logger.warning("Could not access eigenpair nb " + str(nev+1) +"!")
        
        eps.destroy()
        R.destroy()
        return eigVals, eigVecs

    @staticmethod
    def solveEquationSystem(
        A,
        b,
        destroy=False):

        """
        Solve a linear system of equations Ax = b using PETSc.

        Parameters
        ----------
        A : PETSc.Mat
            Matrix representing the linear system.
        b : PETSc.Vec
            Right-hand side of the equation.
        destroy : bool, optional
            Whether to destroy the matrix and vector after solving, by default False.

        Returns
        -------
        x : numpy.ndarray
            Solution vector.
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
        
        solver.getPC().destroy()
        solver.destroy()
        dummy.destroy()
        solution.destroy()

        if destroy:
            A.destroy()
            b.destroy()
        
        return x 


    @staticmethod
    def solveTransposeEquationSystem(
        A,
        b,
        destroy=False):

        """
        Solve the transpose of a linear system A^T x = b using PETSc.

        Parameters
        ----------
        A : PETSc.Mat
            Matrix representing the linear system.
        b : PETSc.Vec
            Right-hand side of the equation.
        destroy : bool, optional
            Whether to destroy the matrix and vector after solving, by default False.

        Returns
        -------
        x : numpy.ndarray
            Solution vector.
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
        
        solver.getPC().destroy()
        solver.destroy()
        dummy.destroy()
        solution.destroy()

        if destroy:
            A.destroy()
            b.destroy()
        
        return x 


    @staticmethod
    def createEquationSystemSolver(
        A):

        """
        Create a KSP PETSc solver to solve a linear equation system.

        This is useful if several linear equation systems with the same matrix are solved, since it stores the preconditioner
        and the calculation time is significantly reduced. To solve the equation system use the method
        "solveEquationSystemWithPredefinedSolver".

        Parameters
        ----------
        A : PETSc.Mat
            Matrix of the linear system.

        Returns
        -------
        solver : PETSc.KSP
            Preconfigured solver for the given matrix.
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
        b,
        destroy=False):

        """
        Solve a linear equation system Ax=b using PETSc and a predefined solver.

        Parameters
        ----------
        solver : PETSc.KSP
            Preconfigured solver, can be created with the method "createEquationSystemSolver".
        b : PETSc.Vec
            Right-hand side of the equation.
        destroy : bool, optional
            Whether to destroy the vector after solving, by default False. Can be useful by repetitive computations to avoid memory leaks.
        
        Returns
        -------
        x : numpy.ndarray
            Solution vector of the linear equation system.
        """
        from petsc4py import PETSc
        solution = b.copy()
        
        solver.solve(b, solution)
        
        x = solution.getArray()

        solution.destroy()

        if destroy:
            b.destroy()
        
        return x 


    @staticmethod
    def solveTransposeEquationSystemWithPredefinedSolver(
        solver,
        b,
        destroy=False):

        """
        Solve a linear equation system A^Tx=b using PETSc and a predefined solver (defined by matrix A).

        Parameters
        ----------
        solver : PETSc.KSP
            Preconfigured solver, can be created with the method "createEquationSystemSolver".
        b : PETSc.Vec
            Right-hand side of the equation.
        destroy : bool, optional
            Whether to destroy the vector after solving, by default False. Can be useful by repetitive computations to avoid memory leaks.
        
        Returns
        -------
        x : numpy.ndarray
            Solution vector of the linear equation system.
        """
        from petsc4py import PETSc
        solution = b.copy()
        
        solver.solveTranspose(b, solution)
        
        x = solution.getArray()

        solution.destroy()

        if destroy:
            b.destroy()
        
        return x 



class ResolventOperator(object):         
    """
    Matrix-free representation of the Resolvent operator multiplied with its conjugate transpose.
    Used to conduct the singular value decomposition of the system for resolvent analysis. Contains a method called "mult",
    which is called by the eigenvalue solver and returns a matrix-vector product of the represented matrix. The resolvent operator
    requires additional full-size quadratic matrices to calculate the forcing and response norms, as well as (possibly rectangular)
    restrictor matrices that limit the spatial domain and variable dimensions.

    **Initialize the ResolventOperator object**

    Parameters
    ----------
    ResolventOperator : PETSc.Mat
        Matrix representing the linear system.
    FEMWeightMatrix_fullSystem : PETSc.Mat
        Weight matrix for the full FEM system.
    FEMWeightMatrix_forcingNorm : PETSc.Mat
        Weight matrix for the forcing norm. Has default size of full system (surplus DOFs will be ignored).
    FEMWeightMatrix_responseNorm : PETSc.Mat
        Weight matrix for the response norm. Has default size of full system (surplus DOFs will be ignored).
    RestrictorMatrix_forcing : PETSc.Mat
        Restrictor matrix for forcing. Rectangular matrix of appropriate size without FEM weights. Spatial restrictor values can be between 0 and 1.
    RestrictorMatrix_response : PETSc.Mat
        Restrictor matrix for response. Rectangular matrix of appropriate size without FEM weights. Spatial restrictor values can be between 0 and 1.
    """

    def __init__(self,                      
            ResolventOperator, #A-i*omega*B                      
            FEMWeightMatrix_fullSystem,                      
            FEMWeightMatrix_forcingNorm,                      
            FEMWeightMatrix_responseNorm,  
            RestrictorMatrix_forcing,                      
            RestrictorMatrix_response):    

        """
        Initialize the ResolventOperator object.

        Parameters
        ----------
        ResolventOperator : PETSc.Mat
            Matrix representing the linear system.
        FEMWeightMatrix_fullSystem : PETSc.Mat
            Weight matrix for the full FEM system.
        FEMWeightMatrix_forcingNorm : PETSc.Mat
            Weight matrix for the forcing norm.
        FEMWeightMatrix_responseNorm : PETSc.Mat
            Weight matrix for the response norm.
        RestrictorMatrix_forcing : PETSc.Mat
            Restrictor matrix for forcing.
        RestrictorMatrix_response : PETSc.Mat
            Restrictor matrix for response.
        """               

        from petsc4py import PETSc                 

        self._size = FEMWeightMatrix_forcingNorm.getSize()    #: :meta private:            

        self._P_forcing  = RestrictorMatrix_forcing                 
        self._P_response = RestrictorMatrix_response                 
        self._W_FEM      = FEMWeightMatrix_fullSystem                 
        self._W_response = FEMWeightMatrix_responseNorm                 
        self._W_forcing  = FEMWeightMatrix_forcingNorm             

        self._O1, self._O2  = ResolventOperator.getVecs()                 
        self._R1, self._R2  = self._W_response.getVecs()                 
        self._F1, self._F2  = self._W_forcing.getVecs()                 

        # create KSP1: This is a solver for the System Operator*x=y.                 
        self._ksp1 = PETSc.KSP().create()                 
        self._ksp1.setOperators(ResolventOperator)                 
        self._ksp1.setType(PETSc.KSP.Type.PREONLY)                 
        self._ksp1.getPC().setType(PETSc.PC.Type.LU)                 
        self._ksp1.getPC().setFactorSolverType('mumps')                 
        self._ksp1.setUp()                 
        # create KSP2: This is a solver for the System conj(OP)*x=y. It will later be used to solve the transposed system, thus 
        #effectively solving OP^H *x=y, which is the Hermitian transpose of the system.                 
        # TODO Sophie: unfortuantely there is no "solveHermitianTranspose" in the petsc4py (yet?). 
        # Thus we have to do an additional LU decomposistion.... Change as soon as this is included in the petsc4py!                 
        OP_H = ResolventOperator.copy()   #create a new matrix, s.t. the original one will not be overwritten                 
        OP_H.conjugate()                 
        OP_H.assemble()                 
        self._ksp2 = PETSc.KSP().create()                 
        self._ksp2.setOperators(OP_H)                 
        self._ksp2.setType(PETSc.KSP.Type.PREONLY)                 
        self._ksp2.getPC().setType(PETSc.PC.Type.LU)                 
        self._ksp2.getPC().setFactorSolverType('mumps')                 
        self._ksp2.setUp()                 
        ## create KSP3: This is a solver for the System Q_f*x=y (will later be used to solve the transposed system).                 
        ##Qf.conjugate() #=> is this needed?                 
        self._ksp3 = PETSc.KSP().create()                 
        self._ksp3.setOperators(FEMWeightMatrix_forcingNorm)                 
        self._ksp3.setType(PETSc.KSP.Type.PREONLY)                 
        self._ksp3.getPC().setType(PETSc.PC.Type.LU)                 
        self._ksp3.getPC().setFactorSolverType('mumps')                 
        self._ksp3.setUp()         

    def getSize(self):
        """
        Get the size of the operator.

        Returns
        -------
        tuple
            The size of the operator (rows, columns).
        """
        return self._size

    def getVecs(self):
        """
        Get PETSc vectors for the operator.

        Returns
        -------
        tuple
            A tuple of PETSc.Vec objects used internally by the operator.
        """                 
        return self._F1, self._F2         

    def mult(self, mat, X, Y):
        """
        Compute the matrix-vector product Y = mat * X.

        Parameters
        ----------
        mat : PETSc.Mat
            The matrix represented by the operator.
        X : PETSc.Vec
            Input vector.
        Y : PETSc.Vec
            Output vector.
        
        Returns
        -------
        PETSc.Vec
            The result of the matrix-vector product.
        """                 
        # returns Y=mat*X                 
        # mat = (W_forcing)^-1 * P_forcing^T * W_FEM^T * (R^H)^-1 * P_response^T * W_response * P_response * R^-1 * W_FEM * P_forcing                 
        self._P_forcing.mult            (X,        self._O1)  #O1 = P_f * X                 
        self._W_FEM.mult                (self._O1, self._O2)  #O2 = W_FEM * O1                 
        self._ksp1.solve                (self._O2, self._O1)  #O1 = OP^-1 * O2                    #

        self._P_response.mult           (self._O1, self._R1)  #R1 = P_r * O1                 
        self._W_response.mult           (self._R1, self._R2)  #R2 = W_r * R1                 
        self._P_response.multTranspose  (self._R2, self._O1)  #O1 = P_r^T * R2                
        self._ksp2.solveTranspose       (self._O1, self._O2)  #O2 = (OP^H)^-1 * O1                 

        self._W_FEM.multTranspose       (self._O2, self._O1)  #O1 = W_FEM^T * O2                 
        self._P_forcing.multTranspose   (self._O1, self._F1)  #F1  = P_f^T * O1                 
        self._ksp3.solve                (self._F1, Y)         #Y  = (W_f^T)^-1 * F1                 

        return Y         

    def getKSP(self):   
        """
        Get the KSP solver for the operator.

        Returns
        -------
        PETSc.KSP
            KSP solver instance that solves the resolvent equation (without its conjugate transpose). Can be used to get the response to a given forcing.
        """
        return self._ksp1

    def destroySelf(self):                 
        """
        Clean up resources to avoid memory leaks.

        This method should be called if multiple resolvent SVDs are performed consecutively.
        """
        self._ksp1.getPC().destroy()                 
        self._ksp2.getPC().destroy()                 
        self._ksp3.getPC().destroy()                 
        self._ksp1.destroy()                 
        self._ksp2.destroy()                 
        self._ksp3.destroy() 


