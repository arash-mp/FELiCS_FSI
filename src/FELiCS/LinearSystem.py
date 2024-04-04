#Standard libraries
import time
import os
import sys
import multiprocessing

# Third party libraries
import scipy.sparse.linalg as splin
import numpy as np
from scipy.sparse import (
	csr_matrix,
	csc_matrix,
	)
from scipy.sparse.linalg import norm as sparse_norm
from scipy.io import savemat

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
from FELiCS.functions import (
	printError,
	printWarning,
	printDebug,
	)
from FELiCS.fluctuationClass import fluctuationSolutions




class linearSystem:
	"""
	This Class provides methods to solve the Discretized System with the
	specified AnalysisMode. Therefore it gets the matrix_dict, which contains
	the discrete problem. Depending on the desired Analysis Mode the public
	methods solve_GEVP(), solveInputOutput() or solveResolvent()
	can be called.

	Private attributes:
	----------
	__matrix_dict : dict
		A dictionary containing the discretized-problem as csr_matrices
	__FEMSpaces : FEMSpaces
		A Object of the FEMSpaces-Class, containing the FE-Spaces
	__Preconditioner: String
		defining the kind of used preconditioning of the matrices
	__nSolut : Desfines the number of Solutions
	__omegas : Given Eigenvalues (?)
	__n_omegas : Number of given Eigenvalues
	__nCPU : Defines the number of used CPU-cores to solve the linear-system
	__results_GEVP : A dictionary in which the solution of the GEVP is saved
	it contains complex-numpy arrays
	__results : A dictionary containing the solution of the InputOutput
	 Analysis
	__gains : The calculated Gains for the given forcing for the
	Resolvent-Analysis
	__forcing : forcing-vector for the resolvent-analysis
	__response :
	__results_resolvent : Solution of the Resolvent-Analysis
	__tol : Given tolerance for the calculation of the eigenvalues

	Protected attributes:
	----------


	Public attributes:
	----------

	
	"""
	
	@staticmethod
	def runInParallel(
			args
			):
		"""runInParallel Wrapper for parallel execution

		This Function acts as a Wrapper to make it possible to use the Pickle
		methods inside the Pool.map() method, which creates several instances
		of the Pickle-Method which one omega respectively. Each function is
		then computed in its own subprocess, which makes it possible to use
		several CPU-Cores.

		Parameters
		----------
		args : list
			List of arguments for the Pickle-Methods

		Returns
		-------
		Parallel Pickle Object
		"""
		if args[1] == 'Resolvent':
			# The content of this warning is due to the fact, that properties necessary to calculate the
			# response cannot be pickled... A solution to this should be found...
			return args[0].__ParallelResolventPickle(
										args[0],
										args[2],
										args[3],
										args[4],
										)

		elif args[1] == 'GEVP':
			return args[0].__ParallelGEVPPythonPickle(
										args[0],
										args[2],
										args[3],
										args[4],
										args[5],
										args[6],
										)

		elif args[1] == 'InputOutput':
			return args[0].__ParallelInputOutputPickle(
										args[0],
										args[2],
										args[3],
										)


	def __init__(
		self,
		matrix_dict,
		matrix_dict_petsc,
		FEMSpaces,
		param,
		meanFlow,
		):
		"""
		The Constructor of the LinearSystem-Class. It constructs the
		LinearSystem Object and initializes the attributes of the object.
		The results attributes of the 3 possible
		AnalysisModes are initalized.

		Function Arguments:
		- matrix_dict: Dictionary containing the discretized linear system
		- FEMspaces: Object of the FEMSpacesClass containing the FE-Spaces
		- omegas:
		- Preconditioner:
		- nSolut: Number of Solutions
		- nCPU: Number of used Cores for solving the linear system
		"""
		self.__matrix_dict = matrix_dict
		self.__param = param
		self.__FEMSpaces = FEMSpaces
		self.__n_omegas = len(self.__param.IOResolvent.Omegas)
		self.__meanFlow = meanFlow
		self.__n_dof=np.shape(self.__matrix_dict['A'])[0]

		self.__matrix_dict_petsc = matrix_dict_petsc #store petsc matrices for PETSc/SLEPc; keep the others in as long as implementation is not finished

		self.__tol = 1e-12



	def solveResolvent(
			self,
			WeakFormulationClass,
			):
		"""
		This function solves the resolvent problem's eigenvalue problem.

		Function arguments:
		- WeakFormulationClass

		Function returns:
		- fluctSolutObjList: List of fluctuation Objects, where each Object
		contains one forcing or response
		"""

		fluctSolutObjList = []

		# Allocate matrices/vectors of the solutions
		gains = np.zeros((self.__param.Numerics.nSolut,self.__n_omegas),\
			'complex')
		forcings = np.zeros((self.__n_dof,self.__param.Numerics.nSolut,\
			self.__n_omegas),'complex')
		responses = np.zeros((self.__n_dof,self.__param.Numerics.nSolut,\
			self.__n_omegas),'complex')

		printDebug(True, '--------------------------------')
		printDebug(True, '-- Getting matrices limiting response/forcing...')
		self.__matrix_dict_petsc['Pu'] = WeakFormulationClass.getPMat()
		self.__matrix_dict_petsc['Cr'] = WeakFormulationClass.getCrMat()
  		# SD: it would make more sense to move this to WeakFormulationCollection!
		# similarly to what we do with forcing_vf 

		printDebug(True, '-- Getting weighting matrix for forcing...')
  		# Modified the weight matrix: 
    	#	(i) the multiplication with Pu^T and Pu is necessary to 
		#	apply the weigths only in domain of the restrictor. 
		# 	(ii) TEMPORARY: we use "B_femWeight" instead of "B" to be
		# 	able to force on any variable (e.g. in compressible case).
		#	But a more consistent implementation must be done!
                #       Q = Pu^T * B_femWeight * Pu
		self.__matrix_dict_petsc['Q'] = \
      			self.__matrix_dict_petsc['Pu'].transposeMatMult(\
                        self.__matrix_dict_petsc['B_femWeight'].matMult(\
                        self.__matrix_dict_petsc['Pu']))

		printDebug(True, '-- Done.')

		tic_res = time.perf_counter()

		if self.__param.Numerics.nCPU > 1:
			printError('Parallel computation with "nCPU" > 1 is currently not working. (Will be added in the restructuring process.)') 			
			exit()
			## Parallel computation of the resolvent --------------------------------
			#printDebug(True, '-- Parallel computation of forcing and gains...')
			#pool = multiprocessing.Pool(self.__param.Numerics.nCPU)

			#args_map = [(linearSystem, 'Resolvent', self.__matrix_dict, \
			#	self.__param.Numerics.nSolut, arg) for arg in \
			#	self.__param.IOResolvent.Omegas]
			##pool.map(run_in_parallel, args_map)
			#resultsPool = pool.map(self.runInParallel, args_map)

			## The gains, forcing, and responses are obtained in the parallel loop
			#for i in range(len(resultsPool)):
			#	gains[:, i] = resultsPool[i][0]
			#	forcings[:, :, i] = resultsPool[i][1]
			#	responses[:, :, i] = resultsPool[i][2]

		else:
			# Serial computation of the resolvent --------------------------------
			printDebug(True, '--------------------------------')
			printDebug(True, '-- Serial computation of forcing, gains, and responses...')
   

			for i in range(self.__n_omegas):
				
				# Frequency 
				omega = self.__param.IOResolvent.Omegas[i]
				printDebug(True, "-- Performing resolvent analysis for omega = " + str(omega))

				OP_petsc          = self.__matrix_dict_petsc['A']-omega*self.__matrix_dict_petsc['B']
				resolventOperator = ResolventOperator(
                                        OP_petsc, 
                                        self.__matrix_dict_petsc['Q'], 
                                        self.__matrix_dict_petsc['Pu'],
                                        self.__matrix_dict_petsc['Cr'],
                                        self.__matrix_dict_petsc['B_forcing'],
                                        self.__matrix_dict_petsc['B_response'])
	
				# Perform eigenvalue decomposition of the linear operator defined in the class "ResolventOperator"
				# via the matrix vector multiplation "mult"
				gains[:,i],eigenvectors_c = self.__solveSVDOfResolvent(
                                                                        resolventOperator,
									nev=self.__param.Numerics.nSolut,
									tol=1.e-13,
									max_it=200,
									)



				# Write gains to results dictionary
				gains[:, i] = np.real(gains[:, i])

				# Iterate through the first nSolut gains
				# Compute the respetive forcing and responses with the solution of the SVD ("eigenvetors_c")
                                # TODO Sophie: do this more elegantly
				for k in range(self.__param.Numerics.nSolut):

                                        # get petsc vectors from petsc matrices (=get petsc vectors with correct sizes)
					X1, X2 = self.__matrix_dict_petsc['Pu'].getVecs()
					X1.setValues(range(0,len(eigenvectors_c[:,k])),eigenvectors_c[:,k])
					Y1, Y2 = self.__matrix_dict_petsc['B_femWeight'].getVecs()

					# forcings = Pu*eigenVectors
					self.__matrix_dict_petsc['Pu'].mult(X1,X2)
					forcings[:,k,i] = X2.getValues(range(0,X2.getSize()))

					# Y1 = -1j * B_femWeight * forcings
					self.__matrix_dict_petsc['B_femWeight'].mult(X2,Y1)
					Y1.scale(-1j)

					# solve (A-omega*B)*responses = Y1
					resolventOperator.getKSP().solve(Y1,X2)
					responses[:, k, i] = X2.getValues(range(0, X2.getSize()))


				resolventOperator.destroySelf()

		# construct a fluctuationSolutions Object for the Response, Forcing of
		# every Frequency
		for i, omega in enumerate(self.__param.IOResolvent.Omegas):
			for gainNumb in range(responses.shape[1]):
				fluctSolutForcing = fluctuationSolutions(
									self.__param,
									self.__meanFlow,
									self.__FEMSpaces,
									self.__param.IOResolvent.Omegas[i],
									forcings[:,gainNumb,i],
									False,
									gainNumb,
									gains[gainNumb, i],
									)
				fluctSolutResponse = fluctuationSolutions(
									self.__param,
									self.__meanFlow,
									self.__FEMSpaces,
									self.__param.IOResolvent.Omegas[i],
									responses[:,gainNumb, i],
									True,
									gainNumb,
									gains[gainNumb, i],
									)

				fluctSolutObjList.append(fluctSolutForcing)
				fluctSolutObjList.append(fluctSolutResponse)

		toc_res = time.perf_counter() - tic_res
		printDebug(True, f"-- Solving resolvent took: {toc_res:0.4f} seconds")

		return fluctSolutObjList

	def solveInputOutput(
		self
		):
		"""solveInputOutput 

		Computes the Input-Output Analysis

		Returns
		-------
		fluctSolutObjList

			List of fluctuation Objects, where each Object contains one eigenvector and eigenvalue.
		"""
		fluctSolutObjList = []

		gains = np.zeros((1,self.__n_omegas))
		responses = np.zeros((self.__n_dof,1,self.__n_omegas),'complex')

		n_omegas=len(self.__param.IOResolvent.Omegas)


		if self.__param.Numerics.nCPU > 1:
			printError('Parallel computation with "nCPU" > 1 is currently not working. (Will be added in the restructuring process.)') 			
			exit()
			#pool=multiprocessing.Pool(processes=self.__param.Numerics.nCPU)
			##func= partial(self.__ParallelInputOutputPickle,self.__matrix_dict)
			#args_map = [[linearSystem, 'InputOutput', self.__matrix_dict, \
			#	arg] for arg in self.__param.IOResolvent.Omegas]

			#resultsPool = pool.map(self.runInParallel, args_map)
			#for i in range(len(resultsPool)):
			#	gains[:,i]	= 1
			#	responses[:,0,i] = resultsPool[i]

		else:
			## Repeat analysis for every omega
			for i in range(n_omegas):
				omega=self.__param.IOResolvent.Omegas[i]
				print("Performing input-output analysis for omega="+str(omega))
				### Get linear operator
				# Define OP
				OP               = self.__matrix_dict_petsc['A']-omega*self.__matrix_dict_petsc['B']
				# Use PETSC to solve linear system (with LU decomposition)
				eigenvectors_c   = self.__solveEquationSystem(OP, self.__matrix_dict_petsc['b_forcing'])
				gains[:,i]       = 1
				responses[:,0,i] = eigenvectors_c

		for i in range(len(self.__param.IOResolvent.Omegas)):
			# construct for each EVal and EVec a fluctuationSolution
			for gainNumb in range(responses.shape[1]):
				fluctSolutObjList.append(fluctuationSolutions(
						self.__param,
						self.__meanFlow,
						self.__FEMSpaces,
						self.__param.IOResolvent.Omegas[i],
						responses[:,gainNumb,i],
						True,
						gainNumb,
						gains[gainNumb,i]
						))

		return fluctSolutObjList



	def solveGEVP(
				self,
				adjointFlag=False,
				):
		"""
		solves the General Eigenvalue Problem

		Function Attributes:
		- adjointFlag:  (bool, optional) if set to "True", the adjoint solution is computed alongside the direct solution. If "False" (default), only the direct solution is computed. 

		Function returns:
		- __results_GEVP: Dictionary containing the direct- and (optionally) the adjoint eigenvalues and eigenvectors
		"""


		start= time.time()

		guesses  = self.__param.Numerics.EigenValueGuess
		nSol     = self.__param.Numerics.nSolut

		# allocate solution arrays
		EVal    = np.zeros((nSol*len(guesses)),'complex')  
		EVec    = np.zeros((self.__n_dof,nSol*len(guesses)),'complex') 		
		if adjointFlag: 			
			EVecAdj = np.zeros((self.__n_dof,nSol*len(guesses)),'complex')
		residui = np.zeros((nSol*len(guesses))) 


		# Possibility to run GEVP of different guesses in parallel
		if self.__param.Numerics.nCPU > 1:
			printError('Parallel computation with "nCPU" > 1 is currently not working. (Will be added in the restructuring process.)') 			
			exit()
			#print("-- Entering parallel loop for GEVP")
			#pool=multiprocessing.Pool(processes=self.__param.Numerics.nCPU)
			#args_map = [(linearSystem, 'GEVP', self.__matrix_dict_petsc, \
			#	nSol, guess, adjointFlag, self.__param.Numerics) \
                        #        for guess in guesses]
			#results_pool = pool.map(self.runInParallel, args_map)

			#for i in range(len(results_pool)):
			#	index = list(range(i*nSol, (i+1)*nSol))
			#	EVal[index]    = results_pool[i][0]
			#	for j in range(nSol):
			#		EVec[:,i*nSol+j] = results_pool[i][1][j,:]
			#		if adjointFlag:
			#			EVecAdj[:,i*nSol+j] = results_pool[i][2][j,:]
			#			residui[index] = resoults_pool[i][3]
			#		else:
			#			residui[index] = resoults_pool[i][2]

			#print("-- Assembled results from all guesses.")

		else:
			for i in range(len(guesses)):
				if adjointFlag:
					EigValTemp, EigVecTemp, EigVecAdjTemp, error = self.__solveGEVP_with_SLEPc(self.__matrix_dict_petsc,
                        	                                                                 nSol,
                        	                                                                 guesses[i],
                        	                                                                 adjointFlag=True)
				else:
					EigValTemp, EigVecTemp, error                = self.__solveGEVP_with_SLEPc(self.__matrix_dict_petsc,
                        	                                                                 nSol,
                        	                                                                 guesses[i],
        			                                                                 adjointFlag=False)
				index = list(range(i*nSol,(i+1)*nSol))
				EVal[index]    = EigValTemp
				for j in range(nSol):
					EVec[:,i*nSol+j] = EigVecTemp[j,:]
					if adjointFlag:
						EVecAdj[:,i*nSol+j] = EigVecAdjTemp[j,:]
				residui[index] = error 


		# normalize the solution
		for i in range(self.__param.Numerics.nSolut):
			EVec[:,i] = EVec[:,i]/np.linalg.norm(EVec[:,i])


		fluctSolutObjList    = []
		fluctSolutObjListAdj = []

		for i in range(0,len(EVal)):

			# construct for each EVal and EVec a fluctuationSolution
			fluctSolutObjList.append(fluctuationSolutions(
									self.__param,
									self.__meanFlow,
									self.__FEMSpaces,
									EVal[i],
									EVec[:, i],
								        True,
									 )
			)

			if adjointFlag:
				fluctSolutObjListAdj.append(fluctuationSolutions(
										self.__param,
										self.__meanFlow,
										self.__FEMSpaces,
										EVal[i],
										EVecAdj[:, i],
										False,
										 )
				)


		end = time.time() - start

		residuum = np.sqrt(np.sum(residui[:]**2.))

		printDebug(True, '-- Solving the GEVP took %4g s' % end)
		printDebug(True, '-- Residuum of solutions (Euclidian norm of all residui that are calculated by SLEPc):  %12g' % (residuum))

		if adjointFlag:
			fluctSolutObjList.extend(fluctSolutObjListAdj)

		return fluctSolutObjList




	def __solveGEVP_with_SLEPc(
				self,
                                matrix_dict_petsc, 
                                nSol, 
                                eigenValueGuess, 
                                adjointFlag=False, 
                                Numerics=None,
				):
		"""
		Solves the GEVP with the a SLEPc eigenvalue solver

		Function arguments:
                - adjointFlag       (bool, optional) If set to "True", the adjoint eigenvectors will be calculated additionally. If set to "False", only the direct eigenvectors will be returned.

		Function returns:
                - EigValTot         numpy array; returns the computed eigenvalues 
                - EigVecTot         numpy array; returns the computed right-hand-side eigenvectors in one big array ("stacked" onto each other)
                - EigVecAdjTot      numpy array, only returned if "adjointFlag=True"; returns the computed left-hand-side eigenvectors in one big array ("stacked" onto each other)
                - error             numpy array; returns, for each eigenproblem solution, the corresponding relative residuum of the eigenproblem
		"""

		# solve GEVP using SLEPc

		printDebug(True, "-- Solving for guess: %4a" % eigenValueGuess)

		if adjointFlag:
			EigVal, EigVec, EigVecAdj, error = self.__solveGeneralEigenproblem(\
                                     self.__matrix_dict_petsc['A'], \
                                     self.__matrix_dict_petsc['B'], \
                                     sigma=eigenValueGuess, \
                                     nev=nSol, \
                                     tol=self.__tol,\
                                     max_it=200,\
                                     adjoint=True, \
                                     isForEigenProblem=True)
			return EigVal, EigVec, EigVecAdj, error

		else:
			EigVal, EigVec, error = self.__solveGeneralEigenproblem(\
                                     self.__matrix_dict_petsc['A'], \
                                     self.__matrix_dict_petsc['B'], \
                                     sigma=eigenValueGuess, \
                                     nev=nSol, \
                                     tol=self.__tol,\
                                     max_it=200,\
                                     adjoint=False, \
                                     isForEigenProblem=True)
			return EigVal, EigVec, error



	@staticmethod
	def runInParallel(
			args
			):
		"""
		This Function acts as a Wrapper to make it possible to use the Pickle
		methods inside the Pool.map() method, which creates several instances
		of the Pickle-Method which one omega respectively. Each function is
		then computed in its own subprocess, which makes it possible to use
		several CPU-Cores.

		Function arguments:
		- args: List of arguments for the Pickle-Methods

		Function returns:

		"""
		if args[1] == 'Resolvent':
			# The content of this warning is due to the fact, that properties necessary to calculate the
			# response cannot be pickled... A solution to this should be found...
			return args[0].__ParallelResolventPickle(
										args[0],
										args[2],
										args[3],
										args[4],
										)

		elif args[1] == 'GEVP':
			return args[0].__ParallelGEVPPythonPickle(
										args[0],
										args[2],
										args[3],
										args[4],
										args[5],
										args[6],
										)

		elif args[1] == 'InputOutput':
			return args[0].__ParallelInputOutputPickle(
										args[0],
										args[2],
										args[3],
										)




	def __ParallelGEVPPickle(
				self,
				matrix_dict,
				nSolut,
				EigGuess,
				adjointFlag,
				Numerics,
				):
		"""
		This function is used for distributing the GEVP associated
		with different eigenvalue guesses to various CPUs.

		"""
		current = multiprocessing.current_process()
		print("-- " + current.name + " running for omega = " + str(EigGuess))


		#eigenvalues, eigenvectors = splin.eigs(
		#								A,
		#								k=nSolut,
		#								M=B,
		#								sigma=EigGuess,
		#								ncv=200,
		#								maxiter=200,
		#								tol=1e-12,
		#								return_eigenvectors=True,
		#								)

		return self.__solveGEVP_withSLEPc(matrix_dict, nSolut, EigGuess, adjointFlag, Numerics)




	def __ParallelInputOutputPickle(
				self,
				matrix_dict,
				omega,
				):
		"""
		 This function is used for distributing the input output analysis to
		 various CPUs.

		 Function arguments:
		 - omega:

		 Function returns:
		 -eigenvectors_c:
		 """
		print("Performing input-output analysis for omega="+str(omega))
		### Get linear operator
		# Define OP
		OP=matrix_dict['A']-omega*matrix_dict['B']
		# Make OP sparse vector
		OP=OP.tocsc()
		# Perform Lower-Upper decomposition
		LU = splin.splu(OP,permc_spec=3)
		# Use LU decomposition to solve the linear system
		eigenvectors_c=LU.solve(matrix_dict['b_forcing'])

		return     eigenvectors_c


	def __ParallelResolventPickle(
					self,
					matrix_dict,
					nSolut,
					omega,):
		"""
		This function is used for distributing the resolvent's
		eigenvalueproblems to various CPUs.

		"""
		LUQ = splin.splu(matrix_dict['Q'],permc_spec=3)
		nu = min(np.shape(matrix_dict['Pu']))
		current = multiprocessing.current_process()
		print("-- " + current.name + " running for omega = " + str(omega))
		def op(x):
			y = matrix_dict['Pu'] * x
			z = matrix_dict['B_forcing'] * y
			y = LU.solve(z)
			y1 = matrix_dict['Cr'] * y
			z = matrix_dict['B_response'] * y1
			z1 = matrix_dict['Cr'].transpose() * z
			y = LU.solve(z1, trans='H')
			z = matrix_dict['B_forcing'].transpose() * y
			y = matrix_dict['Pu'].transpose() * z
			w = LUQ.solve(y, trans='H')
			return w

		OP = matrix_dict['A'] - omega * matrix_dict['B']
		OP = OP.tocsc()
		LU = splin.splu(OP,permc_spec=3)
		SOP = splin.LinearOperator((nu,nu),matvec = op,dtype = 'complex')

		gains,eigenvectors_c = splin.eigs(SOP,
								k = nSolut,
								M = None,
								sigma = None,
								which = 'LM',
								maxiter = 100,
								tol = 10-12,
								return_eigenvectors = True,
								)

		# Write gains to results dictionary
		gains = np.real(gains)

		# Iterate through the first nSolut gains
		nDOF = max(np.shape(matrix_dict['B']))
		forcings = np.zeros((nDOF, nSolut), 'complex')
		responses = np.zeros((nDOF, nSolut), 'complex')
		for k in range(nSolut):
			# Write the respective forcing to results dictionary
			forcings[:, k] = matrix_dict['Pu'] * eigenvectors_c[:, k]
			f = -1j * matrix_dict['B'] * forcings[:, k]
			responses[:, k] = LU.solve(f)

		return gains, forcings, responses



	def __solveGeneralEigenproblem(
                self, 
                A, 
                B, 
                sigma, 
                nev, 
                tol=1.e-16, 
                max_it=200, 
                adjoint=False, 
                isForEigenProblem=True, 
                isIncompressible=True):

		#TODO Sophie: create better system to identify matrix kinds

		"""
		Solves the generalized eigenvalue problem (GEVP) 
  		using the SLEPc and PETSc libraries.

		Function arguments:
		- A, B		Matrices defining the GEVP (A-wB)q = 0
		- sigma 	Eigenvalue guess
		- nev 		Number of eigenvalues to compute
		- tol	        (optional) precision of GEVP
		- max_it	(optional) maximum number of iterations 
		- adjoint	(bool, optional) if set to "True", compute the adjoint GEVP alongside the direct one
                - isForEigenProblem, isIncompressible 	    (optional) flags that can be set to determine the kind of matrices  

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


		# create eigenproblem solver 
		eps = SLEPc.EPS().create()
		eps.setOperators(A,B)
		if isForEigenProblem and isIncompressible:
			eps.setProblemType(SLEPc.EPS.ProblemType.PGNHEP) # general non-Hermitian eigenproblem with positive semi-definite M
		
		#calculate adjoint vectors
		if adjoint:
			eps.setTwoSided(True)
		
		eps.setTolerances(tol=tol,max_it=max_it)
		
		eps.setType(SLEPc.EPS.Type.KRYLOVSCHUR) # is standard, does not need to be set
		eps.getST().setType(SLEPc.ST.Type.SINVERT)
		eps.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)
		eps.setTarget(sigma)
		eps.setDimensions(nev=nev)
		
		eps.getST().getKSP().getPC().setType('lu')
		eps.getST().getKSP().getPC().setFactorSolverType('mumps')
		
		eps.setFromOptions()
		eps.setUp()
		
		eps.solve()
		
		dim = A.getSize()[0]
		eigVals,  eigVecs  = np.empty(nev,complex), np.empty([nev,dim],complex)
		if adjoint:
			eigVecs_adjoint  =  np.empty([nev,dim],complex)
		vec_real, vec_imag = A.getVecs()
		error = np.empty(nev)
		for i in range(nev):
			try:
				eigVals[i]   = eps.getEigenpair(i,vec_real,vec_imag)
				eigVecs[i,:] = vec_real.getArray() + 1j * vec_imag.getArray()
				if adjoint:
					# get adjoint solution
					eps.getLeftEigenvector(i,vec_real,vec_imag)
					eigVecs_adjoint[i,:] = vec_real.getArray() + 1j * vec_imag.getArray()
				error[i] = eps.computeError(i, SLEPc.EPS.ErrorType.RELATIVE)
				#printDebug(True,f"SLEPc error relative: {error[i]}")
			
			except:
				printWarning("Could not access eigenpair nb ", nev+1, "!")
		
		eps.getST().getKSP().getPC().destroy()
		eps.getST().getKSP().destroy()
		eps.getST().destroy()
		eps.destroy()
		if adjoint:
		        return eigVals, eigVecs, eigVecs_adjoint, error
		else:
		        return eigVals, eigVecs, error
	


	def __solveSVDOfResolvent(
                self, 
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



	def __solveEquationSystem(
		self,
		A,
		b):
	
		"""
		Solves a linear equation system Ax=b, using the PETSc libraries.

		Function arguments:
		- A     PETSc matrix,  matrix of the linear system
		- b     PETSc vector,  rhs of the linear system

		Function returns:
		- x     numpy array, solution of the linear equation system
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
		A.destroy()
		b.destroy()
		dummy.destroy()
		
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




        ###########################
        ##### DEPRECATED ##########
        ###########################
#	def checkMatrix(self, eq='ux', searchPoint=[1, 2]):
#		import csv
#		"""
#        This function is used for debugging. It can be used to examine matrix A
#        at a given search point for a given equation.
#        Inputs:
#        - eq: equation - ux, uy, ut, p
#        - searchPoint: search point as tuple [x,y]
#        It prints
#        - the diagonal element at calculated index
#        - the matrix row corresponding to coordinate and variable
#        - the amount of diagonal elements in the whole matrix that are > 10^20
#        - the amount of diagonal elements belonging to the examined subspace
#        that are > 10^20 and the coordinates they are located at
#        Output:
#        - returns index of matrix that corresponds to equation and closest 
#        coordinate
#        """
#
#		checkMat = self.__matrix_dict['A']
#
#		# GENERAL ATTRIBUTES
#		mesh = self.__FEMSpaces.P2.mesh
#		gdim = mesh.geometry.dim
#		if eq == 'ux':
#			dofs = self.__FEMSpaces.VMixed.sub(0).sub(0).collapse()[1]
#		elif eq == 'uy':
#			dofs = self.__FEMSpaces.VMixed.sub(0).sub(1).collapse()[1]
#		elif eq == 'ut':
#			dofs = self.__FEMSpaces.VMixed.sub(1).collapse()[1]
#		elif eq == 'p' :
#			dofs = self.__FEMSpaces.VMixed.sub(2).collapse()[1]
#		else:
#			print('equation not found', file=sys.stderr)
#			return
#
#		# DOF COORDINATES OF MIXED SPACE
#		nDofsVmixed = Function(self.__FEMSpaces.VMixed).vector[:].shape[0]
#		dofs_coord = np.zeros((nDofsVmixed, gdim))
#		for i in range(self.__FEMSpaces.VMixed.num_sub_spaces):
#			if self.__FEMSpaces.VMixed.sub(i).num_sub_spaces > 0:
#				for j in range(self.__FEMSpaces.VMixed.sub(i).num_sub_spaces):
#					dofCoordsOfSubSpace = \
#						self.__FEMSpaces.VMixed.sub(i).sub(j).collapse()[
#							0].tabulate_dof_coordinates()[:, 0:gdim]
#					indices = self.__FEMSpaces.VMixed.sub(i).sub(j).collapse()[
#						1]
#					dofs_coord[indices] = dofCoordsOfSubSpace
#			else:
#				dofCoordsOfSubSpace = self.__FEMSpaces.VMixed.sub(i).collapse()[
#										  0].tabulate_dof_coordinates()[:,
#									  0:gdim]
#				indices = self.__FEMSpaces.VMixed.sub(i).collapse()[1]
#				dofs_coord[indices] = dofCoordsOfSubSpace
#		dofs_coord = dofs_coord.reshape((-1, gdim))
#
#		# FIND THE CLOSEST COORDINATE TO SEARCH POINT
#		# FINDS ALL FOUR INDICES FOR GIVEN SEARCH POINT
#		nearestIndex = np.where(list(
#			map(lambda x: np.linalg.norm(x - searchPoint),
#				dofs_coord)) == min(list(
#			map(lambda x: np.linalg.norm(x - searchPoint),
#				dofs_coord))))
#		#print("closest coordinate: ", dofs_coord[nearestIndex[0][0]], file=sys.stderr)
#		# CHECKS WHICH INDEX BELONGS TO EXAMINED SUBSPACE
#		matchingIndices = np.empty(shape=(0, 0))
#		for index in nearestIndex[0]:
#			if index in dofs:
#				matchingIndices = np.append(matchingIndices, index)
#
#		# PRINT MATRIX ROW CORRESPONDING TO COORDINATE
#		for index in matchingIndices:
#			#print("line ", index, " in A: ", np.shape(checkMat[int(index)]), file=sys.stderr)
#			nonzeroList = checkMat[int(index)].nonzero()
#			row = []
#			for jndex in nonzeroList[1]:
#				row.append(checkMat[int(index), jndex])
#
#		# PRINT EVERY MATRIX ELEMENT IN EVERY ROW CORRESPONDING TO EQUATION
#		maxNonZeroNumb = 0
#		for kindex in dofs:
#			line = []
#			line.append(dofs_coord[kindex][0])
#			line.append(dofs_coord[kindex][1])
#			nonzeroList = checkMat[kindex].nonzero()
#			if len(nonzeroList[1]) > maxNonZeroNumb:
#				maxNonZeroNumb = len(nonzeroList[1])
#			row = []
#			for lindex in nonzeroList[1]:
#				row.append(checkMat[kindex, lindex])
#			row.sort()
#			line.extend(row)
#
#		# HOW MANY DIAGONAL ELEMENTS ARE > 10^20 IN THE WHOLE MATRIX?
#		diagonal = checkMat.diagonal()
#		homogDiric = sum(x > 10e20 for x in diagonal)
#		#print("Amount of diagonal elements in whole matrix > 10^20: ", homogDiric, file=sys.stderr)
#
#		# HOW MANY DIAGONAL ELEMENTS BELONGING TO THE EXAMINED SUBSPACE ARE
#		# > 10^20 AND AT WHAT COORDINATES ARE THEY LOCATED?
#		counter = 0
#		coordinates = []
#		for index in dofs:
#			if checkMat[index, index] > 10e20:
#				counter += 1
#				coordinates.append(tuple(dofs_coord[index, :]))
#		#print("Amount of diagonal elements in examined subspace > 10^20: ", counter, file=sys.stderr)
#		# print("Coordinates corresponding to values > 10^20 ", coordinates, file=sys.stderr)
#
#		return
        ###########################
        ###########################




