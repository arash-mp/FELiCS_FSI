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
from functions import (
	printError,
	printWarning,
	printDebug,
	)
from fluctuationClass import fluctuationSolutions

import pdb

class linearSystem:
	"""
	This Class provides methods to solve the Discretized System with the
	specified AnalysisMode. Therefore it gets the matrix_dict, which contains
	the discrete problem. Depending on the desired Analysis Mode the public
	methods solve_GEVP(), solveInputOutput() or solveResolvent()
	can be called.

	Private attributes:
	- __matrix_dict: A dictionary containing the discretized-problem as
	csr_matrices
	- __FEMSpaces: A Object of the FEMSpaces-Class, containing the FE-Spaces
	- __Preconditioner: String, defining the kind of used preconditioning of
	the matrices
	- __nSolut: Desfines the number of Solutions
	- __omegas: Given Eigenvalues (?)
	- __n_omegas: Number of given Eigenvalues
	- __nCPU: Defines the number of used CPU-cores to solve the linear-system
	- __results_GEVP: A dictionary in which the solution of the GEVP is saved
	it contains complex-numpy arrays
	- __results: A dictionary containing the solution of the InputOutput
	 Analysis
	- __gains: The calculated Gains for the given forcing for the
	Resolvent-Analysis
	- __forcing: forcing-vector for the resolvent-analysis
	- __response:
	- __results_resolvent: Solution of the Resolvent-Analysis
	- __tol: Given tolerance for the calculation of the eigenvalues

	Protected attributes:

	Public attributes:

	"""
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


	def __init__(
		self,
		matrix_dict,
		FEMSpaces,
		param,
		meanFlow,
		):
		"""
		The Constructor of the LinearSystem-Class. It constructs the
		LinearSystem Object and initalizes the attributes of the object.
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

		self.__tol = 1e-12

	def __IntegrateForcingVector(
				self,
				eigenvectors_imag,
				eigenvectors_real
				):
		"""
		Computes the Forcing-Vector for the Resolvent Analysis

		Function arguments:
		- eigenvectors_imag: imaginatinal part of the eigenvector
		- eigenvectors_real: real part of the eigenvector

		Function returns:
		- f: forcing vector
		"""
		# Define the variational formulations for the forcings real
		# and imaginary part
		f_real_vf = 0
		f_imag_vf = 0

		loop_count = 0
		# Loop through all solutions
		for i_name, name in enumerate(self.__param.SolutionList):

			# if nSolut = 1, it shouldnt iterate over all solutions u, ut and p.
			# It should just iterate over u
			# so a switch is implemented, which leaves the outer loop abortive
			loop_count = loop_count + 1
			if loop_count > self.__param.Numerics.nSolut:
				break
			# If it is the velocity, all velocity components must be considered
			if name == 'u':
			# Loop through the velocity components
				for (i_component, component) in enumerate(\
				self.__param.VelocityComponents):

					# Add the imaginary part and the real part for the
					# respective component to the respective equation
					# Here, the multiplication by '-i' is taken into account.
					# Could also be added later, when the
					# linear system is solved
					f_real_vf += self.__R*eigenvectors_imag.\
					split()[i_name].split()[i_component] \
					*self.__X[i_name][i_component] *dx

					f_imag_vf += -self.__R*eigenvectors_real.split()[i_name].\
					split()[i_component] \
					*self.__X[i_name][i_component] *dx

			else:
				# Add the imaginary part and the real part for
				# the respective component to the respective equation
				# Here, the multiplication by '-i' is taken into account.
				# Could also be added later, when the
				# linear system is solved
				f_real_vf += self.__R*eigenvectors_imag.split()[i_name] \
													  *self.__X[i_name] *dx
				f_imag_vf += -self.__R*eigenvectors_real.split()[i_name] \
													   *self.__X[i_name] *dx

			# Define the vectors, and assemple them using the
			# variational formulation defined above

		f_real =assemble_vector(form(f_real_vf))
		f_real.assemble()

		f_imag =assemble_vector(form(f_imag_vf))
		f_imag.assemble()

		#pdb.set_trace()
		return f_real.array + 1j*f_imag.array

	def __solve_with_python(
				self,
				adjointFlag=False,
				):
		"""
		Solves the GEVP with the python solver

		Function arguments:

		Function returns:

		"""

		# EVal = np.zeros((self.__param.Numerics.nSolut), 'complex')
		# EVec = np.zeros((self.__n_dof, self.__param.Numerics.nSolut), 'complex')

		if adjointFlag:
			printDebug(True,'--------------------------------')
			printDebug(True, "-- Solving adjoint GEVP")
			A, B, f = self.__preconditionMatrices(
										self.__matrix_dict['A'].getH(),
										self.__matrix_dict['B'],
										self.__param.Numerics.Preconditioner,
												)
		else:
			printDebug(True,'--------------------------------')
			printDebug(True, "-- Solving direct GEVP")
			A, B, f = self.__preconditionMatrices(
						self.__matrix_dict['A'],
						self.__matrix_dict['B'],
						self.__param.Numerics.Preconditioner
												)

		# Allocate space for complete matrices
		nGuesses = len(self.__param.Numerics.EigenValueGuess)
		nSol = self.__param.Numerics.nSolut
		EigValTot = np.zeros((nSol*nGuesses),'complex')
		EigVecTot = np.zeros((self.__n_dof,nSol*nGuesses),'complex')

		# solve GEVP using eigs for each guess
		for i in range(nGuesses):
			eigenValueGuess = self.__param.Numerics.EigenValueGuess[i]
			printDebug(True, "-- Solving for guess: %4a" % eigenValueGuess)
			EigValTemp, EigVecTemp = splin.eigs(
			 	A,
			 	k=self.__param.Numerics.nSolut,
			 	M=B,
			 	sigma=eigenValueGuess,
			 	ncv=200,
			 	maxiter=100,
			 	tol=self.__tol,
			 	return_eigenvectors=True,
			 	)
			index = list(range(i*nSol,(i+1)*nSol))
			EigValTot[index] = EigValTemp
			EigVecTot[:,index] = EigVecTemp

		return EigValTot, EigVecTot

	def __solve_with_matlab(
				self,
				adjointFlag=False,
				export=False
				):
		"""
		solves the GEVP using the matlab solver

		Function Arguments:
		- export: Specifies if the solution should be exported.
		Default Value is False

		Function returns:

		"""
		import matlab.engine
		eng = matlab.engine.start_matlab()
		eng.addpath (__file__.rsplit('/',1)[0]+'/Matlab', nargout= 0 )

		if adjointFlag:
			print("-- Solving adjoint GEVP")
			# precondition matrices
			A, B, f = self.__preconditionMatrices(
								self.__matrix_dict['A'].getH(),
								self.__matrix_dict['B'],
								self.__param.Numerics.Preconditioner,
								)

			tol = 1e-12


		else:
			print("-- Solving direct GEVP")
			A, B, f = self.__preconditionMatrices(
								self.__matrix_dict['A'],
								self.__matrix_dict['B'],
								self.__param.Numerics.Preconditioner,
								)

			tol = 10e-11

		# save matrices to file
		savemat('A.mat',{'A':A})
		savemat('B.mat',{'B':B})

		if export == True:
			print('Matrix exported. Ending program...')
			exit()

		#run the GEVP solver in the matlab script
		if len(self.__param.Numerics.EigenValueGuess) == 1:
			eigenValueGuess = self.__param.Numerics.EigenValueGuess[0]
		else: printError('matlab GEVP solver only takes one eigenvalue guess')

		temp = eng.MatlabEigs(self.__param.Numerics.nSolut,
			eigenValueGuess,
			200,
			100,
			tol,
			False,
			nargout=2)

		#remove matrices
		os.remove("A.mat")
		os.remove("B.mat")

		return np.array(temp[0]), np.array(temp[1])



	def __checkGEVP(
			self,
			eigVal,
			eigVec,
			A,
			B,
			):
		"""
		Calculating the frombenius norm and the infinite norm of the
		residuum of the solutions to a GEVP

		Function arguments:
		- eigVal:
		- eigVec:
		- A:
		- B:

		"""

		# First normalize the vector (additionally with the norm of matrix A)
		normA = sparse_norm(A)
		eigVecNorm = eigVec/np.linalg.norm(eigVec)/normA
		# Calculate residuum of Ax = omega Bx
		residuum_vec=A@eigVecNorm-eigVal*B@eigVecNorm
		residuum_euc = np.linalg.norm(residuum_vec)
		return residuum_euc

	def __ParallelGEVPPythonPickle(
				self,
				matrix_dict,
				nSolut,
				EigGuess,
				FlagAdjoint,
				Numerics,
				):
		"""
		This function is used for distributing the GEVP associated
		with different eigenvalue guesses to various CPUs.

		"""
		current = multiprocessing.current_process()
		print("-- " + current.name + " running for omega = " + str(EigGuess))

		if FlagAdjoint:

			A, B, f = self.__preconditionMatrices(
										self,
										matrix_dict['A'].getH(),
										matrix_dict['B'],
										'None',
												)
		else:

			A, B, f = self.__preconditionMatrices(
										self,
										matrix_dict['A'],
										matrix_dict['B'],
										'None',
												)

		eigenvalues, eigenvectors = splin.eigs(
										A,
										k=nSolut,
										M=B,
										sigma=EigGuess,
										ncv=200,
										maxiter=100,
										tol=1e-12,
										return_eigenvectors=True,
										)

		return     eigenvalues,eigenvectors


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


	def __preconditionMatrices(
				self,
				A,
				B,
				Preconditioner,
				f=False,
				):
		"""
		Performing a preconditioning of the matrices A,B

		Function arguments:
		- A:
		- B:
		- Preconditioner:
		- f:

		Function returns:
		- A:
		- B:
		- f:

		"""
		if not Preconditioner == 'None':
			if Preconditioner == 'Sum of row':
				S = csr_matrix(1/abs(A).sum(1))
				A = A.multiply(S)
				B = B.multiply(S)

				if not type(f) == bool:
					f = np.multiply(np.squeeze(np.array(S.todense())),f)
			else:
				printError('Preconditioner '+ Preconditioner + ' not known.')
		return A,B,f

	def solveResolvent(
			self,
			WeakFormulationClass,
			):
		"""
		This function solves the resolvent problem's eigenvalue problem.
		Matrices are taken from the global_variables module...

		Function arguments:
		-
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

		# Get the matrices defining the forcing/response norms (compressed)
		self.__matrix_dict['B_forcing'] = \
			self.__matrix_dict['B_forcing'].tocsc()
		self.__matrix_dict['B_response'] = \
			self.__matrix_dict['B_response'].tocsc()

		# Matrix containing the FEM-integration weights
		self.__matrix_dict['B_femWeight'] = \
      		self.__matrix_dict['B_femWeight'].tocsc()
  
		printDebug(True, '--------------------------------')
		printDebug(True, '-- Getting matrices limiting response/forcing...')
		self.__matrix_dict['Pu'] = WeakFormulationClass.getPMat()
		self.__matrix_dict['Pu'] = self.__matrix_dict['Pu'].tocsc()
		self.__matrix_dict['Cr'] = WeakFormulationClass.getCrMat()
		self.__matrix_dict['Cr'] = self.__matrix_dict['Cr'].tocsc()
  		# SD: it would make more sense to move this to WeakFormulationCollection!
		# similarly to what we do with forcing_vf 

		printDebug(True, '-- Getting weighting matrix for forcing...')
  		# Modified the weight matrix: 
    	#	(i) the multiplication with Pu^T and Pu is necessary to 
		#	apply the weigths only in domain of the restrictor. 
		# 	(ii) TEMPORARY: we use "B_femWeight" instead of "B" to be
		# 	able to force on any variable (e.g. in compressible case).
		#	But a more consistent implementation must be done!
		self.__matrix_dict['Q'] = \
      			self.__matrix_dict['Pu'].transpose()*\
      			self.__matrix_dict['B_femWeight']*\
             	self.__matrix_dict['Pu']
        # self.__matrix_dict['Q'] = \
		# 		self.__matrix_dict['Pu'].transpose()*\
        # 		self.__matrix_dict['B']*\
		# 		self.__matrix_dict['Pu']
		self.__matrix_dict['Q'] = self.__matrix_dict['Q'].tocsc()
		nu = min(np.shape(self.__matrix_dict['Pu']))
		printDebug(True, '-- Done.')

		tic_res = time.perf_counter()

		if self.__param.Numerics.nCPU > 1:
			# Parallel computation of the resolvent --------------------------------
			printDebug(True, '-- Parallel computation of forcing and gains...')
			pool = multiprocessing.Pool(self.__param.Numerics.nCPU)

			args_map = [(linearSystem, 'Resolvent', self.__matrix_dict, \
				self.__param.Numerics.nSolut, arg) for arg in \
				self.__param.IOResolvent.Omegas]
			#pool.map(run_in_parallel, args_map)
			resultsPool = pool.map(self.runInParallel, args_map)

			# The gains, forcing, and responses are obtained in the parallel loop
			for i in range(len(resultsPool)):
				gains[:, i] = resultsPool[i][0]
				forcings[:, :, i] = resultsPool[i][1]
				responses[:, :, i] = resultsPool[i][2]

		else:
			# Serial computation of the resolvent --------------------------------
			printDebug(True, '--------------------------------')
			printDebug(True, '-- Serial computation of forcing, gains, and responses...')
   
			# LUQ is obtained here as it cannot be pickled
			self.__matrix_dict['LUQ'] = \
				splin.splu(self.__matrix_dict['Q'], permc_spec=3)

			for i in range(self.__n_omegas):
				
				# Frequency 
				omega = self.__param.IOResolvent.Omegas[i]
				printDebug(True, "-- Performing resolvent analysis for omega = " + str(omega))

				def op(x):
					'''
					Build the linear system which is used to get the 
					eigenvalue problem matrix of the resolvent:

     				'''
					y = self.__matrix_dict['Pu']*x
					z = self.__matrix_dict['B_forcing']*y
					y = LU.solve(z)
					y1 = self.__matrix_dict['Cr']*y 			# added line debug response limitation
					z = self.__matrix_dict['B_response']*y1
					z1 = self.__matrix_dict['Cr'].transpose()*z	# added line debug response limitation
					y = LU.solve(z1, trans='H')
					z = self.__matrix_dict['B_forcing'].transpose()*y
					y = self.__matrix_dict['Pu'].transpose()*z
					w = self.__matrix_dict['LUQ'].solve(y, trans='H')
					return w
				
				OP = self.__matrix_dict['A']-omega*self.__matrix_dict['B']
				OP = OP.tocsc()
				
				# Get lower upper decomposition of OP (used in function op())
				LU = splin.splu(OP,permc_spec=3)
				
				# Create handle for the linear operator defined by function op()
				SOP = splin.LinearOperator((nu,nu),matvec=op,dtype='complex')
				
				# Perform eigenvalue decomposition of the linear
				# operator defined in op() using the handl SOP
				gains[:,i],eigenvectors_c = splin.eigs(SOP,
									k=self.__param.Numerics.nSolut,
									M=None,
									sigma=None,
									which='LM',
									maxiter=100,
									tol=10-12,
									return_eigenvectors=True,
									)

				# Write gains to results dictionary
				gains[:, i] = np.real(gains[:, i])

				# Iterate through the first nSolut gains
				# Maybe we don't need the loop here...
				for k in range(self.__param.Numerics.nSolut):
					# Write the respective forcing to results dictionary
					forcings[:, k, i] = self.__matrix_dict['Pu'] * eigenvectors_c[:, k]
					f = -1j * self.__matrix_dict['B_femWeight'] * forcings[:, k, i]
					# f = -1j * forcings[:, k, i]
					responses[:, k, i] = LU.solve(f)

		#pdb.set_trace()
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
		"""
		Computes the Input-Output Analysis

		Function arguments:

		Function returns:
		- fluctSolutObjList: List of fluctuation Objects, where each Object
		contains one eigenvector and eigenvalue.

		"""
		fluctSolutObjList = []

		gains = np.zeros((1,self.__n_omegas))
		responses = np.zeros((self.__n_dof,1,self.__n_omegas),'complex')

		n_omegas=len(self.__param.IOResolvent.Omegas)

		self.__matrix_dict['A'], self.__matrix_dict['B'], \
		self.__matrix_dict['b_forcing'] = self.__preconditionMatrices(
			self.__matrix_dict['A'],
			self.__matrix_dict['B'],
			 self.__param.Numerics.Preconditioner,
			f=self.__matrix_dict['b_forcing'])

		if self.__param.Numerics.nCPU > 1:

			pool=multiprocessing.Pool(processes=self.__param.Numerics.nCPU)
			#func= partial(self.__ParallelInputOutputPickle,self.__matrix_dict)
			args_map = [[linearSystem, 'InputOutput', self.__matrix_dict, \
				arg] for arg in self.__param.IOResolvent.Omegas]

			resultsPool = pool.map(self.runInParallel, args_map)
			for i in range(len(resultsPool)):
				gains[:,i]	= 1
				responses[:,0,i] = resultsPool[i]

		else:
			## Repeat analysis for every omega
			for i in range(n_omegas):
				omega=self.__param.IOResolvent.Omegas[i]
				print("Performing input-output analysis for omega="+str(omega))
				### Get linear operator
				# Define OP
				OP=self.__matrix_dict['A']-omega*self.__matrix_dict['B']

				# Make OP sparse vector
				OP=OP.tocsc()
				# Perform Lower-Upper decomposition
				LU = splin.splu(OP,permc_spec=3)
				# Use LU decomposition to solve the linear system
				eigenvectors_c=LU.solve(self.__matrix_dict['b_forcing'])
				gains[:,i] = 1
				responses[:,0,i] = eigenvectors_c

				#z = results['forcings'][:,0,i]
				#results['responses'][:,0,i]=-LU.solve(z)
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
		- adjoint: Bool-Variable, speficing if the Direct- or the
		adjoint Eigenvalue-Problem is solved. The Default Value is false,
		which solves the Direct-problem

		Function returns:
		- __results_GEVP: Dictionary containing the direct- and the adjoint
		eigenvalues and eigenvectors
		"""


		start= time.time()

		EVal = np.zeros((self.__param.Numerics.nSolut*len(self.__param.Numerics.EigenValueGuess)), 'complex')
		EVec = np.zeros((self.__n_dof, self.__param.Numerics.nSolut*len(self.__param.Numerics.EigenValueGuess)), 'complex')


		if self.__param.Numerics.LinearAlgebraSolver=='python':
			# Added possibility to run GEVP of different guesses in parallel
			if self.__param.Numerics.nCPU > 1:
				print("-- Entering parallel loop for GEVP")
				pool=multiprocessing.Pool(processes=self.__param.Numerics.nCPU)
				args_map = [(linearSystem, 'GEVP', self.__matrix_dict, \
					self.__param.Numerics.nSolut, arg, adjointFlag, self.__param.Numerics) for arg in \
					self.__param.Numerics.EigenValueGuess]
				results_pool = pool.map(self.runInParallel, args_map)

				for i in range(len(results_pool)):
					index = list(range(i*self.__param.Numerics.nSolut, (i+1)*self.__param.Numerics.nSolut))
					EVal[index] = results_pool[i][0]
					EVec[:, index] = results_pool[i][1]

				print("-- Assembled results from all guesses.")

			else:

				EVal, EVec = self.__solve_with_python(adjointFlag)

		elif self.__param.Numerics.LinearAlgebraSolver in \
		['matlab', 'matrix export']:

			if self.__param.Numerics.LinearAlgebraSolver == 'matrix export':
				EVal, EVec = self.__solve_with_matlab(adjointFlag, True)
			else:
				EVal, EVec = self.__solve_with_matlab(adjointFlag, False)

		else:
			printError(self.__param.Numerics.LinearAlgebraSolver+ \
			' is not a valid value for the linearAlgebraSolver flag!')
			exit()

		# normalize the solution
		for i in range(self.__param.Numerics.nSolut):
			EVec[:,i] = EVec[:,i]/np.linalg.norm(EVec[:,i])

		# check the residuum
		maxRes = 0

		fluctSolutObjList = []


		if adjointFlag:
			method = "Adjoint"
		else:
			method = "Direct"

		for i in range(0,len(EVal)):

			# construct for each EVal and EVec a fluctuationSolution
			fluctSolutObjList.append(fluctuationSolutions(
									self.__param,
									self.__meanFlow,
									self.__FEMSpaces,
									EVal[i],
									EVec[:, i],
									not adjointFlag,
									 )
			)

			resTemp = self.__checkGEVP(
					EVal[i],
					EVec[:,i],
					self.__matrix_dict['A'],
					self.__matrix_dict['B'],
					)
			if resTemp > maxRes:
				maxRes = resTemp

		end = time.time() - start

		printDebug(True, '-- Solving the GEVP took %4g s' % end)
		printDebug(True, '-- Max residuum of %s solutions (EUCLIDIAN norm): %12g' % (method,maxRes))

		return fluctSolutObjList


	def checkMatrix(self, eq='ux', searchPoint=[1, 2]):
		import csv
		"""
        This function is used for debugging. It can be used to examine matrix A
        at a given search point for a given equation.
        Inputs:
        - eq: equation - ux, uy, ut, p
        - searchPoint: search point as tuple [x,y]
        It prints
        - the diagonal element at calculated index
        - the matrix row corresponding to coordinate and variable
        - the amount of diagonal elements in the whole matrix that are > 10^20
        - the amount of diagonal elements belonging to the examined subspace
        that are > 10^20 and the coordinates they are located at
        Output:
        - returns index of matrix that corresponds to equation and closest 
        coordinate
        """

		checkMat = self.__matrix_dict['A']

		# GENERAL ATTRIBUTES
		mesh = self.__FEMSpaces.P2.mesh
		gdim = mesh.geometry.dim
		if eq == 'ux':
			dofs = self.__FEMSpaces.VMixed.sub(0).sub(0).collapse()[1]
		elif eq == 'uy':
			dofs = self.__FEMSpaces.VMixed.sub(0).sub(1).collapse()[1]
		elif eq == 'ut':
			dofs = self.__FEMSpaces.VMixed.sub(1).collapse()[1]
		elif eq == 'p' :
			dofs = self.__FEMSpaces.VMixed.sub(2).collapse()[1]
		else:
			print('equation not found', file=sys.stderr)
			return

		# DOF COORDINATES OF MIXED SPACE
		nDofsVmixed = Function(self.__FEMSpaces.VMixed).vector[:].shape[0]
		dofs_coord = np.zeros((nDofsVmixed, gdim))
		for i in range(self.__FEMSpaces.VMixed.num_sub_spaces):
			if self.__FEMSpaces.VMixed.sub(i).num_sub_spaces > 0:
				for j in range(self.__FEMSpaces.VMixed.sub(i).num_sub_spaces):
					dofCoordsOfSubSpace = \
						self.__FEMSpaces.VMixed.sub(i).sub(j).collapse()[
							0].tabulate_dof_coordinates()[:, 0:gdim]
					indices = self.__FEMSpaces.VMixed.sub(i).sub(j).collapse()[
						1]
					dofs_coord[indices] = dofCoordsOfSubSpace
			else:
				dofCoordsOfSubSpace = self.__FEMSpaces.VMixed.sub(i).collapse()[
										  0].tabulate_dof_coordinates()[:,
									  0:gdim]
				indices = self.__FEMSpaces.VMixed.sub(i).collapse()[1]
				dofs_coord[indices] = dofCoordsOfSubSpace
		dofs_coord = dofs_coord.reshape((-1, gdim))

		# FIND THE CLOSEST COORDINATE TO SEARCH POINT
		# FINDS ALL FOUR INDICES FOR GIVEN SEARCH POINT
		nearestIndex = np.where(list(
			map(lambda x: np.linalg.norm(x - searchPoint),
				dofs_coord)) == min(list(
			map(lambda x: np.linalg.norm(x - searchPoint),
				dofs_coord))))
		#print("closest coordinate: ", dofs_coord[nearestIndex[0][0]], file=sys.stderr)
		# CHECKS WHICH INDEX BELONGS TO EXAMINED SUBSPACE
		matchingIndices = np.empty(shape=(0, 0))
		for index in nearestIndex[0]:
			if index in dofs:
				matchingIndices = np.append(matchingIndices, index)

		# PRINT MATRIX ROW CORRESPONDING TO COORDINATE
		for index in matchingIndices:
			#print("line ", index, " in A: ", np.shape(checkMat[int(index)]), file=sys.stderr)
			nonzeroList = checkMat[int(index)].nonzero()
			row = []
			for jndex in nonzeroList[1]:
				row.append(checkMat[int(index), jndex])

		# PRINT EVERY MATRIX ELEMENT IN EVERY ROW CORRESPONDING TO EQUATION
		maxNonZeroNumb = 0
		for kindex in dofs:
			line = []
			line.append(dofs_coord[kindex][0])
			line.append(dofs_coord[kindex][1])
			nonzeroList = checkMat[kindex].nonzero()
			if len(nonzeroList[1]) > maxNonZeroNumb:
				maxNonZeroNumb = len(nonzeroList[1])
			row = []
			for lindex in nonzeroList[1]:
				row.append(checkMat[kindex, lindex])
			row.sort()
			line.extend(row)

		# HOW MANY DIAGONAL ELEMENTS ARE > 10^20 IN THE WHOLE MATRIX?
		diagonal = checkMat.diagonal()
		homogDiric = sum(x > 10e20 for x in diagonal)
		#print("Amount of diagonal elements in whole matrix > 10^20: ", homogDiric, file=sys.stderr)

		# HOW MANY DIAGONAL ELEMENTS BELONGING TO THE EXAMINED SUBSPACE ARE
		# > 10^20 AND AT WHAT COORDINATES ARE THEY LOCATED?
		counter = 0
		coordinates = []
		for index in dofs:
			if checkMat[index, index] > 10e20:
				counter += 1
				coordinates.append(tuple(dofs_coord[index, :]))
		#print("Amount of diagonal elements in examined subspace > 10^20: ", counter, file=sys.stderr)
		# print("Coordinates corresponding to values > 10^20 ", coordinates, file=sys.stderr)

		return
