#/* Copyright (C) 2019 FLOW group TU Berlin - All Rights Reserved
# * You may NOT use, distribute or modify this code without explicit
# * permission of the copyright owner, the FLOW group at TU Berlin!
# * However, permissions to use and modify the code are generally
# * granted when asked for.
# * To ask for permission please contact t.kaiser@tu-berlin.de.
# */
'''
# **********************************************************************
# * This file deals with all FEM Discretization and solution
# *
# * This file was created by Thomas L. Kaiser. Significant contributions
# * were made by
# * -
# *
# *
# ********************
'''



from scipy.sparse import csr_matrix,csc_matrix
from fenics import dx,ds,assemble, PETScMatrix,TestFunctions,TrialFunction,as_backend_type,SpatialCoordinate,Function,project,MeshFunction,CompiledSubDomain,Measure,FunctionAssigner,inner,TrialFunctions,PETScVector,rhs,lhs,DirichletBC
import WeakFormulationCollection
import scipy.sparse.linalg as splin
import numpy as np
import multiprocessing
from functions import *
import os





def solveGEVP(param, FEMSpaces):
	''' This function solves the modal problem's generalized eigenvalue problem ( A x = B x lambda ). Matrices are taken from the global_variables module... '''
	import time
	import global_variables as glob
	# build results dictionary
	n_dof=np.shape(glob.matrix_dict['A'])[0]
	results_GEVP={\
		"EValDirect":	np.zeros((param.Numerics.nSolut),'complex'),\
		"EValAdjoint":	np.zeros((param.Numerics.nSolut),'complex'),\
		"EVecDirect":	np.zeros((n_dof,param.Numerics.nSolut),'complex'),\
		"EVecAdjoint":	np.zeros((n_dof,param.Numerics.nSolut),'complex'),\
		}
	#the tolerance (tol) for the solver is hard coded
	tol=1e-12
	start= time.time()

	if param.Numerics.LinearAlgebraSolver=='python':

		print("Solving direct GEVP")

		#precondition direct matrices
		A, B, f = preconditionMatrices(glob.matrix_dict['A'],
			glob.matrix_dict['B'],
			param.Numerics.Preconditioner)

		# solve direct GEVP using eigs
		results_GEVP["EValDirect"][:], results_GEVP["EVecDirect"][:,:] = splin.eigs(A,
			k=param.Numerics.nSolut,
			M=B,
			sigma=param.Numerics.EigenValueGuess,
			ncv=200,
			maxiter=100,
			tol=10e-12,
			return_eigenvectors=True)

		print("Solving adjoint GEVP")

		#precondition adjoint matrices
		A, B, f = preconditionMatrices(glob.matrix_dict['A'].getH(),
			glob.matrix_dict['B'],
			param.Numerics.Preconditioner)

		#solve adjoint GEVP using eigs
		results_GEVP["EValAdjoint"][:], results_GEVP["EVecAdjoint"][:,:] = splin.eigs(A,
			k=param.Numerics.nSolut,
			M=B,
			sigma=param.Numerics.EigenValueGuess,
			ncv=200,
			maxiter=100,
			tol=tol,
			return_eigenvectors=True)

	elif param.Numerics.LinearAlgebraSolver in ['matlab', 'matrix export']:
		import matlab.engine
		eng = matlab.engine.start_matlab()
		eng.addpath (__file__.rsplit('/',1)[0]+'/Matlab', nargout= 0 )
		from scipy.io import savemat

		#precondition the direct matrices
		A, B, f = preconditionMatrices(glob.matrix_dict['A'],
			glob.matrix_dict['B'],
			param.Numerics.Preconditioner)

		# save matrices to file
		savemat('A.mat',{'A':A})
		savemat('B.mat',{'B':B})

		if param.Numerics.LinearAlgebraSolver in ['matrix export']:
			print('Matrix exported. Ending program...')
			exit()

		#run the GEVP solver in the matlab script
		temp = eng.MatlabEigs(param.Numerics.nSolut,
			param.Numerics.EigenValueGuess,
			200,
			100,
			10e-11,
			False,
			nargout=2)

		#remove matrices
		os.remove("A.mat")
		os.remove("B.mat")

		#write solutions to the results dictionary
		results_GEVP["EValDirect"][:]=np.array(temp[0])
		results_GEVP["EVecDirect"][:,:]=np.array(temp[1])

		print("Solving adjoint GEVP")
		# precondition matrices
		A, B, f = preconditionMatrices(glob.matrix_dict['A'].getH(),
			glob.matrix_dict['B'],
			param.Numerics.Preconditioner)

		#save matrices to file
		savemat('A.mat',{'A':A})
		savemat('B.mat',{'B':B})

		#run the GEVP solver in the matlab script
		temp = eng.MatlabEigs(param.Numerics.nSolut,
			param.Numerics.EigenValueGuess,
			200,
			100,
			1e-12,False,nargout=2)

		#write solutions to the results dictionary
		results_GEVP["EValAdjoint"][:]=np.array(temp[0])
		results_GEVP["EVecAdjoint"][:,:]=np.array(temp[1])

		# Delete files
		os.remove("A.mat")
		os.remove("B.mat")
	else:
		printError(param.Numerics.LinearAlgebraSolver+ ' is not a valid value for the linearAlgebraSolver flag!')
		exit()
	for i in range(param.Numerics.nSolut):
		results_GEVP["EVecDirect"][:,i]=results_GEVP["EVecDirect"][:,i]/np.linalg.norm(results_GEVP["EVecDirect"][:,i])
		results_GEVP["EVecAdjoint"][:,i]=results_GEVP["EVecAdjoint"][:,i]/np.linalg.norm(results_GEVP["EVecAdjoint"][:,i])
	# Checking for correctness of results by calculating residuals
	#Normalize eigen vectors
	for i in range(param.Numerics.nSolut):
		results_GEVP["EVecDirect"][:,i]=results_GEVP["EVecDirect"][:,i]/np.linalg.norm(results_GEVP["EVecDirect"][:,i])
		results_GEVP["EVecAdjoint"][:,i]=results_GEVP["EVecAdjoint"][:,i]/np.linalg.norm(results_GEVP["EVecAdjoint"][:,i])
	maxResDir=0
	maxResAdj=0
	for i in range(0,len(results_GEVP["EValDirect"])):
		resTemp=checkGEVP(results_GEVP["EValDirect"][i],results_GEVP["EVecDirect"][:,i],glob.matrix_dict['A'],glob.matrix_dict['B'])
		if resTemp>maxResDir:
			maxResDir=resTemp
		resTemp=checkGEVP(results_GEVP["EValDirect"][i],results_GEVP["EVecDirect"][:,i],glob.matrix_dict['A'],glob.matrix_dict['B'])
		if resTemp>maxResAdj:
			maxResAdj=resTemp
	end=time.time()
	print('Solving the GEVP took '+str(end)+'s')
	print('The eigensolvers tolerance is set to '+str(tol)+'...')
	print('Maximum residuum of DIRECT solutions measured in EUCLIDIAN norm: ' +str(maxResDir))
	print('Maximum residuum of ADJOINT solutions measured in EUCLIDIAN norm: ' +str(maxResAdj))


	return results_GEVP



def solveInputOutput(param,FEMSpaces):
	import global_variables as glob
	''' This function solves the Input-Output Analysis. Matrices are taken from the global_variables module... '''
	# Get number of omegas
	n_omegas=len(param.IOResolvent.Omegas)
	#the size of the matrix, A, correspond to the number of degrees of freedom, n_dof
	n_dof=np.shape(glob.matrix_dict['A'])[0]
	#response     = np.zeros([n_dof,1],'complex')
	# Define results
	results={\
		"gains":		np.zeros((1,n_omegas)),\
		"responses":	np.zeros((n_dof,1,n_omegas),'complex'),\
		}
	### Build forcing vector
	# Initialize forcing (fMixed) in mixed function space


	#fMixed=Function(FEMSpaces.VMixed)
	#fMixed_r=Function(FEMSpaces.VMixed)
	#fMixed_i=Function(FEMSpaces.VMixed)





	glob.matrix_dict['A'], glob.matrix_dict['B'], glob.matrix_dict['b_forcing'] = preconditionMatrices(glob.matrix_dict['A'],
		glob.matrix_dict['B'],
		param.Numerics.Preconditioner,
		f=glob.matrix_dict['b_forcing'])

	if param.Numerics.nCPU>1:
		from functools import partial
		pool=multiprocessing.Pool(processes=param.Numerics.nCPU)
		#func= partial(ParallelInputOutput)
		func= partial(ParallelInputOutputPickle,glob.matrix_dict)
		results_pool=pool.map(func, param.IOResolvent.Omegas)
		for i in range(len(results_pool)):
			results['gains'][:,i]=1
			results['responses'][:,0,i]=results_pool[i]

	else:
		## Repeat analysis for every omega
		for i in range(n_omegas):
			omega=param.IOResolvent.Omegas[i]
			print("Performing input-output analysis for omega="+str(omega))
			### Get linear operator
			# Define OP
			OP=glob.matrix_dict['A']-omega*glob.matrix_dict['B']
			from scipy.io import savemat
			savemat('OP.mat',{'OP':OP})
			savemat('f.mat',{'f':glob.matrix_dict['b_forcing']})
			# Make OP sparse vector
			OP=OP.tocsc()
			# Perform Lower-Upper decomposition
			LU = splin.splu(OP,permc_spec=3)
			# Use LU decomposition to solve the linear system
			eigenvectors_c=LU.solve(glob.matrix_dict['b_forcing'])
			results['gains'][:,i]=1
			#results_resolvent['omegas'][:,i]=omega
			results['responses'][:,0,i]= eigenvectors_c
			#z = results['forcings'][:,0,i]
			#results['responses'][:,0,i]=-LU.solve(z)
	return results

def ParallelInputOutputPickle(matrix_dict,omega):
	''' This function is used for distributing the input output analysis to various CPUs.'''
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

	return 	eigenvectors_c
#def ParallelInputOutput(omega):
#	''' This function is used for distributing the input output analysis to various CPUs.'''
#	import global_variables as glob
#	print("Performing input-output analysis for omega="+str(omega))
#	### Get linear operator
#	# Define OP
#	OP=glob.matrix_dict['A']-omega*glob.matrix_dict['B']
#	# Make OP sparse vector
#	OP=OP.tocsc()
#	# Perform Lower-Upper decomposition
#	LU = splin.splu(OP,permc_spec=3)
#	# Use LU decomposition to solve the linear system
#	eigenvectors_c=LU.solve(glob.matrix_dict['b_forcing'])
#
#	return 	eigenvectors_c

def solveResolvent(param,FEMSpaces):
	import time
	import global_variables as glob
	''' This function solves the resolvent prompblem's eigenvalue problem. Matrices are taken from the global_variables module... '''
	glob.matrix_dict['B_forcing']=glob.matrix_dict['B_forcing'].tocsc()
	glob.matrix_dict['B_response']=glob.matrix_dict['B_response'].tocsc()
	glob.matrix_dict['B']= glob.matrix_dict['B'].tocsc()
	glob.matrix_dict['Pu']=getPMat(param,FEMSpaces)
	glob.matrix_dict['Q']= glob.matrix_dict['Pu'].transpose()*glob.matrix_dict['B_response']*glob.matrix_dict['Pu']
	glob.matrix_dict['Q']= glob.matrix_dict['Q'].tocsc()
	#glob.matrix_dict['LUQ']=splin.splu(glob.matrix_dict['Q'],permc_spec=3) # Shifted to a later point
	nu=min(np.shape(glob.matrix_dict['Pu']))

	n_omegas=len(param.IOResolvent.Omegas)
	gains=np.zeros((param.Numerics.nSolut,n_omegas),'complex')
	n_dof=np.shape(glob.matrix_dict['A'])[0]
	forcing=np.zeros((n_dof,param.Numerics.nSolut,n_omegas))
	response=np.zeros((n_dof,param.Numerics.nSolut,n_omegas))
	forcing     = np.zeros([n_dof,param.Numerics.nSolut],'complex')
	response     = np.zeros([n_dof,param.Numerics.nSolut],'complex')
	results_resolvent={\
	    "gains":        np.zeros((param.Numerics.nSolut,n_omegas)),\
	    "forcings": np.zeros((n_dof,param.Numerics.nSolut,n_omegas),'complex'),\
	    "responses":    np.zeros((n_dof,param.Numerics.nSolut,n_omegas),'complex'),\
	    }
	if param.Numerics.nCPU>1:
		from functools import partial
		pool=multiprocessing.Pool(param.Numerics.nCPU)
		#func= partial(ParallelResolvent,param.Numerics.nSolut)
		func= partial(ParallelResolventPickle,param.Numerics.nSolut,glob.matrix_dict)
		results_pool=pool.map(func, param.IOResolvent.Omegas)
		for i in range(len(results_pool)):
			results_resolvent['gains'][:,i]=results_pool[i][0]
			results_resolvent['forcings'][:,:,i]=results_pool[i][1]
			results_resolvent['responses'][:,:,i]=results_pool[i][2]
	else:
		glob.matrix_dict['LUQ']=splin.splu(glob.matrix_dict['Q'],permc_spec=3)
		for i in range(n_omegas):
			def op(x):
				y=glob.matrix_dict['Pu']*x
				z=glob.matrix_dict['B_forcing']*y
				y=LU.solve(z)
				z=glob.matrix_dict['B_response']*y
				y=LU.solve(z,trans='H')
				z=glob.matrix_dict['B_forcing'].transpose()*y
				y=glob.matrix_dict['Pu'].transpose()*z
				w=glob.matrix_dict['LUQ'].solve(y,trans='H')
				return w
			omega=param.IOResolvent.Omegas[i]
			print("Performing resolvent analysis for omega="+str(omega))

			# Built resolvent matrix OP
			OP=glob.matrix_dict['A']-omega*glob.matrix_dict['B']
			OP=OP.tocsc()

			# Get lower upper decomposition of OP (used in function op())
			LU = splin.splu(OP,permc_spec=3)

			# Create handle for the linear operator defined by function op()
			SOP       = splin.LinearOperator((nu,nu),matvec=op,dtype='complex')

			# Perform eigenvalue decomposition of the linear operator defined in op() using the handl SOP
			gains[:,i],eigenvectors_c      = splin.eigs(SOP, k=param.Numerics.nSolut, M=None, sigma=None, which='LM', maxiter=100, tol=10-12, return_eigenvectors=True)

			# Write gains to results dictionary
			results_resolvent['gains'][:,i]=np.real(gains[:,i])

			# If Cylindrical coordinates are chosen, the radial coordinate
			# must be taken into account, if not its set to 1
			if param.Case.CoordinateSystem == 'Cylindrical':
				R=SpatialCoordinate(FEMSpaces.P2.mesh())[1]
			else:
				R=1

			# Define test functions
			X=TestFunctions(FEMSpaces.VMixed)
			# Iterate through the first nSolut gains
			for k in range(param.Numerics.nSolut):
				# Write the respective forcing to tesults dictionary
				results_resolvent['forcings'][:,k,i]= glob.matrix_dict['Pu']*eigenvectors_c[:,k]

				# Calculate the responses to the forcings by applying the resolvent operator
				# To do so the rhs forcing vector must be built
				eigenvectors_real = Function(FEMSpaces.VMixed)
				eigenvectors_real.vector()[:] = results_resolvent['forcings'][:,k,i].real.astype(float)
				eigenvectors_imag = Function(FEMSpaces.VMixed)
				eigenvectors_imag.vector()[:] = results_resolvent['forcings'][:,k,i].imag.astype(float)

				# Define the variational formulations for the forcings real and imaginary part
				f_real_vf=0
				f_imag_vf=0

				# Loop through all solutions
				for name,i_name in zip(param.SolutionList,range(param.Numerics.nSolut)):
					# If it is the velocity, all velocity components must be considered
					if name == 'u':
						# Loop through the velocity components
						for (component,i_component) in zip(param.VelocityComponents,range(param.nVelocityComponents)):
							# Add the imaginary part and the real part for the respective component to the respective equation
							# Here, the multiplication by '-i' is taken into account. Could also be added later, when the
							# linear system is solved
							f_real_vf += R*eigenvectors_imag.split()[i_name].split()[i_component] \
									             *X[i_name][i_component] *dx
							f_imag_vf += -R*eigenvectors_real.split()[i_name].split()[i_component] \
									             *X[i_name][i_component] *dx
					else:
						# Add the imaginary part and the real part for the respective component to the respective equation
						# Here, the multiplication by '-i' is taken into account. Could also be added later, when the
						# linear system is solved
						f_real_vf += R*eigenvectors_imag.split()[i_name] \
									              *X[i_name] *dx
						f_imag_vf += -R*eigenvectors_real.split()[i_name] \
									               *X[i_name] *dx

				# Define the vectors, and assemple them using the variational formulation defined above
				f_real = PETScVector()
				f_imag = PETScVector()
				assemble(f_real_vf,tensor=f_real)
				assemble(f_imag_vf,tensor=f_imag)

				# get the local values of the vector
				f = f_real.get_local() + 1j*f_imag.get_local()

				# Solve the linear system, with the forcing vector on the rhs
				results_resolvent['responses'][:,k,i]=LU.solve(f)
	return results_resolvent

def ParallelResolventPickle(nSolut,matrix_dict,omega):
    ''' This function is used for distributing the resolvent's eigenvalueproblems to various CPUs.'''
    LUQ=splin.splu(matrix_dict['Q'],permc_spec=3)
    nu=min(np.shape(matrix_dict['Pu']))
    print (multiprocessing.current_process())
    def op(x):
        y=matrix_dict['Pu']*x
        z=matrix_dict['B_forcing']*y
        y=LU.solve(z)
        z=matrix_dict['B_response']*y
        y=LU.solve(z,trans='H')
        z=matrix_dict['B_forcing'].transpose()*y
        y=matrix_dict['Pu'].transpose()*z
        w=LUQ.solve(y,trans='H')
        return w
    print("CPU "  + " performing resolvent analysis for omega="+str(omega))
    OP=matrix_dict['A']-omega*matrix_dict['B']
    OP=OP.tocsc()
    LU = splin.splu(OP,permc_spec=3)
    SOP       = splin.LinearOperator((nu,nu),matvec=op,dtype='complex')

    gains,eigenvectors_c      = splin.eigs(SOP, k=nSolut, M=None, sigma=None, which='LM', maxiter=100, tol=10-12, return_eigenvectors=True)
    for k in range(nSolut):
        forcings= matrix_dict['Pu']*eigenvectors_c
        z = forcings
        responses=-LU.solve(z)
    return gains,forcings,responses

#def ParallelResolvent(nSolut,omega):
#    ''' This function is used for distributing the resolvent's eigenvalueproblems to various CPUs.'''
#    import global_variables as glob
#    nu=min(np.shape(glob.matrix_dict['Pu']))
#    print (multiprocessing.current_process())
#    def op(x):
#        y=glob.matrix_dict['Pu']*x
#        z=glob.matrix_dict['B_forcing']*y
#        y=LU.solve(z)
#        z=glob.matrix_dict['B_response']*y
#        y=LU.solve(z,trans='H')
#        z=glob.matrix_dict['B_forcing'].transpose()*y
#        y=glob.matrix_dict['Pu'].transpose()*z
#        w=glob.matrix_dict['LUQ'].solve(y,trans='H')
#        return w
#    print("CPU "  + " performing resolvent analysis for omega="+str(omega))
#    OP=glob.matrix_dict['A']+omega*glob.matrix_dict['B']
#    OP=OP.tocsc()
#    LU = splin.splu(OP,permc_spec=3)
#    SOP       = splin.LinearOperator((nu,nu),matvec=op,dtype='complex')
#
#    gains,eigenvectors_c      = splin.eigs(SOP, k=nSolut, M=None, sigma=None, which='LM', maxiter=100, tol=10-12, return_eigenvectors=True)
#    for k in range(nSolut):
#        forcings= glob.matrix_dict['Pu']*eigenvectors_c
#        z = forcings
#        responses=-LU.solve(z)
#    return gains,forcings,responses


def DiscretizeFlow(param,MF,FEMSpaces,mesh):

	''' This function chooses the right weak formulation and then discretizes it on the given mesh'''
	from copy import deepcopy
	# Check for type of case
	WeakForm=WeakFormulationCollection.WeakFormulationCollectionClass(param,FEMSpaces,MF)
	# Definie matrices
	A_real = PETScMatrix()
	A_imag = PETScMatrix()
	B_real = PETScMatrix()
	B_imag = PETScMatrix()
	BC_Diriclet = PETScMatrix()
	# Assamble the weak formulation to the matrices
	assemble(WeakForm.B_real_vf.lhs,tensor=BC_Diriclet)
	BC_Diriclet.zero()
	n_dof=BC_Diriclet.size(0)
	if not WeakForm.A_imag_vf.lhsIsZero():
		assemble(WeakForm.A_imag_vf.lhs,tensor=A_imag)
	else:
		A_imag=0*BC_Diriclet
	if not WeakForm.A_real_vf.lhsIsZero():
		assemble(WeakForm.A_real_vf.lhs,tensor=A_real)
	else:
		A_real=0*BC_Diriclet
	if not WeakForm.B_real_vf.lhsIsZero():
		assemble(WeakForm.B_real_vf.lhs,tensor=B_real)
	else:
		B_real=0*BC_Diriclet
	if not WeakForm.B_imag_vf.lhsIsZero():
		assemble(WeakForm.B_imag_vf.lhs,tensor=B_imag)
	else:
		B_imag=0*BC_Diriclet
	matrix_dict={}
	#Get the Dirichlet BCs set by the user in a list...
	bcs= getListOfDirichletBCs(param.BCs.getBCsDict(),
		FEMSpaces,
		param.BCs.getBoundaries(),
		param.Case,
		param.debug)
	if param.Case.AnalysisMode in ['Input-Output']:
		if not WeakForm.A_imag_vf.rhsIsZero():
			forcing_vec_i_petsc=assemble(WeakForm.A_imag_vf.rhs)
		else:
			forcing_vec_i_petsc=PETScVector()
			forcing_vec_i_petsc.init(n_dof)
		if not WeakForm.A_real_vf.rhsIsZero():
			forcing_vec_r_petsc=assemble(WeakForm.A_real_vf.rhs)
		else:
			forcing_vec_r_petsc=PETScVector()
			forcing_vec_r_petsc.init(n_dof)
		for bc in bcs:
			bc.apply(A_imag,forcing_vec_i_petsc)
			bc.apply(A_real,forcing_vec_r_petsc)


		forcing_vec_i=as_backend_type(forcing_vec_i_petsc).vec().array
		forcing_vec_r=as_backend_type(forcing_vec_r_petsc).vec().array
		b_forcing =  1j*forcing_vec_r  - forcing_vec_i
		matrix_dict['b_forcing'] =  b_forcing
		del b_forcing, forcing_vec_r, forcing_vec_i

	# Get the BCs provided by the user
	# Aplly the BCs
	## Boundary conditions are applied via penalisation method manually. This is in order to keep the rhs matrix invertible.
	else:
		bcs= getListOfDirichletBCs(param.BCs.getBCsDict(),
			FEMSpaces,
			param.BCs.getBoundaries(),
			param.Case,
			param.debug)
		for bc in bcs:
	#		printDebug(param.debug,"Applying following BC: "+str(bc))
			bc.apply(BC_Diriclet)
	# In the next 12 lines the Imaginary and real parts of both the lhs and rhs matrix are combined to the
	# sparse matrix A and B, respectively. Not needed matrices are deleted
	A_real_mat = as_backend_type(A_real).mat()
	matrix_dict['A'] = csr_matrix(A_real_mat.getValuesCSR()[::-1], shape = A_real_mat.size,dtype=complex)
	del A_real, A_real_mat
	A_imag_mat = as_backend_type(A_imag).mat()
	matrix_dict['A'] = matrix_dict['A'] + 1j * csr_matrix(A_imag_mat.getValuesCSR()[::-1], shape = A_imag_mat.size,dtype=complex)
	del A_imag, A_imag_mat
	if not param.Case.AnalysisMode in ['Input-Output']:
		BC_Diriclet_mat = as_backend_type(BC_Diriclet).mat()
		matrix_dict['A'] = matrix_dict['A'] + 10**30*(1+1j) * csr_matrix(BC_Diriclet_mat.getValuesCSR()[::-1], shape = BC_Diriclet_mat.size,dtype=complex)
		del BC_Diriclet, BC_Diriclet_mat
	B_real_mat = as_backend_type(B_real).mat()
	matrix_dict['B'] = csr_matrix(B_real_mat.getValuesCSR()[::-1], shape = B_real_mat.size,dtype=complex)
	del B_real, B_real_mat
	B_imag_mat = as_backend_type(B_imag).mat()
	matrix_dict['B'] = matrix_dict['B'] + 1j * csr_matrix(B_imag_mat.getValuesCSR()[::-1], shape = B_imag_mat.size,dtype=complex)
	del B_imag, B_imag_mat
	#tempMat=1*matrix_dict['B'].transpose()
	#tempMat[1,1]=100
	#input((matrix_dict['B'] != tempMat).nnz==0)
	#input(matrix_dict['B'])
	# must be constructed. So far only the L2 norm is implemented (To be extended!)
	if param.Case.AnalysisMode in ['Resolvent']:
		#forcing_norm,response_norm=getResolventNorms(param,MF,FEMSpaces)
		B_forcing_real = PETScMatrix()
		assemble(WeakForm.forcing_vf,tensor=B_forcing_real)
		#for bc in bcs:
		#	bc.zero(B_forcing_real)
		B_forcing_real_mat = as_backend_type(B_forcing_real).mat()
		matrix_dict['B_forcing'] = csr_matrix(B_forcing_real_mat.getValuesCSR()[::-1], shape = B_forcing_real_mat.size,dtype=complex)
		del B_forcing_real, B_forcing_real_mat
		#B_forcing_imag = PETScMatrix()
		#assemble(rhs_imag,tensor=B_forcing_imag)
		##for bc in bcs:
		##	bc.zero(B_forcing_imag)
		#B_forcing_imag_mat = as_backend_type(B_forcing_imag).mat()
		#matrix_dict['B_forcing'] = matrix_dict['B_forcing']+ 1j * csr_matrix(B_forcing_imag_mat.getValuesCSR()[::-1], shape = B_forcing_imag_mat.size)
		#del B_forcing_imag, B_forcing_imag_mat
		B_response_real = PETScMatrix()
		assemble(WeakForm.response_vf,tensor=B_response_real)
		#for bc in bcs:
		#	bc.zero(B_response_real)
		B_response_real_mat = as_backend_type(B_response_real).mat()
		matrix_dict['B_response'] = csr_matrix(B_response_real_mat.getValuesCSR()[::-1], shape = B_response_real_mat.size,dtype=complex)
		del B_response_real, B_response_real_mat
		#B_response_imag = PETScMatrix()
		#assemble(rhs_imag,tensor=B_response_imag)
		##for bc in bcs:
		##	bc.zero(B_response_imag)
		#B_response_imag_mat = as_backend_type(B_response_imag).mat()
		#matrix_dict['B_response'] = matrix_dict['B_response'] + 1j * csr_matrix(B_response_imag_mat.getValuesCSR()[::-1], shape = B_response_imag_mat.size)
		#del B_response_imag, B_response_imag_mat


	return matrix_dict

def getPMat(param,FEMSpaces):
	''' This function provides the P matrix, which restricts the forcing '''
	index=np.empty(shape=(0,0))
	mesh=FEMSpaces.P2.mesh()
	gdim = mesh.geometry().dim()
	for i in param.IOResolvent.ForcingCoeff:
		if i < param.nVelocityComponents:
			dofs=FEMSpaces.VMixed.sub(0).sub(i).dofmap()
			dofs_coord=FEMSpaces.VMixed.tabulate_dof_coordinates().reshape((-1, gdim))
			for index_local in dofs.dofs():
				if dofs_coord[index_local,0]< 999999999999999:## In order to restrict the forcing to a certain range, pass a posibility to the user to have access to this if clause and formulate it acordingly
					index=np.append(index,index_local)
		else:
			dofs=FEMSpaces.VMixed.sub(i-param.nVelocityComponents+1).dofmap()
			dofs_coord=FEMSpaces.VMixed.tabulate_dof_coordinates().reshape((-1, gdim))
			for index_local in dofs.dofs():
				if dofs_coord[index_local,0]< 999999999999999:## In order to restrict the forcing to a certain range, pass a posibility to the user to have access to this if clause and formulate it acordingly
					index=np.append(index,index_local)
		#index = np.append(index,dofs.dofs())
	index.sort()
	m=len(FEMSpaces.VMixed.dofmap().dofs())
	n=len(index)
	row_ind=index
	col_ind=np.arange(n)
	P=csr_matrix((np.ones(n),(row_ind,col_ind)),(m,n))
	return P


def getListOfDirichletBCs(bcDict,FEMSpaces,boundaries,Case,debug):
	BClist = []
	mesh=FEMSpaces.P2.mesh()

	VelocityComponents=Case.getInternalVelocityComponents()
	SolutionList=Case.getTransportedQuantityList()
	for k,m in zip(list(bcDict.keys()),range(0,len(bcDict.keys()))):
		# Get index of equation/variable i_eqn and if needed the index of the velocity component
		if k[0]=='u' and k[1] in VelocityComponents:
			i_eqn=0
			i_component=VelocityComponents.index(k[1])
		elif k in SolutionList:
			i_eqn = SolutionList.index(k)
		else:
			continue

		for mm in range(0,len(bcDict[k])) :
			if bcDict[k][mm]['type']=='Dirichlet':


				printDebug(debug,"Adding Dirichlet BC for "+str(k)+ " in equation "+str(i_eqn)+" with value "+str(bcDict[k][mm]['value'])+" on boundary with index "+str(bcDict[k][mm]['ID']))
				if k in ['u'+ component for component in VelocityComponents]:
					BClist.append( DirichletBC(FEMSpaces.VMixed.sub(i_eqn).sub(i_component), bcDict[k][mm]['value'], boundaries, bcDict[k][mm]['ID']) )


				else:
					if len(Case.getTransportedQuantityList()) == 1:
						BClist.append( DirichletBC(FEMSpaces.FunctionSpaceList[i_eqn], bcDict[k][mm]['value'], boundaries, bcDict[k][mm]['ID']) )
					else:
						BClist.append( DirichletBC(FEMSpaces.VMixed.sub(i_eqn), bcDict[k][mm]['value'], boundaries, bcDict[k][mm]['ID']) )

	return BClist


def checkGEVP(eigVal,eigVec, A, B):
	''' Calculating the frombenius norm and the infinite norm of the residuum of the solutions to a GEVP '''
	from scipy.sparse.linalg import norm as sparse_norm
	# First normalize the vector (additionally with the norm of matrix A)
	normA = sparse_norm(A)
	eigVecNorm = eigVec/np.linalg.norm(eigVec)/normA
	# Calculate residuum of Ax = omega Bx
	residuum_vec=A@eigVecNorm-eigVal*B@eigVecNorm
	residuum_euc = np.linalg.norm(residuum_vec)
	return residuum_euc

def preconditionMatrices(A,B,Preconditioner,f=False):
	''' Performing a preconditioning of the matrices A,B'''
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
