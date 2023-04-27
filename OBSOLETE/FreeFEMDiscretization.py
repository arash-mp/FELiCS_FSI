import os
import numpy as np
import time 
import Fields as Fields
try:
    import readmat
except ImportError:
	os.system('f2py3 -m readmat -c readmat.f90 --fcompiler=gnu95')
	import readmat
import scipy.sparse as sp
import scipy.sparse.linalg as splin

class FreeFEMDiscretization3():
	def __init__(self,param):
		self.LoadMatrices(param.MatricesDirectory)

	def LoadMatrices(self,directory):
		"""Here the matrices are read"""
		# Get File list of matrix directory
		filelist=os.listdir(directory)
		# If L.dat file exists in matrix Directory import L matrix
		if "L.dat" in filelist:
			self.L = self.LoadMatrix(directory+"/L.dat")	
		else:
			raise IOError("Cannot find L.dat in "+repr(directory))
		# If B.dat file exists in matrix Directory import B matrix
		if "B.dat" in filelist:
			self.B = self.LoadMatrix(directory+"/B.dat")	
		else:
			raise IOError("Cannot find B.dat in "+repr(directory))
			
	def LoadMatrix(self,name):
		print('Reading file',name)
		f = open(name, 'r')
		rr = f.readlines()
		f.close()
		
		# Read the matrix size
		str=rr[3]
		w=str.split(' ')
		#remove blanks
		data=[]
		for i in range(np.size(w)):
		    if len(w[i]) != 0:
		        data.append(w[i])
		
		n  =int(data[0])
		m  =int(data[1])
		nnz=int(data[3])
		print("  Matrix size:",n,'*',m,', nnz',nnz)
		
		# Determine if the matrix is real or complex
		str=rr[4]
		w=str.split(' ')
		data=[]
		for j in range(np.size(w)):
		    if len(w[j]) != 0:
		        data.append(w[j])
		tmp = data[2].split(',')
		if len(tmp) == 2:
		    print('  Type : complex')
		    icoo,jcoo,dcoo = readmat.readcomplexmat(name,nnz)
		else:
		    print('  Type : real')
		    icoo,jcoo,dcoo = readmat.readrealmat(name,nnz)
		icoo = icoo - 1
		jcoo = jcoo - 1
		
		ijcoo = [icoo,jcoo]
		# Create COO matrix
		mat=sp.coo_matrix((dcoo,ijcoo),shape=(n,m))
		del dcoo,ijcoo,icoo,jcoo
		# Convert to CSR format
		mat=mat.tocsc()
		
		return mat
	def solveGEVP(self,param,mesh):
		print("Finding solutions to generalized eigenvalue problem...")
		tic=time.clock()
		eigenvalues, eigenvectors = splin.eigs(self.L, k=param.nSolut, M=self.B, sigma=0.7, ncv=100, maxiter=100, tol=10e-10, return_eigenvectors=True)
		toc=time.clock()
		print(param.nSolut, "solutions to generalize eigenvalue problem found in",toc-tic,"seconds")
		print("The eigenvalues are:")
		for i in range(param.nSolut):
			print(i,": ",eigenvalues[i])
		Fluc = [Fields.Fluctuation(param,mesh) for i in range(param.nSolut+1)]
		print(np.shape(Fluc[0].U[0]))
		print(	range(2*mesh.np2+mesh.np1))
		for j in range(param.nSolut):
			print (j)
			indexP=0
			indexU1=0
			indexU2=0
			indexU3=0
			for i in range(2*mesh.np2+mesh.np1):
				if (int(mesh.dofs[i,0]) == 1):
					Fluc[j].U[0,indexU1]=eigenvectors[i,j]
					indexU1 += 1
				if (int(mesh.dofs[i,0]) == 0):
					Fluc[j].U[1,indexU2]=eigenvectors[i,j]
					indexU2 += 1
				if (int(mesh.dofs[i,0]) == 2):
					Fluc[j].P[indexP]=eigenvectors[i,j]
					indexP += 1
		return Fluc
