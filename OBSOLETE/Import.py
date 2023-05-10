#/* Copyright (C) 2019 FLOW group TU Berlin - All Rights Reserved
# * You may NOT use, distribute or modify this code without explicit
# * permission of the copyright owner, the FLOW group at TU Berlin!
# * However, permissions to use and modify the code are generally 
# * granted when asked for.
# * To ask for permission please contact t.kaiser@tu-berlin.de.
# */
# **********************************************************************
'''
# **********************************************************************
# * This file deals with reading flow data from various sources and
# * Interpolating it on the computational grid. 
# * This file was created by Thomas L. Kaiser. Significant contributions 
# * were done by...
# * - Jacopo Canton, currently (2019) a PostDoc at ETH Zurich, thanks for 
# *   allowing us to use your reader for Nekton files and saving us a lot
# *   of time
# **********************************************************************
'''
import numpy as np
from fenics import project, Function,SpatialCoordinate, HDF5File,plot,Expression,FiniteElement,MixedElement,FunctionSpace,triangle,split,FunctionAssigner
from fenics import interpolate as FenicsInterpolate
import matplotlib.pyplot as plt
from colorama import Fore, Style
from scipy import interpolate
from functions import *
import scipy.io as spio
import h5py
 
def ReadMeanFlow(param,FEMSpaces):
	''' Function reading the base/mean flow set in param.FlowInput.MeanFlowFilePath
	Input:
	\t param: parameter object (see parameters.py) defining the Case
	\t FEMSpaces: FEMSpaces object (see DefineFEMSpaces.py) defining all FEM spaces needed
	Output:
	\t MeanFlowDict: A dictionary containing fenics function fields of all needed quantities
	'''
	filePath=param.FlowInput.MeanFlowFilePath
	
	## If hdf5 file is used (so far implemented for fenics constructed hdf5 files)

	if filePath.split('.')[-1]=='fel':
		MeanFlowDict, notInFileList=importFelicsFile(param,FEMSpaces)
	if filePath.split('.')[-1]=='hdf5':
		MeanFlowDict, notInFileList=importHDF5File(param,FEMSpaces)
	elif filePath.split('.')[-1]=='mat':
		MeanFlowDict, notInFileList=importMatFile(param,FEMSpaces)
	elif filePath.split('.')[-1]=='h5':
		MeanFlowDict,notInFileList=importAVBPFile(param,FEMSpaces)
	elif filePath.split('.')[-1]=='cgns':
		MeanFlowDict, notInFileList=importCGNSFile(param,FEMSpaces)
	elif filePath.split('.')[-1]=='f00001':
		MeanFlowDict, notInFileList=importNektar5000File(param,FEMSpaces)
	elif filePath.split('.')[-1]=='vtk':
		MeanFlowDict, notInFileList=importVTKFile(param,FEMSpaces)

	#If file was interpolated on mesh (all cases except hdf5), then export the interpolated data
	if not filePath.split('.')[-1]=='hdf5':
		exportBaseFlowAsHDF5(param,MeanFlowDict,FEMSpaces)
	# If the cylindrical coordinate system is used, azimuthal and radial velocities are corrected at the axis:
	if param.CoordinateSystem in ['Cylindrical']:
		MeanFlowDict=CorrectMeanAtBoundaries(MeanFlowDict)
		
	if (param.ReynoldsNumberSource=='Constant'):
		MeanFlowDict['nuturb']=Function(FEMSpaces.P2)
	elif (param.ReynoldsNumberSource=='File'):
		## In this case the work is already done...
		print('Importing Reynolds number from file')
	elif (param.ReynoldsNumberSource=='Sutherland'):
		MeanFlowDict['nuMol'] = Function(FEMSpaces.P1)
		MeanFlowDict['nuMol'].vector()[:] = SutherlandLaw(param,MeanFlowDict)
		
	elif (param.ReynoldsNumberSource=='PowerLaw_AVBP'):
		M=28.949
		p=101300
		R=8314.4598/M
		T=p/R/MeanFlowDict['rho'].vector()[:]
		mu0=1.8405e-5
		Tref=300
		exponent=0.6759
		mu=mu0*(T/Tref)**(exponent) #see avbp handbook
		#nu=mu/MeanFlowDict['rho'].vector()[:]
		MeanFlowDict['Re'].vector()[:]=1/mu[:]
	elif (param.ReynoldsNumberSource=='EddyViscosity'):
		print("Calculating Reynolds number based on eddy viscosity...")
		MeanFlowDict['Re'].vector()[:]=1/MeanFlowDict['nut'].vector()[:]
		#print(MeanFlowDict['Re'].vector()[:])
	elif (param.ReynoldsNumberSource=='TKE-based'):
		print("Calculating Reynolds number based on TKE...")
		MeanFlowDict = nutTKEnut(MeanFlowDict, param, FEMSpaces) 
		MeanFlowDict['Re']=Function(FEMSpaces.P2)
		# this might be case specififc, however, if Temperature is atmospheric and constant Sutherland's law ...
		# ... does not alter the viscosity since T=293 (see implementation in functions.py)
		MeanFlowDict['nuMol_P1']=Function(FEMSpaces.P1)
		MeanFlowDict['nuMol']=Function(FEMSpaces.P2)
		MeanFlowDict['nuMol_P1'].vector()[:]  = SutherlandLaw(param,MeanFlowDict)
		MeanFlowDict['nuMol'] =  FenicsInterpolate(MeanFlowDict['nuMol_P1'], FEMSpaces.P2)
		del MeanFlowDict['nuMol_P1']
		# since nut and nuMol are kinematic viscosity, density is not accounted for yet
		# if density field exists, use that one, if not set to fluid specific density (e.g. air->1.12)
		if 'rho' in param.MeanList:
			MeanFlowDict['rhoP2'] = Function(FEMSpaces.P2)
			MeanFlowDict['rhoP2'] =  FenicsInterpolate(MeanFlowDict['rho'], FEMSpaces.P2)
			rho = MeanFlowDict['rhoP2'].vector()[:]
			del MeanFlowDict['rhoP2']
		else:
			rho = 1.12
		
		# this is only here to avoid nans and infs
		tmp = MeanFlowDict['Re'].vector()[:]
		tmp[np.isnan(tmp)] = param.Re
		tmp[np.isinf(tmp)] = param.Re
		
	elif 'Boussinesq' in param.ReynoldsNumberSource:
		if 'xr' in param.ReynoldsNumberSource:
			print("Calculating Reynolds number based on reduced Boussinesq (x-r-components) Ansatz...")
			MeanFlowDict = nutFromReducedBoussinesq_xr(MeanFlowDict, param, FEMSpaces)
		else:
			MeanFlowDict = nutFromBoussinesq(MeanFlowDict, param, FEMSpaces)  # computes turbulent viscosity field from Boussniesq-Ansatz

	for name in notInFileList:
		print(Fore.YELLOW+"WARNING: Field " + name + " not found in file "+ filePath+"! Assuming Field is zero... "+Style.RESET_ALL)
	if param.FlowMode=='Reacting':	
		## Check if mass fractions are larger or smaller than zero
		for specie in getSpeciesListMean(param):
			if np.max(MeanFlowDict[specie].vector()[:])>1 or np.min(MeanFlowDict[specie].vector()[:])<0:
				printWarning('The maximum of mass fraction of species '+specie+' is outside of the boundaries 0 to 1. The respective values are cut. This can result in significant errors in the sum of mass fractions!')
				MeanFlowDict[specie].vector()[:]=np.clip(MeanFlowDict[specie].vector()[:],0.0000001,1)
		## Check sum of species mass ratio
		sumOfMass=MeanFlowDict[getSpeciesListMean(param)[0]].vector()[:]*0
		for name in getSpeciesListMean(param):
			sumOfMass=sumOfMass+MeanFlowDict[name].vector()[:]
		print('Checking sum of mass fractions...')
		tol=1e-10
		error=np.max(np.abs(sumOfMass-1))
		if error>1e-10:
			printWarning('Largest error in sum of mass fractions is '+str(error))
		else:
			print('Sum of mass fractions is equal to one! Tolerance is '+str(tol))		
		#Add reaction
		ReactionList=getReactionList(param)
		print(ReactionList)

		#Check reaction space type P1 or P2
		ReactionSpaceTest=Function(FEMSpaces.P2)
		if np.size(ReactionSpaceTest.vector()[:])==np.size(MeanFlowDict['T'].vector()[:]):
			ReactionSpaceType=FEMSpaces.P2
		else:
			ReactionSpaceType=FEMSpaces.P1
		print(ReactionSpaceType)
		del ReactionSpaceTest
		
		for Reaction in ReactionList:
			#The temperature guess must be fixed...
			#M=28.949
			#p=101300
			#R=8314.4598/M
			#T=p/R/MeanFlowDict['rho'].vector()[:]
			ReactionName=getReactionName(Reaction)

			#ReactionName=''
			#for specie in Reaction['educts']:
			#	ReactionName=ReactionName+specie
			#print(MeanFlowDict['T'].vector())
			MeanFlowDict[ReactionName]=Function(ReactionSpaceType)
			MolecularMass=getMolecularMass()
			MeanFlowDict[ReactionName].vector()[:]=Reaction['A']\
				*MeanFlowDict['T'].vector()[:]**(Reaction['beta'])\
				*np.exp(-Reaction['E_a']/MeanFlowDict['T'].vector()[:]/Reaction['R'])\
				*MeanFlowDict['rho'].vector()[:]**(Reaction['a']+Reaction['b'])\
				*MeanFlowDict[Reaction['educts'][0]].vector()[:]**Reaction['a']\
				*MeanFlowDict[Reaction['educts'][1]].vector()[:]**Reaction['b']\
				/MolecularMass[Reaction['educts'][0]]**Reaction['a']\
				/MolecularMass[Reaction['educts'][1]]**Reaction['b']

		ReverseReactionList=getReverseReactionList(param)
		print(ReverseReactionList)
		for Reaction in ReverseReactionList:
			ReverseReactionName=getReactionName(Reaction)
			ForwardReactionName=''
			for specie in Reaction['products']:
				ForwardReactionName=ForwardReactionName+specie
			
			#Forwaard Reaction rate kf: Q=A*T**beta*exp(-Ea/R/T) 
			ForwardReactionProgress=Function(ReactionSpaceType)
			ForwardReactionProgress.vector()[:]=MeanFlowDict[ForwardReactionName].vector()[:]\
							*MolecularMass[Reaction['products'][0]]**Reaction['fa']\
							*MolecularMass[Reaction['products'][1]]**Reaction['fb']\
							/MeanFlowDict[Reaction['products'][0]].vector()[:]**Reaction['fa']\
							/MeanFlowDict[Reaction['products'][1]].vector()[:]**Reaction['fb']\
							/MeanFlowDict['rho'].vector()[:]**(Reaction['fa']+Reaction['fb'])
			
			
			#Equilibrium constant Keq
			EquiConstant=Function(ReactionSpaceType)
			p0=101300
			R=8.3144598
			
			nu_j=0
			for i in Reaction['educt_mole_factor']:
				nu_j=nu_j+i
			for i in Reaction['product_mole_factor']:
				nu_j=nu_j-i
			#print('nu_j:'+str(nu_j))
			EquiConstant.vector()[:]=(p0/R/MeanFlowDict['T'].vector()[:])**nu_j
			
			janaf_table=getJANAFTable()
			ExpInEquiConstant=Function(ReactionSpaceType)
			ExpInEquiConstant.vector()[:]=0
			for educt_name in Reaction['educts']:
				janaf_coef=janaf_table[educt_name]
				enthalpy0=Function(ReactionSpaceType)
				entropy0=Function(ReactionSpaceType)
				
				for i in range(0,len(MeanFlowDict['T'].vector())):
					tem=MeanFlowDict['T'].vector()[i]
					if tem>1000:
						enthalpy0.vector()[i]=R*(janaf_coef[0]*tem+janaf_coef[1]*tem**2/2+janaf_coef[2]*tem**3/3\
						+janaf_coef[3]*tem**4/4+janaf_coef[4]*tem**5/5+janaf_coef[5])
						entropy0.vector()[i]=R*(janaf_coef[0]*np.log(tem)+janaf_coef[1]*tem+janaf_coef[2]*tem**2/2\
						+janaf_coef[3]*tem**3/3+janaf_coef[4]*tem**4/4+janaf_coef[6])
					else:
						enthalpy0.vector()[i]=R*(janaf_coef[7]*tem+janaf_coef[8]*tem**2/2+janaf_coef[9]*tem**3/3\
						+janaf_coef[10]*tem**4/4+janaf_coef[11]*tem**5/5+janaf_coef[12])
						entropy0.vector()[i]=R*(janaf_coef[7]*np.log(tem)+janaf_coef[8]*tem+janaf_coef[9]*tem**2/2\
						+janaf_coef[10]*tem**3/3+janaf_coef[11]*tem**4/4+janaf_coef[13])
				ExpInEquiConstant.vector()[:]=ExpInEquiConstant.vector()[:]+(entropy0.vector()[:]/R\
				-enthalpy0.vector()[:]/R/MeanFlowDict['T'].vector())*Reaction['educt_mole_factor'][Reaction['educts'].index(educt_name)]
			
			for product_name in Reaction['products']:
				janaf_coef=janaf_table[product_name]
				enthalpy0=Function(ReactionSpaceType)
				entropy0=Function(ReactionSpaceType)
				
				for i in range(0,len(MeanFlowDict['T'].vector())):
					tem=MeanFlowDict['T'].vector()[i]
					if tem>1000:
						enthalpy0.vector()[i]=R*(janaf_coef[0]*tem+janaf_coef[1]*tem**2/2+janaf_coef[2]*tem**3/3\
						+janaf_coef[3]*tem**4/4+janaf_coef[4]*tem**5/5+janaf_coef[5])
						entropy0.vector()[i]=R*(janaf_coef[0]*np.log(tem)+janaf_coef[1]*tem+janaf_coef[2]*tem**2/2\
						+janaf_coef[3]*tem**3/3+janaf_coef[4]*tem**4/4+janaf_coef[6])
					else:
						enthalpy0.vector()[i]=R*(janaf_coef[7]*tem+janaf_coef[8]*tem**2/2+janaf_coef[9]*tem**3/3\
						+janaf_coef[10]*tem**4/4+janaf_coef[11]*tem**5/5+janaf_coef[12])
						entropy0.vector()[i]=R*(janaf_coef[7]*np.log(tem)+janaf_coef[8]*tem+janaf_coef[9]*tem**2/2\
						+janaf_coef[10]*tem**3/3+janaf_coef[11]*tem**4/4+janaf_coef[13])
				ExpInEquiConstant.vector()[:]=ExpInEquiConstant.vector()[:]-(entropy0.vector()[:]/R\
				-enthalpy0.vector()[:]/R/MeanFlowDict['T'].vector())*Reaction['product_mole_factor'][Reaction['products'].index(product_name)]
			EquiConstant.vector()[:]=EquiConstant.vector()[:]*np.exp(ExpInEquiConstant.vector()[:])
			MeanFlowDict['Reverse'+str(ReverseReactionName)]=Function(FEMSpaces.P1)
			MeanFlowDict['Reverse'+str(ReverseReactionName)].vector()[:]=ForwardReactionProgress.vector()[:]/EquiConstant.vector()[:]\
			*MeanFlowDict[educt_name].vector()[:]**Reaction['rc']/MolecularMass[educt_name]**Reaction['rc']\
			*MeanFlowDict['rho'].vector()[:]**(Reaction['rc'])
	return MeanFlowDict #, notInFileList#, InputFlowDict

def importHDF5File(param,FEMSpaces):
	## Define Dictionaries
	MeanFlowDict={}
	notInFileList=[]
	## Get mesh data
	mesh=FEMSpaces.P2.mesh()
	# If hdf5 file is used (so far implemented for fenics constructed hdf5 files)
	hdf5file = HDF5File(mesh.mpi_comm(),param.FlowInput.MeanFlowFilePath, 'r')        
	for name in paramFunction(FEMSpaces.P2).MeanList:
		if name[:-1]=='u_forcing_':
			# If the field is the velocity it must be read in component by component... These are collected in u_compDict
			u_compDict={}
			for component in param.VelocityComponents:
				u_compDict['u'+component+name[1:]]=Function(FEMSpaces.P2)
				try:
					hdf5file.read(u_compDict['u'+component+name[1:]], 'u'+component+name[1:])
				except Exception as errors:
					#print(errors)
					u_compDict['u'+component+name[1:]].vector()[:]=0.0
					printDebug(param.debug,'u'+component+name[1:]+' set to 0.0')
					notInFileList.append('u'+component+name[1:])
			#MeanFlowDict[name]=Function(FEMSpaces.FunctionSpaceVectorVelocity)
			if param.nVelocityComponents == 2:
				MeanFlowDict[name] = project(Expression(("u1", "u2"),\
	 				u1=u_compDict['u'+param.VelocityComponents[0]+name[1:]],\
					u2=u_compDict['u'+param.VelocityComponents[1]+name[1:]],\
					degree=2), FEMSpaces.FunctionSpaceVectorVelocity)	
# REMARK: there might occur problems concerning memory consumption by the fenics-project-function:
# if you run into such problems try this: 
#				MeanFlowDict[name] = project(Expression(("u1", "u2"),\
#	 				u1=u_compDict['u'+param.VelocityComponents[0]+name[1:]],\
#					u2=u_compDict['u'+param.VelocityComponents[1]+name[1:]],\
#					degree=2), FEMSpaces.FunctionSpaceVectorVelocity, solver_type="cg", preconditioner_type="amg")	
			elif param.nVelocityComponents == 3:
				MeanFlowDict[name] = project(Expression(("u1", "u2", "u3"),\
	 				u1=u_compDict['u'+param.VelocityComponents[0]+name[1:]],\
					u2=u_compDict['u'+param.VelocityComponents[1]+name[1:]],\
					u3=u_compDict['u'+param.VelocityComponents[2]+name[1:]],\
					degree=2), FEMSpaces.FunctionSpaceVectorVelocity)	
# REMARK: there might occur problems concerning memory consumption by the fenics-project-function:
# if you run into such problems try this: 
#				MeanFlowDict[name] = project(Expression(("u1", "u2"),\
#	 				u1=u_compDict['u'+param.VelocityComponents[0]+name[1:]],\
#					u2=u_compDict['u'+param.VelocityComponents[1]+name[1:]],\
#					degree=2), FEMSpaces.FunctionSpaceVectorVelocity, solver_type="cg", preconditioner_type="amg")	
		elif name=='u':
			# If the field is the velocity it must be read in component by component... These are collected in u_compDict
			u_compDict={}
			for component in param.VelocityComponents:
				u_compDict[name+component]=Function(FEMSpaces.P2)
				try:
				    hdf5file.read(u_compDict[name+component], name+component)
				except:
				    notInFileList.append(name+component)
			#MeanFlowDict[name]=Function(FEMSpaces.FunctionSpaceVectorVelocity)
			if param.nVelocityComponents == 2:
				MeanFlowDict[name] = project(Expression(("u1", "u2"),\
	 				u1=u_compDict[name+param.VelocityComponents[0]],\
					u2=u_compDict[name+param.VelocityComponents[1]],\
					degree=2), FEMSpaces.FunctionSpaceVectorVelocity)	
			elif param.nVelocityComponents == 3:
				print(np.shape(u_compDict[name+param.VelocityComponents[0]].vector()[:]))
				print(np.shape(u_compDict[name+param.VelocityComponents[1]].vector()[:]))
				print(np.shape(u_compDict[name+param.VelocityComponents[2]].vector()[:]))

				MeanFlowDict[name] = project(Expression(("u1", "u2", "u3"),\
	 				u1=u_compDict[name+param.VelocityComponents[0]],\
					u2=u_compDict[name+param.VelocityComponents[1]],\
					u3=u_compDict[name+param.VelocityComponents[2]],\
					degree=2), FEMSpaces.FunctionSpaceVectorVelocity)	
		elif name[0:3]=='rst':
			MeanFlowDict[name]=Function(FEMSpaces.P2)
			try:
				hdf5file.read(MeanFlowDict[name], name)
			except:
				notInFileList.append(name)
		else:
			MeanFlowDict[name]=Function(FEMSpaces.P2)
			try:
				hdf5file.read(MeanFlowDict[name], name)
			except:
				notInFileList.append(name)
	printDebug(param.debug,"Read in mean flow fields are " + str(MeanFlowDict.keys()))
	printDebug(param.debug,"Not in mean flow fields are " + str(notInFileList))
	return MeanFlowDict, notInFileList
 
def importMatFile(param,FEMSpaces):
	## Define Dictionaries
	MeanFlowDict={}
	notInFileList=[]
	## Get mesh data
	mesh=FEMSpaces.P2.mesh()
	### Get Mean flow names
		
	filePath=param.FlowInput.MeanFlowFilePath
	mat = spio.loadmat(filePath[:-4])   
	printDebug(param.debug,'The fields in the Matlab file are ' +str(mat.keys()))
	if param.CoordinateSystem=='Cartesian':
		x1_mean=mat['x']
		x2_mean=mat['y']
	elif param.CoordinateSystem=='Cylindrical':
		x1_mean=mat['x']
		x2_mean=mat['r']
	#hdf5file = HDF5File(mesh.mpi_comm(),filePath[0:-4]+"__"+param.MeshPath[0:-3].split('/')[-1]+"hdf5", 'w')
	for name in param.MeanList:
		# In case the data to read is the velocity, we must iterate through all velocity components and
		# use P2 elements
		if name[0]=='u':
			# init MeanFlowDict (velocity): has be built from a VectorFunctionSpace
			MeanFlowDict[name] = Function(FEMSpaces.FunctionSpaceVectorVelocity)
			# 'component' is the actual component key; compCNT is the corresponding index to adress the correct sub-space
			for component, compCNT in zip(param.VelocityComponents, range(0,len(param.VelocityComponents))):
				# get the correct component key to read in from
				nameComponent=name[:1]+component+name[1:]
				
				# Vectorspace is subdivided into sub spaces
					# -> get the dofIDX is a index list including the indices pointing to the dof_coordinates of the currently considered sub space
				dofIDX = FEMSpaces.FunctionSpaceVectorVelocity.sub(compCNT).dofmap().dofs()
					# get dof_coordiantes for the currently considered sub space (dofIDX-depending)
				dof_coordinates = FEMSpaces.FunctionSpaceVectorVelocity.tabulate_dof_coordinates()[dofIDX]                   
				try:
					# MeanFlowDict.vector includes all components in one vector, being sorted with respect to dofIDX
					print(1)
					print(np.shape(np.concatenate([x1_mean,x2_mean],axis=1)))
					print(2)
					print(np.shape(mat[nameComponent]))
					print(3)
					
					MeanFlowDict[name].vector()[dofIDX]=np.squeeze(interpolate.griddata(np.concatenate([x1_mean,x2_mean],axis=1), mat[nameComponent],dof_coordinates, method='linear'))	
				except:
					notInFileList.append(nameComponent)
				if (max(np.isnan(MeanFlowDict[name].sub(compCNT).vector()[:]))):
					print(Fore.RED+"Error: After interpolation NaN values occur in Field " +nameComponent+"! Probably the Felics mesh extends beyond the base flow mesh! Ending program..."+Style.RESET_ALL)
					exit()
		else:
			if name == 'p':
				MeanFlowDict[name]=Function(FEMSpaces.P1)
				dof_coordinates = FEMSpaces.P1.tabulate_dof_coordinates()                      
			else:
				MeanFlowDict[name]=Function(FEMSpaces.P2)
				dof_coordinates = FEMSpaces.P2.tabulate_dof_coordinates()                      
			try:
				MeanFlowDict[name].vector()[:]=np.squeeze(interpolate.griddata(np.concatenate([x1_mean,x2_mean],axis=1), mat[name],dof_coordinates, method='linear'))	
			except:
				notInFileList.append(name)
			if (max(np.isnan(MeanFlowDict[name].vector()[:]))):
				print(Fore.RED+"Error: After interpolation NaN values occur in Field " +name+"! Probably the Felics mesh extends beyond the base flow mesh! Ending program..."+Style.RESET_ALL)
				exit()
		#hdf5file.write(MeanFlowDict[name], name)
#	if 'nut' in nameListMean:
#		MeanFlowDict['nut'].vector()[MeanFlowDict['nut'].vector()[:]==0] = np.nan

	return MeanFlowDict, notInFileList
def importNektar5000File(param,FEMSpaces):
    import neksuite as ns
    import os.path
    import pickle
    ## Define Dictionaries
    MeanFlowDict={}
        #The RawFlowDict contains the raw data in the following format: x_raw,y_raw,z_raw,ux_raw,uy_raw,uz_raw...
        #These must later be transformed to the coordinates used in the Felics solver
    RawFlowDict={}  
    notInFileList=[]
    ## Get mesh data
    mesh=FEMSpaces.P2.mesh()
    ## Get Mean flow names
    #nameListMean=getMeanFlowFieldNames(param)
    nameListMean=param.Case.getMeanFlowFieldNames()
    ## Get AVBP mesh file path
    filePath=param.FlowInput.MeanFlowFilePath
    # Check if intermediate file was written, which can be read much faster than the Nekton files...
    if os.path.isfile('TempExport.pkl'): 
        with open('TempExport.pkl','rb') as f:  # Python 3: open(..., 'rb')
            RawFlowDict = pickle.load(f)[0]
        print(type(RawFlowDict))
        print(RawFlowDict)
    else:
        print('Reading Nektar file. This may take a while...')
        NektarFile = ns.readnek(filePath)
        # polynomial order for each direction (usually the same)
        nx = NektarFile.lr1[0]
        ny = NektarFile.lr1[1]
        nz = NektarFile.lr1[2]
 
        dim = NektarFile.ndim
         
        nel= NektarFile.nel
        #Loop over spectral elements
        n_grid=nel*nx*ny*nz
        index=0
        print (nel)
        RawFlowDict={}
        RawFlowDict['x_raw']=np.zeros(n_grid)
        RawFlowDict['z_raw']=np.zeros(n_grid)
        RawFlowDict['y_raw']=np.zeros(n_grid)
        RawFlowDict['ux_raw']=np.zeros(n_grid)
        RawFlowDict['uz_raw']=np.zeros(n_grid)
        RawFlowDict['uy_raw']=np.zeros(n_grid)
        print ('Extracting flow from Nektar data structure...')
        for e in range(nel):
            # loop for comp. nodes
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        RawFlowDict['x_raw'][index]=NektarFile.elem[e].pos[0][i][j][k]
                        RawFlowDict['y_raw'][index]=NektarFile.elem[e].pos[1][i][j][k]
                        RawFlowDict['z_raw'][index]=NektarFile.elem[e].pos[2][i][j][k]
 
                        RawFlowDict['ux_raw'][index]=NektarFile.elem[e].vel[0][i][j][k]
                        RawFlowDict['uy_raw'][index]=NektarFile.elem[e].vel[1][i][j][k]
                        RawFlowDict['uz_raw'][index]=NektarFile.elem[e].vel[2][i][j][k]
                        #increment index
                        index=index+1
        # Write an export, which can be read much quicker than the Nekton files
        print ("Writing intermediate file 'TempExport.pkl 'for faster input next time...")
        with open('TempExport.pkl', 'wb') as f: 
            pickle.dump([RawFlowDict], f)
    # Perform Coordinate Transform to go from carthesian to cylindrical coordinates and rotate into the right direction
    InputFlowDict=CoordinateTransformation('Cart2Cyl',RawFlowDict, param)
     
 
    x1_mean=np.array(InputFlowDict['x'])
    x2_mean=np.array(InputFlowDict['y'])
    x3_mean=np.array(InputFlowDict['z'])
    x1_mean=x1_mean.tolist()
    x2_mean=x2_mean.tolist()
    x3_mean=x3_mean.tolist()
    #FieldsToDelete=[]
    keys=list(InputFlowDict.keys())
    for name in keys:
        if name[0:3]=='rho' and len(name)>3:
            InputFlowDict[name]=InputFlowDict[name]/InputFlowDict['rho']
            InputFlowDict[name[3:]] = InputFlowDict[name]
            del InputFlowDict[name]
        if name in getSpeciesListMean(param):
            InputFlowDict[name]=InputFlowDict[name]/InputFlowDict['rho']
    for name in nameListMean:
        if name[0]=='u':
            MeanFlowDict[name]=Function(FEMSpaces.P2)
            dof_coordinates = FEMSpaces.P2.tabulate_dof_coordinates()                      
        else:
            print(name)
            MeanFlowDict[name]=Function(FEMSpaces.P1)
            dof_coordinates = FEMSpaces.P1.tabulate_dof_coordinates()                      
        try:
            # Expanding into azimuthal direction
            dof_expanded=ExpandForAzimuthalAverage(dof_coordinates,100)
            print('Performing interpolation for ' + name)
            temp_vec=np.squeeze(interpolate.griddata(list(map(list, zip(*[x1_mean,x2_mean,x3_mean]))), InputFlowDict[name],dof_expanded, method='nearest'))
            # Collapsing back to initial mesh
            MeanFlowDict[name].vector()[:]=ContractFromAximuthalAverage(dof_coordinates,100,temp_vec)
        except KeyError as e:
            notInFileList.append(name)
            print( 'I got a KeyError - reason ' + str(e))
     
    return MeanFlowDict, notInFileList
 
 
def importVTKFile(param, FEMSpaces):
	import meshio         # meshio is able to import vtk files
	from scipy.interpolate import LinearNDInterpolator
	import time
	from sklearn.neighbors import KNeighborsRegressor, RadiusNeighborsRegressor, NearestNeighbors
	## Define Dictionaries
	MeanFlowDict={}
        #The RawFlowDict contains the raw data in the following format: x_raw,y_raw,z_raw,ux_raw,uy_raw,uz_raw...
        #These must later be transformed to the coordinates used in the Felics solver
	RawFlowDict={}  
	notInFileList=[]
    ## Get mesh data
	mesh=FEMSpaces.P2.mesh()
	
    ## Get Mean flow names
	#nameListMean=getMeanFlowFieldNames(param)
	nameListMean=param.Case.getMeanFlowFieldNames()

	print('Reading VTK file. This may take a while...')
	vtkDat = meshio.read(param.FlowInput.MeanFlowFilePath)
	nPoints = vtkDat.points.shape[0]
	
	method = 'nearest griddata' # options are: 'nearest griddata', 'linear griddata' , 'parallel'
	nTheta = 100 # the amount of planes in theta direction which results from azimuthal expansion of Felics grid
	
	RawFlowDict['x'] = np.array([vtkDat.points[k][0] for k in range(0,nPoints)])
	RawFlowDict['y'] = np.array([vtkDat.points[k][1] for k in range(0,nPoints)])
	RawFlowDict['z'] = np.array([vtkDat.points[k][2] for k in range(0,nPoints)])
	RawFlowDict['ux'] = np.array([vtkDat.point_data['UMean'][k][0] for k in range(0,nPoints)])
	RawFlowDict['uy'] = np.array([vtkDat.point_data['UMean'][k][1] for k in range(0,nPoints)])
	RawFlowDict['uz'] = np.array([vtkDat.point_data['UMean'][k][2] for k in range(0,nPoints)])
	

	if param.ReynoldsNumberSource in ['TKE-based','Boussinesq','Boussinesq(xr)']:
		# MARIOTHESIS-RELEVANT ONLY....numbering in KIT vtk-Files differed between cold and hot (Mario 21.06.2020)
		# if this if-statement preserves for longer, than boil it down to the second if-clause
		# it was necessary since the ordering in the vtk -file with name 'berlinResolvantePhi' was changed
		# it was copied into this branch to test "validate" the branch against MarioThesis-branch
		if 'berlinResolvantePhi' in param.FlowInput.MeanFlowFilePath:
			RawFlowDict['rstxx_raw'] = np.array([vtkDat.point_data['UPrime2Mean'][k][0] for k in range(0,nPoints)]) #np.array(h5file['Average']['u2'])-np.array(h5file['Average']['u'])**2
			RawFlowDict['rstyy_raw'] = np.array([vtkDat.point_data['UPrime2Mean'][k][3] for k in range(0,nPoints)])
			RawFlowDict['rstzz_raw'] = np.array([vtkDat.point_data['UPrime2Mean'][k][5] for k in range(0,nPoints)])
			RawFlowDict['rstxy_raw'] = np.array([vtkDat.point_data['UPrime2Mean'][k][1] for k in range(0,nPoints)])
			RawFlowDict['rstxz_raw'] = np.array([vtkDat.point_data['UPrime2Mean'][k][2] for k in range(0,nPoints)])
			RawFlowDict['rstyz_raw'] = np.array([vtkDat.point_data['UPrime2Mean'][k][4] for k in range(0,nPoints)])
		else:
			RawFlowDict['rstxx_raw'] = np.array([vtkDat.point_data['UPrime2Mean'][k][0] for k in range(0,nPoints)]) #np.array(h5file['Average']['u2'])-np.array(h5file['Average']['u'])**2
			RawFlowDict['rstyy_raw'] = np.array([vtkDat.point_data['UPrime2Mean'][k][1] for k in range(0,nPoints)])
			RawFlowDict['rstzz_raw'] = np.array([vtkDat.point_data['UPrime2Mean'][k][2] for k in range(0,nPoints)])
			RawFlowDict['rstxy_raw'] = np.array([vtkDat.point_data['UPrime2Mean'][k][3] for k in range(0,nPoints)])
			RawFlowDict['rstxz_raw'] = np.array([vtkDat.point_data['UPrime2Mean'][k][5] for k in range(0,nPoints)])
			RawFlowDict['rstyz_raw'] = np.array([vtkDat.point_data['UPrime2Mean'][k][4] for k in range(0,nPoints)])
	
	if 'rho' in nameListMean:
#		RawFlowDict['rho'] = np.array([vtkDat.point_data['rhoMean'] for k in range(0,nPoints)])
		RawFlowDict['rho'] = np.array([vtkDat.point_data['rhoMean'][k] for k in range(0,nPoints)])
		RawFlowDict['T'] = np.array([vtkDat.point_data['TMean'][k] for k in range(0,nPoints)])
	
	InputFlowDict=CoordinateTransformation('Cart2Cyl',RawFlowDict, param)
	
	# get the Felics domain bounding box plus tolerance (hmax is the maximum Felics cell size) : tolerance makes the bounding box slighlty larger
	# it has to be larger because the imported grid points (here from vtk) has to be cover the Felics domain to ensure a propper interpolation
	Xcoords = [k[0] for k in mesh.coordinates()]
	Ycoords = [k[1] for k in mesh.coordinates()]
	boundingBox = [min(Xcoords)-mesh.hmax(), max(Xcoords)+mesh.hmax(),\
				min(Ycoords)-mesh.hmax(), max(Ycoords)+mesh.hmax(),\
				-mesh.hmax(), mesh.hmax()]
	
	# coordinates of axial and radial component. cannot bound in azimuthal direction since interpolation is in 3 dimensions
	X = InputFlowDict['x']
	R = InputFlowDict['r']
	print('size of InputFlowDict before applying box restriction:' + str(InputFlowDict['r'].shape[0]))
	# find indices fitting into bounding bounding box
	Xidx = np.argwhere( (X>boundingBox[0]) & (X<boundingBox[1])).flatten()
	Ridx = np.argwhere( ((R<boundingBox[3]))).flatten()
	
	# intersection of indices yields indices in bounding box with respect to X AND R
	IDX = np.intersect1d(Xidx, Ridx)
	
	p1CNT=0
	p2CNT=0
	for k in nameListMean:
		if not k=='Re':
			# cut InputFlowDict by indices IDX
			print(k)
			InputFlowDict[k] = InputFlowDict[k][IDX]
			# some field manipulation @ Thomas ?
			if k[0:3]=='rho' and len(k)>3:
				InputFlowDict[k]=InputFlowDict[k]/InputFlowDict['rho']
				InputFlowDict[k[3:]] = InputFlowDict[k]
				del InputFlowDict[k]
			if k in getSpeciesListMean(param):
				InputFlowDict[k]=InputFlowDict[k]/InputFlowDict['rho']
			# count P1 and P2 fields to initialize stacked value array for interpolation
			if k[0] == 'u' :
				p2CNT+=len(param.VelocityComponents)
			elif k[0:3] == 'rst':
				p2CNT+=1
			else:
				p1CNT+=1
	
	InputFlowDict['x'] = InputFlowDict['x'][IDX]
	InputFlowDict['r'] = InputFlowDict['r'][IDX]
	InputFlowDict['theta'] = InputFlowDict['theta'][IDX]
	
	InputFlowDict['y'] = InputFlowDict['y'][IDX]
	InputFlowDict['z'] = InputFlowDict['z'][IDX]
	# concatenate coordinates to suit a [Npoints,3] = [rows, cols] shape
	points = np.vstack((InputFlowDict['x'],InputFlowDict['y'],InputFlowDict['z'])).T
	# value arrays for interpolation are split up here, since we have to differ between interpolation onto P1 and onto P2 grid
	if p1CNT>0: valsP1 = np.zeros((points.shape[0], p1CNT ))
	if p2CNT>0: valsP2 = np.zeros((points.shape[0], p2CNT ))
	# compute dof_coordinates and their azimuthally expanded pendant for each, P1 and P2 fields
	dof_coordinatesP1 = FEMSpaces.P1.tabulate_dof_coordinates()
	dof_expandedP1=ExpandForAzimuthalAverage(dof_coordinatesP1,nTheta)
	dof_coordinatesP2 = FEMSpaces.P2.tabulate_dof_coordinates()
	dof_expandedP2=ExpandForAzimuthalAverage(dof_coordinatesP2,nTheta) 
	 
	print('size of InputFlowDict after applying box restriction:' + str(InputFlowDict['r'].shape[0]))
	print('size azimuthally expanded dofs:' + str(dof_expandedP2.shape[0]))
	
	m = 0 # counter for P2 arrays and lsit
	mm = 0 # counter for P1 arrays and lsit
	# create two name lists representing existing P1 and p2 fields
	namesP2 = []
	namesP1 = []
	# loop over all mean field names -> initialize FEM-Space, fill value array for interpolation, 
	print(nameListMean)
	
	additionalVelCompsTMP = ['u'+k for k in param.AdditionalVelocityComponents]
	for name in param.MeanList:
		if name == 'u':
			MeanFlowDict[name[0]] = Function(FEMSpaces.FunctionSpaceVectorVelocity)
			# 'component' is the actual ocomponent key; compCNT is the corresponding index to adress the correct sub-space
			for component, compCNT in zip(param.VelocityComponents, range(0,len(param.VelocityComponents))):
				valsP2[:,m] = InputFlowDict[name+component]
				namesP2.append(name)
				m+=1
		elif name[0:3] == 'rst' or name in additionalVelCompsTMP:
			MeanFlowDict[name]=Function(FEMSpaces.P2)
			valsP2[:,m] = InputFlowDict[name]
			namesP2.append(name)
			m+=1
		else:
			print(name)
			MeanFlowDict[name]=Function(FEMSpaces.P1) 
			if not name == 'Re':
				valsP1[:,mm] = InputFlowDict[name]
				mm+=1
			namesP1.append(name)
	
	# three interpolation methods available now: KNeighborsRegressor and griddatawith nearest flag and griddata with linear flag
	print('Starting with the interpolation using Method: '+str(method))
	p1time = time.time()	 # for time measure
	
	
	if p1CNT>0: # if p1 fields are not existing skip interpolation
		if method == 'nearest griddata':
			temp_vecP1 = interpolate.griddata(points, valsP1,dof_expandedP1, method='nearest')
		elif method == 'linear griddata':
			temp_vecP1 = interpolate.griddata(points, valsP1,dof_expandedP1, method='linear', fill_value=0)	
			nansIDX = np.isnan(temp_vecP1)
			print(nansIDX)
			print(nansIDX.shape)
			temp_vecP1[nansIDX]=0
			nansIDX=np.isnan(temp_vecP1)
			print(nansIDX)
		elif method == 'parallel':
			# strategy: 
			# build a convex 3d-cube around the input data and divide it into sub cubes
			# divide the destination grid into the same sub-cubes 
			# there is a tolerance band between the sub cubes such that some points are interpolate twice, however, none is forgitten
			# the tolerance is set to mesh.hmax() = maximum grid spacing of Felics grid, CURRENTLY
			# divisionFactor tells how many cuts are performed for one dimesion (3D and divionsFactor = 2 results in 8 cubes) 
			dof_expanded = dof_expandedP1
			# origin data
			vals = valsP1 
			#imagine a cubic box around the data points and divide it into sub-cubes (subCubesInY = amount if cubes in y direction) 
			# HARD CODED.... division like this is very naive for high refinement gradients and non-cubic domains
			# HARD CODED.... the amount of points is NOT equally distributed onto subcubes....might result in inefficient paralliza.
			subCubesInX = 2
			subCubesInY = 3
			subCubesInZ = 2
			
			# get the bounding coordinates of the original cube in which the origin data fits
			dofMin_x = min(dof_expanded[:,0])
			dofMax_x = max(dof_expanded[:,0])
			dofMin_y = min(dof_expanded[:,1])
			dofMax_y = max(dof_expanded[:,1])
			dofMin_z = min(dof_expanded[:,2])
			dofMax_z = max(dof_expanded[:,2])
			# create arrays with coordinates of cuts through the origin cube
			xTMP = np.linspace(dofMin_x,dofMax_x,subCubesInX+1)
			yTMP = np.linspace(dofMin_y,dofMax_y,subCubesInY+1)
			zTMP = np.linspace(dofMin_z,dofMax_z,subCubesInZ+1)
			
			# bottomBounds* and topBounds* -arrays include the boundaries of each subCube
			# image you move along the x-axis: the boundaries of the n-th cube are in the n-th entry of bottomBoundsX and topBoundsX
			bottomBoundsX = xTMP[0:-1]
			topBoundsX = xTMP[1:]
			
			bottomBoundsY = yTMP[0:-1]
			topBoundsY = yTMP[1:]
			
			bottomBoundsZ = zTMP[0:-1]
			topBoundsZ = zTMP[1:]
			
			# cubeBoundary-Array includes the boundaries of each rectangual sub-cube, ....
			# one entry reads: [xMin,xMax,yMin,yMax,zMax,zMin]
			cubeBoundy = np.zeros((subCubesInX*subCubesInY*subCubesInZ,6))
			
			# iterate over the amount of sub cubes in each domain direction
			subCubeCNT = 0 # sub cube counter counts one up for each sub-cube and sets up the second dimension of cubeBoundary-Array
			for n in range(0,subCubesInZ): # z-direction
				for m in range(0,subCubesInY): # y-direction
					for k in range(0,subCubesInX): # x-direction
						cubeBoundy[subCubeCNT,0] = bottomBoundsX[k]
						cubeBoundy[subCubeCNT,1] = topBoundsX[k]
						
						cubeBoundy[subCubeCNT,2] = bottomBoundsY[m]
						cubeBoundy[subCubeCNT,3] = topBoundsY[m]
						
						cubeBoundy[subCubeCNT,4] = bottomBoundsZ[n]
						cubeBoundy[subCubeCNT,5] = topBoundsZ[n]
						subCubeCNT+=1
			
			# parallelInterpolation is a sub-function defined functions.py: further setup and call of pool-processes for interpolation
			temp_vecP1 = parallelInterpolation(cubeBoundy, dof_expanded, points, vals, mesh.hmax(), param.SolutionDirectory)
	
		p1InterpolationTime = time.time()-p1time
		print('stacked interpolation using '+method+ ' for P1-elements took: '+str(p1InterpolationTime))
		for k,kk in zip(namesP1,range(0,len(namesP1))):  # loop over all P1 fields and assign the values to the previously initialzed FEM-Spaces
				if k == 'Re':
					MeanFlowDict[k].vector()[:] = param.Re
				else:
					MeanFlowDict[k].vector()[:] = ContractFromAximuthalAverage(dof_coordinatesP1,nTheta,temp_vecP1[:,kk])	
					
	p2time = time.time()
	# similar procedure as aboves if statement(p1CNT>0)	
	if p2CNT>0:
		if method == 'nearest griddata':
			temp_vecP2 = interpolate.griddata(points, valsP2,dof_expandedP2, method='nearest')
		elif method == 'linear griddata':
			################## Same routine as above for the P1 interpoaltion #######################
			temp_vecP2 = interpolate.griddata(points, valsP2,dof_expandedP2, method='linear', fill_value=0)	
			#nansIDX = np.argwhere(temp_vecP2==np.nan)		
			nansIDX = np.isnan(temp_vecP2)  
			print(nansIDX[nansIDX==True])
			print(nansIDX[nansIDX==True].shape)
			temp_vecP2[nansIDX]=0
			nansIDX=np.isnan(temp_vecP2)
			print(nansIDX)
		elif method == 'parallel':
			# divisionFactor tells how many cuts are performed for one dimesion (3D and divionsFactor = 2 results in 8 cubes)
			dof_expanded = dof_expandedP2
			vals = valsP2
			subCubesInX = 2
			subCubesInY = 3
			subCubesInZ = 2
			
			dofMin_x = min(dof_expanded[:,0])
			dofMax_x = max(dof_expanded[:,0])
			dofMin_y = min(dof_expanded[:,1])
			dofMax_y = max(dof_expanded[:,1])
			dofMin_z = min(dof_expanded[:,2])
			dofMax_z = max(dof_expanded[:,2])
			xTMP = np.linspace(dofMin_x,dofMax_x,subCubesInX+1)
			yTMP = np.linspace(dofMin_y,dofMax_y,subCubesInY+1)
			zTMP = np.linspace(dofMin_z,dofMax_z,subCubesInZ+1)
			
			bottomBoundsX = xTMP[0:-1]
			topBoundsX = xTMP[1:]
			
			bottomBoundsY = yTMP[0:-1]
			topBoundsY = yTMP[1:]
			
			bottomBoundsZ = zTMP[0:-1]
			topBoundsZ = zTMP[1:]
			
			cubeBoundy = np.zeros((subCubesInX*subCubesInY*subCubesInZ,6))
			subCubeCNT = 0
			for n in range(0,subCubesInZ):
				for m in range(0,subCubesInY):
					for k in range(0,subCubesInX):
						cubeBoundy[subCubeCNT,0] = bottomBoundsX[k]
						cubeBoundy[subCubeCNT,1] = topBoundsX[k]
						
						cubeBoundy[subCubeCNT,2] = bottomBoundsY[m]
						cubeBoundy[subCubeCNT,3] = topBoundsY[m]
						
						cubeBoundy[subCubeCNT,4] = bottomBoundsZ[n]
						cubeBoundy[subCubeCNT,5] = topBoundsZ[n]
						subCubeCNT+=1
					
			temp_vecP2 = parallelInterpolation(cubeBoundy, dof_expanded, points, vals, mesh.hmax(), param.SolutionDirectory)


		p2InterpolationTime = time.time()-p2time
		print('stacked interpolation using '+method+ ' for P2-elements took: '+str(p2InterpolationTime))
		compCNT = 0
		# the interpolation above was performed for stacked quantities to avoid the triangulation for each quantity separately
		# split up the stacked quantties here again and potin to the corerct label in MeanFlowDict
		for k,kk in zip(namesP2,range(0,len(namesP2))):
			if k == 'u':
				# Vectorspace is subdivided into sub spaces
				# -> get the dofIDX is a index list including the indices pointing to the dof_coordinates of the currently considered sub space
				dofIDX = FEMSpaces.FunctionSpaceVectorVelocity.sub(compCNT).dofmap().dofs()           
				# setting the radial component to zero on the axis
				if compCNT==7: # to activate this feature set number to 1
					idxR0 = np.argwhere(dof_coordinatesP2[:,1]==0)
					ur = ContractFromAximuthalAverage(dof_coordinatesP2,nTheta,temp_vecP2[:,kk])
					ur[idxR0] = 0
					MeanFlowDict[k].vector()[dofIDX]=ur
				else:
					MeanFlowDict[k].vector()[dofIDX]=ContractFromAximuthalAverage(dof_coordinatesP2,nTheta,temp_vecP2[:,kk])
				
				compCNT+=1
				print(k+str(compCNT))
			else:
				print(k)
				MeanFlowDict[k].vector()[:] = ContractFromAximuthalAverage(dof_coordinatesP2,nTheta,temp_vecP2[:,kk])
	
		
	print('whole process, including contractingfrom azmithally expaned grid, took: '+ str(time.time()-p1time))	
	return MeanFlowDict, notInFileList
	
def nutTKEnut(MeanFlowDict, param, FEMSpaces):
	mesh=FEMSpaces.P2.mesh()
	from sklearn.neighbors import KNeighborsRegressor, RadiusNeighborsRegressor
	# find points on axis of Felics grid (axisPoints) -> variable xx is are x-coordinates of axisPoints
	dof_coordinates = FEMSpaces.P2.tabulate_dof_coordinates()
	# get dof-Indices of the coordinates with zero radius , i.e. of the points on the axis
	axisIDX = dof_coordinates[:,1]==0
	# get coordaintes of these points explicitely 
	axisPoints = dof_coordinates[axisIDX,:]
	# xx is the x-coordiantes only of the axis coordinates
	xx = axisPoints[:,0]
	
	#compute turbulent kinetic energy (TKE) and initialize function space for length scale 
	MeanFlowDict['TKE']=Function(FEMSpaces.P2)
	MeanFlowDict['TKE'].vector()[:] = 0.5*(MeanFlowDict['rstxx'].vector()[:] + MeanFlowDict['rstrr'].vector()[:] + MeanFlowDict['rsttt'].vector()[:])
	# init the P2-FEM field for the characterlistic length scale lm
	# computation of lm is similar to flame thickness in Poinsot Book
	MeanFlowDict['lm']=Function(FEMSpaces.P2)
	MeanFlowDict['lm'].vector()[:] = 1
	
	# UxNum is the numerical derivate of ux_mean wrt to radius - NOT REQUIRED AS P2-Space (only for plotting in paraview)
	MeanFlowDict['UxNum']=Function(FEMSpaces.P2)
	# UxNumMax is the maximum of numerical derivate of ux_mean wrt to radius at the corresponding streamwise position - NOT REQUIRED AS P2-Space (only for plotting in paraview)
	MeanFlowDict['UxNumMax']=Function(FEMSpaces.P2)
	
	# since velocity is vectorFEMspace, require to find the indices where the ux-component is located
	dofIDX_ux = FEMSpaces.FunctionSpaceVectorVelocity.sub(0).dofmap().dofs()    
	# init a temporarily dictionary entry only for ux-component 
	MeanFlowDict['ux'] = Function(FEMSpaces.P2)
	# and fill it with ux values from vector element space
	MeanFlowDict['ux'].vector()[:] = MeanFlowDict['u'].vector()[dofIDX_ux]
	
	
		
	tol = 20*mesh.hmin() # tolerance to build the band in which we look for local y-max on Felics grid
	
	# in the following the derivatives and lm will be computed
	# to avoid noise due to derivates on nearest-interpolated grids, everything will be a little smoothened (smoothFieldWithKernel)
	# smooth the ux_mean field a little -based on this, the derivative will be computed
#	MeanFlowDict['ux_smooth'] = Function(FEMSpaces.P2)
#	MeanFlowDict['ux_smooth'] = smoothFieldWithKernel(MeanFlowDict['ux'],FEMSpaces,1)
	# derivative of ux wrt. radius
	MeanFlowDict['dux_dr'] = Function(FEMSpaces.P2)
	MeanFlowDict['dux_dr'].vector()[:] = project(MeanFlowDict['ux'].dx(1),FEMSpaces.P2).vector()[:]
	# and subsequent smoothin to avoid significant noise 
#	MeanFlowDict['dux_dr_smooth'] = Function(FEMSpaces.P2)
#	MeanFlowDict['dux_dr_smooth'] = smoothFieldWithKernel(MeanFlowDict['dux_dr'],FEMSpaces,10)
	# since we are on unstructured grid its not trivial to find the max on each streamwise pos.
	# strategy: each xx-point (streamwise point) will get a vertical line with equidistant spacing 
	# and a corresponding tolerance band where points wihtin the tolerance are taken into account for interpoaltion
	# in complex geoemtries (non-convex geometries) this might result in issues
	# the relevant quantities will be interpolated onto this line
	# then the maximum will be determined
	# initialize arrays to store the x-coordinates(pointsX), y-coordinates(pointsY), 
	# length-scales(LM), and length-scale(lm) of all the vertical lines within the loop
	# with each loop iteration these lists will be expanded
	pointsX = [] 
	pointsY = []
	LM = []
	UxNum = []
	UxNumMax = []
	lm = np.ones(xx.shape[0])
	# k: x-coordiante itself; kk: iteration counter of loop
	for k,kk in zip(xx, range(0,xx.shape[0])):
		idx = np.argwhere( (dof_coordinates>k-tol) & (dof_coordinates<k+tol)).flatten() # find dof-coordinates within tolerance band
		ymax = dof_coordinates[idx,1].max() # find local y-max of Felics grid within the tolerance band
		Ny = 500 # number of points distributes equally spaced on vertical line
		lineY = np.linspace(0,ymax,Ny) # create vertical line
		pointsY = np.concatenate([pointsY,lineY]) # append y-coordinates vector with the y-coordaintes of the new line
		pointsX = np.concatenate([pointsX,k*np.ones(lineY.shape[0])])  # append x-coordinates vector with the y-coordaintes of the new line
		lineCoords = np.vstack([k*np.ones(lineY.shape[0]), lineY]).T # build coordiante-vector ([x1,y1],[x2,y2],...) of current line to interpolate onto
		originCoords = dof_coordinates[idx,:] # cut dof-coordiantes of Felics grid by tolerance band (not nessecary, but inscreases speed)
		originValsUx = MeanFlowDict['ux'].vector()[idx]
		originVals_dux_dr = MeanFlowDict['dux_dr'].vector()[idx]
		
		neighP1 = KNeighborsRegressor(n_neighbors=1, weights='distance') # initialize KneighborRegressor (find only the 1 closest neighbor, weighted by distance)
		neighP1.fit(originCoords,originValsUx) # build interpolator
		Ux = neighP1.predict(lineCoords) # interpolate onto the eqully spaced points on line
		
		neighP1 = KNeighborsRegressor(n_neighbors=1, weights='distance') # initialize KneighborRegressor (find only the 1 closest neighbor, weighted by distance)
		neighP1.fit(originCoords,originVals_dux_dr) # build interpolator
		UxDiff = neighP1.predict(lineCoords) # interpolate onto the eqully spaced points on line

		maxUxDiff = max(abs(UxDiff))
		
		UxNum = np.concatenate([UxNum,UxDiff])
		UxNumMax = np.concatenate([UxNumMax,maxUxDiff*np.ones(lineY.shape[0])]) 
		
		UxMax = Ux[0]
		UxMin = min(Ux)
#		UxMin = 0
		
		lm[kk] = (UxMax-UxMin)/maxUxDiff
		LM = np.concatenate([LM,lm[kk]*np.ones(lineY.shape[0])])

	
	# set up coordinate vector, MeanFlowdict(lm) and nut to interpolate onto	
	equiPoints = np.vstack([pointsX,pointsY]).T
	
	# interpoalte the line data of the lists filled inside the loop into unstructured grid
	MeanFlowDict['UxNum'].vector()[:] = interpolate.griddata(equiPoints, UxNum, dof_coordinates, method='nearest')
	MeanFlowDict['UxNumMax'].vector()[:] = interpolate.griddata(equiPoints, UxNumMax, dof_coordinates, method='nearest')
	
	# smooth results field 
	MeanFlowDict['lm'].vector()[:] = interpolate.griddata(equiPoints, LM, dof_coordinates, method='nearest')
#	MeanFlowDict['lm'] = smoothFieldWithKernel(MeanFlowDict['lm'],FEMSpaces,10)
	
	MeanFlowDict['nu_t']=Function(FEMSpaces.P2)
	MeanFlowDict['nuturb']=Function(FEMSpaces.P1)
	
	# negative turbulent kinetic energy is unphysical and might occur due to interpolation -> set to zero
	IDX = MeanFlowDict['TKE'].vector()[:]<0
	tmp = MeanFlowDict['TKE'].vector()[:]
	idx = tmp<0
	tmp[idx]=0
	MeanFlowDict['TKE'].vector()[:] = tmp
	#MeanFlowDict['TKE'].vector()[IDX] = 1234

	# interpolate nu_t onto P1 fcuntion space to match the Re-number function space
	MeanFlowDict['nu_t'].vector()[:] = param.TKEnutConst*MeanFlowDict['lm'].vector()[:]*MeanFlowDict['TKE'].vector()[:]**0.5
	MeanFlowDict['nuturb'] = FenicsInterpolate(MeanFlowDict['nu_t'], FEMSpaces.P2)
#	MeanFlowDict['nuturb'] = smoothFieldWithKernel(MeanFlowDict['nuturb'],FEMSpaces,6)
	# turbulent viscosity has to be opositiv
	tmp = MeanFlowDict['nuturb'].vector()[:]
	tmp = abs(tmp)

	return MeanFlowDict
 
def nutFromBoussinesq(MeanFlowDict, param, FEMSpaces):
# it takes the Boussinesq ansatz an rearranges the equation such that nu_t is computed
# the six different nu_t resulting from this ansatz is accounted for by a least square ansatz and already included in the underneath implementation
	
	# since velcoities come from vector-valued FE-spaces, first, the correct indices have to be found
	dofIDX_ux = FEMSpaces.FunctionSpaceVectorVelocity.sub(0).dofmap().dofs()   
	dofIDX_ur = FEMSpaces.FunctionSpaceVectorVelocity.sub(1).dofmap().dofs()    
	
	MeanFlowDict['ux'] = Function(FEMSpaces.P2)
	MeanFlowDict['ur'] = Function(FEMSpaces.P2)
	# setup fields consisting from components
	MeanFlowDict['ux'].vector()[:] = MeanFlowDict['u'].vector()[dofIDX_ux]
	MeanFlowDict['ur'].vector()[:] = MeanFlowDict['u'].vector()[dofIDX_ur]
	  
	# raidus
	R = FEMSpaces.P2.tabulate_dof_coordinates()[:,1]
	# to avoid division by zero, set R(R=0)=NaN
	# this results in NaN later on in the results which will be set to the molecular viscosity
	R[R==0] = np.nan
	
	# HARD_CODED: all derivatives are performed in cylindrical coordinates
	# derivatives are smoothened to avoid noise due to derivative of nearest interpolated data !!!!!'''''''' ############
	# also do not need all intermediate results such as rstij and dui_dj a FEMspace in MeanflowDict, however it is handy to observe them later in paraview
	MeanFlowDict['dux_dx'] = Function(FEMSpaces.P2)
	MeanFlowDict['dux_dx'].vector()[:] = project(MeanFlowDict['ux'].dx(0),FEMSpaces.P2).vector()[:]
#	MeanFlowDict['dux_dx'] = smoothFieldWithKernel(MeanFlowDict['dux_dx'],FEMSpaces,15)
	
	MeanFlowDict['dux_dr'] = Function(FEMSpaces.P2)
	MeanFlowDict['dux_dr'].vector()[:] = project(MeanFlowDict['ux'].dx(1),FEMSpaces.P2).vector()[:]
#	MeanFlowDict['dux_dr'] = smoothFieldWithKernel(MeanFlowDict['dux_dr'],FEMSpaces,15)
	
	MeanFlowDict['dux_dt'] = Function(FEMSpaces.P2)
	MeanFlowDict['dux_dt'].vector()[:] = MeanFlowDict['ux'].vector()[:]*0
#	MeanFlowDict['dux_dt'] = smoothFieldWithKernel(MeanFlowDict['dux_dt'],FEMSpaces,15)
	
	MeanFlowDict['dur_dx'] = Function(FEMSpaces.P2)
	MeanFlowDict['dur_dx'].vector()[:] = project(MeanFlowDict['ur'].dx(0),FEMSpaces.P2).vector()[:]
#	MeanFlowDict['dur_dx'] = smoothFieldWithKernel(MeanFlowDict['dur_dx'],FEMSpaces,15)
	
	MeanFlowDict['dur_dr'] = Function(FEMSpaces.P2)
	MeanFlowDict['dur_dr'].vector()[:] = project(MeanFlowDict['ur'].dx(1),FEMSpaces.P2).vector()[:]
#	MeanFlowDict['dur_dr'] = smoothFieldWithKernel(MeanFlowDict['dur_dr'],FEMSpaces,15)
	
	MeanFlowDict['dut_dt'] = Function(FEMSpaces.P2)
	tmp = MeanFlowDict['ur'].vector()[:]/R
	tmp2 = project(MeanFlowDict['ur'].dx(1),FEMSpaces.P2).vector()[:]
	idx = np.isnan(R)
#	tmp[idx] = tmp2[idx]
	tmp[idx] = 0
	MeanFlowDict['dut_dt'].vector()[:] = tmp 
#	MeanFlowDict['dut_dt'] = smoothFieldWithKernel(MeanFlowDict['dut_dt'],FEMSpaces,15)
	
	
	MeanFlowDict['dut_dx'] = Function(FEMSpaces.P2)
	MeanFlowDict['dut_dx'].vector()[:] = project(MeanFlowDict['ut'].dx(0),FEMSpaces.P2).vector()[:]
#	MeanFlowDict['dut_dx'] = smoothFieldWithKernel(MeanFlowDict['dut_dx'],FEMSpaces,15)
	
	MeanFlowDict['dut_dr'] = Function(FEMSpaces.P2)
	MeanFlowDict['dut_dr'].vector()[:] = project(MeanFlowDict['ut'].dx(1),FEMSpaces.P2).vector()[:]
#	MeanFlowDict['dut_dr'] = smoothFieldWithKernel(MeanFlowDict['dut_dr'],FEMSpaces,15)
	
	MeanFlowDict['dur_dt'] = Function(FEMSpaces.P2)
	tmp = -MeanFlowDict['ut'].vector()[:]/R
	tmp2 = project(MeanFlowDict['ut'].dx(1),FEMSpaces.P2).vector()[:]
	idx = np.isnan(R)
	tmp[idx] = tmp2[idx]
	MeanFlowDict['dur_dt'].vector()[:] = tmp 
#	MeanFlowDict['dur_dt'] = smoothFieldWithKernel(MeanFlowDict['dur_dt'],FEMSpaces,15)
	
	MeanFlowDict['TKE']=Function(FEMSpaces.P2)
	MeanFlowDict['TKE'].vector()[:] = 0.5*(MeanFlowDict['rstxx'].vector()[:] + MeanFlowDict['rstrr'].vector()[:] + MeanFlowDict['rsttt'].vector()[:])
	 
	numerator = (-MeanFlowDict['rstxx'].vector()[:]+2.0/3.0*MeanFlowDict['TKE'].vector()[:])\
	*(MeanFlowDict['dux_dx'].vector()[:]+MeanFlowDict['dux_dx'].vector()[:])\
	+(-MeanFlowDict['rstxr'].vector()[:])*(MeanFlowDict['dux_dr'].vector()[:]+MeanFlowDict['dur_dx'].vector()[:])\
	+(-MeanFlowDict['rstxt'].vector()[:])*(MeanFlowDict['dux_dt'].vector()[:]+MeanFlowDict['dut_dx'].vector()[:])\
	+(-MeanFlowDict['rstxr'].vector()[:])*(MeanFlowDict['dur_dx'].vector()[:]+MeanFlowDict['dux_dr'].vector()[:])\
	+(-MeanFlowDict['rstrr'].vector()[:]+2.0/3.0*MeanFlowDict['TKE'].vector()[:])\
	*(MeanFlowDict['dur_dr'].vector()[:]+MeanFlowDict['dur_dr'].vector()[:])\
	+(-MeanFlowDict['rstrt'].vector()[:])*(MeanFlowDict['dur_dt'].vector()[:]+MeanFlowDict['dut_dr'].vector()[:])\
	+(-MeanFlowDict['rstxt'].vector()[:])*(MeanFlowDict['dut_dx'].vector()[:]+MeanFlowDict['dux_dt'].vector()[:])\
	+(-MeanFlowDict['rstrt'].vector()[:])*(MeanFlowDict['dut_dr'].vector()[:]+MeanFlowDict['dur_dt'].vector()[:])\
	+(-MeanFlowDict['rsttt'].vector()[:]+2.0/3.0*MeanFlowDict['TKE'].vector()[:])\
	*(MeanFlowDict['dut_dt'].vector()[:]+MeanFlowDict['dut_dt'].vector()[:])
	 
	denominator=+(MeanFlowDict['dux_dx'].vector()[:]+MeanFlowDict['dux_dx'].vector()[:])**2\
	+(MeanFlowDict['dux_dr'].vector()[:]+MeanFlowDict['dur_dx'].vector()[:])**2\
	+(MeanFlowDict['dux_dt'].vector()[:]+MeanFlowDict['dut_dx'].vector()[:])**2\
	+(MeanFlowDict['dur_dx'].vector()[:]+MeanFlowDict['dux_dr'].vector()[:])**2\
	+(MeanFlowDict['dur_dr'].vector()[:]+MeanFlowDict['dur_dr'].vector()[:])**2\
	+(MeanFlowDict['dur_dt'].vector()[:]+MeanFlowDict['dut_dr'].vector()[:])**2\
	+(MeanFlowDict['dut_dx'].vector()[:]+MeanFlowDict['dux_dt'].vector()[:])**2\
	+(MeanFlowDict['dut_dr'].vector()[:]+MeanFlowDict['dur_dt'].vector()[:])**2\
	+(MeanFlowDict['dut_dt'].vector()[:]+MeanFlowDict['dut_dt'].vector()[:])**2
	 
	MeanFlowDict['nuturb']=Function(FEMSpaces.P2)
	MeanFlowDict['nuturb'].vector()[:] = numerator / denominator
	#MeanFlowDict['nuturb'] = FenicsInterpolate(MeanFlowDict['nu_t'], FunctionSpace(FEMSpaces.P2.mesh(), 'P', 1))
	tmp = abs(MeanFlowDict['nuturb'].vector()[:])
#	tmp = MeanFlowDict['nuturb'].vector()[:]
#	tmp[tmp<0] = 0
	MeanFlowDict['nuturb'].vector()[:] = tmp
	del MeanFlowDict['ux']
	del MeanFlowDict['ur']
	return MeanFlowDict

def nutFromReducedBoussinesq_xr(MeanFlowDict, param, FEMSpaces):
	MeanFlowDict['nuturb']=Function(FEMSpaces.P2)
	# instead of computing all six different nu_t arising from Boussinesq ansatz
	# only a reduced ansatz is computed such as in Semeraro-Paper 2015 (International Journal of Heat and Fluid Flow)
	# that is: nut = -rstxr/(dux_dr+epsilon) with epsilon to avoid a singularity in the denominator
	# HOWEVER: applying the epsilon is not trivial since just adding epsilon will shift slighlty negative values to zero
	# HENCE: slightly negative values (-1e-10) are shifted with -epsilon and slightly positive (1e-10) with +epsilon
	# HENCE: denominator = dux_dr+(sign(dux_dr)*epsilon)
	
	# since velcoities come from vector-valued FE-spaces, first, the correct indices have to be found
	dofIDX_ux = FEMSpaces.FunctionSpaceVectorVelocity.sub(0).dofmap().dofs()    
	
	MeanFlowDict['ux'] = Function(FEMSpaces.P2)
	# setup fields consisting from components
	MeanFlowDict['ux'].vector()[:] = MeanFlowDict['u'].vector()[dofIDX_ux]
	  
	# derivative
	MeanFlowDict['dux_dr'] = Function(FEMSpaces.P2)
	#MeanFlowDict['dux_dr'].vector()[:] = project(MeanFlowDict['ux'].dx(1),FEMSpaces.P2).vector()[:]
	MeanFlowDict['dux_dr'] = MeanFlowDict['ux'].dx(1)
	epsilon= 0.001 *np.max(np.abs(MeanFlowDict['dux_dr'].vector()[:]))
	# get signs of derivative
	signs = np.sign(MeanFlowDict['dux_dr'].vector()[:])
	#set up numerator and denominator as the procedure description above
	denominator = MeanFlowDict['dux_dr'].vector()[:]+(signs*epsilon)
	numerator=-MeanFlowDict['rstxr']	.vector()[:]
	# execute
	MeanFlowDict['nuturb'].vector()[:] = numerator / denominator
	#project onto P1, since Reynoldsnumber also lives on P1
	#MeanFlowDict['nuturb'] = FenicsInterpolate(MeanFlowDict['nu_t'], FunctionSpace(FEMSpaces.P2.mesh(), 'P', 1))
	# negative eddy-viscosity must not exist --> set everything smaller than zero to zero
	# later there will be nu_t+nu_mol, sich that there is no dange of dividing by zero late when Re~1/nu_eff
	tmp = MeanFlowDict['nuturb'].vector()[:]
	idx = np.argwhere(tmp<0)
	tmp[idx] = 0
	MeanFlowDict['nuturb'].vector()[:] = tmp 
	del MeanFlowDict['ux']
	return MeanFlowDict
	
 
def importAVBPFile(param,FEMSpaces):
    ## Define Dictionaries
	MeanFlowDict={}
	notInFileList=[]
	if param.CoordinateSystem in ['Cylindrical']: #for Mario's avbp data
		## Get mesh data
	#   mesh=FEMSpaces.P2.mesh()
	    ## Get Mean flow names
		#nameListMean=getMeanFlowFieldNames(param)
		nameListMean=param.Case.getMeanFlowFieldNames()
		## Get AVBP mesh file path
		filePath=param.FlowInput.MeanFlowFilePath
		# Open hdf5 file
		h5file = h5py.File(param.BaseFlowMesh, 'r')
	    # Extract coordinates from mesh file and make them a list
		RawFlowDict={}
		RawFlowDict['x_raw']=np.array(h5file['Coordinates']['x'])
		RawFlowDict['z_raw']=np.array(h5file['Coordinates']['y'])
		RawFlowDict['y_raw']=np.array(h5file['Coordinates']['z'])
		del h5file
	    #Now open AVBP solution file
	    # read the time averaged velocities
		h5file = h5py.File(filePath, 'r')
		RawFlowDict['ux_raw'] = np.array(h5file['Average']['u'])
		RawFlowDict['uz_raw'] = np.array(h5file['Average']['v'])
		RawFlowDict['uy_raw'] = np.array(h5file['Average']['w'])
		if param.ReynoldsNumberSource=='Boussinesq':
	        # read the time averaged reynolds stresses, u2 means average(uu)
			RawFlowDict['rstxx_raw'] = np.array(h5file['Average']['u2'])-np.array(h5file['Average']['u'])**2 #np.array(h5file['Average']['u2'])-np.array(h5file['Average']['u'])**2
			RawFlowDict['rstyy_raw'] = np.array(h5file['Average']['v2'])-np.array(h5file['Average']['v'])**2
			RawFlowDict['rstzz_raw'] = np.array(h5file['Average']['w2'])-np.array(h5file['Average']['w'])**2
			RawFlowDict['rstxy_raw'] = np.array(h5file['Average']['uv'])-np.array(h5file['Average']['u'])*np.array(h5file['Average']['v'])
			RawFlowDict['rstxz_raw'] = np.array(h5file['Average']['uw'])-np.array(h5file['Average']['u'])*np.array(h5file['Average']['w'])
			RawFlowDict['rstyz_raw'] = np.array(h5file['Average']['vw'])-np.array(h5file['Average']['v'])*np.array(h5file['Average']['w'])
	#   RawFlowDict['T_raw']  = np.array(h5file['Average']['T'])
	     
	     
		InputFlowDict=CoordinateTransformation('Cart2Cyl',RawFlowDict, param)
	     
		x1_mean=np.array(InputFlowDict['x'])
		x2_mean=np.array(InputFlowDict['y'])
		x3_mean=np.array(InputFlowDict['z'])
		x1_mean=x1_mean.tolist()
		x2_mean=x2_mean.tolist()
		x3_mean=x3_mean.tolist()
		#FieldsToDelete=[]
		keys=list(InputFlowDict.keys())
		del h5file
	     
		for name in keys:
			if name[0:3]=='rho' and len(name)>3:
				InputFlowDict[name]=InputFlowDict[name]/InputFlowDict['rho']
				InputFlowDict[name[3:]] = InputFlowDict[name]
				del InputFlowDict[name]
			if name in getSpeciesListMean(param):
				InputFlowDict[name]=InputFlowDict[name]/InputFlowDict['rho']
			if name[0:3] == 'rst':
				MeanFlowDict[name]=Function(FEMSpaces.P2)
				dof_coordinates = FEMSpaces.P2.tabulate_dof_coordinates() 
				try:
	                # Expanding into azimuthal direction
					dof_expanded=ExpandForAzimuthalAverage(dof_coordinates,100)
					print('Performing interpolation for ' + name)
					temp_vec=np.squeeze(interpolate.griddata(list(map(list, zip(*[x1_mean,x2_mean,x3_mean]))), InputFlowDict[name],dof_expanded, method='nearest'))
	                # Collapsing back to initial mesh
					MeanFlowDict[name].vector()[:]=ContractFromAximuthalAverage(dof_coordinates,100,temp_vec)
				except KeyError as e:
					notInFileList.append(name)
					print( 'I got a KeyError - reason ' + str(e))
		for name in nameListMean:
			if name[0]=='u':
				MeanFlowDict[name]=Function(FEMSpaces.P2)
				dof_coordinates = FEMSpaces.P2.tabulate_dof_coordinates()                      
			else:
				print(name)
				MeanFlowDict[name]=Function(FEMSpaces.P1)
				dof_coordinates = FEMSpaces.P1.tabulate_dof_coordinates()                      
			try:
	            # Expanding into azimuthal direction
				dof_expanded=ExpandForAzimuthalAverage(dof_coordinates,100)
				print('Performing interpolation for ' + name)
				temp_vec=np.squeeze(interpolate.griddata(list(map(list, zip(*[x1_mean,x2_mean,x3_mean]))), InputFlowDict[name],dof_expanded, method='nearest'))
				# Collapsing back to initial mesh
				MeanFlowDict[name].vector()[:]=ContractFromAximuthalAverage(dof_coordinates,100,temp_vec)
			except KeyError as e:
				notInFileList.append(name)
				print( 'I got a KeyError - reason ' + str(e))
	elif ('LowMach' in param.FlowMode or 'Reacting' in param.FlowMode) and param.nVelocityComponents == 2 and 'Cartesian' in param.CoordinateSystem: #for slot flame avbp data
		mesh=FEMSpaces.P2.mesh()
		## Get Mean flow names
		#nameListMean=getMeanFlowFieldNames(param)
		nameListMean=param.Case.getMeanFlowFieldNames()
		## Get AVBP mesh file path
		filePath=param.FlowInput.MeanFlowFilePath
		# Open hdf5 file
		h5file = h5py.File(param.BaseFlowMesh, 'r')
		# Extract coordinates from mesh file and make them a list
		x1_mean=np.array(h5file['Coordinates']['x'])
		x2_mean=np.array(h5file['Coordinates']['y'])
		x1_mean=x1_mean.tolist()
		x2_mean=x2_mean.tolist()
		del h5file
		#Now open AVBP solution file
		h5file = h5py.File(filePath, 'r')
		
		InputFlowDict={}
		translatorAVBP2Felics={'ux':'u','uy':'v','T':'temperature'}
		for element in h5file.items():
			if element[0]=='GaseousPhase' or element[0]=='RhoSpecies' or element[0]=='Additionals':
				for name in element[1].keys():
					InputFlowDict[name]=np.array(element[1][name])
		del h5file
		FieldsToDelete=[]
		keys=list(InputFlowDict.keys())
		#hdf5file = HDF5File(mesh.mpi_comm(),filePath[0:-4]+"__"+param.MeshPath[0:-3].split('/')[-1]+"hdf5", 'w')
		print(keys)
		for name in keys:
			if name[0:3]=='rho' and len(name)>3:
				print(name)
				InputFlowDict[name]=InputFlowDict[name]/InputFlowDict['rho']
				InputFlowDict[name[3:]] = InputFlowDict[name]
				del InputFlowDict[name]
			if name in getSpeciesListMean(param):
				InputFlowDict[name]=InputFlowDict[name]/InputFlowDict['rho']
		for name in nameListMean:
			if name[0]=='u':
				MeanFlowDict[name]=Function(FEMSpaces.P2)
				dof_coordinates = FEMSpaces.P2.tabulate_dof_coordinates()                      
			else:
				#print(name)
				MeanFlowDict[name]=Function(FEMSpaces.P1)
				dof_coordinates = FEMSpaces.P1.tabulate_dof_coordinates()                      
			try:
				printDebug(param.debug,"Performing flame interpolation for " + name)
				if name in translatorAVBP2Felics.keys():
					MeanFlowDict[name].vector()[:]=np.squeeze(interpolate.griddata(list(map(list, zip(*[x1_mean,x2_mean]))), InputFlowDict[translatorAVBP2Felics[name]],dof_coordinates, method='linear'))	
				else:
					MeanFlowDict[name].vector()[:]=np.squeeze(interpolate.griddata(list(map(list, zip(*[x1_mean,x2_mean]))), InputFlowDict[name],dof_coordinates, method='linear'))
			except KeyError as e:
				notInFileList.append(name)
				if name=='Re':
					print('Re does not need to be interpolated in this case.')
				else:
					print( 'I got a KeyError - reason ' + str(e))

		if param.nVelocityComponents == 2:
			printDebug(param.debug,"Constructing vector u")
			MeanFlowDict['u'] = project(Expression(("u1", "u2"),\
 				u1=MeanFlowDict['ux'],\
				u2=MeanFlowDict['uy'],\
				degree=2), FEMSpaces.FunctionSpaceVectorVelocity)
		elif param.nVelocityComponents == 3:
			print('Constructing vector u not implemented!')
	return MeanFlowDict, notInFileList
 
def importFelicsFile(param, FEMSpaces):
	import time
    ## Define Dictionaries
	MeanFlowDict={}
	notInFileList=[]
	RawFlowDict={}
	## Get mesh data
	mesh=FEMSpaces.P2.mesh()
	## Get Mean flow names
	nameListMean=param.MeanList
	## Get AVBP mesh file path
	filePath=param.FlowInput.MeanFlowFilePath
	# Open hdf5 file
	h5file = h5py.File(filePath, 'r')

	
	print('nameListMean')
	print(nameListMean)
	#Analyze dimension and coordinate system of input data based on the available coordinates
	# Copy the respective coordinates to the RawFlowDict at the same time
	ndim_rawData=0
	CoordinateSystemInputData='Unknown'
	if 'x' in list(h5file['MeanFlow'].keys()):
		RawFlowDict['x']=np.array(h5file['MeanFlow']['x'])
		ndim_rawData += 1
	if 'y' in list(h5file['MeanFlow'].keys()):
		RawFlowDict['y']=np.array(h5file['MeanFlow']['y'])
		ndim_rawData += 1
		CoordinateSystemInputData='Cartesian'
	if 'z' in list(h5file['MeanFlow'].keys()):
		RawFlowDict['z']=np.array(h5file['MeanFlow']['z'])
		ndim_rawData += 1
		CoordinateSystemInputData='Cartesian'
	if 't' in list(h5file['MeanFlow'].keys()):
		RawFlowDict['t']=np.array(h5file['MeanFlow']['t'])

		ndim_rawData += 1
		if CoordinateSystemInputData ==  'Cartesian':
			printError('Ambigous input data: both cylindrical coordinates and cartesian coordinates present. Check input data!')
		else:
			CoordinateSystemInputData = 'Cylindrical'
	if 'r' in list(h5file['MeanFlow'].keys()):
		RawFlowDict['r']=np.array(h5file['MeanFlow']['r'])
		ndim_rawData += 1
		if CoordinateSystemInputData ==  'Cartesian':
			printError('Ambigous input data: both cylindrical coordinates and cartesian coordinates present. Check input data!')
		else:
			CoordinateSystemInputData = 'Cylindrical'
		
	# Now copy the velocity components depending on the respective coordinate system to the RawFlowDict
	if CoordinateSystemInputData == 'Cartesian':
		RawFlowDict['ux']=np.array(h5file['MeanFlow']['ux'])
		RawFlowDict['uy']=np.array(h5file['MeanFlow']['uy'])
		if ndim_rawData > 2:
			RawFlowDict['uz']=np.array(h5file['MeanFlow']['uz'])
	elif CoordinateSystemInputData =='Cylindrical':
		RawFlowDict['ux']=np.array(h5file['MeanFlow']['ux'])
		RawFlowDict['ur']=np.array(h5file['MeanFlow']['ur'])
		if 'ut' in list(h5file['MeanFlow'].keys()):
			RawFlowDict['ut']=np.array(h5file['MeanFlow']['ut'])
		else:
			notInFileList.append('ut')

	# Copy all remaining fields to the RawFlowDict
	for name in nameListMean: 
		# velocities were already copied, so only copy the rest
		if not name in ['u','ut']:
			# In case of Input-Output analysis, forcings should be read. If so, then this must be done for all components
			if  name[0] == 'u' and name not in ['ut','ut_forcing_r','ut_forcing_i']:
				for component in param.VelocityComponents:
					nameComponent=name[:1]+component+name[1:]
					if nameComponent in list(h5file['MeanFlow'].keys()):
						RawFlowDict[nameComponent]=np.array(h5file['MeanFlow'][nameComponent])
					else:
						notInFileList.append(nameComponent)
			# For all not velocity related fields...
			else:
				if name in list(h5file['MeanFlow'].keys()):
					RawFlowDict[name]=np.array(h5file['MeanFlow'][name])
				else:
					notInFileList.append(name)
	# Die mittelung muss in die Parameters und in der GUI abgefragt werden(Name in GUI: Averaging in homogeneous direction), anschliessend folgende Zeile loeschen
	param.SpatialBaseFlowAverageing=True
	linear = False
	
	if ndim_rawData == 3 and param.CoordinateSystem=='Cylindrical' and CoordinateSystemInputData=='Cartesian':
		InputFlowDict=CoordinateTransformation('Cart2Cyl',RawFlowDict, param)
	else:
		#InputFlowDict=CoordinateTransformation('None',RawFlowDict,param)
		InputFlowDict=RawFlowDict

	# For three dimensional dara bases restrict domain to reduce the number of basis points and accelerate the interpolation
	if ndim_rawData>2:	
		Bound = True
	# For 2D flows this most often is not necessary
	else:
		Bound=False
		
	if Bound:
		Xcoords = [k[0] for k in mesh.coordinates()]
		Ycoords = [k[1] for k in mesh.coordinates()]
		boundingBox = [min(Xcoords), max(Xcoords),\
					min(Ycoords), max(Ycoords),\
					-mesh.hmax(), mesh.hmax()]
		
		X = InputFlowDict['x']
		R = InputFlowDict['r']
		
		Xidx = np.argwhere( (X>boundingBox[0]) & (X<boundingBox[1])).flatten()
		Ridx = np.argwhere( ((R<boundingBox[3]))).flatten()
		
		IDX = np.intersect1d(Xidx, Ridx)
		x1_mean=np.array(InputFlowDict['x'][IDX])
		x2_mean=np.array(InputFlowDict['y'][IDX])
		x3_mean=np.array(InputFlowDict['z'][IDX])
		points = np.vstack((x1_mean,x2_mean,x3_mean)).T
		for k in InputFlowDict.keys():
			InputFlowDict[k] = InputFlowDict[k][IDX]
	else:
		x1_mean=np.array(InputFlowDict['x'])
		x1_mean=x1_mean.tolist()
		if param.CoordinateSystem in ['Cartesian']:
			x2_mean=np.array(InputFlowDict['y'])
		elif param.CoordinateSystem in ['Cylindrical']:
			x2_mean=np.array(InputFlowDict['r'])
		x2_mean=x2_mean.tolist()
		if ndim_rawData >2:
			x3_mean=np.array(InputFlowDict['z'])
			x3_mean=x3_mean.tolist()
			points = np.vstack((x1_mean,x2_mean,x3_mean)).T
		else:
			points = np.vstack((x1_mean,x2_mean)).T
	


	# get coordinates of FELiCS mesh
	dof_coordinatesP2 = FEMSpaces.P2.tabulate_dof_coordinates()	
	# In case of a cylindrical coordinate system in the FELiCS calculation and a 3D input flow, azimuthal averaging must be performed 
	if ndim_rawData == 3 and param.CoordinateSystem=='Cylindrical' and CoordinateSystemInputData=='Cartesian':	
		dof_InterpolationFELiCSMeshP2=ExpandForAzimuthalAverage(dof_coordinatesP2,100) 
	# For other cases the interplation points are identical to the FELiCS mesh coordinates
	else:
		dof_InterpolationFELiCSMeshP2=dof_coordinatesP2
	namesP2 = []
	
	
	# Get number of fields to interpolate
	nFieldsToInterpolate=len(nameListMean)+param.nVelocityComponents-1
	if param.AnalysisMode in ['Input-Output']:
		nFieldsToInterpolate+=2*(param.nVelocityComponents-1)
	nFieldsToInterpolate -= len(notInFileList)
	# Define vmatrix for interpolation basis values, valsP2
	valsP2=np.zeros((len(InputFlowDict['ux']),nFieldsToInterpolate+6))	

	m = 0	
	#Fill valsP2 with raw data
	print('nameListMean: '+str(nameListMean))
	print('InputFlowDict: '+str(InputFlowDict.keys()))
	for name in nameListMean:
		if name[0]=='u' and name not in ['ut','ut_forcing_r','ut_forcing_i']:
			MeanFlowDict[name]=Function(FEMSpaces.FunctionSpaceVectorVelocity)
			# All inplane velocity components are defined as vectors. Therefore for these, iterate through the components
			for component in param.VelocityComponents:
				#Get the name in plus component
				nameComponent=name[:1]+component+name[1:]
				if nameComponent in list(InputFlowDict.keys()):
					valsP2[:,m] = np.array(InputFlowDict[nameComponent])
					namesP2.append(nameComponent)
					m += 1
		# Do the same as above, for all scalars. Here no iteration through components is necessary
		else:
			MeanFlowDict[name]=Function(FEMSpaces.P2)
			if name in list(InputFlowDict.keys()):
				valsP2[:,m] = np.array(InputFlowDict[name])
				namesP2.append(name)
				# increment m
				m += 1			 

	#print("points")
	#print(points[:,0])
	#print(points[:,1])
	#print("valsP2")
	#print(valsP2[:,0])
	#print(valsP2[:,1])
	#print(valsP2[:,2])
	#print("dof_InterpolationFELiCSMeshP2")
	#print(dof_InterpolationFELiCSMeshP2[:,0])
	#print(dof_InterpolationFELiCSMeshP2[:,1])
	#print("shapes")
	#print("points")
	#print(np.shape(points))
	#print("valsP2")
	#print(np.shape(valsP2))
	#print("dof_InterpolationFELiCSMeshP2")
	#print(np.shape(dof_InterpolationFELiCSMeshP2))
	# Perform interpolation nearest (unprecise)
	temp_vecP2nearest = interpolate.griddata(points, valsP2,dof_InterpolationFELiCSMeshP2, method='nearest')
	# Perfrm interpolation nearest (more precise)
	try:
		temp_vecP2 = interpolate.griddata(points, valsP2,dof_InterpolationFELiCSMeshP2, method='linear')
		temp_vecP2[np.isnan(temp_vecP2)]=temp_vecP2nearest[np.isnan(temp_vecP2)]
	except:
		printWarning("Linear interpolation failed... Continue with nearest interplation. This may cause strong inaccuracies!")
		temp_vecP2=np.array(temp_vecP2nearest)
	# In case of nan values in the linear inerpolation result, fill the nan values with the results from nearest interpolation
	temp_vecP2[np.isnan(temp_vecP2)]=temp_vecP2nearest[np.isnan(temp_vecP2)]
	# Now write the interpolated Values to the MeanFlowDict
	m=0
	# Do that for every entry in nameListMean
	for name in nameListMean:
		if name[0]=='u' and name not in ['ut','ut_forcing_r','ut_forcing_i']:
			# All inplane velocity components are defined as vectors. Therefore for these, iterate through the components
			for idx,component in enumerate(param.VelocityComponents):
				nameComponent=name[:1]+component+name[1:]
				if not nameComponent in notInFileList:
					# Get the indices of the components entries in the vector
					dofIDX = FEMSpaces.FunctionSpaceVectorVelocity.sub(idx).dofmap().dofs()
					if ndim_rawData == 3 and param.CoordinateSystem=='Cylindrical' and CoordinateSystemInputData=='Cartesian':
						# For this case a azimuthal average is performed by using the function ContractFromAzimuthalAverage
						MeanFlowDict[name].sub(idx).vector()[dofIDX] = ContractFromAximuthalAverage(dof_coordinatesP2,100,temp_vecP2[:,m])
					else:
						# In this case a simple copy of the interolation results is sufficient
						MeanFlowDict[name].sub(idx).vector()[dofIDX]=np.array(temp_vecP2[:,m])
				# Increment m
					m+=1
		else:
			if (not name in notInFileList) or (name in ['rstxx','rstrr','rsttt','rstxr','rstxt','rstrt']) :
				# Do the same as above, for all scalars. Here no iteration through components is necessary
				if ndim_rawData == 3 and param.CoordinateSystem=='Cylindrical' and CoordinateSystemInputData=='Cartesian':
					MeanFlowDict[name].vector()[:] = ContractFromAximuthalAverage(dof_coordinatesP2,100,temp_vecP2[:,m])
				else:
					MeanFlowDict[name].vector()[:]=np.array(temp_vecP2[:,m])
				m+=1
	return MeanFlowDict, notInFileList

def importCGNSFile(param, FEMSpaces):
	import time
    ## Define Dictionaries
	MeanFlowDict={}
	notInFileList=[]
	## Get mesh data
	mesh=FEMSpaces.P2.mesh()
	## Get Mean flow names
	#nameListMean=getMeanFlowFieldNames(param)
	nameListMean=param.Case.getMeanFlowFieldNames()

	## Get AVBP mesh file path
	filePath=param.FlowInput.MeanFlowFilePath
	# Open hdf5 file
	h5file = h5py.File(filePath, 'r')
	
#	DictNames = assignNames(h5file)
    
	RawFlowDict={}
	
#	if 'Default Domain' in h5file['Base'].keys():
#		domain_string='fuid_mixture Domain'
#	elif 'Default Domain Modified' in h5file['Base'].keys():
#		domain_string='fuid_mixture'
	
	RawFlowDict['x_raw']=np.array(h5file['Base']['Default Domain']['GridCoordinates']['CoordinateX'][' data'])
	RawFlowDict['y_raw']=np.array(h5file['Base']['Default Domain']['GridCoordinates']['CoordinateY'][' data'])
	RawFlowDict['z_raw']=np.array(h5file['Base']['Default Domain']['GridCoordinates']['CoordinateZ'][' data'])
	RawFlowDict['ux_raw']=np.array(h5file['Base']['Default Domain']['FlowSolution']['Velocity.TrnavgX'][' data'])
	RawFlowDict['uy_raw']=np.array(h5file['Base']['Default Domain']['FlowSolution']['Velocity.TrnavgY'][' data'])
	RawFlowDict['uz_raw']=np.array(h5file['Base']['Default Domain']['FlowSolution']['Velocity.TrnavgZ'][' data'])
	if param.ReynoldsNumberSource == 'EddyViscosity':
		RawFlowDict['nut']=np.array(h5file['Base']['fluid_mixture']['FlowSolution']['ViscosityEddy'][' data'])
	if param.ReynoldsNumberSource == 'Boussinesq':
		nameListMean.extend(['rstxx', 'rstrr', 'rsttt', 'rstxr', 'rstxt', 'rstrt'])
		
		RawFlowDict['rstxx_raw'] = np.array(h5file['Base']['Default Domain']['FlowSolution']['StatisticalReynoldsStressuu'][' data'])
		RawFlowDict['rstyy_raw'] = np.array(h5file['Base']['Default Domain']['FlowSolution']['StatisticalReynoldsStressvv'][' data'])
		RawFlowDict['rstzz_raw'] = np.array(h5file['Base']['Default Domain']['FlowSolution']['StatisticalReynoldsStressww'][' data'])
		RawFlowDict['rstxy_raw'] = np.array(h5file['Base']['Default Domain']['FlowSolution']['StatisticalReynoldsStressuv'][' data'])
		RawFlowDict['rstxz_raw'] = np.array(h5file['Base']['Default Domain']['FlowSolution']['StatisticalReynoldsStressuw'][' data'])
		RawFlowDict['rstyz_raw'] = np.array(h5file['Base']['Default Domain']['FlowSolution']['StatisticalReynoldsStressvw'][' data'])
	if 'rho' in nameListMean:
		RawFlowDict['rho']=np.array(h5file['Base']['Default Domain']['FlowSolution']['SpeciesDensity'][' data'])
		
	
	# An mario: ndim muss bestimmt werden. Entweder aus dem input file auslesen (besser), oder in der GUI abfragn (Notloesung)
	ndim=3
	# Die mittelung muss in die Parameters und in der GUI abgefragt werden(Name in GUI: Averaging in homogeneous direction), anschliessend folgende Zeile loeschen
	param.SpatialBaseFlowAverageing=True
	Bound = 1
	linear = False
	
	if param.CoordinateSystem=='Cylindrical':
		InputFlowDict=CoordinateTransformation('Cart2Cyl',RawFlowDict, param)
	else:
		InputFlowDict=CoordinateTransformation('None',RawFlowDict,param)
    
	x1_mean=np.array(InputFlowDict['x'])
	x1_mean=x1_mean.tolist()
	if ndim>1:
		x2_mean=np.array(InputFlowDict['y'])
		x2_mean=x2_mean.tolist()
	if ndim>2:
		x3_mean=np.array(InputFlowDict['z'])
		x3_mean=x3_mean.tolist()
	if ndim==2:
		mean_coordinates=[x1_mean,x2_mean]
	if ndim==3:
		mean_coordinates=[x1_mean,x2_mean,x3_mean]
	
	if Bound == 1:
		Xcoords = [k[0] for k in mesh.coordinates()]
		Ycoords = [k[1] for k in mesh.coordinates()]
		boundingBox = [min(Xcoords), max(Xcoords),\
					min(Ycoords), max(Ycoords),\
					-mesh.hmax(), mesh.hmax()]
		
		X = InputFlowDict['x']
		R = InputFlowDict['r']
		
		Xidx = np.argwhere( (X>boundingBox[0]) & (X<boundingBox[1])).flatten()
		Ridx = np.argwhere( ((R<boundingBox[3]))).flatten()
		
		IDX = np.intersect1d(Xidx, Ridx)
		x1_mean=np.array(InputFlowDict['x'][IDX])
		x2_mean=np.array(InputFlowDict['y'][IDX])
		x3_mean=np.array(InputFlowDict['z'][IDX])
		points = np.vstack((x1_mean,x2_mean,x3_mean)).T
		for k in InputFlowDict.keys():
			InputFlowDict[k] = InputFlowDict[k][IDX]
	else:
		x1_mean=np.array(InputFlowDict['x'])
		x2_mean=np.array(InputFlowDict['y'])
		x3_mean=np.array(InputFlowDict['z'])
		x1_mean=x1_mean.tolist()
		x2_mean=x2_mean.tolist()
		x3_mean=x3_mean.tolist()
	
	
	p1CNT = 0
	p2CNT = 0
	for name in nameListMean:
		if name[0]=='u' or name[0:3] == 'rst':
			p2CNT += 1
		elif not name == 'Re':
			p1CNT += 1
	valsP1 = np.zeros((InputFlowDict['ux'].shape[0], p1CNT ))
	valsP2 = np.zeros((InputFlowDict['ux'].shape[0], p2CNT ))
	dof_coordinatesP1 = FEMSpaces.P1.tabulate_dof_coordinates()
	dof_expandedP1=ExpandForAzimuthalAverage(dof_coordinatesP1,100)
	dof_coordinatesP2 = FEMSpaces.P2.tabulate_dof_coordinates()
	dof_expandedP2=ExpandForAzimuthalAverage(dof_coordinatesP2,100) 
 
	#FieldsToDelete=[]
	keys=list(InputFlowDict.keys())
	namesP2 = []
	namesP1 = []
	for name in keys:
		if name[0:3]=='rho' and len(name)>3:
			InputFlowDict[name]=InputFlowDict[name]/InputFlowDict['rho']
			InputFlowDict[name[3:]] = InputFlowDict[name]
			del InputFlowDict[name]
		if name in getSpeciesListMean(param):
			InputFlowDict[name]=InputFlowDict[name]/InputFlowDict['rho']
	
	m = 0
	mm = 0
	for name in nameListMean:
		if name[0]=='u' or name[0:3] == 'rst':
			MeanFlowDict[name]=Function(FEMSpaces.P2)
			valsP2[:,mm] = InputFlowDict[name]
			namesP2.append(name)
			mm += 1
		else:
			MeanFlowDict[name]=Function(FEMSpaces.P1) 
			if not name == 'Re':
				valsP1[:,m] = InputFlowDict[name]
				m += 1
				namesP1.append(name)

	p1time = time.time()	
	if linear == True and p1CNT>0:
		print('started griddata for P1, stacked, linear')
		temp_vecP1 = interpolate.griddata(points, valsP1,dof_expandedP1, method='linear')
	elif linear == False and p1CNT>0:
		print('started griddata for P1, stacked, nearest')
		temp_vecP1 = interpolate.griddata(points, valsP1,dof_expandedP1, method='nearest')

	if p1CNT>0: print('time for griddata stacked command: ' +str(time.time()-p1time))
	for k,kk in zip(namesP1,range(0,len(namesP1))):
		if k == 'Re':
			MeanFlowDict[k].vector()[:] = param.Re
		else:
			MeanFlowDict[k].vector()[:] = ContractFromAximuthalAverage(dof_coordinatesP1,100,temp_vecP1[:,kk])
	print('P1 interpolation: it took: '+str(time.time()-p1time))
		
	p2time = time.time()	
	if linear == True and p2CNT>0:
		print('started griddata for P2, stacked, linear')
		temp_vecP2 = interpolate.griddata(points, valsP2,dof_expandedP2, method='linear')
	elif linear == False and p2CNT>0:
		print('started griddata for P2, stacked, nearest')
		temp_vecP2 = interpolate.griddata(points, valsP2,dof_expandedP2, method='nearest')
	
	print('time for interpolation command: ' +str(time.time()-p2time))
	for k,kk in zip(namesP2,range(0,len(namesP2))):
		MeanFlowDict[k].vector()[:] = ContractFromAximuthalAverage(dof_coordinatesP2,100,temp_vecP2[:,kk])
	print('P2 interpolation: it took: '+str(time.time()-p2time))

	return MeanFlowDict, notInFileList

def exportBaseFlowAsHDF5(param,MeanFlowDict,FEMSpaces):
	mesh = MeanFlowDict[list(MeanFlowDict.keys())[0]].function_space().mesh()
	hdf5file = HDF5File(mesh.mpi_comm(),param.FlowInput.MeanFlowFilePath[0:-4]+"_"+param.MeshPath[0:-3].split('/')[-1]+"hdf5", 'w')
	for name in MeanFlowDict.keys():
		if name == 'u':
			for comp, cnt in zip(param.VelocityComponents,range(0,len(param.VelocityComponents))):
				print(comp)
				dofIDX = FEMSpaces.FunctionSpaceVectorVelocity.sub(cnt).dofmap().dofs()
				uCompDict = {}
				uCompDict['u'+comp] = Function(FEMSpaces.P2)
				uCompDict['u'+comp].vector()[:] = MeanFlowDict['u'].sub(cnt).vector()[dofIDX]
				hdf5file.write(uCompDict['u'+comp], name+comp)
		else:
			hdf5file.write(MeanFlowDict[name], name)
 
 
def CoordinateTransformation(transformType,RawFlowDict, param):
	# This function performs coordinate transformations and in doing so creates InputFlowDicts from RawFlowDicts
	# For all cases first the basis coordinate system is shifted to the axis of the Felics mesh (which is by defeult the x-axis)
	# For the moment only Cart2Cyl is implemented and only calculates the cylindrical velocity components from the carthesian ones
	
	#First calculate the rotation matrix R based on the two vectors A and B
	#A (the vector of origin) is given by the user and corresponds to the axis 
	#of rotational symmetry or to the symmetry plane in 2D carthesian meshes, which must intersect the origin, currently
	#B is the vector of destination and is currently by default the x axis (The 
	#axis or plane of symmetry for the linear analysis is always the x axis)
    
	# Edit:
	# User is now required to set the axis mentioned aboved within the gui settings
	# A dropdown menu is provided were the user can choose from x,y,z for both, 
	# the destination(B) and the origin axis(A)
    
	# if origin axis is y or z the rotation matrix R is formed according to choice
	A=np.zeros(3)
	if 	 param.AveragingFromAxis == 'x': A[0]=1 
	elif param.AveragingFromAxis == 'y': A[1]=1
	elif param.AveragingFromAxis == 'z': A[2]=1
	B=np.zeros(3)
	B[0]=1
		   
	if np.all(A-B ==0):
		print("Origin vector and target vector are identical. Skipping rotation...")
		[RawFlowDict['x'],RawFlowDict['y'],RawFlowDict['z']]= [RawFlowDict['x'],RawFlowDict['y'],RawFlowDict['z']]
	else:
		print("Rotating the coordinate system to align with the Felics mesh...")
		v=np.cross(A,B)
		s=np.linalg.norm(v)
		c=np.dot(A,B)
		v_matrix=np.zeros((3,3))
		v_matrix[1,0]=v[2]
		v_matrix[2,0]=-v[1]
		v_matrix[2,1]=v[0]
		v_matrix[0,1]=-v[2]
		v_matrix[0,2]=v[1]
		v_matrix[1,2]=-v[0]
		print("Rotating base flow mesh to felics mesh...")
		print("Rotation from ("+str(A[0])+","+str(A[1])+","+str(A[2])+") to ("+str(B[0])+","+str(B[1])+","+str(B[2])+")")
		R=np.identity(3)+v_matrix+np.dot(v_matrix,v_matrix)*(1-c)/s**2
		# Multiply R on the raw coordinates (x,y,z) to obtain the final coordinates 
		[RawFlowDict['x'],RawFlowDict['y'],RawFlowDict['z']]= list(np.dot(R,np.array([RawFlowDict['x'],RawFlowDict['y'],RawFlowDict['z']])))
	print("Performing Coordinate Transform from cartesian to cylindrical coordinates...")
	# Rotate velocities (ux,uy,uz) to the new coordinate system via rotation matrix R so they become 
	if np.all(A-B ==0):
		[RawFlowDict['ux'],RawFlowDict['uy'],RawFlowDict['uz']]= [RawFlowDict['ux'],RawFlowDict['uy'],RawFlowDict['uz']]
	else:
		[RawFlowDict['ux'],RawFlowDict['uy'],RawFlowDict['uz']] = list(np.dot(R,np.array([RawFlowDict['ux'],RawFlowDict['uy'],RawFlowDict['uz']])))
		
	if transformType == 'Cart2Cyl':

		# Calculate the radial and azimuthal coordinate
		RawFlowDict['r']=(RawFlowDict['z']**2+RawFlowDict['y']**2)**0.5
		#RawFlowDict['x']=RawFlowDict['x']
		RawFlowDict['theta']=np.arctan2(RawFlowDict['z'],RawFlowDict['y'])
		# Obtain radial and tangential velocity (cylindrical coordinates) from the angle theta and the velocities ux and uy  (carthesian coordinates)
		RawFlowDict['ur']=np.cos(RawFlowDict['theta'])* RawFlowDict['uy'] +np.sin(RawFlowDict['theta'])*RawFlowDict['uz']
		RawFlowDict['ut']=-np.sin(RawFlowDict['theta'])* RawFlowDict['uy'] +np.cos(RawFlowDict['theta'])*RawFlowDict['uz']
#		axisIDX = np.argwhere(RawFlowDict['r']<1e-10)
#		RawFlowDict['ur'][axisIDX] = 0
#		RawFlowDict['ut'][axisIDX] = 0
		
		
		
		if 'rstxx' in list(RawFlowDict.keys()):
			#RawFlowDict['rstxx']= RawFlowDict['rstxx_raw']
			RawFlowDict['rstrr']= np.cos(RawFlowDict['theta'])**2*RawFlowDict['rstyy'] + np.sin(RawFlowDict['theta'])**2* RawFlowDict['rstzz'] \
								+2*np.sin(RawFlowDict['theta'])*np.cos(RawFlowDict['theta'])* RawFlowDict['rstyz']
			RawFlowDict['rsttt']= np.cos(RawFlowDict['theta'])**2*RawFlowDict['rstzz'] + np.sin(RawFlowDict['theta'])**2* RawFlowDict['rstyy'] \
								-2*np.sin(RawFlowDict['theta'])*np.cos(RawFlowDict['theta'])* RawFlowDict['rstyz']
								
			RawFlowDict['rstxr']= np.cos(RawFlowDict['theta'])*RawFlowDict['rstxy'] + np.sin(RawFlowDict['theta'])* RawFlowDict['rstxz'] 
			RawFlowDict['rstxt']= -np.sin(RawFlowDict['theta'])* RawFlowDict['rstxy'] + np.cos(RawFlowDict['theta'])*RawFlowDict['rstxz']
			RawFlowDict['rstrt']= np.cos(RawFlowDict['theta'])**2*RawFlowDict['rstyz'] - np.sin(RawFlowDict['theta'])**2* RawFlowDict['rstyz'] \
								+np.sin(RawFlowDict['theta'])*np.cos(RawFlowDict['theta'])* RawFlowDict['rstzz']\
								-np.sin(RawFlowDict['theta'])*np.cos(RawFlowDict['theta'])* RawFlowDict['rstyy']
								
			del  RawFlowDict['rstyy'], RawFlowDict['rstzz'],\
			RawFlowDict['rstxy'], RawFlowDict['rstxz'], RawFlowDict['rstyz']
	
	return RawFlowDict 	

def ExpandForAzimuthalAverage(coordinates,n_cuts):
    ## This function extends the fenics grid to three dimensions to prepare for interpolation in 3D grids
    size=np.shape(coordinates)[0]
    b=np.zeros((size,3))
    b[:,0]=coordinates[:,0]
    b[:,1]=coordinates[:,1]
    coordinates=b
     
    temp_coordinates=np.zeros(np.shape(coordinates))
    coordinates_out=np.zeros((0,3))
    for angle in np.linspace(2*np.pi/n_cuts,2*np.pi,n_cuts):
        temp_coordinates[:,0]=coordinates[:,0]
        temp_coordinates[:,1]=coordinates[:,1]*np.cos(angle)
        temp_coordinates[:,2]=coordinates[:,1]*np.sin(angle)
        coordinates_out=np.concatenate((coordinates_out,temp_coordinates))
    return coordinates_out
 
def ContractFromAximuthalAverage(coordinates,n_cuts,vec_in):
    ## This function contracts the expanded interpolated data from the expaned grid (see function ExpandForZaimuthalAverage) back
    ## to the 2D base grid
    length=np.shape(coordinates)[0]
    vec_out=np.zeros(length)
    for i in range(n_cuts):
        vec_out=vec_out+vec_in[i*length:i*length+length]
    vec_out=vec_out/n_cuts  
    return vec_out
     
def assignNames(h5file):
    DictNames = {}
    for k in list(h5file['Base'].keys()):
        if 'omain' in k:
            DictNames['domainN'] = k
            for kk in list(h5file['Base'][DictNames['DomainN']]['GridCoordinates'].keys()):
                if 'x' in kk or 'X' in kk:
                     DictNames['x']
                elif 'y' in kk or 'Y' in kk:
                      DictNames['y']
                elif 'z' in kk or 'Z' in kk:
                     DictNames['z']
                      
            for kk in list(h5file['Base'][DictNames['DomainN']]['FlowSolution'].keys()):
                if 'elocity' in kk:
                    if 'X' in kk or 'x' in kk:
                        DictNames['ux'] = kk
                        DictNames['ux_raw'] = 'ux_raw'
                         
                    elif 'Y' in kk or 'y' in kk:
                        DictNames['uy'] = kk
                        DictNames['uy_raw'] = 'uy_raw'
                         
                    elif 'Z' in kk or 'z' in kk:
                        DictNames['uz'] = kk
                        DictNames['uz_raw'] = 'uz_raw'
                         
                    else:
                        DictNames[kk] = kk
          
             
     
     
    return DictNames

def CorrectMeanAtBoundaries(MeanFlowDict):
	
	# Get dimension of mesh
	gdim = MeanFlowDict[list(MeanFlowDict.keys())[0]].function_space().mesh().geometry().dim()
	if 'ut' in list(MeanFlowDict.keys()):
		# get coordinates of nodes
		dofs_coord=MeanFlowDict['ut'].function_space().tabulate_dof_coordinates().reshape((-1,gdim))
		# Get indices of coordinates at axis
		indices=[i for i, val in enumerate(dofs_coord[:,1]==0.00) if val]
		# Set azimuthal velocity of these nodes to zero
		MeanFlowDict['ut'].vector()[indices] = 0
	if 'u' in list(MeanFlowDict.keys()):
		#get the indices of ur in the velocity vector
		dofs_ur=MeanFlowDict['u'].function_space().sub(1).dofmap().dofs()
		# Split velocity into components
		ux, ur = MeanFlowDict['u'].split(deepcopy=True)

		# Get the coordinates of the sub function
		dofs_coord=ur.function_space().tabulate_dof_coordinates().reshape((-1,gdim))
		# Get the indexes at the axis
		indices=[i for i, val in enumerate(dofs_coord[:,1]==0.00) if val]
		# set the component to 0 at axis
		ur.vector()[indices] = 0
		#Copy the cmponent back to the vector
		MeanFlowDict['u'].vector()[dofs_ur]=ur.vector()[:]
	return MeanFlowDict



def ExpandForAverage(coordinates,param):
	''' This function extends the fenics grid to three dimensions to prepare for interpolation in 3D grids. Afet interpolation use the function ContractAfterAverage to project the 3D data again on the 2D FELiCS mesh.
	    \t Input:
	    \t\t -coordinates: Coordinates of the FELiCS mesh
	    \t\t -param: FELiCS parameter object
	    \t Ouput:
	    \t\t -coordinates_out: FELiCS mesh coordinates expanded to 3 dimensions
	    '''
	n_cuts=100
	size=np.shape(coordinates)[0]
	b=np.zeros((size,3))
	b[:,0]=coordinates[:,0]
	b[:,1]=coordinates[:,1]
	coordinates=b
	 
	temp_coordinates=np.zeros(np.shape(coordinates))
	coordinates_out=np.zeros((0,3))
	if param.CoordinateSystem=='Cylindrical'
		for angle in np.linspace(2*np.pi/n_cuts,2*np.pi,n_cuts):
			temp_coordinates[:,0]=coordinates[:,0]
			temp_coordinates[:,1]=coordinates[:,1]*np.cos(angle)
			temp_coordinates[:,2]=coordinates[:,1]*np.sin(angle)
			coordinates_out=np.concatenate((coordinates_out,temp_coordinates))
	elif param.CoordinateSystem=='Carthesian'
		zmin=0.01
		zmax=0.01
		for angle in np.linspace(zmin,zmax,n_cuts):
			temp_coordinates[:,0]=coordinates[:,0]
			temp_coordinates[:,1]=coordinates[:,1]
			temp_coordinates[:,2]=coordinates[:,1]*0+zmin
			coordinates_out=np.concatenate((coordinates_out,temp_coordinates))
	return coordinates_out
 
def ContractAfterAverage(coordinates,n_cuts,vec_in):
    ## This function contracts the expanded interpolated data from the expaned grid (see function ExpandForZaimuthalAverage) back
    ## to the 2D base grid
    length=np.shape(coordinates)[0]
    vec_out=np.zeros(length)
    for i in range(n_cuts):
        vec_out=vec_out+vec_in[i*length:i*length+length]
    vec_out=vec_out/n_cuts  
    return vec_out
