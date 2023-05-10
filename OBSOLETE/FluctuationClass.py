from fenics import TrialFunctions, Constant, Function
class FluctuationClass():
	def __init__(self,param,mean,FEMSpaces,postProcessing=False,solution=None):
		self.__param=param
		self.__ZeroField=Function(FEMSpaces.P2)
		self.__flucDict={}
		self.__FEMSpaces=FEMSpaces
		self.__SolutionList=param.Case.getTransportedQuantityList()
		transportedQuantities = self.__param.Case.getTransportedQuantityList() 
		if not postProcessing:
			self.__fluc = TrialFunctions(FEMSpaces.VMixed)
			for field in transportedQuantities:
				self.__flucDict[field] = self.__fluc[self.__SolutionList.index(field)]
		else:
			self.__fluc = Function(FEMSpaces.VMixed)
			self.__fluc.vector()[:] = solution.astype(float)
			
			for field in transportedQuantities:
				if len(transportedQuantities)>1:
					temp = self.__fluc.split()[self.__SolutionList.index(field)]
				else:
					temp = self.__fluc
				self.__flucDict[field] = temp
		
		if (self.__param.Case.HeatTransfer or self.__param.Case.Compressible or self.__param.Case.Reaction) and not self.__param.Case.Mixture.getReactionMechanism()['type'] in ['NOx']:
			self.init_T(mean)
			self.init_h(mean)
		self.initLamDiff(mean,param)
		if self.__param.Case.Reaction and self.__param.Case.Mixture.getReactionMechanism()['type']=='2S-SM2':
			self.__calculateHeatRelease(mean.reaction)
	
	def init_T(self,mean):
		from fenics import project,Expression
		if not self.__param.Case.Mixture.getReactionMechanism()['type'] in ['2S-SM2']:
			if mean.molarMass.vector()[0]==0:
				M=28.949
				print('using constant molarMass')
			else:
				M=mean.molarMass*1e3
				print('using calculated molarMass field')
			
			p=101300
			R=8314.4598/M
			Tm=p/R/mean.rho
			
			flucT=-self.rho/mean.rho*mean.T # CA: Using mean.T instead of above calculated Tm
			
		else:
			flucT = -self.rho/mean.rho*mean.T
			#self.__T=flucT
			if self.__param.Case.Compressible:
				flucT += self.p / mean.p * mean.T
			
		self.__flucDict['T']=flucT

		if self.__param.Case.Mixture.ReactionMechanism['type']=='TwoStep':
			print('Powerlaw law for visco perturbations')
			power_expo=0.6759
			self.__flucDict['viscoLaw']=-power_expo*self.rho/mean.rho #visco fluctuation factor(PowerLaw)
		#elif self.__param.Case.MolViscPerturbModel: #for Max heat conduction and 1 Step case #self.__param.Case.Mixture.getReactionMechanism()['type']=='Global':
		#	
		#	print('Sutherland law for visco perturbations')
		##elif self.__param.Case.MolViscPerturbModel: #for Max heat conduction and 1 Step case #self.__param.Case.Mixture.getReactionMechanism()['type']=='Global':
		#	print('Sutherland law for visco perturbations')
		#	Ts=170.672
		#	self.__flucDict['viscoLaw']=(mean.T+3*Ts)/(2*(mean.T+Ts))*(-self.rho/mean.rho) #visco fluctuation factor(Sutherland)
	
	def initLamDiff(self,mean,param):
		if (param.Case.MolViscPerturbModel == 'Sutherland') and ('T' in list(self.__flucDict.keys())): 
			As       =       param.Case.Mixture.Viscosity['Constants']['As']
			Ts       =       param.Case.Mixture.Viscosity['Constants']['Ts']
			self.__flucDict['nulam'] = As*self.T*(0.5*mean.T**(3/2)+1.5*Ts*mean.T**(0.5))/((mean.T+Ts)**2)
			Pr=param.Case.Mixture.Pr

			self.__flucDict['alpha'] = self.__flucDict['nulam']/Pr 
			for specie in param.Case.Mixture.getSpeciesList('transported'):
				Sc=param.Case.Mixture.species[specie]['Sc']
				self.__flucDict['D_'+specie]=self.__flucDict['nulam']/Sc
		elif (param.Case.MolViscPerturbModel == 'Sutherland mean'):
			Ts       =       param.Case.Mixture.Viscosity['Constants']['Ts']
			fluct=(mean.T+3*Ts)/(2*(mean.T+Ts))*(-self.rho/mean.rho)
			self.__flucDict['nulam']=mean.nulam*fluct
			self.__flucDict['alpha']=mean.alpha*fluct
			for specie in self.__param.Case.Mixture.getSpeciesList('transported'):
				self.__flucDict['D_'+specie]=mean.D(specie)*fluct
		else:
			self.__flucDict['nulam']=self.__ZeroField	
			self.__flucDict['alpha']=self.__ZeroField
			for specie in self.__param.Case.Mixture.getSpeciesList('transported'):
				self.__flucDict['D_'+specie]=self.__ZeroField

	def init_h(self,mean):
		if not self.__param.Case.Mixture.ReactionMechanism['type']=='NOx':
			self.__flucDict['h']=mean.cp*self.__flucDict['T']


	def __calculateHeatRelease(self,reaction):
		Qtot,QList=reaction.dQ(self)
		self.__flucDict['Q']=Qtot
		for Q in QList:
			self.__flucDict['Q'+str(QList.index(Q))]=Q

	def Y(self,specie):
		return self.__fluc[self.__SolutionList.index(specie)]

	def D(self,specie):
		return self.__flucDict['D_'+specie]

	@property
	def u(self):
		if 'u' in self.__param.Case.getTransportedQuantityList():
			return self.__flucDict['u']
		else:
			u=Function(self.__FEMSpaces.FunctionSpaceVectorVelocity)
			return u


	@property
	def ut(self):
		if 'ut' in self.__param.Case.getTransportedQuantityList():
			return self.__flucDict['ut']
		else:
			return Constant(0)

	@property
	def p(self):
		if 'p' in self.__param.Case.getTransportedQuantityList():
			return self.__flucDict['p']
		else:
			return Constant(0)

	@property
	def rho(self):
		if 'rho' in self.__param.Case.getTransportedQuantityList():
			return self.__flucDict['rho']
		else:
			return self.__ZeroField

	@property
	def T(self):
		return self.__flucDict['T']

	@property
	def Tm(self):
		return self.__flucDict['Tm']

	@property
	def h(self):
		return self.__flucDict['h']

	@property
	def Q(self):
		return self.__flucDict['Q']

# 	@property
# 	def alpha(self):
# 		return self.__alpha

	@property
	def viscoLaw(self):
		if 'viscoLaw' in list(self.__flucDict.keys()):
			return self.__flucDict['viscoLaw']
		else:
			return self.__ZeroField
	
	@property
	def nulam(self):
		return self.__flucDict['nulam']

	@property
	def alpha(self):
		return self.__flucDict['alpha']

	@property
	def fluc(self):
		return self.__fluc

	@property
	def FlucDict(self):
		from fenics import project,Function
		OutputDict={}
		for key in list(self.__flucDict.keys()):
			if not isinstance(self.__flucDict[key],Function):
				if key == 'u':
					FEMSpace=self.__FEMSpaces.FunctionSpaceVectorVelocity
				else:
					FEMSpace=self.__FEMSpaces.P2
				OutputDict[key] = project(self.__flucDict[key], FEMSpace)		
			else:
				OutputDict[key]=self.__flucDict[key]
		return OutputDict

