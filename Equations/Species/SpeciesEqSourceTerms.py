from ufl import div,inner,dx
def addSpeciesEqSourceTerms(self,X,mean,species,param):
	if param.Case.AnalysisMode == 'Input-Output' and param.IOResolvent.ForcingMode == 'Body':
		self.A_vf.add(self.R*mean.forcing_r(species)*X*dx)
		self.A_vf.add(1j * self.R*mean.forcing_i(species)*X*dx)
