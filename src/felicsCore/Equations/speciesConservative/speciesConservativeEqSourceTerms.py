from ufl import div,inner,dx
def addSpeciesConservativeEqSourceTerms(self,X,mean,species,param):
	if param.Case.AnalysisMode == 'Input-Output' and param.IOResolvent.ForcingMode == 'Body':
		self.A_real_vf.add(self.R*mean.forcing_r(species)*X*dx)
		self.A_imag_vf.add(self.R*mean.forcing_i(species)*X*dx)
