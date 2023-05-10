from ufl import inner
def ThirdMomentumEqBCs(self,fluc,X,mean,param):
	for Boundary in param.BCs.getBCsDict()['ut']:
		if (param.Case.AnalysisMode in ['Input-Output']):
			# If forcing is applied at the boundary, and Input-Output mode is on...
			if (Boundary['ID'] in param.IOResolvent.ForcingBoundaryIndices):
				# First subtract the part added in a few lines above...
				self.A_imag_vf.add(+fluc[self.ThirdVelCompIndex]*inner(mean['u'],self.n_BC)*X*self.R*self.ds(Boundary['ID']))
				# Then add the forcing of the respective species given in the mean flow dict at the respective bounary
				self.A_imag_vf.add(-mean['ut_forcing_i']*inner(mean['u'],self.n_BC)*X*self.R*self.ds(Boundary['ID']))
				self.A_real_vf.add(-mean['ut_forcing_r']*inner(mean['u'],self.n_BC)*X*self.R*self.ds(Boundary['ID']))
