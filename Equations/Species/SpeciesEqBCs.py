from ufl import (
	inner,
	Identity,
	i,
	j,
)
def addSpeciesEqBCs(self,fluc,X,mean,species,param):

	I=Identity( fluc.u.geometric_dimension() )
	for Boundary in param.BCs.getBCsDict()[species]:
		# If forcing is applied at the boundary, and Input-Output mode is on...
		if (param.Case.AnalysisMode in ['Input-Output']) and (Boundary['ID'] in param.IOResolvent.ForcingBoundaryIndices) :
			# First subtract the part added in a few lines above...
			self.A_imag_vf.add(self.n_BC[j]*mean.u[i]*self.R*fluc.Y(species)*X*I[i,j]*self.ds(Boundary['ID']))  ### tlk: Why is therer no density in the equation???
			# Then add the forcing of the respective species given in the mean flow dict at the respective bounary
			self.A_imag_vf.add( -self.n_BC[j]*mean.u[i]*self.R*mean.forcing_i(species)*X*I[i,j]*self.ds(Boundary['ID']))
			self.A_real_vf.add( -self.n_BC[j]*mean.u[i]*self.R*mean.forcing_r(species)*X*I[i,j]*self.ds(Boundary['ID']))
		else:
			## Check if Boundary condition is Dirichlet or Neumann
			if Boundary['type'] in ['Dirichlet']:
				# If it is Dirichlet, the BC is applied in the weak formulation
				# First subtract the part added in a few lines above...
				self.A_imag_vf.add(+inner(self.n_BC,mean.u)*self.R*fluc.Y(species)*X*self.ds(Boundary['ID']))
				# Then add the respective Boundary condition
				# If BC value is 0, then the weak formulation throws an error, therefore check if it is zero...
				# ... if the value is zero, a treatment is not necessary anyway
				if not Boundary['value'] in [0.0]:
					self.A_real_vf.add(-inner(self.n_BC,mean.u)*self.R*Boundary['value']*X*self.ds(Boundary['ID']))
			if Boundary['type'] in ['Neumann']:
				if not Boundary['value'] in [0.0]:
					printError('So far only homogeneous Neumann conditions are implemented... Please either change to another BC or - even better -  implement it yourself and upload your well documented implementation to gitlab...')
