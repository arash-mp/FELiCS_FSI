from ufl import (
				inner,
				div,
				grad,
				jump,
				avg,
				dx,
				dS,
				dot,
				i,
				j,
				Identity,
				conj,
				)
def MassEqAdvection(self,fluc,X,mean,param):


	I=Identity( fluc.u.geometric_dimension() )
	self.A_vf.add(1j * conj(X).dx(i)*self.R*fluc.rhou[j]*I[i,j]*dx)
	self.A_vf.add(1j * conj(X).dx(0)*self.R*fluc.rhou[0]*dx)
	self.A_vf.add(1j * conj(X).dx(1)*self.R*fluc.rhou[1]*dx)

	#Integrate boundary term due to integration by parts
	self.A_vf.add(1j * -self.R*fluc.rhou[j]*self.n_BC[i]*conj(X)*I[i,j]*self.all_ds)


	if param.Case.CoordinateSystem in ['Cylindrical']:
		if not param.Case.m == 0:
			self.A_vf.add(conj(X)*mean.rho*self.m*fluc.u[2]*dx)
			self.A_vf.add(conj(X)*fluc.rho*self.m*mean.u[2]*dx)

	# In case of Input-Output analysis, we must adapt the boundary terms...
	if param.Case.AnalysisMode in ['Input-Output']:
		# Iterate through all boundaries, at which forcing is applied
		for boundary_index in param.IOResolvent.ForcingBoundaryIndices:
			# First subtract the part added in a few lines above...
			self.A_vf.add(1j * self.R*inner(fluc.u,self.n_BC)*conj(X)*self.ds(boundary_index))
			# Then add the forcing at the inlet...
			self.A_vf.add(-self.R*inner(mean.u_forcing_r,self.n_BC)*conj(X)*self.ds(boundary_index))
			self.A_vf.add(1j * -self.R*inner(mean.u_forcing_i,self.n_BC)*conj(X)*self.ds(boundary_index))

	if param.NumericalScheme in ['Discontinuous Galerkin']:
		F_rho_LFF = dot(avg(fluc.u),self.n_BC('+'))
		self.A_vf.add(1j * -inner(jump(X),F_rho_LFF)*dS) # possibly conj(X) instead of X
