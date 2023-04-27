from ufl import (
				i,
				j,
				dx,
				Dx,
				inner,
				grad,
				)
def MomentumEqDiffusion(self,fluc,X,mean,param):
	# term I: mu Laplace u (cartesian portion)
	# there are parts which require IbP --> spatially variable viscosity is taken into account
	self.A_imag_vf += -self.R*Dx(fluc['u'][i],j)*Dx(X[i]*mean['MuTot'],j)*dx
	# term I (cylindrical portion)
	if param.Case.CoordinateSystem in ['Cylindrical']:
		self.A_imag_vf += -1*mean['MuTot']*X[1]/self.R*fluc['u'][1]*dx
		if not param.Case.m == 0:
			self.A_imag_vf += -1*mean['MuTot']*self.m**2/self.R*inner( X,fluc['u'] )*dx
			self.A_real_vf += 2*mean['MuTot']*self.m/self.R*X[1]*fluc['ut']*dx

	# IbP of term I yields boudnary integrals
	self.A_imag_vf += 1*mean['MuTot']*X[0]*inner( grad(fluc['u'][0]) , self.n_BC )*self.all_ds
	self.A_imag_vf += 1*mean['MuTot']*X[1]*inner( grad(fluc['u'][1]) , self.n_BC )*self.all_ds

	# term II: (grad u + grad^T u)*(grad mu) - does not require any IbP
	self.A_imag_vf += self.R*X[0]*mean['MuTot'].dx(0)*2*fluc['u'][0].dx(0)*dx
	self.A_imag_vf += self.R*X[0]*mean['MuTot'].dx(1)*(fluc['u'][0].dx(1)+fluc['u'][1].dx(0))*dx

	self.A_imag_vf += self.R*X[1]*mean['MuTot'].dx(0)*(fluc['u'][0].dx(1)+fluc['u'][1].dx(0))*dx
	self.A_imag_vf += self.R*X[1]*mean['MuTot'].dx(1)*2*fluc['u'][1].dx(1)*dx

	# (cartesian portion) (there might be a more sophisticated implementation, suggestion below, but not workin yet)
	#self.A_imag_vf += self.R*( Dx(fluc['u'][i],j)+Dx(fluc['u'][j],i))*Dx(mean['MuTot'],i)*X[i]*dx
	if param.Case.AnalysisMode in ['Input-Output']:
		for boundary_index in param.IOResolvent.ForcingBoundaryIndices:
			# Add physical boundary terms in imaginary part, where the forcing is applied...
			self.A_imag_vf += -self.R*mean['MuTot']*(-X[0] * (self.n_BC[0] * (fluc['u'][0].dx(0)) +
                                                               		self.n_BC[1] *   (fluc['u'][0].dx(1))) -
                                               		X[1] * (self.n_BC[0] *   (fluc['u'][1].dx(0)) +
                                                               		self.n_BC[1] *   (fluc['u'][1].dx(1))) ) *self.ds(boundary_index)
			# Add stabilization terms on imaginary part according to Baumann and Oden JFM 2016 vol 798
			self.A_imag_vf += -self.R*mean['MuTot']*((self.n_BC[0] * (X[0].dx(0)) +
                           				self.n_BC[1] *   (X[0].dx(1))) * (fluc['u'][0]-mean['u_forcing_i'][0]) +#jvs
                          			       (self.n_BC[0] *  (X[1].dx(0)) +
                           				self.n_BC[1] *   (X[1].dx(1))) * (fluc['u'][1]-mean['u_forcing_i'][1])) *self.ds(boundary_index) #jvs

			# Add stabilization terms on imaginary part according to Baumann and Oden JFM 2016 vol 798
			self.A_real_vf += -self.R*mean['MuTot']*((self.n_BC[0] * (X[0].dx(0)) +
                           				self.n_BC[1] * (X[0].dx(1) )) * (-mean['u_forcing_r'][0]) +
                          			       (self.n_BC[0] * (X[1].dx(0) ) +
                           				self.n_BC[1] * (X[1].dx(1))) * (-mean['u_forcing_r'][1])) *self.ds(boundary_index)


	if param.NumericalScheme in ['Discontinuous Galerkin']:
		                #F_vis_Int = ((self.n_BC[0]('+') * avg(4 / 3 * X[0].dx(0) - 2 / 3 * X[1].dx(1)) +
				#	self.n_BC[1]('+') * avg(X[0].dx(1) + X[1].dx(0))) * jump(fluc['u'][0]) +
				#	(self.n_BC[0]('+') * avg(X[0].dx(1) + X[1].dx(0)) +
				#	 self.n_BC[1]('+') * avg(4 / 3 * X[1].dx(1) - 2 / 3 * X[0].dx(0))) * jump(fluc['u'][1]) -
				#	jump(X[0]) * (n[0]('+') * avg(4 / 3 * fluc['u'][0].dx(0) - 2 / 3 * fluc['u'][1].dx(1)) +
				#	                n[1]('+') * avg(fluc['u'][0].dx(1) + fluc['u'][1].dx(0))) -
				#	jump(Xiu[1]) * (n[0]('+') * avg(fluc['u'][0].dx(1) + fluc['u'][1].dx(0)) +
				#	                n[1]('+') * avg(4 / 3 * fluc['u'][1].dx(1) - 2 / 3 * fluc['u'][0].dx(0))))
				#F_vis_Int *= avg(1*mean['MuTot'])
				F_vis_Int_own = ((   self.n_BC[0]('+') * avg(X[0].dx(0)) +
						     self.n_BC[1]('+') * avg(X[0].dx(1) )) * jump(fluc['u'][0]) +
						    (self.n_BC[0]('+') * avg(X[1].dx(0) ) +
						     self.n_BC[1]('+') * avg(X[1].dx(1))) * jump(fluc['u'][1]) -
						jump(X[0]) * (self.n_BC[0]('+') * avg(fluc['u'][0].dx(0)) +
						                      self.n_BC[1]('+') * avg(fluc['u'][0].dx(1))) -
						jump(X[1]) * (self.n_BC[0]('+') * avg(fluc['u'][1].dx(0)) +
						                      self.n_BC[1]('+') * avg(fluc['u'][1].dx(1))))
				F_vis_Int_own *= avg(1*mean['MuTot'])
				self.A_imag_vf += -F_vis_Int_own *dS
