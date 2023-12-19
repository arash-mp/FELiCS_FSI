from ufl import (
				i,
				j,
				Dx,
				dx,
				inner,
				)


def MomentumEqDiffusion(self,fluc,X,mean,param):
	# Integrate dissipative terms in domain (integration by parts is applied)
	self.A_imag_vf += -self.R*mean['MuTot'] * Dx(fluc['u'][i],j)*Dx(X[i],j)*dx
	if param.Case.AnalysisMode in ['Input-Output']:
		for boundary_index in param.IOResolvent.ForcingBoundaryIndices:
			# Add physical boundary terms in imaginary part, where the forcing is applied...
			self.A_imag_vf += -self.R*mean['MuTot']*(-X[0] * (self.n_BC[0] * (fluc['u'][0].dx(0)) +
                                                               		self.n_BC[1] *   (fluc['u'][0].dx(1))) -
                                               		X[1] * (self.n_BC[0] *   (fluc['u'][1].dx(0)) +
                                                               		self.n_BC[1] *   (fluc['u'][1].dx(1))) ) *self.ds(boundary_index)
			# Add stabilization terms on imaginary part according to Baumann and Oden JFM 2016 vol 798
			self.A_imag_vf += -self.R*mean['MuTot']*((self.n_BC[0] * (X[0].dx(0)) +
                           				self.n_BC[1] *   (X[0].dx(1))) * (fluc['u'][0]-mean['u'+param.VelocityComponents[0]+'_forcing_i']) +
                          			       (self.n_BC[0] *  (X[1].dx(0)) +
                           				self.n_BC[1] *   (X[1].dx(1))) * (fluc['u'][1]-mean['u'+param.VelocityComponents[1]+'_forcing_i'])) *self.ds(boundary_index)

			# Add stabilization terms on imaginary part according to Baumann and Oden JFM 2016 vol 798
			self.A_real_vf += -self.R*mean['MuTot']*((self.n_BC[0] * (X[0].dx(0)) +
                           				self.n_BC[1] * (X[0].dx(1) )) * (-mean['u'+param.VelocityComponents[0]+'_forcing_r']) +
                          			       (self.n_BC[0] * (X[1].dx(0) ) +
                           				self.n_BC[1] * (X[1].dx(1))) * (-mean['u'+param.VelocityComponents[1]+'_forcing_r'])) *self.ds(boundary_index)

	if param.Case.CoordinateSystem in ['Cylindrical']:
		# find index of third velocity component:
		ThirdVelCompIndex=param.SolutionList.index('ut')
		# Add m**2 terms in x-r momentum equations
		if not param.Case.m == 0:
			self.A_imag_vf += -1.0*mean['MuTot']/self.R*inner(X, fluc['u'] )*self.m**2*dx
		# Add remaining terms in radial momentum equation
		if not param.Case.m == 0:
			self.A_real_vf += 2.0*mean['MuTot']*X[1]/self.R*self.m*fluc['ut']*dx
		self.A_imag_vf += -1.0*mean['MuTot']*X[1]/self.R*fluc['u'][1]*dx

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
				F_vis_Int_own *= avg(1.0*mean['MuTot'])
				self.A_imag_vf += -F_vis_Int_own *dS



