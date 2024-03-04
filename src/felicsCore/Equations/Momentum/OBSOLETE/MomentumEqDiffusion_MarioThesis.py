from ufl import (
				i,
				j,
				dx,
				Dx,
				inner,
				grad,
				)

def MomentumEqDiffusion(self,fluc,X,mean,param):
		# ************************************************************
		# The odl derivation neglected the transposed gradient tensor within the dissipative terms
		# The change within the linear equations also affects the dilatation terms
		# within Mario-Thesis this was changed, but it is not sophisticated
		self.A_imag_vf += -1/mean['Re']*self.R*(
											  2*X[0].dx(0)*fluc['u'][0].dx(0)
											+ X[0].dx(1)*(fluc['u'][1].dx(0)+fluc['u'][0].dx(1))
											)*dx
		if not param.Case.m == 0:
			self.A_imag_vf += -1/mean['Re']*X[0]*self.m**2/self.R*fluc['u'][0]*dx
			self.A_real_vf += 1/mean['Re']*self.m*X[0]*fluc['ut'].dx(0)*dx



		# r-component
		self.A_imag_vf += -1/mean['Re']*self.R*(
											  X[1].dx(0)*(fluc['u'][1].dx(0)+fluc['u'][0].dx(1))
											+ 2*X[1].dx(1)*fluc['u'][1].dx(1)
											)*dx
		if not param.Case.m == 0:
			self.A_imag_vf += -1/mean['Re']*X[1]*self.m**2/self.R*fluc['u'][1]*dx
			self.A_real_vf += 1/mean['Re']*(
										self.m*X[1]*fluc['ut'].dx(1)
										-self.m*X[1]*fluc['ut']/self.R
										)*dx

#

		if param.NumericalScheme in ['Discontinuous Galerkin']:
			                #F_vis_Int = ((self.n_BC[0]('+') * avg(4 / 3 * X[0].dx(0) - 2 / 3 * X[1].dx(1)) +
					#	self.n_BC[1]('+') * avg(X[0].dx(1) + X[1].dx(0))) * jump(fluc['u'][0]) +
					#	(self.n_BC[0]('+') * avg(X[0].dx(1) + X[1].dx(0)) +
					#	 self.n_BC[1]('+') * avg(4 / 3 * X[1].dx(1) - 2 / 3 * X[0].dx(0))) * jump(fluc['u'][1]) -
					#	jump(X[0]) * (n[0]('+') * avg(4 / 3 * fluc['u'][0].dx(0) - 2 / 3 * fluc['u'][1].dx(1)) +
					#	                n[1]('+') * avg(fluc['u'][0].dx(1) + fluc['u'][1].dx(0))) -
					#	jump(Xiu[1]) * (n[0]('+') * avg(fluc['u'][0].dx(1) + fluc['u'][1].dx(0)) +
					#	                n[1]('+') * avg(4 / 3 * fluc['u'][1].dx(1) - 2 / 3 * fluc['u'][0].dx(0))))
					#F_vis_Int *= avg(1/mean['Re'])
					F_vis_Int_own = ((   self.n_BC[0]('+') * avg(X[0].dx(0)) +
							     self.n_BC[1]('+') * avg(X[0].dx(1) )) * jump(fluc['u'][0]) +
							    (self.n_BC[0]('+') * avg(X[1].dx(0) ) +
							     self.n_BC[1]('+') * avg(X[1].dx(1))) * jump(fluc['u'][1]) -
							jump(X[0]) * (self.n_BC[0]('+') * avg(fluc['u'][0].dx(0)) +
							                      self.n_BC[1]('+') * avg(fluc['u'][0].dx(1))) -
							jump(X[1]) * (self.n_BC[0]('+') * avg(fluc['u'][1].dx(0)) +
							                      self.n_BC[1]('+') * avg(fluc['u'][1].dx(1))))
					F_vis_Int_own *= avg(1/mean['Re'])
					self.A_imag_vf += -F_vis_Int_own *dS
