from ufl import (
				i,
				j,
				dx,
				Dx,
				inner,
				grad,
				div,
				)
				
def tempAllTerms_MarioThesis(self,fluc,X,mean,param):
	Pr = 0.9
	self.A_imag_vf += -mean['rho']**3*X*self.R*div(fluc['u'])*dx
	self.A_imag_vf += -3*mean['rho']**2*fluc['rho']*X*self.R*div(mean['u'])*dx
	if param.Case.CoordinateSystem in ['Cylindrical']: # these terms come from divergence operator in cylinder coords
		 self.A_imag_vf += -mean['rho']**3*X*fluc['u'][1]*dx
		 self.A_imag_vf += -3*mean['rho']**2*fluc['rho']*X*mean['u'][1]*dx
		 if abs(self.m)>0: # in azimuthal direction
			 self.A_real_vf += mean['rho']**3*X*self.m*fluc['ut']*dx
	# term arising from integration by parts of nabla rho_hat
	self.A_imag_vf += 1/mean['Re']/Pr*inner(grad(mean['rho']*self.R*X),grad(fluc['rho']))*dx
	if abs(self.m)>0:
		self.A_imag_vf += 1/mean['Re']/Pr*mean['rho']*X*self.m**2/self.R*fluc['rho']*dx
	# boundary integral from integration by parts:
#		self.A_imag_vf += -1/mean['Re']/Pr*mean['rho']*self.R*X*inner(grad(fluc['rho']),self.n_BC)*self.all_ds
	# term arising from laplace(rho_mean)
	self.A_imag_vf += -self.R*X*fluc['rho']/mean['Re']/Pr*div(grad(mean['rho']))*dx
	if param.Case.CoordinateSystem in ['Cylindrical']:
		self.A_imag_vf += -self.R*X*fluc['rho']/mean['Re']/Pr*mean['rho'].dx(1)*dx
	# cross terms arising from linerization of grad(rho) dot grad(rho):
	self.A_imag_vf += 4*X*self.R/mean['Re']/Pr*inner(grad(fluc['rho']),grad(mean['rho']))*dx
