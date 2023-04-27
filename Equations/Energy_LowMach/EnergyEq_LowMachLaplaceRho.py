from ufl import (
				i,
				j,
				dx,
				Dx,
				inner,
				grad,
				div,
				)

def EnergyEq_LowMachLaplaceRho(self,fluc,X,mean,param):
	# term arising from integration by parts of Laplace rho_hat --> IbP
	self.A_imag_vf += mean['MuTot']/self.Pr*inner(grad(mean['rho']*self.R*X),grad(fluc['rho']))*dx
	if abs(self.m)>0:
		self.A_imag_vf += mean['MuTot']/self.Pr*mean['rho']*X*self.m**2/self.R*fluc['rho']*dx

	# term arising from laplace(rho_mean) --> no IbP
	self.A_imag_vf += -self.R*X*fluc['rho']*mean['MuTot']/self.Pr*div(grad(mean['rho']))*dx
	if param.Case.CoordinateSystem in ['Cylindrical']:
		self.A_imag_vf += -self.R*X*fluc['rho']*mean['MuTot']/self.Pr*mean['rho'].dx(1)*dx
