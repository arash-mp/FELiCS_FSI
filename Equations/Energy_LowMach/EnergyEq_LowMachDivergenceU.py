from ufl import (
				i,
				j,
				dx,
				Dx,
				inner,
				grad,
				div,
				)

def EnergyEq_LowMachDivergenceU(self,fluc,X,mean,param):
	self.A_imag_vf += -mean['rho']**3*X*self.R*div(fluc['u'])*dx
	self.A_imag_vf += -3*mean['rho']**2*fluc['rho']*X*self.R*div(mean['u'])*dx
	if param.Case.CoordinateSystem in ['Cylindrical']: # these terms come from divergence operator in cylinder coords
		 self.A_imag_vf += -mean['rho']**3*X*fluc['u'][1]*dx
		 self.A_imag_vf += -3*mean['rho']**2*fluc['rho']*X*mean['u'][1]*dx
		 if abs(self.m)>0: # in azimuthal direction
			 self.A_real_vf += mean['rho']**3*X*self.m*fluc['ut']*dx
