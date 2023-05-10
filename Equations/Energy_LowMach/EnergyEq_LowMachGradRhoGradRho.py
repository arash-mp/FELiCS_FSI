from ufl import (
				i,
				j,
				dx,
				Dx,
				inner,
				grad,
				div,
				)
				
def EnergyEq_LowMachGradRhoGradRho(self,fluc,X,mean,param):
	# cross terms arising from linerization of grad(rho) dot grad(rho):
	self.A_imag_vf += 4*X*self.R*mean['MuTot']/self.Pr*inner(grad(fluc['rho']),grad(mean['rho']))*dx
