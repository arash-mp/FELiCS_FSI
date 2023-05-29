from ufl import (
				i,
				inner,
				dx,
				)
def MomentumEqTimeDerivative(self,fluc,X,mean):
	self.B_vf.add(mean.rho*self.R*fluc.u[i]*X[i] *dx)
