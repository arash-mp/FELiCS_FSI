from ufl import (
				i,
				inner,
				dx,
				conj,
				)
def MomentumEqTimeDerivative(self,fluc,X,mean):
	self.B_vf.add(mean.rho*self.R*fluc.u[i]*conj(X)[i] *dx)
