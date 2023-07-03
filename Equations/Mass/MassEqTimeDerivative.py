from ufl import (
				dx,
				conj,
				)
def MassEqTimeDerivative(self,fluc,X,mean,param):
	self.B_vf.add( fluc.rho*self.R*conj(X)*dx)
