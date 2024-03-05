from ufl import (
				inner,
				dx,
				grad,
				div,
				conj,
				)
def EnthalpyEqTimeDerivative(self,fluc,X,mean):
	self.B_vf.add( mean.rho*fluc.h*conj(X)*self.R*dx-fluc.p*conj(X)*self.R*dx)
