from ufl import (
				inner,
				dx,
				grad,
				div,
				)
def EnthalpyEqTimeDerivative(self,fluc,X,mean):
	self.B_real_vf.add( mean.rho*fluc.h*X*self.R*dx-fluc.p*X*self.R*dx)
