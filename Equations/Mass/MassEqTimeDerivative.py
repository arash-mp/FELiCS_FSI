from ufl import dx
def MassEqTimeDerivative(self,fluc,X,mean,param):
	self.B_vf.add( fluc.rho*self.R*X*dx)
