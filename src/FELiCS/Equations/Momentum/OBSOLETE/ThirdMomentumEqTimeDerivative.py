from ufl import dx
def ThirdMomentumEqTimeDerivative(self,fluc,X,mean,param):
	self.B_real_vf.add(mean.rho*self.R*fluc.ut*X*dx)
