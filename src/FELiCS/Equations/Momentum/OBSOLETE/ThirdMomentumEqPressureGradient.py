from ufl import dx
def ThirdMomentumEqPressureGradient(self,fluc,X,mean,param):
	if not param.Case.m == 0:
		self.A_real_vf.add(+X*self.m*fluc.p*dx)
