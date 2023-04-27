from ufl import (
				div,
				dx,
				)
def ThirdMomentumEqDilatation(self,fluc,X,mean,param):
	# concerning the coordinate system, the dilatation terms mainly depend on the div u -term (same as in addMomentumEquationDilatation)

	'''
<<<<<<< HEAD
	DivergenceForm_real = div(fluc.u) # cartesian portion of div(uHAT)is strictly real
	if param.CoordinateSystem in ['Cylindrical']: # in cyl coods additional term to real portion
		DivergenceForm_real += fluc.u[1]/self.R
		if not param.m == 0:
			DivergenceForm_imag = self.m/self.R*fluc.ut # if m>0 additional imaginary divergence portion
=======
	'''
	DivergenceForm_real = div(fluc.u[0]) # cartesian portion of div(uHAT)is strictly real
	if param.Case.CoordinateSystem in ['Cylindrical']: # in cyl coods additional term to real portion
		DivergenceForm_real += fluc.u[1]/self.R
		if not param.Case.m == 0:
			DivergenceForm_imag = self.m/self.R*self.hat[1] # if m>0 additional imaginary divergence portion
#>>>> master

	# term III: +1/3 mu grad div u, no IbP required for theta-equation since leading derivative (grad) is given by im/r
	# the underneath line incorporates cartesian as well as cylindrical, because both is covered in above DivergenceForm_i
	# due to the derivative im/r in theta direction, the dilation term in theta-momentum equation only applies of m>0
	if abs(self.m)>0:
		self.A_imag_vf.add(-1./3.*mean.nuTot*self.m*X*DivergenceForm_imag*dx)
		self.A_real_vf.add(-1./3.*mean.nuTot*self.m*X*DivergenceForm_real*dx)
