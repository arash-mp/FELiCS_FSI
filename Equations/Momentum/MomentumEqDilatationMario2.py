from ufl import (
				i,
				j,
				dx,
				Dx,
				inner,
				grad,
				div,
				)
def MomentumEqDilatation(self,fluc,X,mean,param):
	# concerning the coordinate system, the dilatation terms mainly depend on the div u -term
	DivergenceForm_real = div(fluc.u) # cartesian portion of div(uHAT)is strictly real
	if param.Case.CoordinateSystem in ['Cylindrical']: # in cyl coods additional term to real portion
		DivergenceForm_real += fluc.u[1]/self.R
		if not param.Case.m == 0:
			DivergenceForm_imag = self.m/self.R*fluc.ut # if m>0 additional imaginary divergence portion

	# term III: +1/3 mu grad div u, becomes (dx X1 * div u) + (dr X2 * div u) = div(X) * div u due to IbP
	# however, IbP only applies to x/r, since derivative in theta is im/r
	# the underneath line incorporates cartesian as well as cylindrical, because both is covered in above DivergenceForm_i
	self.A_imag_vf.add(-self.R/3*div(X*mean.nuTot)*DivergenceForm_real*dx)
	# but when cylindrical an m > 0, then, DivergenceForm gets an imaginary part for which the same treatment as the above line applies
	if abs(self.m)>0 and param.Case.CoordinateSystem in ['Cylindrical']:
		self.A_real_vf.add(self.R/3*div(X*mean.nuTot)*DivergenceForm_imag*dx)


	# term IV: -2/3 (div u) EYE (grad mu) does not require any IbP and only applies to x/r -equation since d(mu)/d(theta) = 0
	# concerning cartesian and cylindrical, the above holds
	self.A_imag_vf.add(-2/3*self.R*inner( X , grad(mean.nuTot) )*DivergenceForm_real*dx)
	if not param.Case.m == 0:
		self.A_real_vf.add(2/3*self.R*inner( X , grad(mean.nuTot) )*DivergenceForm_imag*dx)
