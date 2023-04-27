from ufl import (
				grad,
				inner,
				dx,
				)
def addSpeciesEqDiffusion(self,fluc,X,mean,specie,param):
	self.A_imag_vf.add(-mean.D(specie)*(self.R*inner(grad(fluc.Y(specie)),grad(X)))*dx) # Shouldn't the r be in the derivative???
	self.A_imag_vf.add(-fluc.D(specie)*(self.R*inner(grad(mean.Y(specie)),grad(X)))*dx)

	if param.Case.CoordinateSystem in ['Cylindrical']:
		if not param.Case.m == 0:
			printWarning("Cylindrical coordinates for m <> 0 are not validated. Treat results with care")
			self.A_imag_vf.add(mean.D*X*param.Case.m**2/self.R*fluc.Y(specie)*dx)
