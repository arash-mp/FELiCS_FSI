from ufl import (
				inner,
				dx,
				grad,
				i,
				)
def EnthalpyEqDiffusion(self,fluc,X,mean,param):
	self.A_imag_vf.add(-inner(mean.alpha*grad(fluc.h),grad(X*self.R))*dx)	#Gleichung angepasst (CA)
	self.A_imag_vf.add(-fluc.alpha*mean.he.dx(i)*X.dx(i)*dx)					#Alte Gleichung
	#self.A_imag_vf.add(-mean.alpha*fluc.viscoLaw*mean.he.dx(i)*X.dx(i)*dx)

	if param.Case.CoordinateSystem in ['Cylindrical']:
		self.A_imag_vf.add(mean.alpha*fluc.h.dx(1)*X*self.R*dx)
		if abs(self.m)>0:
			self.A_imag_vf.add(-mean.alpha*(self.m)*(self.m)/self.R*fluc.h*X*dx)

	#Debug test CA:
	#self.A_imag_vf.add(-inner(mean.alpha*mean.cp*grad(fluc.T),grad(X*self.R))*dx)	#Gleichung angepasst (CA)

	#if param.Case.CoordinateSystem in ['Cylindrical']:
	#	self.A_imag_vf.add(mean.alpha*mean.cp*fluc.T.dx(1)*X*self.R*dx)
	#	if abs(self.m)>0:
	#		self.A_imag_vf.add(-mean.alpha*mean.cp*(self.m)*(self.m)/self.R*fluc.T*X*dx)
