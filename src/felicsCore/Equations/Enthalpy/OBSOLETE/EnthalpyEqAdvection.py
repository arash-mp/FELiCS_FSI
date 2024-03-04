from ufl import (
				inner,
				dx,
				div,
				conj,
				)

def EnthalpyEqAdvection(self,fluc,X,mean,param):

	#Add volume integral of advection terms that remain after partial integration:
	self.A_vf.add(1j * inner(div(conj(X)*mean.rho*mean.u*self.R),conj(fluc.h))*dx) # check whether conj(fluc.h) is right
	self.A_vf.add(1j * inner(div(conj(X)*mean.rho*fluc.u*self.R),mean.he)*dx)
	self.A_vf.add(1j * inner(div(conj(X)*fluc.rho*mean.u*self.R),mean.he)*dx)

	#Add boundary integrals resulting from said partial integration:
	self.A_vf.add(1j * -inner(conj(X)*mean.rho*mean.u*fluc.h*self.R,self.n_BC) * self.all_ds)
	self.A_vf.add(1j * -inner(conj(X)*mean.rho*fluc.u*mean.he*self.R,self.n_BC) * self.all_ds)
	self.A_vf.add(1j * -inner(conj(X)*fluc.rho*mean.u*mean.he*self.R,self.n_BC) * self.all_ds)

	#If coordinates are cylindrical add extra term for m not equal to zero:
	if param.Case.CoordinateSystem in ['Cylindrical'] and abs(self.m)>0:
		self.A_vf.add( X * self.m * mean.rho * mean.ut * fluc.h * dx)