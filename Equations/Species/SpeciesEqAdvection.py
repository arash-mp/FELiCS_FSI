from ufl import (
				i,
				j,
				div,
				inner,
				dx,
				Dx,
				Identity,
				)
def addSpeciesEqAdvection(self,fluc,X,mean,species,param):
	I=Identity( fluc.u.geometric_dimension() )

	self.A_imag_vf.add( fluc.Y(species)*Dx(mean.u[i]*X*self.R*mean.rho,j)*I[i,j]*dx)
	self.A_imag_vf.add( mean.Y(species)*Dx(fluc.u[i]*X*self.R*mean.rho,j)*I[i,j]*dx)
	self.A_imag_vf.add( mean.Y(species)*Dx(mean.u[i]*X*self.R*fluc.rho,j)*I[i,j]*dx)

	self.A_imag_vf.add(-self.n_BC[j]*fluc.u[i]*self.R*mean.Y(species)*X*mean.rho*I[i,j]*self.all_ds)
	self.A_imag_vf.add(-self.n_BC[j]*mean.u[i]*self.R*fluc.Y(species)*X*mean.rho*I[i,j]*self.all_ds)
	self.A_imag_vf.add(-self.n_BC[j]*mean.u[i]*self.R*mean.Y(species)*X*fluc.rho*I[i,j]*self.all_ds)

	if param.Case.CoordinateSystem in ['Cylindrical']:
		if not param.Case.m == 0:
			printWarning("Cylindrical coordinates for m <> 0 are not validated. Treat results with care")
			self.A_real_vf += X*param.Case.m*fluc.Y(i_eqn)*mean.ut*dx
