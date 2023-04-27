
from ufl import (
				i,
				j,
				Dx,
				dx,
				inner,
				div,
				grad,
				Identity,
				)
def ThirdMomentumEqAdvection(self,fluc,X,mean,param):

	I=Identity( fluc.u.geometric_dimension() )
	## Integration by parts is applied
	#self.A_imag_vf.add(fluc.ut*mean.rho*Dx(mean.u[i]*X*self.R,j)* I[i,j]*dx)
	#self.A_imag_vf.add(mean.ut*mean.rho*Dx(fluc.u[i]*X*self.R,j)* I[i,j]*dx)

	## Integrate on boundary terms due to integration by parts
	#self.A_imag_vf.add(fluc.ut*mean.rho*mean.u[i]*X*self.R*self.n_BC[j]* I[i,j]*self.all_ds)
	#self.A_imag_vf.add(mean.ut*mean.rho*fluc.u[i]*X*self.R*self.n_BC[j]* I[i,j]*self.all_ds)

	# Source terms due to cylindrical coordinates
	self.A_imag_vf.add(-X*mean.rho*mean.u[1]*fluc.ut*dx)
	self.A_imag_vf.add(-X*mean.rho*mean.ut*fluc.u[1]*dx)
	if not param.Case.m == 0:
		self.A_real_vf.add(X*mean.rho*param.Case.m*mean.ut*fluc.ut*dx)

		# in case of non-uniform density, additional convective terms occur
		#################### rhoHAT (uMEAN dot div) uMEAN ###############################
	if 'rho' in param.Case.getTransportedQuantityList():
		self.A_imag_vf.add(-fluc.rho*X*self.R*inner(mean.u,grad(mean.ut))*dx)
		self.A_imag_vf.add(-fluc.rho*X*mean.ut*mean.u[1]*dx)
